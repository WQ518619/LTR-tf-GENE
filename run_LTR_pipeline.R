# =========================================================
# LTRnet: End-to-end seven-step pipeline
# File: R/run_LTR_pipeline.R
# =========================================================

#' Run the End-to-End LTR-Associated Regulatory Network Pipeline
#'
#' @description
#' Runs the seven-module LTRnet workflow. Steps 1-4 are shared by
#' forward and reverse modes. STRING PPI (Step 5) is run only in
#' reverse mode. Step 6 retains the complete network; Step 7
#' selects a smaller network for plotting only.
#'
#' @param tefile Path to TEtranscripts differential expression results.
#' @param annotation_file Path to TE annotation (.gtf or .rds). If NULL,
#'   download/load the package-provided annotation for `species`.
#' @param chipfile Path to ChIP-Atlas peaks (.rds or .bed). If NULL,
#'   load/download tissue-specific data.
#' @param tissue Tissue used when `chipfile` is NULL.
#' @param degfile Path to host-gene differential expression results;
#'   required for reverse mode.
#' @param mode Either "forward" or "reverse".
#' @param species Either "human" or "mouse".
#' @param outdir Root output directory.
#' @param custom_config List overriding `step1`, `step2`, `step3`,
#'   `step4`, `step5`, or `step7` default settings.
#' @param cache_dir Optional directory for downloaded reference files.
#'
#' @return Invisibly returns a list containing `nodes`, `edges`,
#'   `plot`, `out_inter`, and `out_final`.
#' @importFrom utils modifyList
#' @export
run_LTR_pipeline <- function(
    tefile,
    annotation_file = NULL,
    chipfile = NULL,
    tissue = NULL,
    degfile = NULL,
    mode = c("forward", "reverse"),
    species = c("human", "mouse"),
    outdir = "LTRnet_Results",
    custom_config = list(),
    cache_dir = NULL
) {
  mode <- match.arg(mode)
  species <- match.arg(species)

  # =====================================================
  # 0. Validate inputs and resolve reference data
  # =====================================================
  if (!is.character(tefile) || length(tefile) != 1L ||
      is.na(tefile) || !nzchar(tefile) || !file.exists(tefile)) {
    stop("Input TEtranscripts file 'tefile' does not exist.")
  }

  if (mode == "reverse" &&
      (is.null(degfile) || length(degfile) != 1L ||
       is.na(degfile) || !nzchar(degfile) || !file.exists(degfile))) {
    stop("Reverse mode requires an existing DEG file ('degfile').")
  }

  if (!is.list(custom_config)) {
    stop("'custom_config' must be a list.")
  }

  if (is.null(annotation_file)) {
    message("  -> No annotation supplied; retrieving annotation for ", species, "...")
    annotation_file <- load_rmsk_gtf(
      species = species,
      local_gtf_path = NULL,
      cache_dir = cache_dir
    )
  }
  if (!is.character(annotation_file) || length(annotation_file) != 1L ||
      is.na(annotation_file) || !file.exists(annotation_file)) {
    stop("LTR annotation file was not found: ", paste(annotation_file, collapse = ", "))
  }

  if (is.null(chipfile)) {
    if (is.null(tissue) || length(tissue) != 1L ||
        is.na(tissue) || !nzchar(tissue)) {
      stop("Provide either 'chipfile' or 'tissue' (for example, 'Liver').")
    }
    message("  -> No ChIP-Atlas file supplied; retrieving ",
            species, " / ", tissue, "...")
    chipfile <- load_chip_atlas_data(
      tissue = tissue,
      species = species,
      cache_dir = cache_dir,
      local_data_dir = NULL
    )
  }
  if (!is.character(chipfile) || length(chipfile) != 1L ||
      is.na(chipfile) || !file.exists(chipfile)) {
    stop("ChIP-Atlas file was not found: ", paste(chipfile, collapse = ", "))
  }

  txdb <- if (species == "human") {
    "TxDb.Hsapiens.UCSC.hg38.knownGene"
  } else {
    "TxDb.Mmusculus.UCSC.mm10.knownGene"
  }
  taxonomy_id <- if (species == "human") 9606 else 10090

  # =====================================================
  # Configure the seven-module pipeline
  # =====================================================
  default_cfg <- list(
    step1 = list(lfc = 0.585, padj_thresh = 0.05),
    step2 = list(qval = 0.05, support = 1,
                 allowed_tiers = c("A", "B")),
    step3 = list(window = 100, lfc = 0.585, padj = 0.05),
    step4 = list(dorothea_conf = c("A", "B", "C", "D")),
    step5 = list(score_threshold = 400),
    step7 = list(width = 14, height = 12,
                 top_n_ltr = 5, top_n_tf = 8,
                 top_n_gene_tf = 8, top_n_genes = 3,
                 label_size = 3.1)
  )
  cfg <- utils::modifyList(default_cfg, custom_config)

  dir_inter <- file.path(outdir, "1_Intermediate_Data", mode)
  dir_final <- file.path(outdir, "2_Final_Outputs", mode)
  dir.create(dir_inter, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir_final, recursive = TRUE, showWarnings = FALSE)

  message("\n=======================================================")
  message("LTRnet pipeline started (", mode, " | ", species, ")")
  message("=======================================================")

  # =====================================================
  # Step 1: Upregulated LTR identification and GRanges
  # =====================================================
  message("\n>>> [Step 1/7] Identifying LTRs and genomic regions...")
  step1 <- get_ltr_regions(
    tefile = tefile,
    annotation_file = annotation_file,
    lfc = cfg$step1$lfc,
    padj_thresh = cfg$step1$padj_thresh,
    output_dir = dir_inter
  )
  if (is.null(step1) || length(step1$ltr_names) == 0L ||
      length(step1$query_gr) == 0L || length(step1$universe_gr) == 0L) {
    stop("[Step 1] No valid LTR regions or universe regions found.")
  }

  # =====================================================
  # Step 2: LOLA enrichment of LTR-associated TFs
  # =====================================================
  message("\n>>> [Step 2/7] Running LOLA enrichment...")
  step2 <- run_lola_enrichment(
    query_gr = step1$query_gr,
    universe_gr = step1$universe_gr,
    chipfile = chipfile,
    loladir = file.path(dir_inter, "LOLA_DB"),
    qval = cfg$step2$qval,
    support = cfg$step2$support,
    allowed_tiers = cfg$step2$allowed_tiers,
    outdir = dir_inter
  )
  if (is.null(step2) || is.null(step2$sig_tf) ||
      nrow(step2$sig_tf) == 0L) {
    stop("[Step 2] No LTR-associated TFs passed the selection criteria.")
  }

  # =====================================================
  # Step 3: Nearby gene identification
  # =====================================================
  message("\n>>> [Step 3/7] Identifying nearby candidate genes...")
  step3 <- get_nearby_genes(
    query_gr = step1$query_gr,
    valid_ltr_names = unique(step2$sig_tf$LTR),
    degfile = degfile,
    mode = mode,
    txdb = txdb,
    window = cfg$step3$window,
    lfc = cfg$step3$lfc,
    padj = cfg$step3$padj,
    outdir = dir_inter
  )
  if (is.null(step3) || is.null(step3$final_candidates) ||
      nrow(step3$final_candidates) == 0L) {
    stop("[Step 3] No candidate genes found.")
  }

  # =====================================================
  # Step 4: Mode-specific DoRothEA analysis
  # =====================================================
  message("\n>>> [Step 4/7] Running DoRothEA (", mode, ")...")
  step4 <- run_dorothea_analysis(
    nearby_genes = step3$final_candidates,
    sig_tf = if (mode == "forward") step2$sig_tf else NULL,
    mode = mode,
    species = species,
    dorothea_conf = cfg$step4$dorothea_conf,
    outdir = dir_inter
  )

  # =====================================================
  # Steps 5-6: STRING (reverse only), complete network
  # =====================================================
  if (mode == "forward") {
    if (is.null(step4) || is.null(step4$forward_edges) ||
        nrow(step4$forward_edges) == 0L) {
      stop("[Step 4] No forward LTR-TF-Gene paths were identified.")
    }

    message("\n>>> [Step 5/7] STRING PPI skipped in forward mode.")
    message("\n>>> [Step 6/7] Assembling complete forward network...")
    step6 <- build_network_tables(
      mode = "forward",
      sig_tf = step2$sig_tf,
      step_res = step4,
      outdir = dir_inter
    )
  } else {
    if (is.null(step4) || is.null(step4$reverse_edges) ||
        nrow(step4$reverse_edges) == 0L ||
        length(step4$sig_gene_tfs) == 0L) {
      stop("[Step 4] No reverse TF-Gene associations were identified.")
    }

    message("\n>>> [Step 5/7] Running STRING cross-group PPI analysis...")
    step5 <- run_string_ppi(
      ltr_tfs = step2$ltr_tfs,
      gene_tfs = step4$sig_gene_tfs,
      species_id = taxonomy_id,
      score_threshold = cfg$step5$score_threshold,
      outdir = dir_inter
    )
    if (is.null(step5) || nrow(step5) == 0L) {
      stop("[Step 5] No cross-group STRING associations found.")
    }

    message("\n>>> [Step 6/7] Assembling complete reverse network...")
    step6 <- build_network_tables(
      mode = "reverse",
      sig_tf = step2$sig_tf,
      step_res = step4,
      cross_ppi = step5,
      outdir = dir_inter
    )
  }

  if (is.null(step6) || is.null(step6$nodes) ||
      is.null(step6$edges) || nrow(step6$nodes) == 0L ||
      nrow(step6$edges) == 0L) {
    stop("[Step 6] No complete candidate network could be assembled.")
  }

  # =====================================================
  # Step 7: Simplified publication-ready visualization
  # =====================================================
  message("\n>>> [Step 7/7] Plotting selected network...")
  p <- plot_ltr_tf_gene_network(
    nodes = step6$nodes,
    edges = step6$edges,
    mode = mode,
    forward_paths = if (mode == "forward") step4$forward_edges else NULL,
    outdir = dir_final,
    width = cfg$step7$width,
    height = cfg$step7$height,
    top_n_ltr = cfg$step7$top_n_ltr,
    top_n_tf = cfg$step7$top_n_tf,
    top_n_gene_tf = cfg$step7$top_n_gene_tf,
    top_n_genes = cfg$step7$top_n_genes,
    label_size = cfg$step7$label_size
  )
  if (is.null(p)) {
    warning("[Step 7] No plot was generated; complete network CSVs are available.")
  }

  message("\n=======================================================")
  message("LTRnet pipeline finished (", mode, ").")
  message("Intermediate data: ", dir_inter)
  message("Figures and selected network: ", dir_final)
  message("=======================================================\n")

  invisible(list(
    nodes = step6$nodes,
    edges = step6$edges,
    plot = p,
    out_inter = dir_inter,
    out_final = dir_final
  ))
}
