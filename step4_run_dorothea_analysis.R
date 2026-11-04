# =========================================================
# Combined DoRothEA Analysis: Forward and Reverse
# File: R/run_dorothea_analysis.R
# =========================================================

#' Run DoRothEA Forward or Reverse Regulatory Analysis
#'
#' @description
#' Use DoRothEA TF-target regulatory interactions for two modes.
#' In forward mode, retain LTR-TF-Gene triples supported by
#' the LTR-TF pairs from LOLA and LTR-nearby genes from Step 3.
#' In reverse mode, identify DoRothEA TF regulators of candidate
#' genes from Step 3; no LTR-family-specific gene matching is made.
#'
#' @param nearby_genes Data frame from Step 3 containing `LTR` and
#'   `target_gene` (forward) or `target_gene` (reverse).
#' @param sig_tf Data frame from Step 3 containing `LTR` and `TF`;
#'   required only for forward mode.
#' @param mode Either `"forward"` or `"reverse"`.
#' @param species Either `"human"` or `"mouse"`.
#' @param dorothea_conf DoRothEA confidence levels to retain.
#' @param outdir Optional output directory for CSV files.
#'
#' @return In forward mode, a list with `forward_edges` and
#'   `target_summary`; in reverse mode, a list with `reverse_edges`
#'   and `sig_gene_tfs`.
#' @importFrom magrittr %>%
#' @export
run_dorothea_analysis <- function(
    nearby_genes,
    sig_tf = NULL,
    mode = c("forward", "reverse"),
    species = c("human", "mouse"),
    dorothea_conf = c("A", "B", "C", "D"),
    outdir = NULL
) {
  mode <- match.arg(mode)
  species <- match.arg(species)

  message("【DoRothEA】 Running ", mode, " regulatory analysis...")

  # =====================================================
  # 1. Define mode-specific empty results and output
  # =====================================================
  empty_forward_edges <- data.frame(
    LTR = character(), TF = character(), Gene = character(),
    confidence = character(), stringsAsFactors = FALSE
  )
  empty_target_summary <- data.frame(
    Gene = character(),
    Regulating_LTR_TFs_Count = integer(),
    Regulating_LTR_TFs = character(),
    stringsAsFactors = FALSE
  )
  empty_reverse_edges <- data.frame(
    Gene = character(), TF = character(), confidence = character(),
    stringsAsFactors = FALSE
  )

  save_forward <- function(forward_edges, target_summary) {
    if (!is.null(outdir)) {
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(
        forward_edges,
        file.path(outdir, "04_forward_network_edges.csv"),
        row.names = FALSE
      )
      utils::write.csv(
        target_summary,
        file.path(outdir, "04_forward_core_targets.csv"),
        row.names = FALSE
      )
    }
    list(forward_edges = forward_edges, target_summary = target_summary)
  }

  save_reverse <- function(reverse_edges, sig_gene_tfs) {
    if (!is.null(outdir)) {
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(
        reverse_edges,
        file.path(outdir, "04_reverse_network_edges.csv"),
        row.names = FALSE
      )
    }
    list(reverse_edges = reverse_edges, sig_gene_tfs = sig_gene_tfs)
  }

  # =====================================================
  # 2. Prepare mode-specific input
  # =====================================================
  if (mode == "forward") {
    if (is.null(sig_tf) || nrow(sig_tf) == 0) {
      warning("No LTR-TF pairs provided.")
      return(save_forward(empty_forward_edges, empty_target_summary))
    }
    if (is.null(nearby_genes) || nrow(nearby_genes) == 0) {
      warning("No nearby genes provided.")
      return(save_forward(empty_forward_edges, empty_target_summary))
    }
    if (!all(c("LTR", "TF") %in% names(sig_tf)) ||
        !all(c("LTR", "target_gene") %in% names(nearby_genes))) {
      stop("Forward mode requires sig_tf: LTR, TF and nearby_genes: LTR, target_gene.")
    }

    ltr_tf_pairs <- sig_tf %>%
      dplyr::select(LTR, TF) %>%
      dplyr::mutate(
        LTR = trimws(as.character(LTR)),
        TF = trimws(as.character(TF))
      ) %>%
      dplyr::filter(!is.na(LTR), !is.na(TF), LTR != "", TF != "") %>%
      dplyr::distinct()

    nearby_pairs <- nearby_genes %>%
      dplyr::select(LTR, target_gene) %>%
      dplyr::mutate(
        LTR = trimws(as.character(LTR)),
        target_gene = trimws(as.character(target_gene))
      ) %>%
      dplyr::filter(
        !is.na(LTR), !is.na(target_gene),
        LTR != "", target_gene != ""
      ) %>%
      dplyr::distinct()

    if (nrow(ltr_tf_pairs) == 0 || nrow(nearby_pairs) == 0) {
      warning("No valid LTR-TF or LTR-nearby-gene pairs remain.")
      return(save_forward(empty_forward_edges, empty_target_summary))
    }
    target_tfs <- unique(ltr_tf_pairs$TF)

  } else {
    if (is.null(nearby_genes) || nrow(nearby_genes) == 0) {
      warning("No nearby genes provided.")
      return(save_reverse(empty_reverse_edges, character(0)))
    }
    if (!"target_gene" %in% names(nearby_genes)) {
      stop("Reverse mode requires nearby_genes column: target_gene.")
    }
    target_genes <- unique(trimws(as.character(nearby_genes$target_gene)))
    target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]
    if (length(target_genes) == 0) {
      return(save_reverse(empty_reverse_edges, character(0)))
    }
  }

  # =====================================================
  # 3. Load DoRothEA once for the selected mode
  # =====================================================
  if (!requireNamespace("dorothea", quietly = TRUE)) {
    stop("Please install the 'dorothea' package.")
  }
  db_name <- if (species == "human") "dorothea_hs" else "dorothea_mm"
  db_env <- new.env(parent = baseenv())
  utils::data(list = db_name, package = "dorothea", envir = db_env)
  if (!exists(db_name, envir = db_env, inherits = FALSE)) {
    stop("DoRothEA dataset not found: ", db_name)
  }
  regulon_db <- get(db_name, envir = db_env, inherits = FALSE)

  high_conf_regulon <- regulon_db %>%
    dplyr::filter(confidence %in% dorothea_conf) %>%
    dplyr::transmute(
      TF = as.character(tf),
      Gene = as.character(target),
      confidence = as.character(confidence)
    ) %>%
    dplyr::distinct()

  message("  -> Using DoRothEA confidence levels: ",
          paste(dorothea_conf, collapse = ", "))

  # =====================================================
  # 4A. Forward: exact LTR-TF and LTR-nearby-gene matching
  # =====================================================
  if (mode == "forward") {
    all_edges <- high_conf_regulon %>%
      dplyr::filter(TF %in% target_tfs) %>%
      dplyr::distinct()

    if (nrow(all_edges) == 0) {
      warning("None of the candidate TFs were found in the selected DoRothEA regulon.")
      return(save_forward(empty_forward_edges, empty_target_summary))
    }

    forward_edges <- ltr_tf_pairs %>%
      dplyr::inner_join(
        all_edges,
        by = "TF"
     ) %>%
      dplyr::filter(
        Gene %in% unique(nearby_pairs$target_gene)
     ) %>%
      dplyr::select(
        LTR, TF, Gene, confidence
     ) %>%
      dplyr::distinct()

    if (nrow(forward_edges) == 0) {
      warning("No LTR-TF-Gene connections passed spatial and regulatory filters.")
      return(save_forward(empty_forward_edges, empty_target_summary))
    }

    target_summary <- forward_edges %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(
        Regulating_LTR_TFs_Count = dplyr::n_distinct(TF),
        Regulating_LTR_TFs = paste(sort(unique(TF)), collapse = "; "),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(Regulating_LTR_TFs_Count))

    message("  -> Found ", nrow(forward_edges), " LTR-TF-Gene connections.")
    message("  -> Involved ", dplyr::n_distinct(forward_edges$Gene),
            " unique target genes.")
    return(save_forward(forward_edges, target_summary))
  }

  # =====================================================
  # 4B. Reverse: find upstream TFs of candidate genes
  # =====================================================
  reverse_edges <- high_conf_regulon %>%
    dplyr::filter(Gene %in% target_genes) %>%
    dplyr::select(Gene, TF, confidence) %>%
    dplyr::distinct()

  if (nrow(reverse_edges) == 0) {
    warning("None of the candidate genes have known regulating TFs in DoRothEA.")
    return(save_reverse(empty_reverse_edges, character(0)))
  }

  sig_gene_tfs <- unique(reverse_edges$TF)
  message("  -> Found ", nrow(reverse_edges), " Gene-TF connections.")
  message("  -> Identified ", length(sig_gene_tfs), " upstream TFs.")
  return(save_reverse(reverse_edges, sig_gene_tfs))
}
