# =========================================================
# Step 1: Identify Upregulated LTRs and Build GRanges
# File: R/get_ltr_regions.R
# =========================================================
#' Identify Upregulated LTRs and Construct Genomic Regions
#'
#' @description
#' Identify significantly upregulated LTRs from TEtranscripts
#' differential expression results, map LTR families to
#' RepeatMasker genomic annotations, and construct query
#' and restricted universe regions for LOLA analysis.
#'
#' @param tefile Path to TEtranscripts differential expression results.
#' @param annotation_file Path to TE annotation (.gtf or .rds).
#' @param lfc Log2 fold-change threshold. Default 0.585.
#' @param padj_thresh Adjusted p-value threshold. Default 0.05.
#' @param output_dir Directory for intermediate results.
#'
#' @return A list containing ltr_df, ltr_names,
#'   query_gr, and universe_gr.
#'
#' @import dplyr stringr GenomicRanges IRanges
#' @importFrom magrittr %>%
#' @export
get_ltr_regions <- function(
    tefile,
    annotation_file,
    lfc = 0.585,
    padj_thresh = 0.05,
    output_dir = NULL
) {
  message("【Step 1】 Identifying LTRs and building genomic regions...")
  # =====================================================
  # 1. 读取 TEtranscripts 差异表达结果
  # =====================================================
  if (!file.exists(tefile)) {
    stop("TEtranscripts file not found: ", tefile)
  }
  if (!file.exists(annotation_file)) {
    stop("Annotation file not found: ", annotation_file)
  }
  if (grepl("\\.xlsx$|\\.xls$", tefile, ignore.case = TRUE)) {
    te_df <- readxl::read_excel(tefile)
  } else {
    te_df <- data.table::fread(
      tefile,
      data.table = FALSE,
      showProgress = FALSE
    )
  }
  # 检查必须列
  required_cols <- c("ID", "log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(te_df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  # =====================================================
  # 2. 筛选显著上调 LTR
  # =====================================================
  ltr_filtered <- te_df %>%
    dplyr::mutate(
      log2FoldChange = suppressWarnings(
        as.numeric(log2FoldChange)
      ),
      padj = suppressWarnings(
        as.numeric(padj)
      )
    ) %>%
    dplyr::filter(
      stringr::str_detect(ID, "LTR"),
      !is.na(log2FoldChange),
      log2FoldChange > lfc,
      !is.na(padj),
      padj < padj_thresh
    ) %>%
    dplyr::mutate(
      name = stringr::str_extract(ID, "^[^:]+")
    )
  ltr_names <- unique(ltr_filtered$name)
  ltr_names <- ltr_names[
    !is.na(ltr_names) & ltr_names != ""
  ]
  if (length(ltr_names) == 0) {
    stop("No significantly upregulated LTRs found.")
  }
  message(
    "  -> Identified ",
    length(ltr_names),
    " upregulated LTR families."
  )
  # =====================================================
  # 3. 读取 RepeatMasker 注释
  # =====================================================
  if (grepl("\\.rds$", annotation_file, ignore.case = TRUE)) {
    te_parsed <- readRDS(annotation_file)
    required_annotation_cols <- c(
      "chr", "start", "end",
      "strand", "gene_id", "class_id"
    )
    missing_annotation_cols <- setdiff(
      required_annotation_cols,
      colnames(te_parsed)
    )
    if (length(missing_annotation_cols) > 0) {
      stop(
        "Annotation RDS missing columns: ",
        paste(missing_annotation_cols, collapse = ", ")
      )
    }
  } else {
    te_raw <- data.table::fread(
      annotation_file,
      sep = "\t",
      header = FALSE,
      data.table = FALSE,
      showProgress = FALSE
    )
    if (ncol(te_raw) < 9) {
      stop("Invalid GTF annotation file.")
    }
    colnames(te_raw) <- c(
      "chr", "source", "feature",
      "start", "end", "score",
      "strand", "frame", "attribute"
    )
    te_parsed <- te_raw %>%
      dplyr::mutate(
        gene_id = stringr::str_match(
          attribute,
          'gene_id "?([^";]+)"?'
        )[, 2],
        class_id = stringr::str_match(
          attribute,
          'class_id "?([^";]+)"?'
        )[, 2]
      )
  }
  # =====================================================
  # 4. 提取目标 LTR 坐标
  # =====================================================
  query_regions <- te_parsed %>%
    dplyr::filter(
      class_id == "LTR",
      gene_id %in% ltr_names,
      !is.na(gene_id)
    )
  if (nrow(query_regions) == 0) {
    stop("No matching LTR regions found in annotation.")
  }
  query_gr <- GenomicRanges::GRanges(
    seqnames = query_regions$chr,
    ranges = IRanges::IRanges(
      start = query_regions$start,
      end = query_regions$end
    ),
    strand = query_regions$strand,
    name = query_regions$gene_id
  )
  GenomeInfoDb::seqlevelsStyle(query_gr) <- "UCSC"
  query_gr <- GenomeInfoDb::keepStandardChromosomes(
    query_gr,
    pruning.mode = "coarse"
  )
  query_gr <- unique(query_gr)
  if (length(query_gr) == 0) {
    stop("No LTR regions remain after chromosome filtering.")
  }
  message(
    "  -> Target regions: ",
    length(query_gr)
  )
  # =====================================================
  # 5. 构建 LOLA restricted universe
  # =====================================================
  user_sets <- split(
    query_gr,
    S4Vectors::mcols(query_gr)$name
  )
  universe_gr <- LOLA::buildRestrictedUniverse(
    user_sets
  )
  universe_gr <- GenomeInfoDb::keepStandardChromosomes(
    universe_gr,
    pruning.mode = "coarse"
  )
  universe_gr <- sort(universe_gr)
  message(
    "  -> Background regions: ",
    length(universe_gr)
  )
  # =====================================================
  # 6. 保存结果
  # =====================================================
  if (!is.null(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    utils::write.csv(
      ltr_filtered,
      file.path(
        output_dir,
        "01_upregulated_LTR_details.csv"
      ),
      row.names = FALSE
    )
    saveRDS(
      ltr_names,
      file.path(
        output_dir,
        "01_upregulated_LTR_names.rds"
      )
    )
    saveRDS(
      query_gr,
      file.path(
        output_dir,
        "01_query_gr.rds"
      )
    )
    saveRDS(
      universe_gr,
      file.path(
        output_dir,
        "01_universe_gr.rds"
      )
    )
    utils::write.csv(
      as.data.frame(query_gr),
      file.path(
        output_dir,
        "01_query_gr_coordinates.csv"
      ),
      row.names = FALSE
    )
  }
  return(list(
    ltr_df = ltr_filtered,
    ltr_names = ltr_names,
    query_gr = query_gr,
    universe_gr = universe_gr
  ))
}