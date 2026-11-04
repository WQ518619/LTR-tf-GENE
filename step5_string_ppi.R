# =========================================================
# Step 5: Run STRINGdb PPI Analysis
# File: R/run_string_ppi.R
# =========================================================
#' Run STRINGdb PPI Analysis
#'
#' @description
#' Identify protein-protein associations between
#' LTR-associated TFs and gene-associated TFs
#' using the STRINGdb database.
#'
#' Only cross-group associations are retained.
#' The PPI associations are undirected.
#'
#' For downstream network visualization:
#'   from = LTR-associated TF (TF1)
#'   to   = Gene-associated TF (TF2)
#'
#' @param ltr_tfs Character vector. LTR-associated TFs
#'   identified in Step 3.
#' @param gene_tfs Character vector. Gene-associated TFs
#'   identified in Step 6.
#' @param species_id Numeric. NCBI taxonomy ID:
#'   9606 for human, 10090 for mouse.
#' @param score_threshold Numeric. STRING combined score
#'   threshold. Default is 400.
#' @param outdir Character. Directory to save results.
#'
#' @return A data frame containing:
#'   from, to, combined_score.
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
run_string_ppi <- function(
    ltr_tfs,
    gene_tfs,
    species_id = 9606,
    score_threshold = 400,
    outdir = NULL
) {
  message("【Step 5】 Running STRINGdb PPI analysis...")
  # =====================================================
  # 0. 定义空结果及保存函数
  # =====================================================
  empty_ppi <- data.frame(
    from = character(),
    to = character(),
    combined_score = numeric(),
    stringsAsFactors = FALSE
  )
  save_results <- function(ppi_df) {
    if (!is.null(outdir)) {
      dir.create(
        outdir,
        recursive = TRUE,
        showWarnings = FALSE
      )
      utils::write.csv(
        ppi_df,
        file.path(
          outdir,
          "05_PPI_network.csv"
        ),
        row.names = FALSE
      )
    }
    return(ppi_df)
  }
  # =====================================================
  # 1. 清洗并合并两侧 TF
  # =====================================================
  ltr_tfs <- unique(
    trimws(as.character(ltr_tfs))
  )
  gene_tfs <- unique(
    trimws(as.character(gene_tfs))
  )
  ltr_tfs <- ltr_tfs[
    !is.na(ltr_tfs) & nzchar(ltr_tfs)
  ]
  gene_tfs <- gene_tfs[
    !is.na(gene_tfs) & nzchar(gene_tfs)
  ]
  # 剔除明显的表观遗传标签及非 TF 项目
  epi_pattern <- paste0(
    "(?i)acetyl|H[234]K[0-9]+|",
    "Epitope|IgG|Input"
  )
  ltr_tfs <- ltr_tfs[
    !grepl(epi_pattern, ltr_tfs, perl = TRUE)
  ]
  gene_tfs <- gene_tfs[
    !grepl(epi_pattern, gene_tfs, perl = TRUE)
  ]
  if (length(ltr_tfs) == 0 ||
      length(gene_tfs) == 0) {
    warning(
      "LTR-side or gene-side TF list is empty."
    )
    return(save_results(empty_ppi))
  }
  # 合并全部 TF，用于 STRINGdb 查询
  all_tfs <- unique(
    c(ltr_tfs, gene_tfs)
  )
  if (length(all_tfs) < 2) {
    warning(
      "Less than two unique TFs. ",
      "No cross-group PPI can be identified."
    )
    return(save_results(empty_ppi))
  }
  message(
    "  -> LTR-side TFs: ",
    length(ltr_tfs)
  )
  message(
    "  -> Gene-side TFs: ",
    length(gene_tfs)
  )
  message(
    "  -> Total unique TFs: ",
    length(all_tfs)
  )
  # =====================================================
  # 2. 初始化 STRINGdb
  # =====================================================
  message(
    "  -> Initializing STRINGdb (Species: ",
    species_id,
    ")..."
  )
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop(
      "Please install the 'STRINGdb' package."
    )
  }
  string_db <- STRINGdb::STRINGdb$new(
    species = species_id,
    score_threshold = score_threshold,
    input_directory = tempdir()
  )
  # =====================================================
  # 3. 基因名映射为 STRING ID
  # =====================================================
  tf_df <- data.frame(
    genes = all_tfs,
    stringsAsFactors = FALSE
  )
  message(
    "  -> Mapping ",
    nrow(tf_df),
    " TFs to STRING IDs..."
  )
  mapped <- string_db$map(
    tf_df,
    "genes",
    removeUnmappedRows = TRUE
  )
  if (nrow(mapped) == 0) {
    warning(
      "STRINGdb mapping failed for all TFs."
    )
    return(save_results(empty_ppi))
  }
  message(
    "  -> Successfully mapped ",
    length(unique(mapped$genes)),
    " TF names."
  )
  # =====================================================
  # 4. 获取 STRING 蛋白关联网络
  # =====================================================
  message(
    "  -> Fetching interactions (score >= ",
    score_threshold,
    ")..."
  )
  ppi_raw <- string_db$get_interactions(
    unique(mapped$STRING_id)
  )
  if (nrow(ppi_raw) == 0) {
    warning(
      "No STRING interactions found ",
      "at the current score threshold."
    )
    return(save_results(empty_ppi))
  }
  # =====================================================
  # 5. 还原 TF 名称并筛选跨组关联
  # =====================================================
  message(
    "  -> Identifying cross-group TF interactions..."
  )
  # 5.1 STRING ID 与原始 TF 名称对应关系
  id_map <- mapped %>%
    dplyr::select(
      STRING_id,
      genes
    ) %>%
    dplyr::distinct()
  # 5.2 将 STRING ID 还原为 TF 名称
  ppi_named <- ppi_raw %>%
    dplyr::inner_join(
      id_map %>%
        dplyr::rename(
          from = STRING_id,
          from_gene = genes
        ),
      by = "from"
    ) %>%
    dplyr::inner_join(
      id_map %>%
        dplyr::rename(
          to = STRING_id,
          to_gene = genes
        ),
      by = "to"
    )
  # 5.3 筛选 LTR-TF 与 Gene-TF 之间的关联
  #
  # from = LTR 侧 TF1
  # to   = Gene 侧 TF2
  # 同时考虑 STRING 返回边的两种端点排列。
  cross_ppi <- dplyr::bind_rows(
    # 情况 1：
    # STRING from 属于 LTR 侧
    # STRING to 属于基因侧
    ppi_named %>%
      dplyr::filter(
        from_gene %in% ltr_tfs,
        to_gene %in% gene_tfs,
        from_gene != to_gene
      ) %>%
      dplyr::transmute(
        from = from_gene,
        to = to_gene,
        combined_score = combined_score
      ),
    # 情况 2：
    # STRING to 属于 LTR 侧
    # STRING from 属于基因侧
    ppi_named %>%
      dplyr::filter(
        to_gene %in% ltr_tfs,
        from_gene %in% gene_tfs,
        from_gene != to_gene
      ) %>%
      dplyr::transmute(
        from = to_gene,
        to = from_gene,
        combined_score = combined_score
      )
  ) %>%
    dplyr::group_by(from, to) %>%
    dplyr::summarise(
      combined_score = max(
        combined_score,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(combined_score),
      from,
      to
    )
  if (nrow(cross_ppi) == 0) {
    warning(
      "No cross-group TF interactions found."
    )
    return(save_results(empty_ppi))
  }
  message(
    "  -> Found ",
    nrow(cross_ppi),
    " cross-group PPI connections."
  )
  message(
    "  -> Involved ",
    length(unique(cross_ppi$from)),
    " LTR-side TFs and ",
    length(unique(cross_ppi$to)),
    " gene-side TFs."
  )
  # =====================================================
  # 6. 保存并返回结果
  # =====================================================
  return(
    save_results(cross_ppi)
  )
}