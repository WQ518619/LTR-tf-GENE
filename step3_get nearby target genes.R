# =========================================================
# Step 3: Get Nearby Target Genes
# File: R/get_nearby_genes.R
# =========================================================

#' Get Nearby Target Genes
#'
#' @description
#' Identify host genes whose transcription start sites (TSS)
#' are located within the specified flanking window of LTRs.
#'
#' Reverse mode: retain only significantly upregulated DEGs
#' among nearby genes.
#'
#' Forward mode: retain all nearby genes without DEG filtering.
#'
#' @param query_gr GRanges. Target LTR regions.
#' @param valid_ltr_names Character vector. LTR families with
#'   selected candidate TFs from Step 3.
#' @param degfile Character. DEG file required in reverse mode.
#' @param mode Character. "reverse" or "forward".
#' @param txdb Character. TxDb package name.
#' @param window Numeric. Flanking window in kb. Default 100.
#' @param lfc Numeric. Minimum positive log2FC. Default 0.585.
#' @param padj Numeric. Adjusted P-value threshold. Default 0.05.
#' @param outdir Character. Output directory.
#'
#' @return A list containing final_candidates and all_nearby.
#'
#' @import dplyr GenomicRanges IRanges GenomicFeatures
#' @importFrom magrittr %>%
#' @importFrom S4Vectors mcols queryHits subjectHits
#' @importFrom readxl read_excel
#' @importFrom data.table fread
#' @export

get_nearby_genes <- function(
    query_gr,
    valid_ltr_names = NULL,
    degfile = NULL,
    mode = c("reverse", "forward"),
    txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    window = 100,
    lfc = 0.585,
    padj = 0.05,
    outdir = NULL
) {

  mode <- match.arg(mode)

  message(
    "【Step 3】 Finding nearby target genes (",
    mode,
    " mode)..."
  )

  # =====================================================
  # 0. 只保留具有候选 TF 的 LTR 家族
  # =====================================================

  if (!is.null(valid_ltr_names)) {

    if (length(valid_ltr_names) == 0) {

      warning(
        "No LTR families have selected candidate TFs."
      )

      query_gr <- query_gr[0]

    } else {

      message(
        "  -> Filtering LTR families with ",
        "selected candidate TFs..."
      )

      valid_mask <-
        S4Vectors::mcols(query_gr)$name %in%
        valid_ltr_names

      query_gr <- query_gr[valid_mask]
    }
  }

  # 统一的空结果结构
  empty_genes <- data.frame(
    chr = character(),
    LTR_start = integer(),
    LTR_end = integer(),
    LTR = character(),
    gene_entrez = character(),
    target_gene = character(),
    stringsAsFactors = FALSE
  )

  # 保存结果的内部函数
  save_results <- function(final_candidates, all_nearby) {

    if (!is.null(outdir)) {

      dir.create(
        outdir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      utils::write.csv(
        all_nearby,
        file.path(
          outdir,
          "03_all_ltr_nearby_genes.csv"
        ),
        row.names = FALSE
      )

      utils::write.csv(
        final_candidates,
        file.path(
          outdir,
          "03_final_candidate_genes.csv"
        ),
        row.names = FALSE
      )
    }

    return(list(
      final_candidates = final_candidates,
      all_nearby = all_nearby
    ))
  }

  if (length(query_gr) == 0) {

    warning(
      "No LTR regions remain after filtering."
    )

    return(
      save_results(empty_genes, empty_genes)
    )
  }

  # =====================================================
  # 1. 构建 LTR ±window kb 侧翼窗口
  # =====================================================

  message(
    "  -> Extending LTR windows by ±",
    window,
    " kb..."
  )

  ltr_window <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(query_gr),

    ranges = IRanges::IRanges(
      start = pmax(
        1,
        GenomicRanges::start(query_gr) -
          window * 1000
      ),

      end = GenomicRanges::end(query_gr) +
        window * 1000
    ),

    strand = GenomicRanges::strand(query_gr)
  )

  S4Vectors::mcols(ltr_window) <-
    S4Vectors::mcols(query_gr)

  # =====================================================
  # 2. 加载 TxDb 基因注释并寻找重叠
  # =====================================================

  message("  -> Loading TxDb (", txdb, ")...")

  if (!requireNamespace(txdb, quietly = TRUE)) {
    stop(
      "Missing TxDb package: ",
      txdb,
      ". Please install it via BiocManager."
    )
  }

  txdb_obj <- getExportedValue(txdb, txdb)

  genes_gr <- GenomicFeatures::genes(txdb_obj)

  GenomeInfoDb::seqlevelsStyle(genes_gr) <- "UCSC"

  genes_gr <- GenomeInfoDb::keepStandardChromosomes(
    genes_gr,
    pruning.mode = "coarse"
  )

  # 提取基因 TSS
  # 正链取 start，负链取 end
  gene_tss_gr <- GenomicRanges::resize(
    genes_gr,
    width = 1,
    fix = "start"
  )

  message(
    "  -> Finding overlaps between ",
    "LTR windows and gene TSS..."
  )

  hits <- GenomicRanges::findOverlaps(
    ltr_window,
    gene_tss_gr,
    ignore.strand = TRUE
  )

  if (length(hits) == 0) {

    warning(
      "No nearby genes found within LTR windows."
    )

    return(
      save_results(empty_genes, empty_genes)
    )
  }

  # =====================================================
  # 3. 构建 LTR-邻近基因对应关系
  # =====================================================

  # 使用原始 LTR 坐标，而不是扩展窗口坐标
  ltr_hits <- query_gr[
    S4Vectors::queryHits(hits)
  ]

  gene_hits <- genes_gr[
    S4Vectors::subjectHits(hits)
  ]

  # 选择物种注释数据库
  if (grepl(
    "Mmusculus|mm10|mm39",
    txdb,
    ignore.case = TRUE
  )) {

    if (!requireNamespace(
      "org.Mm.eg.db",
      quietly = TRUE
    )) {
      stop("Please install org.Mm.eg.db")
    }

    org_db <- org.Mm.eg.db::org.Mm.eg.db

  } else {

    if (!requireNamespace(
      "org.Hs.eg.db",
      quietly = TRUE
    )) {
      stop("Please install org.Hs.eg.db")
    }

    org_db <- org.Hs.eg.db::org.Hs.eg.db
  }

  # TxDb genes() 中的基因 ID
  ltr_gene_rel <- data.frame(
    chr = as.character(
      GenomicRanges::seqnames(ltr_hits)
    ),

    LTR_start = GenomicRanges::start(ltr_hits),

    LTR_end = GenomicRanges::end(ltr_hits),

    LTR = as.character(
      S4Vectors::mcols(ltr_hits)$name
    ),

    gene_entrez = as.character(
      names(gene_hits)
    ),

    stringsAsFactors = FALSE
  ) %>%
    dplyr::distinct()

  # Entrez ID 转换为 Gene Symbol
  gene_ids <- unique(
    ltr_gene_rel$gene_entrez
  )

  gene_ids <- gene_ids[
    !is.na(gene_ids) &
    nzchar(gene_ids)
  ]

  if (length(gene_ids) == 0) {

    warning("No valid Entrez gene IDs found.")

    return(
      save_results(empty_genes, empty_genes)
    )
  }

  gene_mapping <- suppressMessages(
    suppressWarnings(
      clusterProfiler::bitr(
        gene_ids,
        fromType = "ENTREZID",
        toType = "SYMBOL",
        OrgDb = org_db
      )
    )
  )

  if (is.null(gene_mapping) ||
      nrow(gene_mapping) == 0) {

    warning(
      "No Gene Symbol mappings were found."
    )

    return(
      save_results(empty_genes, empty_genes)
    )
  }

  ltr_gene_rel <- ltr_gene_rel %>%
    dplyr::inner_join(
      gene_mapping,
      by = c("gene_entrez" = "ENTREZID")
    ) %>%
    dplyr::rename(
      target_gene = SYMBOL
    ) %>%
    dplyr::distinct()

  # =====================================================
  # 4. 根据 mode 筛选
  # =====================================================

  if (mode == "reverse") {

    # 反向分析只保留显著上调 DEG
    if (is.null(degfile) ||
        !file.exists(degfile)) {

      stop(
        "DEG file is required for reverse mode."
      )
    }

    message(
      "  -> Loading DEG file and filtering ",
      "significantly upregulated genes..."
    )

    # 读取 DEG 文件
    if (grepl(
      "\\.xlsx$|\\.xls$",
      degfile,
      ignore.case = TRUE
    )) {

      deg_df <- readxl::read_excel(degfile)

    } else {

      deg_df <- data.table::fread(
        degfile,
        data.table = FALSE,
        showProgress = FALSE
      )
    }

    # 自动识别基因名列
    possible_cols <- c(
      "Gene_Symbol",
      "SYMBOL",
      "gene_name",
      "GeneName",
      "symbol",
      "gene_symbol"
    )

    gene_col <- intersect(
      possible_cols,
      colnames(deg_df)
    )[1]

    if (is.na(gene_col)) {
      gene_col <- colnames(deg_df)[1]

      message(
        "  -> Gene Symbol column not detected; ",
        "using first column: ",
        gene_col
      )
    }

    # 自动识别 log2FC 与 padj 列
    lfc_cols <- c(
      "log2FoldChange",
      "log2FC",
      "logFC",
      "log2_fold_change"
    )

    padj_cols <- c(
      "padj",
      "adj.P.Val",
      "FDR",
      "p_adj"
    )

    lfc_col <- intersect(
      lfc_cols,
      colnames(deg_df)
    )[1]

    padj_col <- intersect(
      padj_cols,
      colnames(deg_df)
    )[1]

    if (is.na(lfc_col) ||
        is.na(padj_col)) {

      stop(
        "DEG file must contain log2FC ",
        "and adjusted P-value columns."
      )
    }

    # =================================================
    # 只保留显著上调的差异表达基因
    #
    # log2FC > lfc
    # padj < padj
    # =================================================

    up_genes <- deg_df %>%
      dplyr::mutate(
        lfc_val = as.numeric(
          .data[[lfc_col]]
        ),

        padj_val = as.numeric(
          .data[[padj_col]]
        )
      ) %>%
      dplyr::filter(
        !is.na(lfc_val),
        !is.na(padj_val),
        lfc_val > .env$lfc,
        padj_val < .env$padj
      ) %>%
      dplyr::pull(
        .data[[gene_col]]
      ) %>%
      as.character() %>%
      unique()

    # LTR 邻近基因与显著上调 DEG 取交集
    final_candidates <- ltr_gene_rel %>%
      dplyr::filter(
        target_gene %in% up_genes
      )

    message(
      "  -> Reverse mode: retained only ",
      "significantly upregulated nearby genes."
    )

  } else {

    # 正向模式：保留全部邻近基因
    final_candidates <- ltr_gene_rel

    message(
      "  -> Forward mode: retained all nearby genes."
    )
  }

  # =====================================================
  # 5. 整理最终结果
  # =====================================================

  final_candidates <- final_candidates %>%
    dplyr::distinct() %>%
    dplyr::arrange(
      LTR,
      chr,
      LTR_start,
      target_gene
    )

  if (nrow(final_candidates) == 0) {

    warning(
      "No candidate genes passed the ",
      "spatial/upregulated-DEG filters."
    )

  } else {

    message(
      "  -> Found ",
      length(unique(final_candidates$target_gene)),
      " candidate target genes."
    )
  }

  # 保存两份结果并返回
  return(
    save_results(
      final_candidates,
      ltr_gene_rel
    )
  )
}

