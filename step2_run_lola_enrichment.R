
# =========================================================
# Step 2: LOLA Enrichment Analysis 
# File: R/run_lola_enrichment.R
# =========================================================

#' Run LOLA Enrichment Analysis
#'
#' @description
#' Load ChIP-Atlas data, build a LOLA region database,
#' and perform LTR-family enrichment analysis.
#'
#' Tier A: LOLA native qValue < qval.
#' Tier B: P < 0.05, OR > 1, support >= threshold.
#' Tier C: OR > 1, support >= threshold.
#'
#' Tier assignment follows A > B > C.
#' If the selected tiers have no candidates, the function
#' automatically falls back to Tier C.
#'
#' @param query_gr GRanges. Target LTR regions.
#' @param universe_gr GRanges. Background regions.
#' @param chipfile Character. ChIP-Atlas .rds or .bed file.
#' @param loladir Character. LOLA database directory.
#' @param qval Numeric. Native LOLA qValue threshold.
#' @param support Numeric. Minimum support for Tier B/C.
#' @param allowed_tiers Character. Default c("A", "B").
#' @param outdir Character. Output directory.
#'
#' @import dplyr stringr GenomicRanges LOLA
#' @importFrom magrittr %>%
#' @importFrom readr read_tsv
#' @export

run_lola_enrichment <- function(
    query_gr,
    universe_gr,
    chipfile,
    loladir = "LOLA_DB",
    qval = 0.05,
    support = 1,
    allowed_tiers = c("A", "B"),
    outdir = NULL
) {

  message(
    "【Step 2】 Running LOLA Enrichment ",
    "Analysis (Tiered Strategy)..."
  )

  # 统一 Tier 名称
  # 支持 c("A","B") 或 c("Tier A","Tier B")
  allowed_tiers <- unique(
    toupper(trimws(as.character(allowed_tiers)))
  )

  allowed_tiers <- sub(
    "^TIER\\s*",
    "",
    allowed_tiers
  )

  if (length(allowed_tiers) == 0 ||
      anyNA(allowed_tiers) ||
      !all(allowed_tiers %in% c("A", "B", "C"))) {
    stop("allowed_tiers must contain A, B and/or C.")
  }

  selected_tiers <- paste("Tier", allowed_tiers)

  # =====================================================
  # 2.1 读取 ChIP-Atlas 数据
  # =====================================================

  message("  -> Loading ChIP-Atlas data...")

  if (grepl("\\.rds$", chipfile, ignore.case = TRUE)) {

    chip <- readRDS(chipfile)

    if (!"TF" %in% colnames(chip)) {

      chip$TF <- stringr::str_extract(
        chip$info,
        "Name=[^;]+"
      ) %>%
        stringr::str_replace("Name=", "") %>%
        stringr::str_replace("%20.*", "") %>%
        stringr::str_trim()
    }

  } else {

    bed_lines <- readLines(chipfile, warn = FALSE)

    bed_body <- bed_lines[
      nzchar(trimws(bed_lines)) &
        !grepl(
          "^\\s*(track|browser|#)",
          bed_lines
        )
    ]

    chip <- readr::read_tsv(
      I(bed_body),
      col_names = FALSE,
      show_col_types = FALSE
    )

    if (ncol(chip) < 5) {
      stop("ChIP BED must contain at least 5 columns.")
    }

    chip <- chip[, 1:5]

    colnames(chip) <- c(
      "chr", "start", "end", "info", "score"
    )

    chip$TF <- stringr::str_extract(
      chip$info,
      "Name=[^;]+"
    ) %>%
      stringr::str_replace("Name=", "") %>%
      stringr::str_replace("%20.*", "") %>%
      stringr::str_trim()
  }

  # 剔除无效 TF 标签及无效坐标
  invalid_tfs <- c(
    "Unclassified",
    "Input",
    "IgG",
    "Control",
    "Mock"
  )

  chip <- chip %>%
    dplyr::mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      TF = as.character(TF)
    ) %>%
    dplyr::filter(
      !is.na(TF),
      nzchar(TF),
      !TF %in% invalid_tfs,
      !grepl(
        "Unclassified",
        TF,
        ignore.case = TRUE
      ),
      !is.na(start),
      !is.na(end),
      is.finite(start),
      is.finite(end),
      start >= 0,
      end > start
    )

  if (nrow(chip) == 0) {
    stop("No valid ChIP peaks remain.")
  }

  # =====================================================
  # 2.2 构建 LOLA 字典
  # 保留原来的数据库文件布局
  # =====================================================

  message("  -> Building LOLA Region Database...")

  coll_name <- "ChIPAtlas_Enrichment"

  regions_dir <- file.path(
    loladir,
    coll_name,
    "regions"
  )

  dir.create(
    regions_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # 清理旧 BED，防止重复运行时引用过期数据
  old_beds <- list.files(
    regions_dir,
    pattern = "\\.bed$",
    full.names = TRUE
  )

  if (length(old_beds) > 0) {
    removed <- file.remove(old_beds)

    if (!all(removed)) {
      stop("Failed to remove old LOLA BED files.")
    }
  }

  tf_list <- unique(chip$TF)

  for (tf in tf_list) {

    tf_df <- chip %>%
      dplyr::filter(TF == tf) %>%
      dplyr::select(chr, start, end) %>%
      dplyr::distinct()

    utils::write.table(
      tf_df,
      file.path(
        regions_dir,
        paste0(tf, ".bed")
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }

  bed_files <- list.files(
    regions_dir,
    pattern = "\\.bed$",
    full.names = FALSE
  )

  utils::write.table(
    data.frame(
      filename = bed_files,
      description = sub("\\.bed$", "", bed_files)
    ),
    file.path(
      loladir,
      coll_name,
      "index.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      collection = coll_name,
      description = "ChIP-Atlas Data"
    ),
    file.path(loladir, "collection.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      collection = coll_name,
      filename = bed_files,
      description = sub("\\.bed$", "", bed_files)
    ),
    file.path(loladir, "regionSets.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  regionDB <- LOLA::loadRegionDB(loladir)

  # =====================================================
  # 2.3 执行 LOLA 富集
  # =====================================================

  message("  -> Splitting target sets by family...")

  userSets <- split(
    query_gr,
    S4Vectors::mcols(query_gr)$name
  )

  family_counts <- lengths(userSets)

  userSets <- userSets[
    family_counts >= 1
  ]

  if (length(userSets) == 0) {
    stop("No LTR families found in query_gr.")
  }

  message(
    "  -> Running LOLA for ",
    length(userSets),
    " families simultaneously..."
  )

  res <- tryCatch({

    LOLA::runLOLA(
      userSets = userSets,
      userUniverse = universe_gr,
      regionDB = regionDB,
      minOverlap = 1,
      cores = 1,
      redefineUserSets = TRUE
    )

  }, error = function(e) {

    stop(
      "LOLA execution failed: ",
      conditionMessage(e)
    )
  })

  # =====================================================
  # 2.4 Tier 分级
  # =====================================================

  # 无 LOLA 结果时，返回空表
  if (is.null(res) || nrow(res) == 0) {

    warning(
      "LOLA returned no enrichment results."
    )

    res_tiered <- data.frame(
      LTR = character(),
      TF = character(),
      support = numeric(),
      qValue = numeric(),
      oddsRatio = numeric(),
      pValue = numeric(),
      score = numeric(),
      Tier = character(),
      stringsAsFactors = FALSE
    )

  } else {

    res_df <- as.data.frame(res)

    # 保留 LOLA 原生 qValue
    if (!"qValue" %in% colnames(res_df)) {

      warning(
        "LOLA did not return qValue. ",
        "Tier A cannot be evaluated. ",
        "Please check the qvalue package."
      )

      res_df$qValue <- NA_real_
    }

    # 计算原始 P 值
    res_df <- res_df %>%
      dplyr::mutate(
        support = as.numeric(support),
        qValue = as.numeric(qValue),
        oddsRatio = as.numeric(oddsRatio),
        pValue = 10^(-as.numeric(pValueLog))
      )

    # 过滤表观遗传修饰等非 TF 标签
    non_tf <- c(
      "Epitope",
      "H3K27ac",
      "H3K4me3",
      "RNAPol2"
    )

    res_df <- res_df %>%
      dplyr::filter(
        !description %in% non_tf
      )

    # -----------------------------------------------
    # Tier A / B / C
    # -----------------------------------------------

    res_tiered <- res_df %>%
      dplyr::mutate(

        # 已确认 userSet 为 LTR 家族名称
        LTR = as.character(userSet),

        TF = stringr::str_remove(
          description,
          "\\.bed$"
        ),

        score = -log10(qValue + 1e-300),

        Tier = dplyr::case_when(

          # Tier A:
          # 仅使用 LOLA 原生 qValue
          !is.na(qValue) &
            qValue < .env$qval ~ "Tier A",

          # Tier B:
          # 原始 P < 0.05 + OR > 1 + support
          !is.na(pValue) &
            pValue < 0.05 &
            !is.na(oddsRatio) &
            oddsRatio > 1 &
            !is.na(support) &
            support >= .env$support ~ "Tier B",

          # Tier C:
          # 不满足 A/B，但存在正向 OR 与重叠支持
          !is.na(oddsRatio) &
            oddsRatio > 1 &
            !is.na(support) &
            support >= .env$support ~ "Tier C",

          TRUE ~ "None"
        )
      ) %>%
      dplyr::arrange(
        LTR,
        factor(
          Tier,
          levels = c(
            "Tier A",
            "Tier B",
            "Tier C",
            "None"
          )
        ),
        qValue
      )
  }

  # =====================================================
  # 2.5 筛选进入下游分析的 TF
  # =====================================================

  fallback_used <- FALSE

  # 默认使用 Tier A + B
  sig_res <- res_tiered %>%
    dplyr::filter(
      Tier %in% selected_tiers
    )

  # 只有选中等级完全没有候选时才回退至 C
  if (nrow(sig_res) == 0 &&
      !("Tier C" %in% selected_tiers)) {

    message(
      "  -> No candidates in selected tiers. ",
      "Checking Tier C..."
    )

    fallback_res <- res_tiered %>%
      dplyr::filter(
        Tier == "Tier C"
      )

    if (nrow(fallback_res) > 0) {

      sig_res <- fallback_res
      fallback_used <- TRUE

      message(
        "  -> Tier C fallback activated ",
        "(exploratory candidates)."
      )
    }
  }

  if (nrow(sig_res) == 0) {

    message(
      "  -> No candidate TFs available ",
      "for downstream analysis."
    )

  } else {

    message(
      "  -> Found ",
      length(unique(sig_res$TF)),
      " unique candidate TFs."
    )
  }

  # 全部 A/B/C 候选
  candidate_res <- res_tiered %>%
    dplyr::filter(
      Tier %in% c(
        "Tier A",
        "Tier B",
        "Tier C"
      )
    )

  # 仅 Tier A
  tierA_res <- res_tiered %>%
    dplyr::filter(
      Tier == "Tier A"
    )

  # =====================================================
  # 2.6 保存结果
  # =====================================================

  if (!is.null(outdir)) {
    dir.create(
      outdir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    # 1. 保存全部 LOLA 结果（包含 Tier A/B/C/None）
    utils::write.csv(
      res_tiered,
      file.path(
        outdir,
        "02_lola_all_results.csv"
      ),
      row.names = FALSE
    )

    # 2. 保存最终进入下游分析的候选 TF
    utils::write.csv(
      sig_res,
      file.path(
        outdir,
        "02_lola_selected_tfs.csv"
      ),
      row.names = FALSE
    )
  }

  # =====================================================
  # 2.7 返回结果
  # =====================================================

  return(list(
    sig_tf = sig_res,
    ltr_tfs = unique(sig_res$TF),
    fallback_used = fallback_used
  ))
}