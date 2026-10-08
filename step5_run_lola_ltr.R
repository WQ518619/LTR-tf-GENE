# =========================================================
# Step 5: LTR LOLA 富集分析
# File:
# R/step5_run_lola_ltr.R
# =========================================================
# ---------------------------------------------------------
# Function:
# run_lola_ltr
#
# Purpose:
# 使用LOLA包利用富集分析，识别出可能与LTR结合的TF
#
# Input:
#
# query_gr
#   LTR query GRanges（来自 Step 2）
#
# universe_gr
#   Universe GRanges（来自 Step 2）
#
# regionDB
#   LOLA 数据库（来自 Step 4）
#
# minOverlap, pval_cutoff, support_cutoff, or_cutoff
#   富集分析参数（可调节）
#
# output_dir
#   可选，保存结果目录
#
# Return:
# list（res_df, sig_tf, ltr_tfs 等）
# ---------------------------------------------------------
run_lola_ltr <- function(query_gr,
                         universe_gr,
                         regionDB,
                         minOverlap ,
                         qval_cutoff, 
                         support_cutoff ,
						 cores,
                         output_dir = "results/step5/") {
  
  cat("【Step 6】开始运行 LTR LOLA 富集分析...\n")
  
  res <- runLOLA(userSets = query_gr,
                 userUniverse = universe_gr,
                 regionDB = regionDB,
                 minOverlap = minOverlap,
                 cores = cores)
  
  # =====================================================
  # 2. 结果整理，全部转数值型
  # =====================================================
  
  res_df <- as.data.frame(res) %>%
    mutate(
      support   = as.numeric(support),
      b         = as.numeric(b),
      c         = as.numeric(c),
      d         = as.numeric(d),
      pValueLog = as.numeric(pValueLog),
      oddsRatio = as.numeric(oddsRatio),
      qValue    = as.numeric(qValue)
    )
  
  cat("LOLA 原始结果行数:", nrow(res_df), "\n")
  
  # =====================================================
  # 3. 过滤非 TF 条目（在 description 列上过滤）
  # =====================================================
  
  non_tf <- c(
    "Epitope",
    "Input",
    "IgG",
    "H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3", "H3K36me3",
    "H3K4me1", "H3K9ac", "H3K14ac", "H4K20me1",
    "RNAPol2", "Pol2", "POLR2A"
  )
  
  res_df <- res_df %>%
    filter(!description %in% non_tf)
  
  cat("过滤非 TF 后行数:", nrow(res_df), "\n")
  
  # =====================================================
  # 4. 筛选显著富集
  # =====================================================
  
  sig_res <- res_df %>%
    filter(
      !is.na(qValue),
      qValue < qval_cutoff,
      support >= support_cutoff
    ) %>%
    arrange(qValue) %>%
    mutate(
      TF = str_remove(description, "\\.bed$"),
	  score = -log10(qValue + 1e-300)
    ) %>%
    relocate(TF, .before = everything())
  
  ltr_tfs <- unique(sig_res$TF)
  
  cat("LTR LOLA 分析完成！显著富集的 TF 数量:", length(ltr_tfs), "\n")
  
  # =====================================================
  # 5. 输出结果
  # =====================================================
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    write.csv(res_df,  file.path(output_dir, "lola_ltr_all.csv"),         row.names = FALSE)
    write.csv(sig_res, file.path(output_dir, "lola_ltr_significant.csv"), row.names = FALSE)
    saveRDS(res_df,    file.path(output_dir, "lola_ltr_res.rds"))
    
    cat("结果已保存至:", output_dir, "\n")
  }
  
  # =====================================================
  # 6. 返回
  # =====================================================
  
  return(list(
    res_df    = res_df,
    sig_tf    = sig_res,
    ltr_tfs   = ltr_tfs,
    sig_count = length(ltr_tfs)
  ))
}