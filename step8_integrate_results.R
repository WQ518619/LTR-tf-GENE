# =========================================================
# Step 8: 整合最终结果
# File:
# R/step8_integrate_results.R
# =========================================================

# ---------------------------------------------------------
# Function:
# integrate_final_results
#
# Purpose:
# 整合：
#   LTR nearby gene
#   LTR LOLA
#   promoter LOLA
#
# 并识别 Bridge TF
#
# Input:
#
# ltr_tfs
#   LTR端显著TF
#
# promoter_tfs
#   promoter端显著TF
#
# final_candidates
#   Step6输出结果
#
# sig_tf
#   Step5显著TF结果
#
# sig_tf_promoter
#   Step7显著TF结果
#
# Return:
#
# list:
#
# final_output
#   最终整合结果
#
# bridge_tfs
#   Bridge TF列表
# ---------------------------------------------------------

integrate_final_results <- function(
    ltr_tfs,
    gene_tfs,              # 接收 Step7 的 sig_gene_tfs (原 promoter_tfs)
    final_candidates,
    sig_tf,                # LTR 端 LOLA 结果
    tftf_summary,          # 接收 Step7 的 master_summary (原 sig_tf_promoter)
    output_dir = "results/step8/"
){
  cat("【Step 8】整合最终结果 (LTR-LOLA 与 Gene-TFTF)...\n")
  
  # =====================================================
  # 1. 寻找核心 Bridge TF (跨组交集)
  # =====================================================
  bridge_tfs <- intersect(ltr_tfs, gene_tfs)
  cat("Bridge TF数量 (同时结合LTR和靶基因):", length(bridge_tfs), "\n")
  
  # =====================================================
  # 2. LTR TF 统计 (处理 LOLA 的 score/qValue)
  # =====================================================
  if (!"score" %in% colnames(sig_tf)) {
    if ("qValue" %in% colnames(sig_tf)) {
      sig_tf <- sig_tf %>% mutate(score = -log10(as.numeric(qValue) + 1e-300))
    } else {
      stop("sig_tf 中缺少 qValue，无法计算 score")
    }
  }
  
  sig_tf <- sig_tf %>% mutate(score = as.numeric(score), support = as.numeric(support))
  
  ltr_tf_stats <- sig_tf %>%
    group_by(TF) %>%
    summarise(
      ltr_max_score  = max(score, na.rm = TRUE),
      ltr_support    = max(support, na.rm = TRUE),
      .groups = "drop"
    )
  
  # =====================================================
  # 3. Gene TF 统计 (将 TFTF 的 Hit_DB_Count 转化为 score 供下游使用)
  # =====================================================
  # 这一步是为了完美兼容你 Step 11 中需要的 sig_tf_promoter$score 格式
  gene_tf_stats <- tftf_summary %>%
    filter(TF %in% gene_tfs) %>%
    mutate(
      score = as.numeric(Hit_DB_Count),     # 用命中的数据库数量代表可信度 score
      support = as.numeric(Target_Count)    # 用调控的基因数代表 support
    ) %>%
    dplyr::select(TF, score, support)       # 瘦身，对齐原先的数据结构
  
  # =====================================================
  # 4. 整理最终 LTR-Gene 距离结果 (保留你原有的优先级打分逻辑)
  # =====================================================
  final_output <- final_candidates %>%
    mutate(
      abs_dist = abs(distance),
      distance_category = case_when(
        abs_dist <= 10000 ~ "Very Close",
        abs_dist <= 50000 ~ "Close",
        TRUE              ~ "Distal"
      ),
      priority_score = 1 / log10(abs_dist + 1000)
    ) %>%
    arrange(desc(priority_score))
  
  # =====================================================
  # 5. 保存并返回结果
  # =====================================================
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(final_output, file.path(output_dir, "final_integrated_results.csv"), row.names = FALSE)
    write.table(bridge_tfs, file.path(output_dir, "bridge_tfs.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  return(list(
    final_output      = final_output,
    bridge_tfs        = bridge_tfs,
    ltr_tf_stats      = ltr_tf_stats,
    gene_tf_stats     = gene_tf_stats  # 新输出，将传递给 Step 11
  ))
}