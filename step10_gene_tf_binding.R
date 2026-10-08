# =========================================================
# Step 10: Gene-TF Binding 分析
# File:
# R/step10_gene_tf_binding.R
# =========================================================

# ---------------------------------------------------------
# Function:
# build_gene_tf_binding
#
# Purpose:
# 基于 ChIP-seq peak 与启动子 overlap，
# 构建 Gene-TF binding 数据集。
#
# 功能包括：
#
# 1. 构建候选基因 promoter
# 2. 与 TF ChIP-seq peak overlap
# 3. 提取 TF-Gene binding 关系
# 4. 生成 gene_tf_summary
#
# ---------------------------------------------------------
#
# Input:
#
# final_candidates
#   Step6输出结果
#
# genes_gr
#   全基因组 gene GRanges
#
# chip_gr
#   TF ChIP-seq peak GRanges
#
# upstream, downstream
#   promoter 区域范围
#
# output_dir
#   输出目录
#
# ---------------------------------------------------------
#
# Return:
#
# list:
#
# gene_tf_binding
#   TF-Gene binding 明细
#
# gene_tf_summary
#   每个基因对应 TF 汇总
#
# ---------------------------------------------------------


build_gene_tf_binding <- function(
    tftf_full_network,   # 接收 Step7 生成的 full_network
    gene_tfs,            # 接收 Step7 筛选出的核心 Hub TFs
    output_dir = "results/step10/"
) {
  
  cat("\n=================== Gene-TF Binding 提取 ===================\n")
  cat("使用 TFTF 多数据库先验证据，跳过坐标重叠计算...\n")
  
  # =====================================================
  # 1. 提取 TFTF 连线并格式化为下游需要的列名
  # =====================================================
  # TFTF 原始列为: TF, Gene, Target_Count
  # 下游 Step 11 需要: TF, gene_symbol, peak_score
  gene_tf_binding <- tftf_full_network %>%
    filter(TF %in% gene_tfs) %>%           # 仅保留核心 TF 的连线
    mutate(
      gene_symbol = Gene,
      peak_score  = 1000                   # 伪造一个满分的 peak_score 供 Step11 生成线条粗细 (weight = 1)
    ) %>%
    dplyr::select(gene_symbol, TF, peak_score)
  
  # =====================================================
  # 2. TF 靶向汇总统计
  # =====================================================
  gene_tf_summary <- gene_tf_binding %>%
    group_by(gene_symbol) %>%
    summarise(
      n_TF = n_distinct(TF),
      TFs = paste(sort(unique(TF)), collapse = "; "),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_TF))
  
  # =====================================================
  # 3. 输出统计与保存
  # =====================================================
  cat("总 binding 数:", nrow(gene_tf_binding), "\n")
  cat("涉及基因数:", length(unique(gene_tf_binding$gene_symbol)), "\n")
  cat("涉及核心 TF 数:", length(unique(gene_tf_binding$TF)), "\n")
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(gene_tf_binding, file.path(output_dir, "gene_tf_binding_detailed.csv"), row.names = FALSE)
    write.csv(gene_tf_summary, file.path(output_dir, "gene_tf_summary.csv"), row.names = FALSE)
  }
  
  return(list(
    gene_tf_binding = gene_tf_binding,
    gene_tf_summary = gene_tf_summary
  ))
}