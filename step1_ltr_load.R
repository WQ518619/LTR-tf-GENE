# =========================================================
# Step 1: Load Upregulated LTR Families
# File:
#   R/step1_ltr_load.R
# =========================================================

# ---------------------------------------------------------
# Function:
#   load_up_ltr_families
#
# Purpose:
#   Load significantly upregulated LTR families.
#
# Input:
#
#   te_file
#       TEtranscripts result file.
#
#   sheet
#       Excel sheet name.
#
#   lfc_cutoff
#       log2FC cutoff.
#
#   padj_cutoff
#       padj cutoff.
#
# Output:
#   Filtered LTR table.
#
# Return:
#   list
#
# Example:
#   ltr_res <- load_up_ltr_families(
#     te_file = "data/TE.xlsx"
#   )
# ---------------------------------------------------------

load_up_ltr_families <- function(
    te_file,
    lfc_cutoff ,
    padj_cutoff,
    output_dir = "results/step1/"
) {
  message("Loading significant LTR families...")
  
  ltr_df <- read_excel(te_file)
  
  # 将 log2FoldChange 和 padj 转换为数值（避免字符比较错误）
  ltr_df <- ltr_df %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      padj           = as.numeric(padj)
    )
  
  # ① 只保留 ID 中含有 "LTR" 的行
  ltr_df <- ltr_df %>%
    filter(str_detect(ID, "LTR"))
  
  # ② 筛选显著上调的
  ltr_df <- ltr_df %>%
    filter(
      log2FoldChange > lfc_cutoff,
      padj < padj_cutoff
    )
  
  # ③ 提取家族名
  ltr_df$family <- str_extract(ltr_df$ID, "^[^:]+")
  ltr_up <- unique(ltr_df$family)
  
  # ④ 可选：保存结果
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(ltr_df, file.path(output_dir, "upregulated_LTR.csv"), row.names = FALSE)
    saveRDS(ltr_up, file.path(output_dir, "ltr_up.rds"))
  }
  
  return(list(
    ltr_df = ltr_df,
    ltr_up = ltr_up
  ))
}