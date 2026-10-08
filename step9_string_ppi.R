# =========================================================
# Step 9: STRINGdb PPI 分析
# File:
# R/step9_string_ppi.R
# =========================================================

# ---------------------------------------------------------
# Function:
# run_string_ppi
#
# Purpose:
# 使用 STRINGdb 构建 TF 蛋白互作网络（PPI）
#
# 功能包括：
# 1. STRINGdb 基因映射
# 2. 获取 PPI interaction
# 3. 提取 LTR TF 与 Promoter TF 的跨组互作
# 4. 删除孤立节点
# 5. 绘制 STRING 网络图
#
# ---------------------------------------------------------
#
# Input:
#
# ltr_tfs
#   LTR 端显著 TF
#
# promoter_tfs
#   Promoter 端显著 TF
#
# score_threshold
#   STRING 最低置信度
#
# species
#   物种 ID
#   人类 = 9606
#
# remove_isolated
#   是否删除孤立节点
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
# cross_ppi
#   LTR-Promoter 跨组 PPI
#
# ppi_interactions
#   所有 STRING PPI
#
# mapped
#   STRING 映射结果
#
# mapped_clean
#   删除孤立节点后的节点
#
# ---------------------------------------------------------

run_string_ppi <- function(
    ltr_tfs,
    promoter_tfs,
    score_threshold ,
    species = 9606,
    string_version ,
    remove_isolated = TRUE,
    output_dir ="results/step9/",
    string_data_dir = ""
) {
  
  cat("\n================ STRINGdb PPI 分析 ================\n")
  
  # =====================================================
  # Step 1. 初始化 STRINGdb
  # =====================================================
  
  string_db <- STRINGdb$new(
    version = string_version,
    species = species,
    score_threshold = score_threshold,
    input_directory = string_data_dir
  )
  
  # =====================================================
  # 数据清理：剔除表观遗传修饰项（如 Pan-acetyllysine, H3K27ac, H3K4me3 等）
  # =====================================================
  
  epi_pattern <- "(?i)acetyl|H[234]K[0-9]+"
  ltr_tfs <- ltr_tfs[!grepl(epi_pattern, ltr_tfs)]
  promoter_tfs <- promoter_tfs[!grepl(epi_pattern, promoter_tfs)]
  
  # =====================================================
  # Step 2. 合并 TF
  # =====================================================
  
  all_tfs <- unique(c(ltr_tfs, promoter_tfs))
  
  cat("总 TF 数量 (过滤表观修饰后):", length(all_tfs), "\n")
  
  tf_df <- data.frame(
    genes = all_tfs,
    stringsAsFactors = FALSE
  )
  
  # =====================================================
  # Step 3. STRING ID 映射
  # =====================================================
  
  mapped <- string_db$map(
    tf_df,
    "genes",
    removeUnmappedRows = TRUE
  )
  
  if (nrow(mapped) == 0) {
    
    cat("STRINGdb 映射失败。\n")
    
    return(NULL)
  }
  
  cat("成功映射 TF 数量:", nrow(mapped), "\n")
  
  # =====================================================
  # Step 4. 获取 STRING interaction
  # =====================================================
  
  ppi_interactions <- string_db$get_interactions(
    mapped$STRING_id
  )
  
  if (nrow(ppi_interactions) == 0) {
    
    cat("未发现 PPI interaction。\n")
    
    return(NULL)
  }
  
  cat("总 PPI interaction 数量:",
      nrow(ppi_interactions), "\n")
  
  # =====================================================
  # Step 5. 添加基因名
  # =====================================================
  
  ppi_interactions <- ppi_interactions %>%
    
    left_join(
      mapped[, c("STRING_id", "genes")],
      by = c("from" = "STRING_id")
    ) %>%
    rename(from_gene = genes) %>%
    
    left_join(
      mapped[, c("STRING_id", "genes")],
      by = c("to" = "STRING_id")
    ) %>%
    rename(to_gene = genes)
  
  # =====================================================
  # Step 6. 提取跨组 PPI
  # =====================================================
  
  cross_ppi <- ppi_interactions %>%
    
    filter(
      (from_gene %in% ltr_tfs &
         to_gene %in% promoter_tfs) |
        
        (from_gene %in% promoter_tfs &
           to_gene %in% ltr_tfs)
    ) %>%
    
    arrange(desc(combined_score))
  
  cat("跨组 PPI 数量:",
      nrow(cross_ppi), "\n")
  
  # =====================================================
  # Step 7. 删除孤立节点
  # =====================================================
  
  mapped_clean <- mapped
  
  if (remove_isolated) {
    
    cat("正在删除孤立节点...\n")
    
    g_temp <- graph_from_data_frame(
      ppi_interactions[, c("from", "to")],
      directed = FALSE
    )
    
    connected_nodes <- names(
      which(degree(g_temp) >= 1)
    )
    
    mapped_clean <- mapped[
      mapped$STRING_id %in% connected_nodes,
    ]
    
    cat(
      "原始节点:", nrow(mapped),
      " | 删除孤立节点后:",
      nrow(mapped_clean),
      "\n"
    )
  }
  
  # =====================================================
  # Step 8. 输出网络图
  # =====================================================
  
  if (!is.null(output_dir)) {
    
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    pdf_file <- file.path(
      output_dir,
      paste0(
        "PPI_Network_STRINGdb_",
        score_threshold,
        "_Clean.pdf"
      )
    )
    
    pdf(
      pdf_file,
      width = 14,
      height = 11
    )
    
    string_db$plot_network(
      mapped_clean$STRING_id,
      required_score = score_threshold,
      add_link = TRUE
    )
    
    dev.off()
    
    cat("PPI 网络图已保存:\n")
    cat(pdf_file, "\n")
  }
  
  # =====================================================
  # Step 9. R 窗口显示网络
  # =====================================================
  
  cat("正在 R 窗口显示 STRING 网络图...\n")
  
  string_db$plot_network(
    mapped_clean$STRING_id,
    required_score = score_threshold,
    add_link = TRUE
  )
  
  cat("\nSTRINGdb PPI 分析完成。\n")
  
  # =====================================================
  # Step 10. 返回结果
  # =====================================================
  
  return(list(
    
    cross_ppi = cross_ppi,
    
    ppi_interactions = ppi_interactions,
    
    mapped = mapped,
    
    mapped_clean = mapped_clean
  ))
}