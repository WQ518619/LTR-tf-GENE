# =========================================================
# Step 7: 靶基因侧核心转录因子预测 (基于 TFTF)
# File:
# R/step7_gene_tftf.R
# =========================================================

# ---------------------------------------------------------
# Function: run_gene_tftf_predict
#
# Purpose:
# 采用 TFTF 包对候选靶基因进行多数据库交叉预测。
# 纯依赖海量先验证据库
# 从而最大程度保留真实的分子调控机制。
# ---------------------------------------------------------

run_gene_tftf_predict <- function(
    gene_input,       # <--- [核心改动] 参数名改为更通用的 gene_input
    db_list,
    fimo_score,
    pwm_pval,
    cut_log2fc,
    down_only,
    min_db_support,
    output_dir = "results/step7/"
){
  
  cat("【Step 7】开始 TFTF 靶基因核心转录因子预测...\n")
  
  # =====================================================
  # 1. 智能提取与清洗靶基因 (支持多态输入)
  # =====================================================
  cat("智能解析输入的靶基因格式...\n")
  
  # 情景 A：用户直接输入了基因字符向量 (如 c("TP53", "MYC") 或单个基因)
  if (is.character(gene_input)) {
    my_genes <- gene_input
    cat("检测到输入格式为：纯字符向量 (Character Vector)。\n")
    
  # 情景 B：用户输入了来自流水线上游的数据框 (Data Frame)
  } else if (is.data.frame(gene_input)) {
    cat("检测到输入格式为：数据框 (Data Frame)，正在检索基因列...\n")
    if ("target_gene" %in% colnames(gene_input)) {
      my_genes <- gene_input$target_gene
    } else if ("gene_symbol" %in% colnames(gene_input)) {
      my_genes <- gene_input$gene_symbol
    } else if ("SYMBOL" %in% colnames(gene_input)) {
      my_genes <- gene_input$SYMBOL
    } else {
      stop("数据框中未找到基因名列！请确保列名为 'target_gene', 'gene_symbol' 或 'SYMBOL'。")
    }
    
  # 情景 C：输入了不支持的奇葩格式
  } else {
    stop("不支持的 gene_input 格式！请输入基因名向量 (如 c('A','B')) 或包含基因列的数据框。")
  }
  
  # 统一进行清洗：去首尾空格、去空值、去重
  my_genes <- trimws(as.character(my_genes))
  my_genes <- unique(my_genes[my_genes != "" & !is.na(my_genes)])
  
  cat("成功获取", length(my_genes), "个唯一的候选靶基因。\n")
  if (length(my_genes) == 0) stop("清洗后候选基因列表为空，请检查输入数据！")
  
  # =====================================================
  # 2. 批量进行 TFTF 预测 
  # =====================================================
  cat("连接所选的", length(db_list), "个数据库进行逆向推断...\n")
  
  all_edges <- data.frame(TF = character(), Gene = character(), Source_DB = character(), stringsAsFactors = FALSE)
  
  for (gene in my_genes) {
    cat("  -> 预测靶基因:", gene, "...\n")
    
    tryCatch({
      res <- predict_TF(
        datasets    = db_list, 
        target      = gene,
        FIMO.score  = fimo_score,
        PWMEnrich.p = pwm_pval,
        cut.log2FC  = cut_log2fc,
        down.only   = down_only,
		cor_DB      = NULL,
        app         = FALSE
      )
      
      for (db_name in db_list) {
        item <- res[[db_name]]
        if (is.null(item)) next
        
        current_db_tfs <- c()
        if (is.data.frame(item) && "TF" %in% colnames(item)) {
          current_db_tfs <- as.character(item$TF)
        } else if (is.vector(item) || is.factor(item)) {
          current_db_tfs <- as.character(item)
        }
        
        current_db_tfs <- unique(trimws(current_db_tfs[!is.na(current_db_tfs) & current_db_tfs != ""]))
        
        if (length(current_db_tfs) > 0) {
          temp_edges <- data.frame(TF = current_db_tfs, Gene = gene, Source_DB = db_name, stringsAsFactors = FALSE)
          all_edges <- rbind(all_edges, temp_edges)
        }
      }
    }, error = function(e) {
      cat("    [警告] 基因", gene, "报错已跳过:", conditionMessage(e), "\n")
    })
  }
  
  all_edges <- unique(all_edges)
  if (nrow(all_edges) == 0) stop("未找到任何 TF-Gene 预测网络连线！")
  
  # =====================================================
  # 3. 拓扑统计与多数据库交集计算
  # =====================================================
  cat("统计多源数据库支持度...\n")
  
  pure_network <- unique(all_edges[, c("TF", "Gene")])
  
  tf_counts <- as.data.frame(table(pure_network$TF), stringsAsFactors = FALSE)
  colnames(tf_counts) <- c("TF", "Target_Count")
  
  db_support_count <- aggregate(Source_DB ~ TF, data = all_edges, FUN = function(x) length(unique(x)))
  colnames(db_support_count) <- c("TF", "Hit_DB_Count")
  
  db_support_names <- aggregate(Source_DB ~ TF, data = all_edges, FUN = function(x) paste(unique(x), collapse = " | "))
  colnames(db_support_names) <- c("TF", "Hit_DB_Names")
  
  master_tf_table <- merge(tf_counts, db_support_count, by = "TF")
  master_tf_table <- merge(master_tf_table, db_support_names, by = "TF")
  master_tf_table <- master_tf_table[order(-master_tf_table$Hit_DB_Count, -master_tf_table$Target_Count), ]
  
  full_edges <- merge(pure_network, tf_counts, by = "TF", all.x = TRUE)
  full_edges <- full_edges[order(-full_edges$Target_Count, full_edges$TF), ]
  
  # =====================================================
  # 4. 提取高置信度 TF
  # =====================================================
  high_conf_tfs <- master_tf_table$TF[master_tf_table$Hit_DB_Count >= min_db_support]
  if (length(high_conf_tfs) == 0) {
    warning("当前 min_db_support 过高，未筛选到高置信 TF，将返回所有网络节点。")
    high_conf_tfs <- master_tf_table$TF
  }
  
  # =====================================================
  # 5. 保存结果
  # =====================================================
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(master_tf_table, file.path(output_dir, "1_TFTF_All_TFs_Master_Summary.csv"), row.names = FALSE)
    write.csv(full_edges,      file.path(output_dir, "2_TFTF_Full_TF_Gene_Network.csv"), row.names = FALSE)
    cat("结果已保存至:", output_dir, "\n")
  }
  
  # =====================================================
  # 6. 返回结果列表
  # =====================================================
  return(list(
    master_summary = master_tf_table,
    full_network   = full_edges,
    sig_gene_tfs   = high_conf_tfs 
  ))
}