# =========================================================
# Step 6b: 基于 TFTF 的无差异基因正向靶向预测
# File: R/step6b_tftf_forward_predict.R
# =========================================================

run_tftf_forward_predict <- function(
    ltr_tfs,
    db_list,
    fimo_score,
    pwm_pval,
    cut_log2fc,
    down_only,
    min_db_support,
    min_tf_count,
    output_dir = "results/step6b/"
){
  
  cat("【Step 6b】开始 TFTF 无差异基因模式的正向靶向预测...\n")
  
  # =====================================================
  # 1. 提取与清洗输入的 LTR-TFs (防御性编程)
  # =====================================================
  ltr_tfs <- unique(trimws(ltr_tfs[ltr_tfs != "" & !is.na(ltr_tfs)]))
  
  cat("成功获取", length(ltr_tfs), "个独立且显著的 LTR-TF 用于靶基因预测。\n")
  if (length(ltr_tfs) == 0) {
    stop("LTR-TF 列表为空，无法进行预测。请检查 Step 5 的 sig_tf 输出。")
  }
  
  # =====================================================
  # 2. 批量进行 TFTF 全基因组预测靶基因
  # =====================================================
  cat("连接所选的", length(db_list), "个数据库进行靶基因正向挖掘...\n")
  
  all_edges <- data.frame(TF = character(), Gene = character(), Source_DB = character(), stringsAsFactors = FALSE)
  
  for (tf in ltr_tfs) {
    cat("  -> 正在挖掘 LTR-TF:", tf, "的全基因组靶点...\n")
    
    tryCatch({
      # 规范调用：直接指定 TFTF 命名空间，完全不依赖 library()
      res <- predict_target(
        datasets    = db_list, 
        tf          = tf,          
        FIMO.score  = fimo_score,
        PWMEnrich.p = pwm_pval,
        cut.log2FC  = cut_log2fc,
        down.only   = down_only,
		cor_DB      = NULL,
        app         = FALSE
      )
      
      # 遍历外部传入的数据库列表，提取证据连线
      for (db_name in db_list) {
        item <- res[[db_name]]
        if (is.null(item)) next
        
        target_col <- intersect(colnames(item), c("Target", "Gene", "symbol"))
        if (length(target_col) > 0) {
          current_db_genes <- as.character(item[[target_col[1]]])
          current_db_genes <- unique(trimws(current_db_genes[!is.na(current_db_genes) & current_db_genes != ""]))
          
          if (length(current_db_genes) > 0) {
            temp_edges <- data.frame(TF = tf, Gene = current_db_genes, Source_DB = db_name, stringsAsFactors = FALSE)
            all_edges <- rbind(all_edges, temp_edges)
          }
        }
      }
    }, error = function(e) {
      cat("    [警告] TF", tf, "报错已跳过:", conditionMessage(e), "\n")
    })
  }
  
  all_edges <- unique(all_edges)
  if (nrow(all_edges) == 0) stop("非常遗憾，未能在数据库中预测到任何靶基因！")
  
  # =====================================================
  # 3. 统计命中数据库的数量 (Hit_DB_Count) 并过滤
  # =====================================================
  cat("统计多源数据库支持度与构建共调控网络...\n")
  
  db_support_count <- aggregate(Source_DB ~ TF + Gene, data = all_edges, FUN = function(x) length(unique(x)))
  colnames(db_support_count)[3] <- "Hit_DB_Count"
  
  # 基于外部传入的 min_db_support 参数进行严苛过滤
  high_conf_edges <- db_support_count[db_support_count$Hit_DB_Count >= min_db_support, ]
  
  # =====================================================
  # 4. 计算共调控网络 (寻找核心枢纽靶基因)
  # =====================================================
  # 规范引入管道符并全程使用 dplyr 命名空间
  `%>%` <- dplyr::`%>%`
  
  target_summary <- high_conf_edges %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      Regulating_LTR_TFs_Count = dplyr::n_distinct(TF),
      Regulating_LTR_TFs       = paste(sort(unique(TF)), collapse = "; "),
      Avg_DB_Support           = mean(Hit_DB_Count),
      .groups = 'drop'
    ) %>%
    # 基于外部传入的 min_tf_count 参数进行核心节点过滤
    dplyr::filter(Regulating_LTR_TFs_Count >= min_tf_count) %>%
    dplyr::arrange(dplyr::desc(Regulating_LTR_TFs_Count), dplyr::desc(Avg_DB_Support))
  
  cat("预测完成！共发现", nrow(target_summary), "个符合严苛阈值的高置信核心靶基因。\n")
  
  final_network_edges <- high_conf_edges %>% dplyr::filter(Gene %in% target_summary$Gene)
  
  # =====================================================
  # 5. 结果落盘保存
  # =====================================================
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(final_network_edges, file.path(output_dir, "1_Forward_TFTF_HighConf_Edges.csv"), row.names = FALSE)
    write.csv(target_summary,      file.path(output_dir, "2_Forward_TFTF_Core_Targets.csv"), row.names = FALSE)
    cat("预测结果(节点与边表)已成功保存至:", output_dir, "\n")
  }
  
  return(list(
    forward_edges  = final_network_edges,
    target_summary = target_summary
  ))
}