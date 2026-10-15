# ==========================================================
# main.R – LTR-TF 调控网络主控（双模式智能分流版）
# ==========================================================

library(here)
setwd(here())

# 0. 加载配置 --------------------------------------------------------------
source(here("R", "config.R"))

# 自动创建所有输出目录
for (d in unlist(config$output)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# 日志系统 ----------------------------------------------------------------
log_file <- file.path(config$output$logs, paste0("analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con  <- file(log_file, open = "wt")
sink(log_con, append = TRUE, split = TRUE)

# 确保即使出错也关闭日志
on.exit({ sink(); close(log_con) }, add = TRUE)

cat("=== LTR-TF Network Analysis ===\n")
cat("Pipeline Mode:", config$pipeline_mode, "\n")
cat("Start time:", Sys.time(), "\n\n")

# 1. 加载所有函数定义 ------------------------------------------------------
source(here("R", "step0_packages.R"))
source(here("R", "step1_ltr_load.R"))
source(here("R", "step2_ltr_granges.R"))
source(here("R", "step2b_ltr_phylogeny.R"))
source(here("R", "step3_chip_ltr_overlap.R"))
source(here("R", "step4_create_lola_db.R"))
source(here("R", "step5_run_lola_ltr.R"))
source(here("R", "step6_nearby_genes.R"))
source(here("R", "step6b_tftf_forward_predict.R"))  
source(here("R", "step7_gene_tftf.R"))              
source(here("R", "step8_integrate_results.R"))
source(here("R", "step9_string_ppi.R"))
source(here("R", "step10_gene_tf_binding.R"))
source(here("R", "step11_build_network.R"))
source(here("R", "step12_plot_network.R"))
source(here("R", "step13_3d_loop.R"))
cat("All functions loaded.\n\n")
load_required_packages()
start_time <- Sys.time()

# ======================== 公共预处理模块 ========================

# -------- Step 1: 加载上调 LTR 家族 --------
cat("[Step 1] Loading upregulated LTR families...\n")
step1_res <- tryCatch(
  load_up_ltr_families(
    te_file     = config$paths$te_file,
    lfc_cutoff  = config$step1$lfc_cutoff,
    padj_cutoff = config$step1$padj_cutoff,
    output_dir  = config$output$step1
  ),
  error = function(e) stop("Step1 failed: ", e$message)
)
ltr_up <- step1_res$ltr_up
cat("Step 1 done.", length(ltr_up), "LTR families retained.\n\n")

# -------- Step 2: 构建 LTR GRanges --------
cat("[Step 2] Building LTR GRanges...\n")
step2_res <- tryCatch(
  build_ltr_granges(
    rmsk_file          = config$paths$rmsk_file,
    te_file            = config$paths$te_file,
    ltr_up             = ltr_up,
    expressed_baseMean = config$step2$expressed_baseMean,
    output_dir         = config$output$step2
  ),
  error = function(e) stop("Step2 failed: ", e$message)
)
query_gr    <- step2_res$query_gr
rmsk_parsed <- step2_res$rmsk_parsed
universe_gr <- step2_res$universe_gr  
cat("Step 2 done. query_gr:", length(query_gr), ", universe_gr:", length(universe_gr), "\n\n")

# -------- Step 2b: 进化树建树与活跃 LTR 提取 --------
cat("[Step 2b] Running LTR Phylogeny & Active Copy Extraction...\n")
step2b_res <- tryCatch(
  run_ltr_phylogeny(
    query_gr_obj    = query_gr,
    dfam_fasta_file = config$paths$dfam_fasta,
    iqtree_path     = config$paths$iqtree_exe,
    top_n_youngest  = config$step2b$top_n_youngest,
    output_dir      = config$output$step2b,
    run_phylogeny   = config$step2b$run_phylogeny
  ),
  error = function(e) stop("Step 2b failed: ", e$message)
)
if (!is.null(step2b_res) && nrow(step2b_res) > 0) {
  cat("Step 2b completed. Trees and Top20 BED saved.\n")
} else {
  cat("Phylogeny skipped or no active copies extracted.\n")
}
cat("\n")

# -------- Step 3: ChIP-Atlas 处理 & LTR-TF overlap --------
cat("[Step 3] ChIP-Atlas processing...\n")
step3_res <- tryCatch(
  process_chip_atlas(
    chip_file    = config$paths$chip_file,
    score_cutoff = config$step3$chip_score_cutoff,
    output_dir   = config$output$step3
  ),
  error = function(e) stop("Step3 failed: ", e$message)
)
chip_gr <- step3_res$chip_gr
chip_df <- step3_res$chip_df
cat("Step 3 done. Peaks:", step3_res$peak_count, "TFs:", step3_res$tf_count, "\n\n")

cat("[Step 3B] LTR-TF binding overlap...\n")
ltr_tf_binding <- tryCatch(
  get_ltr_tf_binding(
    query_gr       = query_gr,
    chip_gr        = chip_gr,
    min_peak_score = config$step3$min_peak_score,
    output_dir     = config$output$step3
  ),
  error = function(e) stop("Step3B failed: ", e$message)
)
cat("Step 3B done. LTR-TF pairs:", nrow(ltr_tf_binding), "\n\n")

# -------- Step 4: 创建 LOLA 数据库 --------
cat("[Step 4] Creating LOLA database...\n")
step4_res <- tryCatch(
  create_lola_db(
    chip_df           = chip_df,
    base_dir          = config$paths$lola_base_dir,
    coll_name         = config$step4$lola_collection_name,
    min_peaks_per_tf  = config$step4$min_peaks_per_tf,
    output_dir        = config$output$step4
  ),
  error = function(e) stop("Step4 failed: ", e$message)
)
regionDB <- step4_res$regionDB
cat("Step 4 done. LOLA DB with", step4_res$tf_count, "TFs.\n\n")

# -------- Step 5: LTR LOLA 富集 --------
cat("[Step 5] Running LTR LOLA enrichment...\n")
step5_res <- tryCatch(
  run_lola_ltr(
    query_gr       = query_gr,
    universe_gr    = universe_gr,
    regionDB       = regionDB,
    minOverlap     = config$step5$minOverlap,
    qval_cutoff    = config$step5$qval_cutoff,      
    support_cutoff = config$step5$support_cutoff,
    cores          = config$step5$cores,            
    output_dir     = config$output$step5
  ),
  error = function(e) stop("Step5 failed: ", e$message)
)
ltr_tfs <- step5_res$ltr_tfs
sig_tf  <- step5_res$sig_tf
cat("Step 5 done. Significant LTR TFs:", length(ltr_tfs), "\n\n")

# -------- Step 6: 邻近上调基因 --------
  cat("[Step 6] Finding nearby upregulated genes...\n")
  step6_res <- tryCatch(
    find_nearby_up_genes(
      query_gr    = query_gr,
      DEG_file    = config$paths$deg_file,
      gene_col    = config$step6$gene_col,
      lfc_col     = config$step6$lfc_col,
      padj_col    = config$step6$padj_col,
      window_kb   = config$step6$window_kb,
      lfc_cutoff  = config$step6$lfc_cutoff,
      padj_cutoff = config$step6$padj_cutoff,
      output_dir  = config$output$step6,
	  pipeline_mode = config$pipeline_mode
    ),
    error = function(e) stop("Step6 failed: ", e$message)
  )
  final_candidates <- step6_res$final_candidates
  cat("Step 6 done. Candidate genes:", step6_res$up_gene_count, "\n\n")

# ======================== 智能分流模块 ========================

if (config$pipeline_mode == "forward") {
  
  cat("======================================================\n")
  cat("【模式探测】：启动无差异基因正向预测模式 (Forward Mode)\n")
  cat("======================================================\n\n")
  
  # 提取 Step 5 算出的显著 TF
  target_tfs <- step5_res$sig_tf$TF
  
  step6b_res <- tryCatch(
    run_tftf_forward_predict(
      ltr_tfs        = target_tfs,
      db_list        = config$step7$db_list,             
      fimo_score     = config$step7$fimo_score,          
      pwm_pval       = config$step7$pwm_pval,            
      cut_log2fc     = config$step7$cut_log2fc,          
      down_only      = config$step7$down_only,           
      min_db_support = config$step6b$min_db_support, 
      min_tf_count   = config$step6b$min_tf_count,   
      output_dir     = config$output$step6b
    ),
    error = function(e) stop("Step 6b failed: ", e$message)
  )
  
  # 为 Step 12 自动组装正向网络数据
  cat("[组装数据] 正在为正向网络生成标准的节点与边表...\n")
  library(dplyr)
  
  core_tfs <- unique(step6b_res$forward_edges$TF)
  
  # LTR -> TF (Top 3 按峰值打分)
  edges_ltr_tf <- ltr_tf_binding %>%
    filter(TF %in% core_tfs) %>%
    group_by(LTR) %>% slice_max(peak_score, n = 3, with_ties = FALSE) %>% ungroup() %>%
    transmute(from = as.character(LTR), to = as.character(TF), type = "LTR-TF", weight = peak_score / 1000)
  
  # TF -> Target Gene
  edges_tf_gene <- step6b_res$forward_edges %>%
    transmute(from = as.character(TF), to = as.character(Gene), type = "TF-Gene", weight = Hit_DB_Count)
  
  edge_list <- bind_rows(edges_ltr_tf, edges_tf_gene)
  
  nodes_df <- data.frame(name = unique(c(edge_list$from, edge_list$to))) %>%
    mutate(node_type = case_when(
      name %in% edges_ltr_tf$from ~ "LTR",
      name %in% edges_tf_gene$to  ~ "Target_Gene",
      TRUE                        ~ "TF"
    ))
    
  cat("正向模式网络数据组装完毕！直接进入绘图。\n\n")

} else {
  
  cat("======================================================\n")
  cat("【模式探测】：启动有差异基因的逆向验证模式 (Reverse Mode)\n")
  cat("======================================================\n\n")
 

  # ------- Step 7: 靶基因侧核心转录因子预测 (TFTF) -------
  cat("[Step 7] Running gene side TFTF prediction...\n")
  step7_res <- tryCatch(
    run_gene_tftf_predict(
      gene_input       = final_candidates,                
      db_list          = config$step7$db_list,	
      fimo_score       = config$step7$fimo_score,
      pwm_pval         = config$step7$pwm_pval,
      cut_log2fc       = config$step7$cut_log2fc,
      down_only        = config$step7$down_only,
      min_db_support   = config$step7$min_db_support,
      output_dir       = config$output$step7
    ),
    error = function(e) stop("Step7 failed: ", e$message)
  )
  gene_side_tfs <- step7_res$sig_gene_tfs
  cat("Step 7 done. Gene side high-confidence TFs:", length(gene_side_tfs), "\n\n")

  # -------- Step 8: 整合结果 --------
  cat("[Step 8] Integrating final results...\n")
  step8_res <- tryCatch(
    integrate_final_results(
      ltr_tfs          = ltr_tfs,
      gene_tfs         = gene_side_tfs,               
      final_candidates = final_candidates,
      sig_tf           = sig_tf,
      tftf_summary     = step7_res$master_summary,    
      output_dir       = config$output$step8
    ),
    error = function(e) stop("Step8 failed: ", e$message)
  )
  bridge_tfs <- step8_res$bridge_tfs
  cat("Step 8 done. Bridge TFs:", length(bridge_tfs), "\n\n")

  # -------- Step 9: STRING PPI --------
  cat("[Step 9] STRING PPI analysis...\n")
  step9_res <- tryCatch(
    run_string_ppi(
      ltr_tfs          = ltr_tfs,
      promoter_tfs     = gene_side_tfs,               
      score_threshold  = config$step9$score_threshold,
      species          = config$project$taxonomy_id,
      string_version   = config$step9$string_version,
      remove_isolated  = TRUE,
      output_dir       = config$output$step9,
      string_data_dir  = config$step9$string_data_dir
    ),
    error = function(e) stop("Step9 failed: ", e$message)
  )
  cross_ppi <- step9_res$cross_ppi
  cat("Step 9 done. Cross-group PPIs:", nrow(cross_ppi), "\n\n")

  # -------- Step 10: Gene-TF binding --------
  cat("[Step 10] Extracting Gene-TF binding from TFTF...\n")
  step10_res <- tryCatch(
    build_gene_tf_binding(
      tftf_full_network = step7_res$full_network,     
      gene_tfs          = gene_side_tfs,              
      output_dir        = config$output$step10
    ),
    error = function(e) stop("Step10 failed: ", e$message)
  )
  gene_tf_binding <- step10_res$gene_tf_binding
  cat("Step 10 done. Gene-TF pairs:", nrow(gene_tf_binding), "\n\n")

  # -------- Step 11: 构建网络表格 --------
  cat("[Step 11] Assembling network tables...\n")
  step11_res <- tryCatch(
    build_network_tables(
      ltr_tf_binding         = ltr_tf_binding,
      gene_tf_binding        = gene_tf_binding,
      cross_ppi              = cross_ppi,
      sig_tf                 = sig_tf,
      sig_tf_promoter        = step8_res$gene_tf_stats, 
      top_ltr_tf_per_family  = config$step11$top_ltr_tf_per_family,
      min_ltr_peak_score     = config$step11$min_ltr_peak_score,
      output_dir             = config$output$step11
    ),
    error = function(e) stop("Step11 failed: ", e$message)
  )
  nodes_df  <- step11_res$nodes_df
  edge_list <- step11_res$edge_list
  cat("Step 11 done. Nodes:", nrow(nodes_df), "Edges:", nrow(edge_list), "\n\n")

  # -------- Step 13: 3D 染色质环化验证 (仅逆向模式) --------
  cat("[Step 13] Verifying 3D Chromatin Loops...\n")
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  genes_gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  seqlevelsStyle(genes_gr) <- "UCSC"
  genes_gr <- keepStandardChromosomes(genes_gr, pruning.mode = "coarse")

  step13_res <- tryCatch(
    run_3d_loop_verification(
      query_gr          = query_gr,                    
      final_candidates  = final_candidates,  
      genes_gr          = genes_gr,                    
      loop_bedpe_file   = config$paths$hic_bedpe,      
      upstream          = config$step13$upstream,       
      downstream        = config$step13$downstream,
      output_dir        = config$output$step13
    ),
    error = function(e) stop("Step 13 failed: ", e$message)
  )
  cat("Step 13 done. 3D Loop interactions analyzed.\n\n")
}

# ======================== 公共绘图模块 ========================

# -------- Step 12: 绘制桑葚图网络 --------
cat("[Step 12] Plotting network (Sankey)...\n")
dynamic_base_name <- paste0(config$step12$base_name, "_", config$pipeline_mode)
step12_res <- tryCatch(
  plot_ltr_tf_gene_network(
    nodes_df   = nodes_df,
    edges_df   = edge_list,
    pipeline_mode   = config$pipeline_mode,
    output_dir      = config$output$step12,
    base_name       = paste0(config$step12$base_name, "_", config$pipeline_mode),
    width           = config$step12$plot_width,
    height          = config$step12$plot_height,
    dpi             = config$step12$plot_dpi,
    exclude_epitope = config$step12$exclude_epitope,
    strict_ppi_mode = config$step12$strict_ppi_mode,
    max_ltrs        = config$step12$max_plot_ltrs,
    max_tfs         = config$step12$max_plot_tfs,
    max_genes       = config$step12$max_plot_genes
  ),
  error = function(e) stop("Step12 failed: ", e$message)
)
cat("Step 12 done. Network plot saved successfully.\n\n")


# ======================== 收尾 ========================
end_time <- Sys.time()
cat("=== Analysis finished at:", end_time, "===\n")
cat("Total time:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")

# 保存会话信息
cat("\n--- Session Info ---\n")
writeLines(capture.output(sessionInfo()), file.path(config$output$logs, "sessionInfo.txt"))
cat("Session info saved.\n")

sink()
close(log_con)
