# =========================================================
# Step 13: 获取 3D 染色质环化 
# File: R/step13_3d_loop.R
# =========================================================

run_3d_loop_verification <- function(
    query_gr,          # 前序步骤一直使用的 LTR 原生 GRanges
    final_candidates,  # Step 6 输出的靶基因列表
    genes_gr,          # 全基因组注释
    loop_bedpe_file,   # Hi-C 3D 环化数据
    upstream = 2000,   
    downstream = 500,  
    output_dir = "results/step13_3d_loops/"
) {
  
  cat("\n=================== Step 13: 3D 环化分析 ===================\n")
  
  if (nrow(final_candidates) == 0) {
    cat("⚠️ 候选靶基因列表为空，跳过 3D 环化验证。\n")
    return(NULL)
  }
  
  # =====================================================
  # 1. 准备 1D 区域 
  # =====================================================
  cat("1. 正在直接读取内存中的 LTR 与启动子坐标...\n")
  
  # LTR 坐标：直接使用内存中的 query_gr，并提取 name 列
  ltr_gr <- query_gr
  # 统一列名以适配后续代码，如果原列名是 name 则转为 ltr_name
  if("name" %in% colnames(S4Vectors::mcols(ltr_gr))) {
    ltr_gr$ltr_name <- ltr_gr$name 
  }
  
  # 启动子坐标：直接利用 final_candidates 中的基因 ID 极速提取
  target_entrez <- unique(as.character(final_candidates$gene_entrez))
  promoter_gr <- GenomicRanges::promoters(
    genes_gr[as.character(genes_gr$gene_id) %in% target_entrez],
    upstream = upstream,
    downstream = downstream
  )
  
  # 映射基因 Symbol
  match_idx <- match(promoter_gr$gene_id, final_candidates$gene_entrez)
  promoter_gr$gene_name <- final_candidates$target_gene[match_idx]
  
  # =====================================================
  # 2. 读取并解析 3D Loops 数据
  # =====================================================
  cat("2. 正在加载公共 3D 染色质环数据...\n")
  loop_data <- readr::read_delim(loop_bedpe_file, delim = "\t", comment = "#", col_names = FALSE, show_col_types = FALSE)
  
  colnames(loop_data)[1:12] <- c(
    "chr1", "start1", "end1", 
    "chr2", "start2", "end2", 
    "name_dot", "score_dot", "strand1_dot", "strand2_dot", "color_dot", 
    "observed_score"
  )
  
  loop_data <- loop_data %>% 
    dplyr::mutate(
      loop_id = paste0("Loop_", dplyr::row_number()),
      loop_span_bp = abs((start1 + end1)/2 - (start2 + end2)/2),
      observed_score = as.numeric(observed_score) 
    )
  
  anchor1_gr <- GenomicRanges::GRanges(seqnames = loop_data$chr1, ranges = IRanges::IRanges(start = loop_data$start1 + 1, end = loop_data$end1), loop_id = loop_data$loop_id)
  anchor2_gr <- GenomicRanges::GRanges(seqnames = loop_data$chr2, ranges = IRanges::IRanges(start = loop_data$start2 + 1, end = loop_data$end2), loop_id = loop_data$loop_id)
  
  # =====================================================
  # 3. 空间碰撞分析 (双向比对)
  # =====================================================
  cat("3. 正在执行 3D 空间碰撞...\n")
  
  # A 侧碰撞
  hit_LTR_A1 <- GenomicRanges::findOverlaps(ltr_gr, anchor1_gr)
  hit_Pro_A2 <- GenomicRanges::findOverlaps(promoter_gr, anchor2_gr)
  res_A <- dplyr::inner_join(
    data.frame(ltr_name = ltr_gr$ltr_name[S4Vectors::queryHits(hit_LTR_A1)], loop_id = anchor1_gr$loop_id[S4Vectors::subjectHits(hit_LTR_A1)]),
    data.frame(gene_name = promoter_gr$gene_name[S4Vectors::queryHits(hit_Pro_A2)], loop_id = anchor2_gr$loop_id[S4Vectors::subjectHits(hit_Pro_A2)]),
    by = "loop_id",
    relationship = "many-to-many" # 👈 核心修复 1：消除多对多警告
  ) %>% dplyr::mutate(Loop_Direction = "LTR_Left_Pro_Right")
  
  # B 侧碰撞
  hit_LTR_A2 <- GenomicRanges::findOverlaps(ltr_gr, anchor2_gr)
  hit_Pro_A1 <- GenomicRanges::findOverlaps(promoter_gr, anchor1_gr)
  res_B <- dplyr::inner_join(
    data.frame(ltr_name = ltr_gr$ltr_name[S4Vectors::queryHits(hit_LTR_A2)], loop_id = anchor2_gr$loop_id[S4Vectors::subjectHits(hit_LTR_A2)]),
    data.frame(gene_name = promoter_gr$gene_name[S4Vectors::queryHits(hit_Pro_A1)], loop_id = anchor1_gr$loop_id[S4Vectors::subjectHits(hit_Pro_A1)]),
    by = "loop_id",
    relationship = "many-to-many" # 👈 核心修复 2：消除多对多警告
  ) %>% dplyr::mutate(Loop_Direction = "LTR_Right_Pro_Left")
  
  # =====================================================
  # 4. 整合与结果输出
  # =====================================================
  final_edge_table <- dplyr::bind_rows(res_A, res_B) %>%
    dplyr::distinct() %>%
    dplyr::left_join(loop_data %>% dplyr::select(loop_id, observed_score, loop_span_bp, chr1, start1, end1, chr2, start2, end2), by = "loop_id") %>%
    dplyr::arrange(dplyr::desc(observed_score)) %>%
    dplyr::select(LTR_Name = ltr_name, Target_Gene = gene_name, loop_id, Loop_Direction, observed_score, loop_span_bp, dplyr::everything())
  
  if(nrow(final_edge_table) > 0) {
    # 👈 核心修复 3：安全兜底机制，防止外部传进来的 output_dir 为空导致 dir.create 崩溃
    if (is.null(output_dir) || is.na(output_dir) || output_dir == "") {
      output_dir <- "results/step13_3d_loops/"
    }
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    readr::write_csv(final_edge_table, file.path(output_dir, "HepG2_HiC_Verified_Loops.csv"))
    
    bedpe_out <- final_edge_table %>% dplyr::select(chr1, start1, end1, chr2, start2, end2, loop_id, observed_score)
    readr::write_delim(bedpe_out, file.path(output_dir, "HepG2_HiC_Verified_Edges.bedpe"), delim = "\t", col_names = FALSE)
    
    cat("🎉 完美！直接利用前期内存对象，完成 3D 物理验证。\n")
  } else {
    cat("\n⚠️ 未发现显著重叠的 3D 染色质环。\n")
  }
  
  return(final_edge_table)
}