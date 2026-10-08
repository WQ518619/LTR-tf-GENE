# =========================================================
# Step 2b: LTR Phylogeny & Top Active Copy Extraction
# File:
#   R/step2b_ltr_phylogeny.R
# =========================================================

# ---------------------------------------------------------
# Function:
#   run_ltr_phylogeny
#
# Purpose:
#   1. 直接接收 Step 2 输出的 query_gr (GRanges 对象)
#   2. 与 Dfam 祖先文库比对，构建高精度进化树 (IQ-TREE)
#   3. 从每棵树中精确打捞 Top 20 最年轻拷贝
#   4. 生成可视化进化树
# ---------------------------------------------------------

#' 运行 LTR 亚家族定根发育树与活跃拷贝提取
#'
#' @param query_gr_obj Step 2 输出的 query_gr (GRanges对象)
#' @param dfam_fasta_file Dfam 全人类祖先文库 FASTA 文件的绝对路径
#' @param iqtree_path 本地 IQ-TREE 执行程序的绝对路径
#' @param top_n_youngest 整数。每个家族提取的最年轻序列数量，默认为 20
#' @param output_dir 字符串。输出目录名，默认为 "results/step2b_phylogeny/"
#' @param run_phylogeny 逻辑值。总开关，TRUE 为执行分析，FALSE 为跳过。
#'
#' @return 返回包含所有提取出的年轻拷贝的数据框 (Data.frame)
#' @export
run_ltr_phylogeny <- function(
    query_gr_obj,
    dfam_fasta_file,
    iqtree_path,
    top_n_youngest ,
    output_dir = "results/step2b_phylogeny/",
    run_phylogeny = TRUE
) {
  
  if (!run_phylogeny) {
    message("⏭️ [模块跳过] 用户设定 run_phylogeny = FALSE，已跳过 LTR 进化树分析。")
    return(NULL)
  }
  
  message("========================================")
  message("Step 2b: LTR Phylogeny & Active Copy Extraction")
  message("========================================")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # =====================================================
  # 1. 解析传入的 GRanges 对象为内部 BED 格式
  # =====================================================
  message("Converting GRanges to data.frame...")
  
  # 无缝对接：将内存中的 GRanges 转回 data.frame 进行序列提取
  bed <- as.data.frame(query_gr_obj) %>%
    dplyr::rename(chr = seqnames, repName = name) %>%
    dplyr::mutate(score = ".")
  
  # 过滤过短的片段 (>=200bp) 保障比对质量
  bed <- dplyr::filter(bed, end - start >= 200)
  message("Filtered fragments (>= 200bp): ", nrow(bed))
  
  # =====================================================
  # 2. 读取 Dfam 祖先序列文库并做交集过滤
  # =====================================================
  message("Reading Dfam consensus library...")
  all_human_consensus <- Biostrings::readDNAStringSet(dfam_fasta_file)
  names(all_human_consensus) <- sub("^[^ ]+ ", "", names(all_human_consensus))
  
  dfam_valid_names <- names(all_human_consensus)
  valid_families <- intersect(unique(bed$repName), dfam_valid_names)
  
  bed_matched <- bed %>% dplyr::filter(repName %in% valid_families)
  ltr_families <- unique(bed_matched$repName)
  
  if (length(ltr_families) == 0) {
    message("⚠️ 过滤后没有任何合格的 LTR 家族可供建树，分析终止。")
    return(NULL)
  }
  
  message("将针对以下 ", length(ltr_families), " 个上调亚家族建树与提取 Top ", top_n_youngest, " 序列：")
  message(paste(ltr_families, collapse=", "))
  
  # =====================================================
  # 3. 按亚家族分组，高精度建树与提取
  # =====================================================
  all_top20_merged <- data.frame()
  
  for (fam in ltr_families) {
    message("\n========== 处理亚家族： ", fam, " ==========")
    fam_clean <- gsub(":", "_", gsub("/", "_", fam)) 
    
    sub_bed <- bed_matched %>% dplyr::filter(repName == fam)
    
    if(nrow(sub_bed) < 3) {
      message("  ⚠️ 拷贝数不足 3 个，无法建树，跳过")
      next
    }
    
    gr <- GenomicRanges::GRanges(seqnames = sub_bed$chr,
                                 ranges = IRanges::IRanges(start = sub_bed$start + 1, end = sub_bed$end),
                                 strand = sub_bed$strand)
    seqs <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, gr)
    
    strand_clean <- ifelse(sub_bed$strand == "+", "plus", "minus")
    names(seqs) <- paste0(sub_bed$chr, "_", sub_bed$start, "_", sub_bed$end, "_", strand_clean)
    
    current_consensus <- all_human_consensus[names(all_human_consensus) == fam]
    root_name <- paste0("Ancestral_Consensus_", fam_clean) 
    names(current_consensus) <- root_name
    seqs <- c(seqs, current_consensus) 
    
    message("  多序列比对中 (100次迭代优化)...")
    aln <- DECIPHER::AlignSeqs(seqs, iterations = 100, refinements = 100, verbose = FALSE)
    aln_dna <- Biostrings::DNAStringSet(aln)
    
    aln_mat <- as.matrix(aln_dna)
    gap_frac <- apply(aln_mat, 2, function(x) sum(x == "-") / length(x))
    aln_trim <- aln_mat[, gap_frac <= 0.2, drop = FALSE] 
    
    if(ncol(aln_trim) < 10) {
      message("  ⚠️ 警告：修剪太严格导致序列长度不足10bp，跳过。")
      next
    }
    
    aln_trim_dna <- Biostrings::DNAStringSet(apply(aln_trim, 1, paste, collapse = ""))
    trim_file <- file.path(output_dir, paste0(fam_clean, "_trim.fa"))
    Biostrings::writeXStringSet(aln_trim_dna, trim_file)
    
    message("  运行 IQ-TREE 建树 (原装高精度模式)...")
    prefix_path <- file.path(output_dir, fam_clean)
    
    system2(iqtree_path, 
            args = c("-s", trim_file,
                     "-m", "MFP",
                     "-B", "1000",
                     "-alrt", "1000",
                     "-T", "AUTO",
                     "--prefix", prefix_path,
                     "-redo"),
            stdout = FALSE, stderr = FALSE) 
    
    tree_file <- paste0(prefix_path, ".treefile")
    
    if(file.exists(tree_file)) {
      tree <- treeio::read.iqtree(tree_file)
      phylo_tree <- tree@phylo 
      extracted_count <- 0
      
      if(root_name %in% phylo_tree$tip.label) {
        phylo_tree <- ape::root(phylo_tree, outgroup = root_name, resolve.root = TRUE)
        tree@phylo <- phylo_tree
        message("  🌳 已成功使用 Ancestral_Consensus 定根！")
        
        dist_matrix <- stats::cophenetic(phylo_tree)
        dist_from_root <- dist_matrix[root_name, ]
        
        sub_bed$strand_clean <- ifelse(sub_bed$strand == "+", "plus", "minus")
        sub_bed$tip_name <- paste0(sub_bed$chr, "_", sub_bed$start, "_", sub_bed$end, "_", sub_bed$strand_clean)
        sub_bed$divergence <- dist_from_root[sub_bed$tip_name]
        
        sub_top_active <- sub_bed %>% 
          dplyr::arrange(divergence) %>% 
          dplyr::slice_head(n = top_n_youngest)
        
        if(nrow(sub_top_active) > 0) {
          extracted_count <- nrow(sub_top_active)
          final_sub_bed <- sub_top_active %>% dplyr::select(chr, start, end, repName, score, strand)
          
          sub_bed_output <- file.path(output_dir, paste0(fam_clean, "_top", top_n_youngest, "_active.bed"))
          utils::write.table(final_sub_bed, sub_bed_output, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
          message(sprintf("  💾 成功导出 Top %d 条最年轻拷贝 -> %s", extracted_count, sub_bed_output))
          
          all_top20_merged <- rbind(all_top20_merged, final_sub_bed)
        }
      }
      
      num_tips <- length(phylo_tree$tip.label)
      dynamic_height <- max(12, num_tips * 0.25)
      
      p <- ggtree::ggtree(tree, layout = "rectangular", branch.length = 'branch.length') +
        ggtree::geom_tippoint(ggplot2::aes(color = (label == root_name)), size = 1.5) +
        ggtree::geom_tiplab(ggplot2::aes(color = (label == root_name), 
                                         fontface = ifelse(label == root_name, "bold", "plain")),
                            size = 2.5,          
                            align = TRUE,       # 保持对齐
                            linetype = "dotted", 
                            linesize = 0.3,      
                            offset = 0.05) +     
        ggtree::geom_treescale(x = 0, y = -1.5, width = 0.1, fontsize = 3, linesize = 0.5) +
        ggplot2::scale_color_manual(values = c("black", "red")) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.5))) + 
        ggplot2::labs(title = paste("LTR Evolutionary Phylogram:", fam),
                      subtitle = paste0("Left -> Right = Distance. Top ", extracted_count, " youngest copies exported."),
                      x = "Substitution distance")
      
      pdf_file <- file.path(output_dir, paste0(fam_clean, "_rectangular_tree.pdf"))
      ggplot2::ggsave(pdf_file, p, width = 18, height = dynamic_height, limitsize = FALSE)
      message("  ✅ 进化树已保存至：", pdf_file)
    } else {
      message("  ❌ 错误：IQ-TREE未能生成树文件")
    }
  }
  
  # =====================================================
  # 4. 结果汇总与返回数据给 Step 3
  # =====================================================
  if(nrow(all_top20_merged) > 0) {
    master_bed_output <- file.path(output_dir, paste0("All_Families_Top", top_n_youngest, "_Merged.bed"))
    utils::write.table(all_top20_merged, master_bed_output, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    message("\n👑 【Step 2b 成果】上调亚家族的最年轻序列已成功提取！")
    message("📊 共计筛出高质量年轻拷贝总数：", nrow(all_top20_merged), " 个")
    message("📂 终极合并 BED 文件保存在：", master_bed_output)
    message("Step 2b completed.")
    
    return(all_top20_merged)
  } else {
    message("\n⚠️ 运行结束，未提取到任何有效序列。")
    message("Step 2b completed.")
    return(NULL)
  }
}