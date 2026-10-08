# =========================================================
# Step 6: LTR 侧翼窗口 + 邻近上调基因分析
# File:
# R/step6_nearby_genes.R
# =========================================================

# ---------------------------------------------------------
# Function:
# find_nearby_up_genes
#
# Purpose:
#   寻找 LTR 侧翼一定距离内的上调基因。
#
# Input:
#
# query_gr
#   LTR 的 GRanges 对象
#
# DEG_file
#   差异表达结果文件
#   支持:
#     .csv
#     .tsv
#     .txt
#     .xlsx
#
# gene_col
#   基因Symbol列名
#
# lfc_col
#   logFC列名
#
# padj_col
#   adjusted pvalue列名
#
# window_kb
#   LTR侧翼窗口大小（kb）
#
# lfc_cutoff
#   上调基因logFC阈值
#
# padj_cutoff
#   adjusted pvalue阈值
#
# output_dir
#   输出目录
#
# Return:
#
# list:
#
# final_candidates
#   最终LTR-邻近上调基因关系
#
# all_nearby
#   所有邻近基因
#
# up_gene_count
#   上调基因数量
#
# total_pairs
#   LTR-Gene关系数量
# ---------------------------------------------------------

find_nearby_up_genes <- function(
    
  query_gr,
  
  DEG_file,
  
  gene_col = "Gene_Symbol",
  
  lfc_col = "log2FoldChange",
  
  padj_col = "padj",
  
  window_kb ,
  
  lfc_cutoff ,
  
  padj_cutoff ,
  
  output_dir = "results/step6/"
  
){
  
  cat(
    "【Step 7】开始LTR邻近上调基因分析\n"
  )
  
  # =====================================================
  # 1. 读取DEG文件
  # =====================================================
  
  cat("读取DEG文件...\n")
  
  if(grepl("\\.csv$", DEG_file)){
    
    res_annotated <- read.csv(
      DEG_file,
      stringsAsFactors = FALSE
    )
    
  } else if(
    grepl("\\.tsv$", DEG_file) |
    grepl("\\.txt$", DEG_file)
  ){
    
    res_annotated <- read.delim(
      DEG_file,
      stringsAsFactors = FALSE
    )
    
  } else if(
    grepl("\\.xlsx$", DEG_file)
  ){
    
    res_annotated <- readxl::read_excel(
      DEG_file
    )
    
  } else {
    
    stop(
      "不支持的DEG文件格式"
    )
  }
  
  cat(
    "DEG数据读取完成:",
    nrow(res_annotated),
    "行\n"
  )
  
  # =====================================================
  # 2. 构建LTR侧翼窗口
  # =====================================================
  
  cat(
    "构建LTR窗口: ±",
    window_kb,
    "kb\n"
  )
  
  ltr_window <- GRanges(
    
    seqnames = seqnames(query_gr),
    
    ranges = IRanges(
      
      start = pmax(
        1,
        start(query_gr) - window_kb * 1000
      ),
      
      end = end(query_gr) +
        window_kb * 1000
    ),
    
    strand = strand(query_gr)
  )
  
  # 保留LTR metadata
  mcols(ltr_window) <- mcols(query_gr)
  
  # =====================================================
  # 3. 加载基因注释
  # =====================================================
  
  cat(
    "加载TxDb基因注释...\n"
  )
  
  genes_gr <- genes(
    TxDb.Hsapiens.UCSC.hg38.knownGene
  )
  
  seqlevelsStyle(genes_gr) <- "UCSC"
  
  genes_gr <- keepStandardChromosomes(
    genes_gr,
    pruning.mode = "coarse"
  )
  
  # 提取TSS
  gene_tss_gr <- resize(
    genes_gr,
    width = 1,
    fix = "start"
  )
  
  # =====================================================
  # 4. overlap分析
  # =====================================================
  
  cat(
    "寻找LTR邻近基因...\n"
  )
  
  hits <- findOverlaps(
    
    ltr_window,
    
    gene_tss_gr,
    
    ignore.strand = FALSE
  )
  
  if(length(hits) == 0){
    
    cat(
      "警告: 未发现邻近基因\n"
    )
    
    return(NULL)
  }
  
  ltr_hits  <- ltr_window[
    queryHits(hits)
  ]
  
  gene_hits <- genes_gr[
    subjectHits(hits)
  ]
  
  # =====================================================
  # 5. 计算LTR与基因距离
  # =====================================================
  
  cat(
    "计算LTR-TSS距离...\n"
  )
  
  ltr_center <- (
    start(ltr_hits) +
      end(ltr_hits)
  ) / 2
  
  gene_tss_pos <- start(
    
    resize(
      gene_hits,
      width = 1,
      fix = "start"
    )
  )
  
  dist_val <- ifelse(
    
    strand(gene_hits) == "+",
    
    gene_tss_pos - ltr_center,
    
    ltr_center - gene_tss_pos
  )
  
  # =====================================================
  # 6. 构建LTR-Gene关系表
  # =====================================================
  
  cat(
    "构建LTR-Gene关系表...\n"
  )
  
  ltr_gene_rel <- data.frame(
    
    chr = as.character(
      seqnames(ltr_hits)
    ),
    
    LTR_start = start(ltr_hits),
    
    LTR_end = end(ltr_hits),
    
    LTR_subfamily = mcols(
      ltr_hits
    )$name,
    
    gene_entrez = as.character(
      mcols(gene_hits)$gene_id
    ),
    
    distance = round(
      dist_val
    )
    
  ) %>%
    
    inner_join(
      
      bitr(
        
        .$gene_entrez,
        
        fromType = "ENTREZID",
        
        toType = "SYMBOL",
        
        OrgDb = org.Hs.eg.db
        
      ),
      
      by = c(
        "gene_entrez" = "ENTREZID"
      )
      
    ) %>%
    
    rename(
      target_gene = SYMBOL
    ) %>%
    
    distinct()
  
  cat(
    "发现邻近基因:",
    length(
      unique(ltr_gene_rel$target_gene)
    ),
    "\n"
  )
  
  # =====================================================
# 7. 提取上调基因
# =====================================================

cat("筛选上调基因...\n")

# ⭐ 先确保两列是数值型
res_annotated <- res_annotated %>%
  mutate(
    lfc_val  = as.numeric(.data[[lfc_col]]),
    padj_val = as.numeric(.data[[padj_col]])
  )

# ⭐ 同时确保阈值也是数值型
lfc_cutoff  <- as.numeric(lfc_cutoff)
padj_cutoff <- as.numeric(padj_cutoff)

# 检查类型
cat("lfc 类型:", class(res_annotated$lfc_val), "\n")
cat("padj 类型:", class(res_annotated$padj_val), "\n")
cat("lfc_cutoff:", lfc_cutoff, "  padj_cutoff:", padj_cutoff, "\n")

# 筛选
up_genes_list <- res_annotated %>%
  filter(
    !is.na(lfc_val),
    !is.na(padj_val),
    lfc_val > lfc_cutoff,
    padj_val < padj_cutoff
  ) %>%
  pull(.data[[gene_col]]) %>%
  unique()

cat("上调基因数量:", length(up_genes_list), "\n")
  
  # =====================================================
  # 8. 提取邻近上调基因
  # =====================================================
  
  final_candidates <- ltr_gene_rel %>%
    
    filter(
      target_gene %in% up_genes_list
    ) %>%
    
    mutate(
      
      abs_dist = abs(distance),
      
      distance_category = case_when(
        
        abs_dist <= 10000
        ~ "Very Close",
        
        abs_dist <= 50000
        ~ "Close",
        
        TRUE
        ~ "Distal"
      )
    )
  
  cat(
    "最终LTR邻近上调基因数量:",
    length(
      unique(
        final_candidates$target_gene
      )
    ),
    "\n"
  )
  
  cat(
    "LTR-Gene关系总数:",
    nrow(final_candidates),
    "\n"
  )
  
  # =====================================================
  # 9. 保存结果
  # =====================================================
  
  if(!is.null(output_dir)){
    
    cat(
      "保存结果...\n"
    )
    
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    write.csv(
      
      ltr_gene_rel,
      
      file.path(
        output_dir,
        "all_ltr_nearby_genes.csv"
      ),
      
      row.names = FALSE
    )
    
    write.csv(
      
      final_candidates,
      
      file.path(
        output_dir,
        "final_ltr_up_genes.csv"
      ),
      
      row.names = FALSE
    )
  }
  
  # =====================================================
  # 10. 返回结果
  # =====================================================
  
  return(
    
    list(
      
      final_candidates =
        final_candidates,
      
      all_nearby =
        ltr_gene_rel,
      
      up_gene_count =
        length(
          unique(
            final_candidates$target_gene
          )
        ),
      
      total_pairs =
        nrow(
          final_candidates
        )
    )
  )
}