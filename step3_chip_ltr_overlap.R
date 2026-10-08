# =========================================================
# Step 3: ChIP-Atlas 处理 + LTR-TF 位点级 overlap
# File:
# R/step3_chip_ltr_overlap.R
# =========================================================
# ---------------------------------------------------------
# Function:
# process_chip_atlas
#
# Purpose:
# 处理 ChIP-Atlas 文件，提取 TF 信息并构建 GRanges 对象。
#
# Input:
#
# chip_file
#   ChIP-Atlas BED/TSV 文件路径。
#
# score_cutoff
#   峰强度过滤阈值，默认 300。
#
# output_dir
#   可选，保存处理后文件的目录。
#
# Output:
# chip_df 和 chip_gr。
#
# Return:
# list
#
# Example:
# chip_res <- process_chip_atlas(
#     chip_file = "data/chip_atlas.tsv",
#     score_cutoff = 300,
#     output_dir = "results/step3/"
# )
# ---------------------------------------------------------
process_chip_atlas <- function(chip_file,
                               score_cutoff ,
                               output_dir = "results/step3/") {
  
  cat("【Step 3.1】正在处理 ChIP-Atlas 数据...\n")
  
  chip <- read_tsv(chip_file, col_names = FALSE, comment = "track")
  chip <- chip[, 1:5]
  
  colnames(chip) <- c("chr", "start", "end", "info", "score")
  
  chip$TF <- str_extract(chip$info, "Name=[^;]+") %>%
    str_replace("Name=", "") %>%
    str_replace("%20.*", "") %>%
    str_trim()
  
  chip <- chip %>% filter(score >= score_cutoff)
  
  chip_gr <- GRanges(
    seqnames = chip$chr,
    ranges = IRanges(start = chip$start, end = chip$end),
    TF = chip$TF,
    score = chip$score
  )
  
  seqlevelsStyle(chip_gr) <- "UCSC"
  chip_gr <- keepStandardChromosomes(chip_gr, pruning.mode = "coarse")
  chip_gr <- unique(chip_gr)
  
  cat("ChIP-Atlas 处理完成！共", length(chip_gr), "个峰，涉及",
      n_distinct(chip$TF), "个 TF (score ≥", score_cutoff, ")\n")
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(chip_gr, file.path(output_dir, "chip_gr.rds"))
    write.csv(chip, file.path(output_dir, "chip_processed.csv"), row.names = FALSE)
  }
  
  return(list(
    chip_df = chip,
    chip_gr = chip_gr,
    tf_count = n_distinct(chip$TF),
    peak_count = nrow(chip)
  ))
}

# ---------------------------------------------------------
# Function:
# get_ltr_tf_binding
#
# Purpose:
# 执行 LTR 与 ChIP-Atlas TF 峰的位点级 overlap 分析。
#
# Input:
#
# query_gr
#   LTR GRanges（来自 Step 3）。
#
# chip_gr
#   ChIP-Atlas GRanges。
#
# min_peak_score
#   可选，进一步过滤峰强度。
#
# output_dir
#   可选，保存结果目录。
#
# Return:
# data.frame 或 NULL
# ---------------------------------------------------------
get_ltr_tf_binding <- function(query_gr,
                               chip_gr,
                               min_peak_score ,
                               output_dir = "results/step3/") {
  
  cat("【Step 3.2】执行 LTR-TF 位点级 overlap...\n")
  
  if (!is.null(min_peak_score)) {
    chip_gr <- chip_gr[mcols(chip_gr)$score >= min_peak_score]
  }
  
  hits <- findOverlaps(query_gr, chip_gr, ignore.strand = TRUE)
  
  if (length(hits) == 0) {
    cat("警告：未发现任何 LTR-TF 重叠！\n")
    return(NULL)
  }
  
  ltr_tf_binding <- data.frame(
    chr = as.character(seqnames(query_gr))[queryHits(hits)],
    LTR_start = start(query_gr)[queryHits(hits)],
    LTR_end = end(query_gr)[queryHits(hits)],
    LTR = mcols(query_gr)$name[queryHits(hits)],
    TF = mcols(chip_gr)$TF[subjectHits(hits)],
    peak_score = mcols(chip_gr)$score[subjectHits(hits)],
    peak_start = start(chip_gr)[subjectHits(hits)],
    peak_end = end(chip_gr)[subjectHits(hits)]
  ) %>%
    distinct() %>%
    arrange(LTR, desc(peak_score))
  
  cat("LTR-TF binding 关系数量:", nrow(ltr_tf_binding), "\n")
  cat("涉及 LTR 数量:", n_distinct(ltr_tf_binding$LTR), "\n")
  cat("涉及 TF 数量:", n_distinct(ltr_tf_binding$TF), "\n")
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(ltr_tf_binding, file.path(output_dir, "ltr_tf_binding.csv"), row.names = FALSE)
  }
  
  return(ltr_tf_binding)
}