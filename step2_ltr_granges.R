# =========================================================
# Step 2: Build LTR GRanges
# File:
#   R/step2_ltr_granges.R
# =========================================================

# ---------------------------------------------------------
# Function:
#   build_ltr_granges
#
# Purpose:
#   Construct:
#     1. query_gr      -> 上调 LTR
#     2. universe_gr   -> 所有表达的 LTR（优化背景）
#
# Input:
#
#   rmsk_file
#       RepeatMasker annotation file
#
#   te_file
#       TEtranscripts result file
#
#   ltr_up
#       Significant upregulated LTR families
#
#   expressed_baseMean
#       Expression threshold for universe
#
#   output_dir
#       Output directory
#
# Return:
#   list(
#       query_gr,
#       universe_gr,
#       rmsk_parsed
#   )
# ---------------------------------------------------------

# =========================================================
# Step 2: Build LTR GRanges (修正版)
# File: R/step2_ltr_granges.R
# =========================================================

build_ltr_granges <- function(
    rmsk_file,
    te_file,
    ltr_up,
    expressed_baseMean,
    output_dir = "results/step2/"
){
  
  message("========================================")
  message("Building LTR GRanges")
  message("========================================")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # =====================================================
  # 1. Read RepeatMasker annotation
  # =====================================================
  
  message("Reading RepeatMasker annotation...")
  
  rmsk <- read.table(
    rmsk_file,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  colnames(rmsk) <- c(
    "chr", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
  )
  
  # =====================================================
  # 2. Parse RepeatMasker attributes
  # =====================================================
  
  message("Parsing RepeatMasker attributes...")
  
  rmsk_parsed <- rmsk %>%
    dplyr::mutate(
      gene_id   = stringr::str_match(attribute, 'gene_id "?([^";]+)"?')[,2],
      family_id = stringr::str_match(attribute, 'family_id "?([^";]+)"?')[,2],
      class_id  = stringr::str_match(attribute, 'class_id "?([^";]+)"?')[,2]
    )
  
  message("Parsed rows: ", nrow(rmsk_parsed))
  
  # =====================================================
  # 3. Build query_gr
  #    上调 LTR
  # =====================================================
  
  message("Constructing query_gr...")
  
  query_regions <- rmsk_parsed %>%
    dplyr::filter(
      class_id == "LTR",
      gene_id %in% ltr_up
    )
  
  query_gr <- GenomicRanges::GRanges(
    seqnames = query_regions$chr,
    ranges   = IRanges::IRanges(
      start = query_regions$start,
      end   = query_regions$end
    ),
    strand = query_regions$strand,
    name   = query_regions$gene_id
  )
  
  GenomeInfoDb::seqlevelsStyle(query_gr) <- "UCSC"
  
  query_gr <- GenomeInfoDb::keepStandardChromosomes(
    query_gr,
    pruning.mode = "coarse"
  )
  
  query_gr <- unique(query_gr)
  
  message("query_gr size: ", length(query_gr))
  
  # =====================================================
  # 4. Read TEtranscripts result
  #    ⭐ 新增：显式转换为数值型
  # =====================================================
  
  message("Reading TEtranscripts result...")
  
  te_all <- readxl::read_excel(te_file)
  
  # ⭐ 将所有可能参与比较的数值列显式转换为 numeric
  te_all <- te_all %>%
    dplyr::mutate(
      baseMean       = suppressWarnings(as.numeric(baseMean)),
      log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
      padj           = suppressWarnings(as.numeric(padj)),
      pvalue         = if ("pvalue" %in% colnames(.)) {
        suppressWarnings(as.numeric(pvalue))
      } else NA_real_,
      stat           = if ("stat" %in% colnames(.)) {
        suppressWarnings(as.numeric(stat))
      } else NA_real_
    )
  
  # 检查类型确认
  message("baseMean class: ", class(te_all$baseMean))
  message("log2FoldChange class: ", class(te_all$log2FoldChange))
  message("padj class: ", class(te_all$padj))
  
  # =====================================================
  # 5. Extract expressed LTR families
  # =====================================================
  
  message("Extracting expressed LTR families...")
  
  expressed_ltr <- te_all %>%
    dplyr::filter(
      stringr::str_detect(ID, "LTR"),
      !is.na(baseMean),
      baseMean > expressed_baseMean
    ) %>%
    dplyr::mutate(
      LTR_family = stringr::str_extract(ID, "^[^:]+")
    ) %>%
    dplyr::pull(LTR_family) %>%
    unique()
  
  expressed_ltr <- unique(as.character(expressed_ltr))
  expressed_ltr <- expressed_ltr[expressed_ltr != "" & !is.na(expressed_ltr)]
  
  message("Expressed LTR family count: ", length(expressed_ltr))
  
  # =====================================================
  # 6. Build optimized universe_gr
  #    仅使用表达的 LTR
  # =====================================================
  
  message("Constructing optimized universe_gr...")
  
  universe_regions <- rmsk_parsed %>%
    dplyr::filter(
      class_id == "LTR",
      gene_id %in% expressed_ltr,
      !is.na(gene_id)
    )
  
  universe_gr <- GenomicRanges::GRanges(
    seqnames = universe_regions$chr,
    ranges   = IRanges::IRanges(
      start = universe_regions$start,
      end   = universe_regions$end
    ),
    strand = universe_regions$strand,
    name   = universe_regions$gene_id
  )
  
  GenomeInfoDb::seqlevelsStyle(universe_gr) <- "UCSC"
  
  universe_gr <- GenomeInfoDb::keepStandardChromosomes(
    universe_gr,
    pruning.mode = "coarse"
  )
  
  universe_gr <- unique(universe_gr)
  
  message("universe_gr size: ", length(universe_gr))
  
  # =====================================================
  # 7. Save intermediate files
  # =====================================================
  
  saveRDS(query_gr,    file.path(output_dir, "query_gr.rds"))
  saveRDS(universe_gr, file.path(output_dir, "universe_gr.rds"))
  
  write.csv(
    data.frame(expressed_ltr = expressed_ltr),
    file.path(output_dir, "expressed_ltr_families.csv"),
    row.names = FALSE
  )
  
  # =====================================================
  # 8. Return
  # =====================================================
  
  message("Step2 completed.")
  
  return(list(
    query_gr       = query_gr,
    universe_gr    = universe_gr,
    rmsk_parsed    = rmsk_parsed,
    expressed_ltr  = expressed_ltr
  ))
}