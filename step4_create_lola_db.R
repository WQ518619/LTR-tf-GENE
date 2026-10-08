# =========================================================
# Step 4: 创建 LOLA 数据库（ChIP-Atlas Liver TF）
# File:
# R/step4_create_lola_db.R
# =========================================================
# ---------------------------------------------------------
# Function:
# create_lola_db
#
# Purpose:
# 将 ChIP-Atlas 处理结果拆分为按 TF 的 BED 文件，并创建 LOLA 数据库索引。
#
# Input:
#
# chip_df
#   来自 process_chip_atlas 的 chip_df。
#
# base_dir
#   LOLA 数据库根目录。
#
# coll_name
#   集合名称。
#
# min_peaks_per_tf
#   每个 TF 至少保留的峰数量。
#
# output_dir
#   可选，保存备份目录。
#
# Return:
# list
# ---------------------------------------------------------
create_lola_db <- function(chip_df,
                           base_dir = "LOLA_DB",
                           coll_name = "ChIPAtlas_Liver",
                           min_peaks_per_tf ,
                           output_dir = "results/step4/") {
  
  cat("【Step 4】开始创建 LOLA 数据库...\n")
  
  regions_dir <- file.path(base_dir, coll_name, "regions")
  dir.create(regions_dir, recursive = TRUE, showWarnings = FALSE)
  
  tf_list <- unique(chip_df$TF)
  cat("共检测到", length(tf_list), "个 TF，正在按 TF 拆分 BED 文件...\n")
  
  saved_count <- 0
  
  for(tf in tf_list){
    tf_df <- chip_df %>%
      filter(TF == tf) %>%
      dplyr::select(chr, start, end) %>%
      distinct()
    
    if(nrow(tf_df) < min_peaks_per_tf) next
    
    bed_file <- file.path(regions_dir, paste0(tf, ".bed"))
    write.table(tf_df, bed_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    
    saved_count <- saved_count + 1
  }
  
  cat("成功保存", saved_count, "个 TF 的 BED 文件（峰数量 ≥", min_peaks_per_tf, "）\n")
  
  bed_files <- list.files(regions_dir, pattern = "\\.bed$", full.names = FALSE)
  
  index_df <- data.frame(filename = bed_files,
                         description = gsub(".bed", "", bed_files))
  write.table(index_df, file.path(base_dir, coll_name, "index.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  collection_df <- data.frame(
    collection = coll_name,
    description = "Liver TF ChIP-seq from ChIP-Atlas"
  )
  write.table(collection_df, file.path(base_dir, "collection.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  regionSets_df <- data.frame(
    collection = coll_name,
    filename = bed_files,
    description = gsub(".bed", "", bed_files)
  )
  write.table(regionSets_df, file.path(base_dir, "regionSets.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  regionDB <- loadRegionDB(base_dir)
  
  cat("✅ LOLA 数据库创建完成！\n")
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(regionDB, file.path(output_dir, "regionDB.rds"))
    write.csv(chip_df, file.path(output_dir, "chip_for_lola.csv"), row.names = FALSE)
  }
  
  return(list(
    regionDB = regionDB,
    coll_name = coll_name,
    base_dir = base_dir,
    tf_count = length(bed_files),
    regions_dir = regions_dir
  ))
}