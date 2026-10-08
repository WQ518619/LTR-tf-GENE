# =========================================================
# Step 0: Load Required Packages
# =========================================================

# ---------------------------------------------------------
# Function:
#   load_required_packages
#
# Purpose:
#   Load all required R packages.
#
# Input:
#   None
#
# Output:
#   Loaded package environment
#
# Return:
#   NULL
# ---------------------------------------------------------

load_required_packages <- function(){
  
  pkg_list <- c(
    
    "dplyr",
    "stringr",
    "readr",
    "readxl",
    
    "GenomicRanges",
    "GenomeInfoDb",
    "IRanges",
    
    "LOLA",
    "rtracklayer",
    
    "clusterProfiler",
    "org.Hs.eg.db",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    
    "STRINGdb",
    "DESeq2",
    "biomaRt",
	"TFTF",
	# ==========================================
    # 【新增】Step 2b: 进化树与多序列比对依赖包
    # ==========================================
    "BSgenome.Hsapiens.UCSC.hg38",
    "Biostrings",
    "DECIPHER",
    "ape",
    "treeio",
    "ggtree",
    # ==========================================
    
    "igraph",
    "ggplot2",
    "ggforce",
    
    "tidyverse"
    
  )
  
  invisible(
    
    lapply(
      pkg_list,
      library,
      character.only = TRUE
    )
    
  )
  
  options(stringsAsFactors = FALSE)
  
  message("All required packages loaded.")
  
  invisible(pkg_list)
}