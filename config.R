# R/config.R
# =========================================================
# Global configuration for LTR-TF-Gene pipeline
# =========================================================

library(here)
config <- list(
  # =========================================================
  # 全局分析模式控制 (Dual-Mode Switch)
  # =========================================================
  # "reverse" : 有差异基因（DEG）时的逆向验证模式（执行原有的所有分析）
  # "forward" : 无差异基因时的正向探索模式（直接通过 TFTF 预测靶基因）
  pipeline_mode = "reverse", 
  # =====================================================
  # Project Information
  # =====================================================
  project = list(
    
    project_name = "LTR_TF_Gene_Pipeline",
    
    species = "human",
    
    taxonomy_id = 9606,
    
    genome = "hg38",
    
    seed = 1234
    
  ),
  
  # =====================================================
  # Input Files
  # =====================================================
  paths = list(
    
    # TEtranscripts differential TE result
    te_file = here(
      "data",
      "TEtranscripts_out_gene_TE_analysis30.xlsx"
    ),
    
    # RepeatMasker annotation
    rmsk_file = here(
      "data",
      "GRCh38_rmsk_TE.gtf"
    ),
    
    # DEG result
    deg_file = here(
      "data",
      "all_genes_annotated_filtered.xlsx"
    ),
    
    # ChIP-Atlas raw peaks
    chip_file = here(
      "data",
      "Oth.Liv.20.AllAg.AllCell.bed"
    ),
    
    # LOLA database root
    lola_base_dir = here(
      "LOLA_DB"
    ),
	# Dfam 祖先共识序列库 FASTA 文件
    dfam_fasta = here(
      "data",
      "dfam-fasta-download.fasta"
    ),
    
    # 本地 IQ-TREE 执行程序路径
    iqtree_exe = here(
      "data",
      "iqtree3.exe"
    ),
	
	#3D Hi-C 环化数据
	hic_bedpe = here(
      "data",
      "ENCFF463OJZ.bedpe"
    )
  ),
  
  # =====================================================
  # Step 1: Load Upregulated LTR Families
  # =====================================================
  step1 = list(
    
    lfc_cutoff = 0.585,
    
    padj_cutoff = 0.05
  ),
  
  # =====================================================
  # Step 2: Build LTR GRanges
  # =====================================================
  step2 = list(
    
    expressed_baseMean = 5
  ),
  # =====================================================
  # Step 2b: LTR Phylogeny & Active Copy Extraction
  # =====================================================
  step2b = list(
    
    # 模块总开关：TRUE 为执行进化树分析，FALSE 为跳过
    run_phylogeny = FALSE,
    
    # 每个亚家族打捞的最年轻活跃拷贝数量
    top_n_youngest = 20
  ),
  
  # =====================================================
  # Step 3: ChIP-Atlas Processing
  # =====================================================
  step3 = list(
    
    chip_score_cutoff = 300,
    
    min_peak_score = 300	
  ),
  
  # =====================================================
  # Step 4: Create LOLA Database
  # =====================================================
  step4 = list(
    
    lola_collection_name = "ChIPAtlas_Liver",
    
    min_peaks_per_tf = 100
  ),
  
  # =====================================================
  # Step 5: LTR LOLA Enrichment
  # =====================================================
  step5 = list(
    
    minOverlap = 1,
    
    qval_cutoff    = 0.05,
    
    support_cutoff = 10,
    cores = 1
  ),
  
  # =====================================================
  # Step 6: Nearby Upregulated Genes
  # =====================================================
  step6 = list(
    
    window_kb = 100,
    
    gene_col = "Gene_Symbol",
    
    lfc_col = "log2FoldChange",
    
    padj_col = "padj",
    
    lfc_cutoff = 0.585,
    
    padj_cutoff = 0.05
  ),
  
# =========================================================
# Step 6b: 正向探索模式专属参数 (Forward Prediction)
# =========================================================
  step6b = list(
    min_tf_count   = 2,    # 某个靶基因至少需要被几个 LTR-TF 共同调控
    min_db_support = 7    # TF-Gene 的调控连线至少需要命中的数据库数量
), 
# =========================================================
# Step 7: 靶基因侧核心转录因子预测 (基于 TFTF)
# =========================================================
step7 = list(
  db_list        = c("hTFtarget", "KnockTF", "FIMO_JASPAR", "PWMEnrich_JASPAR", 
                     "ENCODE", "CHEA", "TRRUST", "GTRD", "ChIP_Atlas"),
  fimo_score = 10,          
  pwm_pval   = 0.1,         
  cut_log2fc = 1,           
  down_only  = TRUE,  
  min_db_support =8  
  
),
  
  # =====================================================
  # Step 8: Final Integration
  # =====================================================
  step8 = list(
    
    enable_bridge_tf = TRUE
  ),
  
  # =====================================================
  # Step 9: STRINGdb PPI
  # =====================================================
  step9 = list(
    
    string_version = "11.5",
    
    score_threshold = 700,
    
    remove_isolated = TRUE,
	string_data_dir = here("data")
  ),
  
  # =====================================================
  # Step 10: Gene-TF Binding
  # =====================================================
  step10 = list(
    
    upstream = 2000,
    
    downstream = 500
  ),
  
  # =====================================================
  # Step 11: Network Construction
  # =====================================================
  step11 = list(
    
    top_ltr_tf_per_family = 3,
    
    min_ltr_peak_score = 500
  ),
  
  # =====================================================
  # Step 12: Network Visualization
  # =====================================================
  step12 = list(
    
    plot_width = 18,
    
    plot_height = 14,
    
    plot_dpi = 300,
    
    base_name = "LTR_TF_Gene_Network",
    
    save_pdf = TRUE,
    
    save_png = TRUE,
    
    save_svg = TRUE,
	exclude_epitope = TRUE,      # 是否排除 Epitope 等非特异性标签
    strict_ppi_mode = TRUE,     # 正向分析中建议设为 FALSE，逆向有 PPI 时设为 TRUE
    max_plot_ltrs   = 20,        # 限制网络图左侧最多显示几个 LTR
    max_plot_tfs    = Inf,        # 限制网络图中间最多显示几个 TF
    max_plot_genes  = 20         # 限制网络图右侧最多显示几个靶基因
  ),
  
  #  Step 13 
 step13 = list(
   upstream   = 2000,   # 定义靶基因 Promoter 向上游延伸的距离（用于 3D 空间碰撞）
   downstream = 500     # 定义靶基因 Promoter 向下游延伸的距离
),
  
  # =====================================================
  # Parallel Settings
  # =====================================================
  parallel = list(
    
    workers = 4
  ),
  
  # =====================================================
  # Cache Settings
  # =====================================================
  cache = list(
    
    use_cache = TRUE,
    
    overwrite = FALSE
  ),
  
  # =====================================================
  # Output Directories
  # =====================================================
  output = list(
    
    step1 = here(
      "results",
      "step1"
    ),
    
    step2 = here(
      "results",
      "step2"
    ),
	
	step2b = here(
      "results",
      "step2b"
    ),
    
    step3 = here(
      "results",
      "step3"
    ),
    
    step4 = here(
      "results",
      "step4"
    ),
    
    step5 = here(
      "results",
      "step5"
    ),
    
    step6 = here(
      "results",
      "step6"
    ),
	
	step6b = here(
      "results",
      "step6b"
    ),
    
    step7 = here(
      "results",
      "step7"
    ),
    
    step8 = here(
      "results",
      "step8"
    ),
    
    step9 = here(
      "results",
      "step9"
    ),
    
    step10 = here(
      "results",
      "step10"
    ),
    
    step11 = here(
      "results",
      "step11"
    ),
    
    step12 = here(
      "results",
      "step12"
    ),
    
    logs = here(
      "logs"
    ),
    
    cache = here(
      "cache"
    )
  )
)

# =========================================================
# Initialize global settings
# =========================================================

set.seed(
  config$project$seed
)

# =========================================================
# Auto-create output directories
# =========================================================

dir_list <- unlist(
  config$output
)

invisible(
  lapply(
    dir_list,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
)