# =========================================================
# Step 11: 提取网络节点与边
# File:
# R/step11_build_network.R
# =========================================================

# ---------------------------------------------------------
# Function:
# build_network_tables
#
# Purpose:
# 构建 Cytoscape/networkD3 网络：
#
# 1. LTR -> TF
# 2. TF -> TF (PPI)
# 3. TF -> Gene
#
# 并输出：
#
# - node_attributes
# - edge_list
#
# ---------------------------------------------------------

build_network_tables <- function(
    
  ltr_tf_binding,
  gene_tf_binding,
  cross_ppi,
  
  sig_tf,
  sig_tf_promoter,
  
  top_ltr_tf_per_family ,
  min_ltr_peak_score ,
  
  output_dir = "results/step11/"
) {
  
  cat("\n=================== 构建网络节点与边 ===================\n")
  
  # =====================================================
  # Step 1. 显著 TF
  # =====================================================
  
  sig_ltr_tfs <- unique(sig_tf$TF)
  
  sig_promoter_tfs <- unique(
    sig_tf_promoter$TF
  )
  
  # =====================================================
  # Step 2. LTR -> TF
  # =====================================================
  
  edges_ltr_tf <- ltr_tf_binding %>%
    
    filter(
      TF %in% sig_ltr_tfs,
      peak_score >= min_ltr_peak_score
    ) %>%
    
    group_by(LTR) %>%
    
    slice_max(
      order_by = peak_score,
      n = top_ltr_tf_per_family,
      with_ties = FALSE
    ) %>%
    
    ungroup() %>%
    
    transmute(
      
      from = as.character(LTR),
      
      to = as.character(TF),
      
      type = "LTR-TF",
      
      weight = peak_score / 1000
      
    ) %>%
    
    distinct()
  
  cat("LTR-TF 边数:",
      nrow(edges_ltr_tf), "\n")
  
  # =====================================================
  # Step 3. TF -> TF (PPI)
  # =====================================================
  
  edges_ppi <- cross_ppi %>%
    
    transmute(
      
      from = as.character(from_gene),
      
      to = as.character(to_gene),
      
      type = "PPI",
      
      weight = combined_score / 1000
      
    ) %>%
    
    distinct()
  
  cat("PPI 边数:",
      nrow(edges_ppi), "\n")
  
  # =====================================================
  # Step 4. TF -> Gene
  # =====================================================
  
  edges_tf_gene <- gene_tf_binding %>%
    
    filter(TF %in% sig_promoter_tfs) %>%
    
    transmute(
      
      from = as.character(TF),
      
      to = as.character(gene_symbol),
      
      type = "TF-Gene",
      
      weight = peak_score / 1000
      
    ) %>%
    
    distinct()
  
  cat("TF-Gene 边数:",
      nrow(edges_tf_gene), "\n")
  
  # =====================================================
  # Step 5. 合并边
  # =====================================================
  
  all_edges <- bind_rows(
    
    edges_ltr_tf,
    
    edges_ppi,
    
    edges_tf_gene
  )
  
  # =====================================================
  # Step 6. 提取节点
  # =====================================================
  
  all_active_nodes <- unique(
    c(all_edges$from, all_edges$to)
  )
  
  nodes_df <- data.frame(
    name = all_active_nodes
  ) %>%
    
    mutate(
      
      node_type = case_when(
        
        name %in% unique(ltr_tf_binding$LTR)
        ~ "LTR",
        
        name %in%
          unique(gene_tf_binding$gene_symbol)
        ~ "Target_Gene",
        
        (name %in% sig_ltr_tfs) &
          (name %in% sig_promoter_tfs)
        ~ "Bridge_TF",
        
        name %in% sig_ltr_tfs
        ~ "LTR_TF",
        
        name %in% sig_promoter_tfs
        ~ "Promoter_TF",
        
        TRUE ~ "TF"
      )
    ) %>%
    
    left_join(
      
      sig_tf %>%
        dplyr::select(
          TF,
          ltr_logP = score
        ) %>%
        distinct(),
      
      by = c("name" = "TF")
    ) %>%
    
    left_join(
      
      sig_tf_promoter %>%
        dplyr::select(
          TF,
          gene_logP = score
        ) %>%
        distinct(),
      
      by = c("name" = "TF")
    )
  
  # =====================================================
  # Step 7. 最终 edge_list
  # =====================================================
  
  edge_list <- all_edges %>%
    
    filter(
      from %in% nodes_df$name &
        to %in% nodes_df$name
    )
  
  # =====================================================
  # Step 8. 输出结果
  # =====================================================
  
  cat("\n网络统计:\n")
  
  cat("总节点数:",
      nrow(nodes_df), "\n")
  
  cat("总边数:",
      nrow(edge_list), "\n")
  
  # =====================================================
  # Step 9. 保存
  # =====================================================
  
  if (!is.null(output_dir)) {
    
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    write.csv(
      nodes_df,
      file.path(output_dir,
                "node_attributes_clean.csv"),
      row.names = FALSE
    )
    
    write.csv(
      edge_list,
      file.path(output_dir,
                "edge_list_clean.csv"),
      row.names = FALSE
    )
    
    cat("网络文件已保存。\n")
  }
  
  # =====================================================
  # Step 10. 返回结果
  # =====================================================
  
  return(list(
    
    nodes_df = nodes_df,
    
    edge_list = edge_list,
    
    edges_ltr_tf = edges_ltr_tf,
    
    edges_ppi = edges_ppi,
    
    edges_tf_gene = edges_tf_gene
  ))
}