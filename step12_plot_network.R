# =========================================================
# Step 12: 网络可视化（无重叠优化版）
# File:
# R/step12_plot_network.R
# =========================================================

# ---------------------------------------------------------
# Function:
# plot_ltr_tf_gene_network
#
# Purpose:
# 对 Step11 生成的节点表与边表进行网络可视化。
#
# 主要功能：
# 1. 自动进行网络分层
# 2. 横向拉宽网络结构
# 3. 避免文字与连线重叠
# 4. 输出高清 PNG 图
# 5. 输出网络统计信息
#
# Input:
#
# nodes_df
#   Step11 生成的节点属性表
#
# edges_df
#   Step11 生成的边表
#
# output_dir
#   输出目录
#
# base_name
#   图片名称
#
# width,height
#   图片尺寸
#
# dpi
#   分辨率
#
# Return:
# list
#   p            ggplot对象
#   nodes_pos    节点坐标
#   edge_data    边数据
#
# ---------------------------------------------------------

# =========================================================
# Step 12: 核心网络可视化模块 (LTR-TF-Gene 调控网络)
# =========================================================
plot_ltr_tf_gene_network <- function(
    nodes_df,
    edges_df,
    pipeline_mode,     # 自动识别: "forward"(3列) 或 "reverse"(4列)
    output_dir,
    base_name,
    width,
    height,
    dpi,
    exclude_epitope,   # 是否排除 Epitope 标签
    strict_ppi_mode,   # 仅逆向模式生效：严格切断 PPI 捷径
    max_ltrs,          # LTR 显示上限
    max_tfs,           # TF 显示上限
    max_genes          # 靶基因显示上限
){
  
  message(sprintf("\n[Step 12] 检测到 %s 模式，正在按配置组装网络布局...", pipeline_mode))
  
  # ---------------------------------------------------------
  # 1. 基础清洗与节点标记
  # ---------------------------------------------------------
  epi_pattern <- "(?i)acetyl|H[234]K[0-9]+"
  nodes_df <- nodes_df %>% filter(!grepl(epi_pattern, name))
  edges_df <- edges_df %>% filter(from %in% nodes_df$name & to %in% nodes_df$name)
  
  if (exclude_epitope) {
    nodes_df <- nodes_df %>% filter(!grepl("Epitope", name, ignore.case = TRUE))
    edges_df <- edges_df %>% filter(from %in% nodes_df$name & to %in% nodes_df$name)
  }
  
  # 给下游 Target 加上后缀区分身份
  edges_df <- edges_df %>%
    mutate(to = ifelse(type == "TF-Gene", paste0(to, " (Target)"), to))
  
  target_nodes_df <- data.frame(
    name = unique(edges_df$to[edges_df$type == "TF-Gene"]),
    node_type = "Target_Gene"
  )
  
  nodes_df <- nodes_df %>%
    mutate(node_type = ifelse(node_type == "Target_Gene", "TF", node_type)) %>%
    bind_rows(target_nodes_df) %>%
    distinct(name, .keep_all = TRUE)
  
  # ---------------------------------------------------------
  # 2. 网络颜值保护 (动态数量截断)
  # ---------------------------------------------------------
  ltr_nodes <- nodes_df %>% filter(node_type == "LTR") %>% pull(name)
  gene_nodes <- nodes_df %>% filter(node_type == "Target_Gene") %>% pull(name)
  tf_nodes_all <- nodes_df %>% filter(!node_type %in% c("LTR", "Target_Gene")) %>% pull(name)
  
  ltr_edges <- c(edges_df$from, edges_df$to)[c(edges_df$from, edges_df$to) %in% ltr_nodes]
  gene_edges <- c(edges_df$from, edges_df$to)[c(edges_df$from, edges_df$to) %in% gene_nodes]
  tf_edges <- c(edges_df$from, edges_df$to)[c(edges_df$from, edges_df$to) %in% tf_nodes_all]
  
  keep_ltrs <- names(sort(table(ltr_edges), decreasing = TRUE))[1:min(max_ltrs, length(unique(ltr_edges)))]
  keep_genes <- names(sort(table(gene_edges), decreasing = TRUE))[1:min(max_genes, length(unique(gene_edges)))]
  keep_tfs <- names(sort(table(tf_edges), decreasing = TRUE))[1:min(max_tfs, length(unique(tf_edges)))]
  
  keep_nodes <- c(keep_ltrs, keep_genes, keep_tfs)
  
  nodes_df <- nodes_df %>% filter(name %in% keep_nodes)
  edges_df <- edges_df %>% filter(from %in% nodes_df$name & to %in% nodes_df$name)
  
  # ---------------------------------------------------------
  # 3. 逆向 PPI 处理模块 (仅在 reverse 模式下激活)
  # ---------------------------------------------------------
  if (pipeline_mode == "reverse" && strict_ppi_mode) {
    ltr_tfs_current <- edges_df %>% filter(type == "LTR-TF") %>% pull(to)
    prom_tfs_current <- edges_df %>% filter(type == "TF-Gene") %>% pull(from)
    bridge_tfs <- intersect(ltr_tfs_current, prom_tfs_current)
    
    edges_df <- edges_df %>% filter(!(type == "TF-Gene" & from %in% bridge_tfs))
    
    valid_tf1s <- edges_df %>% filter(type == "LTR-TF") %>% pull(to) %>% unique()
    valid_tf2s <- edges_df %>% filter(type == "TF-Gene") %>% pull(from) %>% unique()
    
    edges_df <- edges_df %>%
      mutate(
        temp_from = ifelse(type == "PPI" & from %in% valid_tf2s & to %in% valid_tf1s, to, from),
        temp_to = ifelse(type == "PPI" & from %in% valid_tf2s & to %in% valid_tf1s, from, to)
      ) %>%
      mutate(from = temp_from, to = temp_to) %>%
      dplyr::select(-temp_from, -temp_to) %>%
      filter(type != "PPI" | (from %in% valid_tf1s & to %in% valid_tf2s))
  }
  
  # ---------------------------------------------------------
  # 4. 智能分层与坐标生成核心逻辑
  # ---------------------------------------------------------
  active_nodes <- unique(c(edges_df$from, edges_df$to))
  nodes_df <- nodes_df %>% filter(name %in% active_nodes)
  
  if (pipeline_mode == "forward") {
    # 🌟 正向模式：强行收束为 3 列瀑布流 (LTR -> TF -> Gene)
    nodes_df <- nodes_df %>%
      mutate(
        group = case_when(
          node_type == "LTR" ~ 1,
          node_type == "Target_Gene" ~ 3,
          TRUE ~ 2
        )
      )
    nodes_pos <- nodes_df %>%
      group_by(group) %>% arrange(name) %>% 
      mutate(
        x_base = group * 8,
        x = case_when(group == 1 ~ x_base - 2, group == 3 ~ x_base + 2, TRUE ~ x_base),
        y = (row_number() - (n() + 1) / 2) * 3,
        text_hjust = case_when(group == 1 ~ 1, group == 3 ~ 0, TRUE ~ 0.5)
      ) %>% ungroup()
      
  } else {
    # 🌟 逆向模式：保留 4 列分布式瀑布流 (LTR -> TF1 -> PPI -> TF2 -> Gene)
    nodes_df <- nodes_df %>%
      mutate(
        group = case_when(
          name %in% edges_df$from[edges_df$type == "LTR-TF"] ~ 1,  
          name %in% edges_df$to[edges_df$type == "LTR-TF"] ~ 2,    
          name %in% edges_df$from[edges_df$type == "TF-Gene"] ~ 3, 
          name %in% edges_df$to[edges_df$type == "TF-Gene"] ~ 4,   
          TRUE ~ 2.5 
        )
      )
    nodes_pos <- nodes_df %>%
      group_by(group) %>% arrange(desc(node_type), name) %>% 
      mutate(
        x_base = group * 8,
        x = case_when(group == 1 ~ x_base - 3, group == 4 ~ x_base + 3, TRUE ~ x_base),
        y = (row_number() - (n() + 1) / 2) * 3,
        text_hjust = case_when(group == 1 ~ 1, group == 4 ~ 0, TRUE ~ 0.5)
      ) %>% ungroup()
  }
  
  # ---------------------------------------------------------
  # 5. 颜色映射与出图
  # ---------------------------------------------------------
  edge_data <- edges_df %>%
    inner_join(nodes_pos, by = c("from" = "name")) %>%
    inner_join(nodes_pos, by = c("to" = "name"), suffix = c(".u", ".v")) %>%
    mutate(
      x_start = x.u + (x.v - x.u) * 0.1,
      x_end   = x.v - (x.v - x.u) * 0.1,
      edge_id = row_number(),
      edge_color = case_when(
        type == "PPI" ~ "#d95f02",      
        type == "LTR-TF" ~ "#1b9e77",   
        type == "TF-Gene" ~ "#7570b3",  
        TRUE ~ "grey70"
      )
    )
  
  p <- ggplot() +
    geom_diagonal(
      data = edge_data, aes(x = x_start, y = y.u, xend = x_end, yend = y.v, group = edge_id, color = edge_color),
      strength = 0.6, linewidth = 0.7, alpha = 0.55
    ) +
    geom_text(
      data = nodes_pos, aes(x = x, y = y, label = name, hjust = text_hjust),
      color = "black", size = 3.6, fontface = "bold", family = "serif", position = position_nudge(y = 0.15)
    ) +
    scale_color_identity() + theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(30, 80, 30, 80)) +
    coord_cartesian(clip = "off")
  
  if(!is.null(output_dir)){
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    base_name <- tools::file_path_sans_ext(base_name)
    ggsave(file.path(output_dir, paste0(base_name, ".png")), p, width = width, height = height, dpi = dpi)
    ggsave(file.path(output_dir, paste0(base_name, ".pdf")), p, width = width, height = height)
    message(sprintf("[Step 12] 智能布局网络图绘制完毕 (%s模式)，已保存至: %s", pipeline_mode, output_dir))
  }
  
  return(list(p = p, nodes_pos = nodes_pos, edge_data = edge_data))
}