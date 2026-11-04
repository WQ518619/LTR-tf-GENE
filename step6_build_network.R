# =========================================================
# Step 6: Build Network Tables (Nodes & Edges)
# File: R/build_network_tables.R
# =========================================================
#' Build Network Tables for Plotting
#'
#' @description
#' Assemble network nodes and edges for Sankey or
#' network visualization.
#'
#' Forward mode:
#'   LTR - TF - Gene
#'
#' Reverse mode:
#'   LTR - TF1 - TF2 - Gene
#'
#' In reverse mode, TF1-TF2 represents an undirected
#' STRING protein association. The left-to-right
#' order is used only for visualization.
#'
#' Only edges participating in complete candidate
#' paths are retained.
#'
#' Node IDs are role-specific to prevent identical
#' TF names on different sides from being merged.
#'
#' @param mode Character. "forward" or "reverse".
#' @param sig_tf Data frame. Step 3 output containing
#'   LTR and TF.
#' @param step_res List. Step 5 forward output or
#'   Step 6 reverse output.
#' @param cross_ppi Data frame. Step 7 output
#'   containing from and to. Required in reverse mode.
#' @param outdir Character. Output directory.
#'
#' @return A list containing nodes and edges.
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @export
build_network_tables <- function(
    mode = c("forward", "reverse"),
    sig_tf,
    step_res,
    cross_ppi = NULL,
    outdir = NULL
) {
  mode <- match.arg(mode)
  message(
    "【Step 6】 Assembling network nodes and edges (",
    mode,
    " mode)..."
  )
  # =====================================================
  # 0. 定义空结果及保存函数
  # =====================================================
  empty_nodes <- data.frame(
    id = character(),
    name = character(),
    node_type = character(),
    stringsAsFactors = FALSE
  )
  empty_edges <- data.frame(
    from = character(),
    to = character(),
    type = character(),
    weight = numeric(),
    stringsAsFactors = FALSE
  )
  save_results <- function(nodes_df, edges_df) {
    if (!is.null(outdir)) {
      dir.create(
        outdir,
        recursive = TRUE,
        showWarnings = FALSE
      )
      utils::write.csv(
        nodes_df,
        file.path(outdir, "06_nodes.csv"),
        row.names = FALSE
      )
      utils::write.csv(
        edges_df,
        file.path(outdir, "06_edges.csv"),
        row.names = FALSE
      )
    }
    return(list(
      nodes = nodes_df,
      edges = edges_df
    ))
  }
  # =====================================================
  # 1. 根据模式构建完整候选路径
  # =====================================================
  if (mode == "forward") {
    # =================================================
    # Forward: LTR - TF - Gene
    # =================================================
    message(
      "  -> Building forward LTR-TF-Gene network..."
    )
    if (is.null(step_res) ||
        is.null(step_res$forward_edges) ||
        nrow(step_res$forward_edges) == 0) {
      warning(
        "No forward candidate connections found."
      )
      return(save_results(empty_nodes, empty_edges))
    }
    # 使用 Step 5 已筛选完成的三元组
    # 不重新引入 Step 3 的全部 LTR-TF 边
    complete_paths <- step_res$forward_edges %>%
      dplyr::transmute(
        LTR = trimws(as.character(LTR)),
        TF = trimws(as.character(TF)),
        Gene = trimws(as.character(Gene))
      ) %>%
      dplyr::filter(
        !is.na(LTR),
        !is.na(TF),
        !is.na(Gene),
        LTR != "",
        TF != "",
        Gene != ""
      ) %>%
      dplyr::distinct()
    if (nrow(complete_paths) == 0) {
      warning(
        "No valid forward paths remain."
      )
      return(save_results(empty_nodes, empty_edges))
    }
    # 1.1 LTR - TF
    edges_ltr_tf <- complete_paths %>%
      dplyr::transmute(
        from = paste0("LTR::", LTR),
        to = paste0("LTR_TF::", TF),
        type = "LTR-TF"
      ) %>%
      dplyr::distinct()
    # 1.2 TF - Gene
    edges_tf_gene <- complete_paths %>%
      dplyr::transmute(
        from = paste0("LTR_TF::", TF),
        to = paste0("GENE::", Gene),
        type = "TF-Gene"
      ) %>%
      dplyr::distinct()
    # 1.3 合并正向网络
    all_edges <- dplyr::bind_rows(
      edges_ltr_tf,
      edges_tf_gene
    )
    message(
      "  -> Valid LTR-TF-Gene paths: ",
      nrow(complete_paths)
    )
  } else {
    # =================================================
    # Reverse: LTR - TF1 - TF2 - Gene
    # =================================================
    message(
      "  -> Building reverse LTR-TF1-TF2-Gene network..."
    )
    # 反向模式检查对应数据是否为空
    if (is.null(sig_tf) ||
        nrow(sig_tf) == 0) {
      warning("No LTR-TF pairs available.")
      return(save_results(empty_nodes, empty_edges))
    }
    if (is.null(cross_ppi) ||
        nrow(cross_ppi) == 0) {
      warning("No PPI connections available.")
      return(save_results(empty_nodes, empty_edges))
    }
    if (is.null(step_res) ||
        is.null(step_res$reverse_edges) ||
        nrow(step_res$reverse_edges) == 0) {
      warning(
        "No reverse TF-Gene connections available."
      )
      return(save_results(empty_nodes, empty_edges))
    }
    # =================================================
    # 1.1 Step 3: LTR-TF1 配对
    # =================================================
    ltr_tf_pairs <- sig_tf %>%
      dplyr::transmute(
        LTR = trimws(as.character(LTR)),
        TF1 = trimws(as.character(TF))
      ) %>%
      dplyr::filter(
        !is.na(LTR),
        !is.na(TF1),
        LTR != "",
        TF1 != ""
      ) %>%
      dplyr::distinct()
    # =================================================
    # 1.2 Step 7: TF1-TF2 跨组 PPI
    # =================================================
    # Step 7:
    # from = LTR 侧 TF1
    # to   = 基因侧 TF2
    ppi_pairs <- cross_ppi %>%
      dplyr::transmute(
        TF1 = trimws(as.character(from)),
        TF2 = trimws(as.character(to))
      ) %>%
      dplyr::filter(
        !is.na(TF1),
        !is.na(TF2),
        TF1 != "",
        TF2 != "",
        TF1 != TF2
      ) %>%
      dplyr::distinct()
    # =================================================
    # 1.3 Step 6: TF2-Gene 调控关系
    # =================================================
    tf_gene_pairs <- step_res$reverse_edges %>%
      dplyr::transmute(
        TF2 = trimws(as.character(TF)),
        Gene = trimws(as.character(Gene))
      ) %>%
      dplyr::filter(
        !is.na(TF2),
        !is.na(Gene),
        TF2 != "",
        Gene != ""
      ) %>%
      dplyr::distinct()
    # =================================================
    # 1.4 构建完整 LTR-TF1-TF2-Gene 路径
    # =================================================
    # 只保留三种关系均存在的候选路径
    complete_paths <- ltr_tf_pairs %>%
      dplyr::inner_join(
        ppi_pairs,
        by = "TF1"
      ) %>%
      dplyr::inner_join(
        tf_gene_pairs,
        by = "TF2"
      ) %>%
      dplyr::distinct()
    if (nrow(complete_paths) == 0) {
      warning(
        "No complete reverse network paths found."
      )
      return(save_results(empty_nodes, empty_edges))
    }
    message(
      "  -> Valid LTR-TF1-TF2-Gene paths: ",
      nrow(complete_paths)
    )
    # =================================================
    # 1.5 提取 LTR-TF1 边
    # =================================================
    edges_ltr_tf <- complete_paths %>%
      dplyr::transmute(
        from = paste0("LTR::", LTR),
        to = paste0("LTR_TF::", TF1),
        type = "LTR-TF"
      ) %>%
      dplyr::distinct()
    # =================================================
    # 1.6 提取 TF1-TF2 PPI 边
    # =================================================
    # PPI 没有生物学调控方向
    # from/to 仅表示桑基图节点位置
    edges_ppi <- complete_paths %>%
      dplyr::transmute(
        from = paste0("LTR_TF::", TF1),
        to = paste0("GENE_TF::", TF2),
        type = "PPI"
      ) %>%
      dplyr::distinct()
    # =================================================
    # 1.7 提取 TF2-Gene 边
    # =================================================
    edges_tf_gene <- complete_paths %>%
      dplyr::transmute(
        from = paste0("GENE_TF::", TF2),
        to = paste0("GENE::", Gene),
        type = "TF-Gene"
      ) %>%
      dplyr::distinct()
    # =================================================
    # 1.8 合并反向网络
    # =================================================
    all_edges <- dplyr::bind_rows(
      edges_ltr_tf,
      edges_ppi,
      edges_tf_gene
    )
  }
  # =====================================================
  # 2. 整理边表
  # =====================================================
  all_edges <- all_edges %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      weight = 1
    )
  # =====================================================
  # 3. 构建节点表
  # =====================================================
  message("  -> Building node table...")
  all_node_ids <- unique(
    c(all_edges$from, all_edges$to)
  )
  nodes_df <- data.frame(
    id = all_node_ids,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      # 按节点 ID 前缀确定类型
      node_type = dplyr::case_when(
        startsWith(id, "LTR::") ~ "LTR",
        startsWith(id, "LTR_TF::") ~ "LTR_TF",
        startsWith(id, "GENE_TF::") ~ "Gene_TF",
        startsWith(id, "GENE::") ~ "Target_Gene",
        TRUE ~ "Other"
      ),
      # 提取原始显示名称
      name = sub(
        "^(LTR|LTR_TF|GENE_TF|GENE)::",
        "",
        id
      )
    ) %>%
    dplyr::select(
      id,
      name,
      node_type
    ) %>%
    dplyr::distinct()
  # =====================================================
  # 4. 控制台统计
  # =====================================================
  message(
    "  -> Total nodes: ",
    nrow(nodes_df)
  )
  message(
    "  -> Total edges: ",
    nrow(all_edges)
  )
  node_counts <- table(
    nodes_df$node_type
  )
  for (nt in names(node_counts)) {
    message(
      "     - ",
      nt,
      ": ",
      node_counts[[nt]]
    )
  }
  edge_counts <- table(
    all_edges$type
  )
  for (et in names(edge_counts)) {
    message(
      "     - ",
      et,
      " edges: ",
      edge_counts[[et]]
    )
  }
  # =====================================================
  # 5. 保存并返回结果
  # =====================================================
  return(
    save_results(
      nodes_df,
      all_edges
    )
  )
}