# =========================================================
# Step 7: Publication-ready LTR-associated network plot
# File: R/plot_ltr_tf_gene_network.R
# Forward: preserve LOLA LTR--TF and DoRothEA TF--Gene evidence.
# Reverse: reconstruct supported paths from Step 6 edges.
# Visualization filtering does not modify Step 6 full network.
# =========================================================
#' Plot a simplified LTR-associated regulatory network
#'
#' @description
#' Forward mode uses Step 4 LTR--TF--Gene associations; LTR--TF
#' pairs are retained as supplied. TF and Gene nodes are shared globally.
#' Gene need not be proximal to the displayed LTR family.
#' Reverse mode uses role-specific edges from Step 6; STRING
#' protein associations are undirected (visual order only).
#'
#' @param nodes Full Step 6 network nodes (id, name, node_type).
#' @param edges Full Step 6 network edges (from, to, type, weight).
#' @param mode Either "forward" or "reverse".
#' @param outdir Output directory or NULL.
#' @param width,height Figure dimensions in inches.
#' @param top_n_ltr Number of LTR families to show.
#' @param top_n_tf Number of distinct TFs globally (both modes).
#' @param top_n_gene_tf Number of TF2 nodes (reverse only).
#' @param top_n_genes Forward: genes per TF; Reverse: genes per TF2.
#' @param label_size ggplot text size.
#' @param forward_paths Step 4 forward_edges with LTR, TF, Gene.
#' @return A ggplot with selected nodes, edges and paths as attributes.
#' @importFrom magrittr %>%
#' @export
plot_ltr_tf_gene_network <- function(
    nodes, edges, mode = c("forward", "reverse"), outdir = NULL,
    width = 14, height = 10,
    top_n_ltr = 5, top_n_tf = 8, top_n_gene_tf = 8,
    top_n_genes = 3, label_size = 3.1, forward_paths = NULL
) {
  mode <- match.arg(mode)
  limits <- c(top_n_ltr, top_n_tf, top_n_gene_tf, top_n_genes)
  if (anyNA(limits) || any(!is.finite(limits)) ||
      any(limits < 1) || any(limits != floor(limits))) {
    stop("All top_n_* values must be positive integers.")
  }
  if (!is.numeric(label_size) || length(label_size) != 1L ||
      !is.finite(label_size) || label_size <= 0) {
    stop("'label_size' must be a positive number.")
  }

  # Select only supported input associations, not arbitrary unobserved edges.
  keep_top <- function(tbl, key, partner, n) {
    ranking <- tbl %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
      dplyr::summarise(degree = dplyr::n_distinct(.data[[partner]]),
                       .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(degree), .data[[key]]) %>%
      dplyr::slice_head(n = n)
    dplyr::semi_join(tbl, ranking, by = key)
  }

  if (mode == "forward") {
    if (is.null(forward_paths) ||
        !all(c("LTR", "TF", "Gene") %in% names(forward_paths))) {
      stop("Forward mode requires Step 4 forward_edges with LTR, TF, Gene.")
    }
    paths <- forward_paths %>%
      dplyr::transmute(
        LTR = trimws(as.character(LTR)),
        TF = trimws(as.character(TF)),
        Gene = trimws(as.character(Gene))
      ) %>%
      dplyr::filter(!is.na(LTR), !is.na(TF), !is.na(Gene),
                    nzchar(LTR), nzchar(TF), nzchar(Gene)) %>%
      dplyr::distinct()
    if (!nrow(paths)) {
      warning("No valid forward triples to plot.")
      return(NULL)
    }

    # Select LTR families by number of distinct associated TFs.
    paths <- keep_top(paths, "LTR", "TF", top_n_ltr)

    # Select up to top_n_tf distinct TFs in the ENTIRE Forward plot.
    # Prefer TFs associated with more LTR families, then more genes.
    tf_keep <- paths %>%
      dplyr::group_by(TF) %>%
      dplyr::summarise(
        n_ltr = dplyr::n_distinct(LTR),
        n_genes = dplyr::n_distinct(Gene),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(n_ltr), dplyr::desc(n_genes), TF) %>%
      dplyr::slice_head(n = top_n_tf) %>%
      dplyr::pull(TF)
    paths <- paths %>% dplyr::filter(TF %in% tf_keep)

    # Select up to top_n_genes distinct genes PER TF (not per LTR).
    # Prefer genes with support from more selected TFs.
    gene_rank <- paths %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(tf_support = dplyr::n_distinct(TF),
                       .groups = "drop")
    gene_keep <- paths %>%
      dplyr::distinct(TF, Gene) %>%
      dplyr::left_join(gene_rank, by = "Gene") %>%
      dplyr::group_by(TF) %>%
      dplyr::arrange(dplyr::desc(tf_support), Gene, .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n_genes) %>%
      dplyr::ungroup() %>%
      dplyr::select(TF, Gene)
    paths <- dplyr::semi_join(paths, gene_keep, by = c("TF", "Gene"))

    # ONE plot node per LTR, TF, and Gene. The LTR--TF pairs are
    # faithful to LOLA; TF--Gene pairs are faithful to DoRothEA.
    # Shared nodes do NOT establish LTR-specific cis gene regulation.
    paths <- paths %>%
      dplyr::mutate(
        LTR_ID = paste0("LTR::", LTR),
        TF1_ID = paste0("LTR_TF::", TF),
        Gene_ID = paste0("GENE::", Gene)
      )

    selected_edges <- dplyr::bind_rows(
      paths %>% dplyr::transmute(from = LTR_ID, to = TF1_ID,
                                 type = "LTR-TF", weight = 1),
      paths %>% dplyr::transmute(from = TF1_ID, to = Gene_ID,
                                 type = "TF-Gene", weight = 1)
    ) %>% dplyr::distinct()

    selected_nodes <- dplyr::bind_rows(
      paths %>% dplyr::distinct(id = LTR_ID, name = LTR) %>%
        dplyr::mutate(node_type = "LTR"),
      paths %>% dplyr::distinct(id = TF1_ID, name = TF) %>%
        dplyr::mutate(node_type = "LTR_TF"),
      paths %>% dplyr::distinct(id = Gene_ID, name = Gene) %>%
        dplyr::mutate(node_type = "Target_Gene")
    ) %>% dplyr::distinct(id, .keep_all = TRUE)

    # The global TF/Gene nodes do not belong to one LTR-specific group.
    selected_nodes$group_rank <- 1L

  } else {
    if (is.null(nodes) || is.null(edges) || !nrow(nodes) || !nrow(edges)) {
      warning("Empty Step 6 network.")
      return(NULL)
    }
    if (!all(c("id", "name", "node_type") %in% names(nodes)) ||
        !all(c("from", "to", "type") %in% names(edges))) {
      stop("Step 6 nodes/edges are missing required columns.")
    }
    ltr_tf <- edges %>%
      dplyr::filter(type == "LTR-TF") %>%
      dplyr::transmute(LTR_ID = as.character(from),
                        TF1_ID = as.character(to)) %>% dplyr::distinct()
    ppi <- edges %>%
      dplyr::filter(type == "PPI") %>%
      dplyr::transmute(TF1_ID = as.character(from),
                        TF2_ID = as.character(to)) %>% dplyr::distinct()
    tf_gene <- edges %>%
      dplyr::filter(type == "TF-Gene") %>%
      dplyr::transmute(TF2_ID = as.character(from),
                        Gene_ID = as.character(to)) %>% dplyr::distinct()
    paths <- ltr_tf %>% dplyr::inner_join(ppi, by = "TF1_ID") %>%
      dplyr::inner_join(tf_gene, by = "TF2_ID") %>% dplyr::distinct()
    if (!nrow(paths)) {
      warning("No complete reverse paths to plot.")
      return(NULL)
    }
    paths <- keep_top(paths, "LTR_ID", "TF1_ID", top_n_ltr)
    paths <- keep_top(paths, "TF1_ID", "LTR_ID", top_n_tf)
    paths <- keep_top(paths, "TF2_ID", "TF1_ID", top_n_gene_tf)
    gene_keep <- paths %>%
      dplyr::distinct(TF2_ID, Gene_ID) %>%
      dplyr::group_by(TF2_ID) %>%
      dplyr::arrange(Gene_ID, .by_group = TRUE) %>%
      dplyr::slice_head(n = top_n_genes) %>% dplyr::ungroup()
    paths <- dplyr::semi_join(paths, gene_keep,
                              by = c("TF2_ID", "Gene_ID"))
    selected_edges <- dplyr::bind_rows(
      paths %>% dplyr::transmute(from = LTR_ID, to = TF1_ID, type = "LTR-TF"),
      paths %>% dplyr::transmute(from = TF1_ID, to = TF2_ID, type = "PPI"),
      paths %>% dplyr::transmute(from = TF2_ID, to = Gene_ID, type = "TF-Gene")
    ) %>% dplyr::distinct()
    if ("weight" %in% names(edges)) {
      selected_edges <- selected_edges %>%
        dplyr::left_join(
          edges %>% dplyr::select(from, to, type, weight) %>% dplyr::distinct(),
          by = c("from", "to", "type")
        )
    } else {
      selected_edges$weight <- 1
    }
    ids <- unique(c(selected_edges$from, selected_edges$to))
    selected_nodes <- nodes %>%
      dplyr::filter(id %in% ids) %>% dplyr::distinct(id, .keep_all = TRUE)
    selected_nodes$group_rank <- 1L
  }

  # Arrange nodes by layers and center every layer around y = 0.
  max_layer <- if (mode == "forward") 3L else 4L
  x_positions <- seq(0, by = 5, length.out = max_layer)
  selected_nodes <- selected_nodes %>%
    dplyr::mutate(
      layer = dplyr::case_when(
        node_type == "LTR" ~ 1L,
        node_type == "LTR_TF" ~ 2L,
        node_type == "Gene_TF" ~ 3L,
        node_type == "Target_Gene" & mode == "forward" ~ 3L,
        node_type == "Target_Gene" & mode == "reverse" ~ 4L,
        TRUE ~ NA_integer_
      )
    )
  if (anyNA(selected_nodes$layer)) {
    stop("Unrecognized node_type in selected network.")
  }
  plot_nodes <- selected_nodes %>%
    dplyr::group_by(layer) %>%
    dplyr::arrange(group_rank, name, id, .by_group = TRUE) %>%
    dplyr::mutate(
      y = ((dplyr::n() + 1) / 2 - dplyr::row_number()) * 1.15
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(x = x_positions[layer])

  # Lines originate at node anchors; labels are offset away from anchors.
  plot_edges <- selected_edges %>%
    dplyr::inner_join(
      plot_nodes %>% dplyr::transmute(from = id, x1 = x, y1 = y, from_layer = layer),
      by = "from"
   ) %>%
    dplyr::inner_join(
      plot_nodes %>% dplyr::transmute(to = id, x2 = x, y2 = y, to_layer = layer),
      by = "to"
   ) %>%
    dplyr::mutate(
      x_start = ifelse(from_layer == 1L, x1 + 0.15, x1 + 0.50),
      x_end   = ifelse(to_layer == max_layer, x2 - 0.15, x2 - 0.50),
      edge_id = dplyr::row_number()
   )

  plot_nodes <- plot_nodes %>%
    dplyr::mutate(
      label_x = dplyr::case_when(
        layer == 1L ~ x - 0.12,
        layer == max_layer ~ x + 0.12,
        TRUE ~ x
      ),
      label_y = y,
      label_hjust = dplyr::case_when(
        layer == 1L ~ 1,
        layer == max_layer ~ 0,
        TRUE ~ 0.5
      )
    )

  p <- ggplot2::ggplot() +
    ggforce::geom_diagonal(
      data = plot_edges,
      mapping = ggplot2::aes(
        x = x_start, y = y1, xend = x_end, yend = y2,
        color = type, group = edge_id
      ),
      strength = 0.5, alpha = 0.48, linewidth = 0.65
    ) +
    ggplot2::geom_text(
      data = plot_nodes,
      mapping = ggplot2::aes(
        x = label_x, y = label_y,
        label = name, hjust = label_hjust
      ),
      size = label_size, color = "#202020", fontface = "plain",
      check_overlap = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("LTR-TF" = "#229B85", "PPI" = "#D58B48",
                 "TF-Gene" = "#797BC7"), name = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = NULL, expand = ggplot2::expansion(add = c(1.8, 2.4))
    ) +
    ggplot2::scale_y_continuous(
      breaks = NULL, expand = ggplot2::expansion(add = 1.2)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(20, 32, 20, 32)
    )

  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(paths,
      file.path(outdir, paste0("07_plot_paths_", mode, ".csv")),
      row.names = FALSE)
    utils::write.csv(selected_nodes,
      file.path(outdir, paste0("07_plot_nodes_", mode, ".csv")),
      row.names = FALSE)
    utils::write.csv(selected_edges,
      file.path(outdir, paste0("07_plot_edges_", mode, ".csv")),
      row.names = FALSE)
    ggplot2::ggsave(
      file.path(outdir, paste0("07_network_", mode, ".pdf")),
      plot = p, width = width, height = height,
      device = grDevices::cairo_pdf, limitsize = FALSE
    )
    ggplot2::ggsave(
      file.path(outdir, paste0("07_network_", mode, ".png")),
      plot = p, width = width, height = height,
      dpi = 300, limitsize = FALSE
    )
  }

  attr(p, "plot_nodes") <- selected_nodes
  attr(p, "plot_edges") <- selected_edges
  attr(p, "plot_paths") <- paths
  message("  -> Plotted ", nrow(selected_nodes), " nodes and ",
          nrow(selected_edges), " edges (", mode,"). Full Step 6 unchanged.")
  p
}
