# volcano plot for "conditon"_vs_Mock comparison
create_enhanced_volcano <- function(deg_results, comparison_name) {

  # Define theme for all plots
  large_text_theme <- theme(
    text = element_text(size = 30),           # Increase general text size by 100% (from ~8 to 16)
    axis.text = element_text(size = 24),      # Axis text size
    axis.text.y = element_text(size = 28),
    axis.text.x = element_text(size = 24),
    axis.title = element_text(size = 32),     # Axis title size
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),  # Plot title
    #legend.text = element_text(size = 24),    # Legend text
    #legend.title = element_text(size = 30),    # Legend title
    legend.position = "none"
  )

  # Prepare the data
  volcano_data <- deg_results %>%
    as.data.frame() %>%
    mutate(
      log10FDR = -log10(FDR),
      significance = case_when(
        FDR < 0.05 & logFC > 1.5 ~ "Upregulated",
        FDR < 0.05 & logFC < -1.5 ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      gene_name = rownames(deg_results),
    )

  # Create a combined significance column that includes the pathway information
  volcano_data <- volcano_data %>%
    mutate(
      display_group = case_when(
        significance == "Upregulated" ~ "Upregulated",
        significance == "Downregulated" ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )

  up_count <- sum(volcano_data$significance == "Upregulated")
  down_count <- sum(volcano_data$significance == "Downregulated")

  # Create the base volcano plot
  p <- ggplot(volcano_data, aes(x = logFC, y = log10FDR, color = display_group)) +
    geom_point(aes(alpha = ifelse(display_group == "Not Significant", 0.5, 1)), size = 3) +
    scale_alpha_identity() +
    scale_color_manual(
      values = c(
        "Upregulated" = "red",
        "Downregulated" = "blue",
        "Not Significant" = "grey70"
      ),
      breaks = c("Upregulated", "Downregulated", "Not Significant")
    ) +

    # Add vertical and horizontal lines
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +

    annotate("text",
             x = min(volcano_data$logFC, na.rm = TRUE) + 1,
             y = max(volcano_data$log10FDR, na.rm = TRUE) - 0.5,
             label = paste0("Down: ", down_count, "   Up: ", up_count ),
             size = 12,
             hjust = 0,
             vjust = 18,
             color = "black",
             fontface = "bold",
             lineheight = 1) +

    # Theme and labels
    labs(
      x = "log2(Fold Change)",
      y = "-log10(FDR)",
      color = "Significance"
    ) + large_text_theme

    legend_data <- data.frame(
      Group = factor(c("Upregulated", "Downregulated", "Not Significant"),
                     levels = c("Upregulated", "Downregulated", "Not Significant")),
      Color = c("red", "blue", "grey70")
    )

    # Legend 
    create_legend_grob <- function() {
      vp <- viewport(x = 0.5, y = 0.5, width = 0.9, height = 0.9)
      
      legend_grob <- gTree(vp = vp, children = gList(
        rectGrob(gp = gpar(fill = "white", col = "white")),
        
        pointsGrob(x = seq(0, 0.7, length.out = 3), y = rep(0.9, 3),
                   pch = 21, size = unit(1.5, "char"),
                   gp = gpar(fill = legend_data$Color, col = "black")),
        
        # Legend text
        textGrob(x = seq(0.025, 0.725, length.out = 3), y = rep(0.9, 3),
                 label = legend_data$Group,
                 just = "left", gp = gpar(fontsize = 22))
      ))
      
      return(legend_grob)
    }
    
    # Create the legend grob
    legend_grob <- create_legend_grob()

  # Combine volcano plot and legend
  combined_plot <- grid.arrange(
    p, legend_grob,
    nrow = 2,
    heights = c(4, 1)
  )

  # Save the volcano plot with white background
  volcano_filename <- paste0("Volcano_Plot_", comparison_name, ".png")
  ggsave(volcano_filename, plot = combined_plot, width = 10, height = 10, dpi = 300, bg = "white")
  cat(paste("Volcano plot saved:", volcano_filename, "\n"))

  return(list(volcano_plot = combined_plot))
}
