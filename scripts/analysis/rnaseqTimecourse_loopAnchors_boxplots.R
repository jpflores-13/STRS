## Integrated Timecourse Analysis: Visualization with Statistical Testing
## Combines clean boxplot visualization with rigorous consecutive timepoint analysis

library(InteractionSet)
library(mariner)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(RColorBrewer)
library(cowplot)
library(rstatix)
library(stringr)

## Load and process data
dge <- readRDS("data/processed/rna/timecourse/output/deseqObjs/LRTtimecourse.rds")
dge_gr <- results(dge, format = "GRanges") 
mcols(dge_gr) <- rowData(dge)
dge_gr_prom <- promoters(dge_gr) |> 
  keepStandardChromosomes(pruning.mode = "coarse")
seqlevelsStyle(dge_gr_prom) <- "UCSC"

## Process loop data
loops <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions() |> 
  as.data.frame() |> 
  as_ginteractions()

## Separate loops into gained, lost, static categories
gainedLoops <- loops |> 
  subset(padj < 0.1 & log2FoldChange > 0)

lostLoops <- loops |> 
  subset(padj < 0.1 & log2FoldChange < 0)

staticLoops <- loops |> 
  subset(padj > 0.1)

## Find overlaps and create combined dataset
gained_hits <- subsetByOverlaps(dge_gr_prom, gainedLoops)
lost_hits <- subsetByOverlaps(dge_gr_prom, lostLoops)
static_hits <- subsetByOverlaps(dge_gr_prom, staticLoops)

gained_df <- as.data.frame(gained_hits)
gained_df$type <- "gained" 

lost_df <- as.data.frame(lost_hits)
lost_df$type <- "lost" 

static_df <- as.data.frame(static_hits)
static_df$type <- "static" 

combined <- bind_rows(gained_df, lost_df, static_df) |> 
  select(Time_1h_vs_0h,
         Time_3h_vs_0h,
         Time_6h_vs_0h,
         Time_9h_vs_0h,
         Time_12h_vs_0h,
         Time_24h_vs_0h,
         type,
         symbol,
         gene_id) |> 
  pivot_longer(cols = starts_with("Time"),
               values_to = "log2FoldChange", 
               names_to = "timepoint") |> 
  mutate(timepoint = case_when(
    timepoint == "Time_1h_vs_0h" ~ "1",
    timepoint == "Time_3h_vs_0h" ~ "3",
    timepoint == "Time_6h_vs_0h" ~ "6",
    timepoint == "Time_9h_vs_0h" ~ "9",
    timepoint == "Time_12h_vs_0h" ~ "12",
    timepoint == "Time_24h_vs_0h" ~ "24"))

# Set proper factor ordering for consistent analysis and visualization
combined$timepoint <- factor(combined$timepoint, 
                             levels = c("1", "3", "6", "9", "12", "24"))
combined$type <- factor(combined$type, 
                        levels = c("static", "gained", "lost"))

# Define consistent color palette for all visualizations
loop_colors <- c(
  "static" = "#999999",   # Neutral gray for static loops
  "gained" = "#F8766D",   # Red for gained loops  
  "lost" = "#619CFF"      # Blue for lost loops
)

## CONSECUTIVE TIMEPOINT ANALYSIS FUNCTIONS

# Function to prepare paired consecutive timepoint data for statistical testing
# This transformation creates expression change values between adjacent timepoints
# which we then test against zero to identify periods of systematic regulation
prepare_consecutive_data <- function(data) {
  
  # Transform data to capture changes between consecutive timepoints
  # Each row represents one gene's expression change between two adjacent timepoints
  consecutive_pairs <- data |>
    arrange(gene_id, type, timepoint) |>
    group_by(gene_id, type) |>
    mutate(
      # Calculate expression change from the previous timepoint
      # This is the core measurement for our statistical tests
      previous_expression = lag(log2FoldChange),
      previous_timepoint = lag(timepoint),
      expression_change = log2FoldChange - previous_expression,
      
      # Create descriptive labels for each transition period
      transition = paste0(previous_timepoint, "→", timepoint)
    ) |>
    filter(!is.na(expression_change)) |>  # Remove first timepoint (no previous comparison)
    ungroup()
  
  return(consecutive_pairs)
}

# Function to conduct statistical testing using Wilcoxon signed-rank tests
# Tests whether median expression changes between consecutive timepoints 
# are significantly different from zero (no change)
conduct_consecutive_tests <- function(consecutive_data) {
  
  # Storage for results from all three loop types
  all_test_results <- list()
  
  # Analyze each loop type separately to respect biological grouping
  loop_types <- c("static", "gained", "lost")
  
  for(loop_type in loop_types) {
    
    # Extract data for this specific loop type
    type_data <- consecutive_data |> filter(type == loop_type)
    
    # For each transition period, test whether changes are systematically different from zero
    # Wilcoxon signed-rank test is ideal because it tests whether median of difference
    # values differs from a specified value (zero) without requiring normal distribution
    transition_results <- type_data |>
      group_by(transition) |>
      rstatix::wilcox_test(expression_change ~ 1, mu = 0) |>  # Test against zero change
      mutate(
        # Calculate median change for biological interpretation
        median_change = type_data |> group_by(transition) |> 
          summarise(med = median(expression_change, na.rm = TRUE)) |> pull(med),
        
        # Count genes contributing to each test for context
        n_genes = type_data |> group_by(transition) |> 
          summarise(n = dplyr::n()) |> pull(n)
      ) |>
      ungroup()
    
    # Apply Benjamini-Hochberg correction within each loop type
    # Controls false discovery rate among the five transitions tested per loop type
    transition_results$p_adjusted <- p.adjust(transition_results$p, method = "BH")
    
    # Convert adjusted p-values to significance indicators using standard conventions
    transition_results$significance_stars <- case_when(
      transition_results$p_adjusted < 0.001 ~ "***",
      transition_results$p_adjusted < 0.01 ~ "**",
      transition_results$p_adjusted < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    
    # Store results for visualization
    all_test_results[[loop_type]] <- transition_results
  }
  
  return(all_test_results)
}

# Main plotting function that creates panels with statistical significance indicators
# Integrates clean boxplot visualization with comprehensive statistical testing results
create_integrated_plot <- function(data, test_results) {
  
  # Function to create individual panels with statistical annotations
  create_panel_with_stats <- function(loop_type, show_y_axis = TRUE, show_y_title = TRUE, show_x_title = TRUE) {
    
    # Filter data for this specific loop type
    plot_data <- data |> filter(type == loop_type)
    
    # Get statistical test results for this loop type
    loop_test_results <- test_results[[loop_type]]
    
    # Build base plot with clean boxplot styling
    p <- ggplot(plot_data, aes(x = timepoint, y = log2FoldChange)) +
      
      # Boxplots with white fill and colored borders for clean appearance
      geom_boxplot(
        fill = "white",
        color = loop_colors[loop_type],
        outlier.shape = NA,  # Remove outliers for cleaner visualization
        width = 0.7,
        linewidth = 0.5
      ) +
      
      # Reference line at zero - critical for interpreting fold changes
      geom_hline(yintercept = 0, color = "gray30", linetype = "dashed", linewidth = 0.4) +
      
      # Connected median points to visualize temporal expression patterns
      stat_summary(
        fun = median,
        geom = "line",
        aes(group = 1),
        color = loop_colors[loop_type],
        linewidth = 0.8
      ) +
      
      # Highlight individual median points for precise reading
      stat_summary(
        fun = median,
        geom = "point",
        color = loop_colors[loop_type],
        size = 2,
        shape = 18
      ) +
      
      # Loop type label in upper-left corner following scientific conventions
      annotate("text", 
               x = -Inf, y = Inf,
               label = stringr::str_to_title(loop_type),
               hjust = -0.1, vjust = 1.5,
               size = 4, fontface = "bold",
               color = loop_colors[loop_type])
    
    # Add statistical significance indicators for all tested transitions
    # This provides complete transparency about which comparisons were performed
    
    # Define x-axis positions for precise bar placement
    timepoint_positions <- c("1" = 1, "3" = 2, "6" = 3, "9" = 4, "12" = 5, "24" = 6)
    
    # Position significance bars in upper portion of plot area
    y_position <- 0.85
    bar_spacing <- 0.08
    
    # Process all transitions to show complete statistical analysis
    for(i in 1:nrow(loop_test_results)) {
      
      transition <- loop_test_results$transition[i]
      stars <- loop_test_results$significance_stars[i]
      
      # Parse transition string to identify compared timepoints
      timepoints <- str_split(transition, "→")[[1]]
      start_timepoint <- timepoints[1]
      end_timepoint <- timepoints[2]
      
      # Convert to x-axis positions
      x_start <- timepoint_positions[start_timepoint]
      x_end <- timepoint_positions[end_timepoint]
      
      # Create shorter bars for refined appearance (reduce length by 30% total)
      bar_length_reduction <- 0.15
      x_start_adjusted <- x_start + bar_length_reduction
      x_end_adjusted <- x_end - bar_length_reduction
      
      # Calculate vertical position for this significance indicator
      current_y <- y_position - (i - 1) * bar_spacing
      
      # Use consistent loop colors with transparency to indicate significance status
      # This maintains biological category identity while showing statistical results
      bar_alpha <- if(stars != "ns") 1.0 else 0.4  # Transparent for non-significant
      text_alpha <- if(stars != "ns") 1.0 else 0.6
      
      # Add horizontal significance bar
      p <- p + 
        annotate("segment",
                 x = x_start_adjusted, xend = x_end_adjusted,
                 y = current_y, yend = current_y,
                 color = alpha(loop_colors[loop_type], bar_alpha),
                 linewidth = 0.6) +
        
        # Add significance indicator text above bar center
        annotate("text",
                 x = (x_start + x_end) / 2,
                 y = current_y + 0.04,
                 label = stars,
                 size = if(stars != "ns") 4 else 3.5,
                 fontface = if(stars != "ns") "bold" else "plain",
                 color = alpha(loop_colors[loop_type], text_alpha))
    }
    
    # Apply symmetric coordinate system for balanced visual presentation
    p <- p +
      coord_cartesian(ylim = c(-1.5, 1.5)) +
      theme_minimal() +
      theme(
        # Clean panel styling
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        
        # Professional axis styling
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = if (show_y_axis) element_text(size = 9, color = "black") else element_blank(),
        axis.title.x = if (show_x_title) element_text(size = 10, face = "bold") else element_blank(),
        axis.title.y = if (show_y_title) element_text(size = 10, face = "bold") else element_blank(),
        axis.ticks.y = if (show_y_axis) element_line(color = "black", linewidth = 0.5) else element_blank(),
        
        # Appropriate margins for multi-panel layout
        plot.margin = margin(20, 25, 20, 15)
      ) +
      
      # Descriptive axis labels
      labs(
        x = if (show_x_title) "Hours after hyperosmotic stress" else "",
        y = if (show_y_title) "log2FoldChange(treated/untreated)" else ""
      )
    
    return(p)
  }
  
  # Create three panels with appropriate axis labeling
  # Only leftmost panel shows y-axis, only middle panel shows x-axis title
  p_static <- create_panel_with_stats("static", 
                                      show_y_axis = TRUE, 
                                      show_y_title = TRUE, 
                                      show_x_title = FALSE)
  
  p_gained <- create_panel_with_stats("gained", 
                                      show_y_axis = FALSE, 
                                      show_y_title = FALSE, 
                                      show_x_title = TRUE)
  
  p_lost <- create_panel_with_stats("lost", 
                                    show_y_axis = FALSE, 
                                    show_y_title = FALSE, 
                                    show_x_title = FALSE)
  
  # Combine panels into horizontal layout
  combined_plots <- plot_grid(p_static, p_gained, p_lost, ncol = 3, align = "h")
  
  # Add external legend for significance criteria
  # Positions legend outside plot area to avoid visual interference
  final_plot_with_legend <- ggdraw(combined_plots) +
    draw_text(
      "** p < 0.01, *** p < 0.001",
      x = 0.98,  # Far right positioning
      y = 0.02,  # Bottom margin positioning
      hjust = 1,
      vjust = 0,
      size = 9,
      color = "gray60",
      lineheight = 0.9
    )
  
  return(final_plot_with_legend)
}

## EXECUTE INTEGRATED ANALYSIS

# Prepare data for consecutive timepoint analysis
consecutive_data <- prepare_consecutive_data(combined)

# Conduct statistical testing for all loop types
test_results <- conduct_consecutive_tests(consecutive_data)

# Create integrated visualization with statistical annotations
integrated_plot <- create_integrated_plot(combined, test_results)

# Display the final integrated plot
print(integrated_plot)

## SAVE OUTPUTS

# Save the integrated plot with statistical testing
ggsave("plots/rnaseqTimecourse_loopAnchors.pdf", 
       integrated_plot, width = 12, height = 4.5, units = "in")

# Save square version for presentations
ggsave("plots/rnaseqTimecourse_loopAnchors_square.pdf", 
       integrated_plot, width = 8, height = 8, units = "in")

## SUMMARY STATISTICS

message("=== Integrated Analysis Summary ===")
for(loop_type in c("static", "gained", "lost")) {
  type_data <- combined |> filter(type == loop_type)
  n_genes <- length(unique(type_data$gene_id))
  n_observations <- nrow(type_data)
  
  # Statistical results summary
  results <- test_results[[loop_type]]
  significant_transitions <- results$transition[results$significance_stars != "ns"]
  
  message(sprintf("%s loops: %d genes, %d observations", 
                  stringr::str_to_title(loop_type), n_genes, n_observations))
  
  if(length(significant_transitions) > 0) {
    message(sprintf("  Significant transitions: %s", 
                    paste(significant_transitions, collapse = ", ")))
  } else {
    message("  No significant consecutive transitions detected")
  }
}
