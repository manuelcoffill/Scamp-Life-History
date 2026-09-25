library(bayestestR)
library(readxl)
library(ggplot2)


# Read the first sheet
postfreel0 <- read_excel("PosteriorsFreeL0.xlsx")
postknownl0 <- read_excel("PosteriorsKnownL0.xlsx")

head(postfreel0)
head(postknownl0)
colnames(postfreel0)





###############################################################################
# Posterior Difference Comparisons Between Growth Model Parameter Estimates
#
# Purpose:
# Compare posterior distributions of parameter differences between:
#
#   1. Free L0 model:
#      Multiple reader - Consensus reads
#
#   2. Known L0 model:
#      Multiple reader - Consensus reads
#
#   3. Multiple reader estimates:
#      Free L0 - Known L0
#
#   4. Consensus read estimates:
#      Free L0 - Known L0
#
# For each comparison:
#   - Calculate posterior difference distribution
#   - Calculate MAP estimate
#   - Calculate 90% equal-tailed interval (ETI)
#   - Plot posterior density
#
###############################################################################

library(ggplot2)
library(dplyr)
library(bayestestR)


###############################################################################
# Function to compare posterior distributions for one parameter
###############################################################################

plot_parameter_difference <- function(parameter,
                                      postfreel0,
                                      postknownl0){
  
  ###########################################################################
  # Identify posterior columns
  #
  # Example:
  # parameter = "Linf"
  #
  # Creates:
  #   Linf[1] = Multiple reader
  #   Linf[2] = Consensus reads
  ###########################################################################
  
  col1 <- paste0(parameter, "[1]")
  col2 <- paste0(parameter, "[2]")
  
  
  ###########################################################################
  # Calculate posterior distributions of differences
  ###########################################################################
  
  diff_df <- bind_rows(
    
    #-----------------------------------------------------------------------
    # Within Free L0 model
    #-----------------------------------------------------------------------
    
    tibble(
      Difference =
        postfreel0[[col1]] - postfreel0[[col2]],
      
      Comparison =
        "Free L0: Multiple reader - Consensus reads"
    ),
    
    
    #-----------------------------------------------------------------------
    # Within Known L0 model
    #-----------------------------------------------------------------------
    
    tibble(
      Difference =
        postknownl0[[col1]] - postknownl0[[col2]],
      
      Comparison =
        "Known L0: Multiple reader - Consensus reads"
    ),
    
    
    #-----------------------------------------------------------------------
    # Multiple reader comparison between models
    #-----------------------------------------------------------------------
    
    tibble(
      Difference =
        postfreel0[[col1]] - postknownl0[[col1]],
      
      Comparison =
        "Multiple reader: Free L0 - Known L0"
    ),
    
    
    #-----------------------------------------------------------------------
    # Consensus read comparison between models
    #-----------------------------------------------------------------------
    
    tibble(
      Difference =
        postfreel0[[col2]] - postknownl0[[col2]],
      
      Comparison =
        "Consensus reads: Free L0 - Known L0"
    )
    
  )
  
  
  ###########################################################################
  # Calculate MAP and 90% credible intervals
  ###########################################################################
  
  summary_df <- diff_df %>%
    
    group_by(Comparison) %>%
    
    summarise(
      
      MAP =
        point_estimate(
          Difference,
          centrality = "MAP",
          threshold = 0.10
        )$MAP,
      
      
      CI_low =
        ci(
          Difference,
          ci = 0.90
        )$CI_low,
      
      
      CI_high =
        ci(
          Difference,
          ci = 0.90
        )$CI_high,
      
      
      .groups = "drop"
      
    )
  
  
  ###########################################################################
  # Create annotation text
  ###########################################################################
  
  annotation_text <- paste0(
    
    summary_df$Comparison,
    
    "\nMAP = ",
    round(summary_df$MAP, 2),
    
    "\n90% ETI = [",
    round(summary_df$CI_low, 2),
    
    ", ",
    round(summary_df$CI_high, 2),
    
    "]"
    
  )
  
  
  annotation_text <- paste(
    annotation_text,
    collapse = "\n\n"
  )
  
  
  ###########################################################################
  # Define colors for each comparison
  ###########################################################################
  
  comparison_colors <- c(
    
    "Free L0: Multiple reader - Consensus reads" =
      "#0072B2",
    
    "Known L0: Multiple reader - Consensus reads" =
      "#D55E00",
    
    "Multiple reader: Free L0 - Known L0" =
      "#009E73",
    
    "Consensus reads: Free L0 - Known L0" =
      "#CC79A7"
    
  )
  
  
  ###########################################################################
  # Generate density plot
  ###########################################################################
  
  p <- ggplot(
    
    diff_df,
    
    aes(
      x = Difference,
      colour = Comparison
    )
    
  ) +
    
    
    # Posterior density curves
    geom_density(
      linewidth = 1.2
    ) +
    
    
    # Zero difference reference line
    geom_vline(
      
      xintercept = 0,
      
      linetype = "dashed",
      
      color = "black",
      
      linewidth = 0.8
      
    ) +
    
    
    # MAP lines
    geom_vline(
      
      data = summary_df,
      
      aes(
        xintercept = MAP,
        colour = Comparison
      ),
      
      linewidth = 0.8,
      
      show.legend = FALSE
      
    ) +
    
    
    # Summary statistics annotation
    annotate(
      
      "text",
      
      x = Inf,
      
      y = Inf,
      
      label = annotation_text,
      
      hjust = 1.05,
      
      vjust = 1.05,
      
      size = 3.5,
      
      family = "mono"
      
    ) +
    
    
    scale_colour_manual(
      
      values = comparison_colors
      
    ) +
    
    
    labs(
      
      title =
        paste(
          parameter
          # "Posterior Difference Comparisons"
        ),
      
      x =
        "Posterior Difference",
      
      y =
        "Posterior Density",
      
      colour =
        NULL
      
    ) +
    
    
    theme_bw(
      
      base_size = 14
      
    ) +
    
    theme(
      
      legend.position = "bottom",
      
      legend.direction = "vertical",
      
      panel.grid = element_blank(),
      
      plot.title =
        element_text(
          face = "bold"
        )
      
    )
  
  
  return(p)
  
}



###############################################################################
# Generate posterior difference plots
###############################################################################

plot_parameter_difference(
  "Linf",
  postfreel0,
  postknownl0
)


plot_parameter_difference(
  "k",
  postfreel0,
  postknownl0
)


plot_parameter_difference(
  "Lknot",
  postfreel0,
  postknownl0
)


plot_parameter_difference(
  "cv_proc",
  postfreel0,
  postknownl0
)


plot_parameter_difference(
  "tknot",
  postfreel0,
  postknownl0
)


plot_parameter_difference(
  "nu_proc",
  postfreel0,
  postknownl0
)























################################################################################
###############################################################################
# Combine all parameter plots into one manuscript figure
###############################################################################

library(patchwork)
library(cowplot)


###############################################################################
# Generate individual parameter plots
###############################################################################

p_Linf <- plot_parameter_difference(
  "Linf",
  postfreel0,
  postknownl0
)

p_k <- plot_parameter_difference(
  "k",
  postfreel0,
  postknownl0
)

p_Lknot <- plot_parameter_difference(
  "Lknot",
  postfreel0,
  postknownl0
)

p_cv_proc <- plot_parameter_difference(
  "cv_proc",
  postfreel0,
  postknownl0
)

p_tknot <- plot_parameter_difference(
  "tknot",
  postfreel0,
  postknownl0
)

p_nu_proc <- plot_parameter_difference(
  "nu_proc",
  postfreel0,
  postknownl0
)


###############################################################################
# Extract a single shared legend
###############################################################################

shared_legend <- get_legend(
  
  p_Linf +
    theme(
      legend.position = "bottom"
    )
  
)


###############################################################################
# Remove legends from all plots
###############################################################################

p_Linf <- p_Linf + theme(legend.position="none")
p_k <- p_k + theme(legend.position="none")
p_Lknot <- p_Lknot + theme(legend.position="none")
p_cv_proc <- p_cv_proc + theme(legend.position="none")
p_tknot <- p_tknot + theme(legend.position="none")
p_nu_proc <- p_nu_proc + theme(legend.position="none")


###############################################################################
# Create final 2-column x 3-row layout
#
# Sixth panel replaced with legend
###############################################################################

final_figure <-
  
  (p_Linf | p_k) /
  
  (p_Lknot | p_cv_proc) /
  
  (p_tknot | wrap_elements(shared_legend))
final_figure

###############################################################################
# Save figure
###############################################################################

ggsave(
  
  filename =
    "Posterior_difference_comparisons_all_parameters.png",
  
  plot =
    final_figure,
  
  width =
    12,
  
  height =
    16,
  
  dpi =
    600
  
)








###############################################################################
# Combine selected parameter plots into one figure
#
# Layout:
#
#              Column 1                 Column 2
#
# Row 1       Linf                      k
#
# Row 2       cv_proc                   Shared legend
#
# Uses a single x-axis and y-axis label for the entire figure
###############################################################################

library(patchwork)
library(cowplot)


###############################################################################
# Generate plots
###############################################################################

p_Linf <- plot_parameter_difference(
  "Linf",
  postfreel0,
  postknownl0
)

p_k <- plot_parameter_difference(
  "k",
  postfreel0,
  postknownl0
)

p_cv_proc <- plot_parameter_difference(
  "cv_proc",
  postfreel0,
  postknownl0
)


###############################################################################
# Extract shared legend
###############################################################################

shared_legend <- get_legend(
  
  p_Linf +
    theme(
      legend.position = "bottom"
    )
  
)


###############################################################################
# Remove legends and individual axis labels
###############################################################################

remove_labels <- function(p){
  
  p +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  
}


p_Linf <- remove_labels(p_Linf)

p_k <- remove_labels(p_k)

p_cv_proc <- remove_labels(p_cv_proc)



###############################################################################
# Create empty space for legend panel
###############################################################################

legend_panel <- wrap_elements(shared_legend)



###############################################################################
# Assemble 2 x 2 figure
###############################################################################

final_figure <-
  
  (p_Linf | p_k) /
  
  (p_cv_proc | legend_panel) +
  
  plot_annotation(
    
    # title =
    #   "Posterior Difference Comparisons",
    
    theme =
      theme(
        
        plot.title =
          element_text(
            face="bold",
            size=16,
            hjust=0.5
          )
        
      )
    
  ) &
  
  theme(
    
    # Shared axis labels
    plot.margin =
      margin(10,10,10,10)
    
  )


###############################################################################
# Add unified axis labels
###############################################################################

final_figure <-
  
  final_figure +
  
  plot_annotation(
    
    theme =
      theme(
        plot.margin = margin(10,10,10,10)
      )
    
  )

final_figure

###############################################################################
# Save figure
###############################################################################

ggsave(
  
  filename =
    "Posterior_difference_comparisons_selected_parameters.png",
  
  plot =
    final_figure,
  
  width =
    12,
  
  height =
    10,
  
  dpi =
    600
  
)








































###############################################################################
# Extract MAP and 90% ETI for all posterior difference comparisons
#
# Output:
#   Parameter
#   Comparison
#   MAP
#   90% ETI lower
#   90% ETI upper
#
# Saves:
#   Posterior_difference_summary.xlsx
###############################################################################

library(dplyr)
library(bayestestR)
library(openxlsx)


###############################################################################
# Function to extract posterior difference summaries
###############################################################################

extract_parameter_summary <- function(parameter,
                                      postfreel0,
                                      postknownl0){
  
  ###########################################################################
  # Define posterior columns
  ###########################################################################
  
  col1 <- paste0(parameter,"[1]")
  col2 <- paste0(parameter,"[2]")
  
  
  ###########################################################################
  # Calculate posterior differences
  ###########################################################################
  
  diff_df <- bind_rows(
    
    tibble(
      Difference =
        postfreel0[[col1]] - postfreel0[[col2]],
      
      Comparison =
        "Free L0: Multiple reader - Consensus reads"
    ),
    
    
    tibble(
      Difference =
        postknownl0[[col1]] - postknownl0[[col2]],
      
      Comparison =
        "Known L0: Multiple reader - Consensus reads"
    ),
    
    
    tibble(
      Difference =
        postfreel0[[col1]] - postknownl0[[col1]],
      
      Comparison =
        "Multiple reader: Free L0 - Known L0"
    ),
    
    
    tibble(
      Difference =
        postfreel0[[col2]] - postknownl0[[col2]],
      
      Comparison =
        "Consensus reads: Free L0 - Known L0"
    )
    
  )
  
  
  ###########################################################################
  # Calculate MAP and 90% ETI
  ###########################################################################
  
  summary_df <- diff_df %>%
    
    group_by(Comparison) %>%
    
    summarise(
      
      MAP =
        point_estimate(
          Difference,
          centrality = "MAP",
          threshold = 0.10
        )$MAP,
      
      
      ETI_90_lower =
        ci(
          Difference,
          ci = 0.90
        )$CI_low,
      
      
      ETI_90_upper =
        ci(
          Difference,
          ci = 0.90
        )$CI_high,
      
      
      .groups = "drop"
      
    )
  
  
  ###########################################################################
  # Add parameter name
  ###########################################################################
  
  summary_df <- summary_df %>%
    
    mutate(
      Parameter = parameter
    ) %>%
    
    select(
      Parameter,
      Comparison,
      MAP,
      ETI_90_lower,
      ETI_90_upper
    )
  
  
  return(summary_df)
  
}



###############################################################################
# Run extraction for all parameters
###############################################################################

parameters <- c(
  "Linf",
  "k",
  "Lknot",
  "cv_proc",
  "tknot",
  "nu_proc"
)


parameter_summary <-
  
  lapply(
    
    parameters,
    
    extract_parameter_summary,
    
    postfreel0 = postfreel0,
    
    postknownl0 = postknownl0
    
  ) %>%
  
  bind_rows()



###############################################################################
# View summary table
###############################################################################

parameter_summary



###############################################################################
# Save as Excel workbook
###############################################################################

write.xlsx(
  
  parameter_summary,
  
  file =
    "Posterior_difference_summary.xlsx",
  
  overwrite = TRUE
  
)
