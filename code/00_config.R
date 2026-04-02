# ==============================================================================
# GLOBAL CONFIGURATION FILE (00_config.R)
# Project: Systemic Sclerosis Disease Activity Multi-Omics
# ==============================================================================
cat("Loading Global Configuration...\n")

# 1. GLOBAL LIBRARIES REQUIRED FOR CONFIG
library(ggplot2)

# 2. GLOBAL SEED
# Guarantees identical permutations and plot jittering forever
set.seed(2026)

# 3. GLOBAL GROUP DEFINITIONS
# Using these constants prevents typos in downstream scripts
GRP_ACTIVE     <- "Active"
GRP_NON_ACTIVE <- "Non-Active"
GRP_CONTROL    <- "Control"

# 4. GLOBAL COLOR PALETTE (Colorblind Safe)
COLOR_ACTIVE     <- "#D55E00" # Vermilion (Red/Orange)
COLOR_NON_ACTIVE <- "#0072B2" # Blue
COLOR_CONTROL    <- "#999999" # Grey

# ------------------------------------------------------------------------------
# DYNAMIC COLOR PALETTE MAPPER
# ------------------------------------------------------------------------------
# Safely maps your specific groups to their global colors, regardless of
# which specific comparisons are running (e.g., if a plot only has 2 groups)
get_group_colors <- function(group_levels) {

  # Master dictionary
  color_map <- setNames(
    c(COLOR_ACTIVE, COLOR_NON_ACTIVE, COLOR_CONTROL),
    c(GRP_ACTIVE, GRP_NON_ACTIVE, GRP_CONTROL)
  )

  # Extract requested colors
  out_colors <- color_map[group_levels]

  # Fallback for any unexpected group names to prevent ggplot crashes
  out_colors[is.na(out_colors)] <- "#333333"
  names(out_colors) <- group_levels

  return(out_colors)
}

# ------------------------------------------------------------------------------
# GLOBAL GGPLOT THEME OVERRIDE
# ------------------------------------------------------------------------------
theme_project_base <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2)),
      plot.subtitle = element_text(color = "grey30", size = rel(0.9)),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = rel(1.1))
    )
}

cat("Configuration Successfully Loaded!\n")
