# ==============================================================================
# GLOBAL CONFIGURATION FILE (00_config.R)
# Project: Systemic Sclerosis Disease Activity Multi-Omics
# ==============================================================================
cat("Loading Global Configuration...\n")

# 1. GLOBAL LIBRARIES REQUIRED FOR CONFIG
library(ggplot2)

# 2. GLOBAL SEED
set.seed(2026)

# 3. GLOBAL GROUP DEFINITIONS
GRP_ACTIVE     <- "Active"
GRP_NON_ACTIVE <- "Non-Active"
GRP_CONTROL    <- "Control"
SEX_MALE       <- "Male"
SEX_FEMALE     <- "Female"

# 4. GLOBAL COLOR PALETTE
COLOR_ACTIVE     <- "#D55E00"
COLOR_NON_ACTIVE <- "#0072B2"
COLOR_CONTROL    <- "#999999"
COLOR_MALE       <- "#56B4E9"
COLOR_FEMALE     <- "#E69F00"
COLOR_HIGHLIGHT <- "#E6007E" # Magenta/Pink for shared trajectory drivers

# Age Gradient (Purple)
COLOR_AGE_LOW <- "white"
COLOR_AGE_HIGH <- "purple4"

# Heatmap Aesthetics
HM_Z_LOW  <- "dodgerblue4"
HM_Z_MID  <- "white"
HM_Z_HIGH <- "firebrick4"

# ------------------------------------------------------------------------------
# NEW: MASTER COLOR DICTIONARY
# ------------------------------------------------------------------------------
# A single, unified function that can map colors for any variable in the project
get_project_colors <- function(requested_levels) {

  color_map <- setNames(
    c(COLOR_ACTIVE, COLOR_NON_ACTIVE, COLOR_CONTROL, COLOR_MALE, COLOR_FEMALE),
    c(GRP_ACTIVE, GRP_NON_ACTIVE, GRP_CONTROL, SEX_MALE, SEX_FEMALE)
  )

  out_colors <- color_map[requested_levels]

  # Fallback for unexpected names
  out_colors[is.na(out_colors)] <- "#333333"
  names(out_colors) <- requested_levels

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
