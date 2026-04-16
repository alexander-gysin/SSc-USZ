# ==============================================================================
# GLOBAL CONFIGURATION FILE (00_config.R)
# Project: Systemic Sclerosis Disease Activity Multi-Omics
# ==============================================================================

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
COLOR_HIGHLIGHT  <- "#E6007E" # Magenta/Pink for shared trajectory drivers

# Age Gradient (Purple)
COLOR_AGE_LOW <- "white"
COLOR_AGE_HIGH <- "purple4"

# Heatmap Aesthetics
HM_Z_LOW  <- "dodgerblue4"
HM_Z_MID  <- "white"
HM_Z_HIGH <- "firebrick4"


# Fallback Palettes for Exploratory Plotting (PCA, etc.)
# Used automatically when a variable doesn't have a strict color mapped above
FALLBACK_CAT_PALETTE  <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#0072B2", "#F0E442", "#999999")
FALLBACK_CONT_PALETTE <- "viridis" # Continuous fallback palette (viridis, plasma, magma, etc.)

# 5. CLINICAL VARIABLES SUBSETS
# Clinical description of SSC:
ssc_clinical_parameters <- c("EULAR/ACR 2013 Criteria",
                             "Abnormal Capillaroscopy",
                             "Raynaud phenomenon",
                             "Modified rodnan skin score",
                             "Subset_Diffuse",
                             "Puffy Fingers",
                             "Digital ulcer",
                             "Telangiectasia (any)",
                             "Joint involvement",
                             "Tendon friction rubs",
                             "Proximal muscle weakness",
                             "Anti-centromere (ACA)",
                             "Anti-Scl70",
                             "Anti-RNA polymerase III",
                             "CRP",
                             "CK value in serum (U/L)",
                             "ILD diagnosed via HRCT",
                             "DLCO (SB) (% predicted)",
                             "FVC %pred",
                             "Diastolic Dysf. (Echo)",
                             "LVEF %",
                             "Esophageal symptom at the time of visit",
                             "Stomach Symptoms",
                             "Intestinal Symptoms",
                             "eustarAI"
)

ssc_clinical_parameters_categorical <- c("EULAR/ACR 2013 Criteria",
                                         "Abnormal Capillaroscopy",
                                         "Raynaud phenomenon",
                                         "Subset_Diffuse",
                                         "Puffy Fingers",
                                         "Digital ulcer",
                                         "Telangiectasia (any)",
                                         "Joint involvement",
                                         "Tendon friction rubs",
                                         "Proximal muscle weakness",
                                         "Anti-centromere (ACA)",
                                         "Anti-Scl70",
                                         "Anti-RNA polymerase III",
                                         "ILD diagnosed via HRCT",
                                         "Diastolic Dysf. (Echo)",
                                         "Esophageal symptom at the time of visit",
                                         "Stomach Symptoms",
                                         "Intestinal Symptoms"
)

ssc_clinical_parameters_continuous <- c("Modified rodnan skin score",
                                        "CRP",
                                        "CK value in serum (U/L)",
                                        "DLCO (SB) (% predicted)",
                                        "FVC %pred",
                                        "LVEF %",
                                        "eustarAI"
)

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

# Custom Box Generators
what_this_does <- function(...) {
  points <- list(...)
  li_items <- paste0("    <li>", unlist(points), "</li>", collapse = "\n")
  html <- sprintf('
<div class="info-box">
  <strong>What this code does:</strong>
  <ul>
%s
  </ul>
</div>', li_items)
  return(html)
}

analysis_decision <- function(...) {
  points <- list(...)
  li_items <- paste0("    <li>", unlist(points), "</li>", collapse = "\n")
  html <- sprintf('
<div class="decision-box">
  <strong>Analysis Decisions:</strong>
  <ul>
%s
  </ul>
</div>', li_items)
  return(html)
}

key_insight <- function(...) {
  points <- list(...)
  li_items <- paste0("    <li>", unlist(points), "</li>", collapse = "\n")
  html <- sprintf('
<div class="insight-box">
  <strong>Key Insight:</strong>
  <ul>
%s
  </ul>
</div>', li_items)
  return(html)
}
