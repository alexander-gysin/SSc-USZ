# GLOBAL CONFIGURATION FILE (00_config.R)
# Project: Systemic Sclerosis Disease Activity Multi-Omics


# 1. SCRIPT GROUPS ----------------------

base_scripts <- c(
  "analysis/index.Rmd",
  "analysis/about.Rmd",
  "analysis/license.Rmd"
)

clinical_scripts <- c(
  "analysis/01a_data_import_and_qc.Rmd",
  "analysis/01b_clinical_feature_engineering.Rmd",
  "analysis/02a_cohort_characterization.Rmd",
  "analysis/02b_clinical_clustering.Rmd"
)

proteomics_scripts <- c(
  "analysis/03a_proteomics_variance.Rmd",
  "analysis/03b_proteomics_differential.Rmd",
  "analysis/03c_proteomics_clustering.Rmd",
  "analysis/03d_proteomics_subgroups.Rmd",
  "analysis/03e_proteomics_correlations.Rmd",
  "analysis/03f_proteomics_machine_learning.Rmd"
)

lipidomics_scripts <- c(
  "analysis/04a_lipidomics_topography.Rmd",
  "analysis/04b_lipidomics_differential.Rmd",
  "analysis/04c_lipidomics_machine_learning.Rmd"
)

multiomics_scripts <- c(
  "analysis/05a_investigator.Rmd"
)

sideproject_scripts <- c(
  "analysis/z01_sideproject_ellen.Rmd",
  "analysis/z02_sideproject_vascular_2.Rmd"
)

# 2. AUTOMATIC MASTER LIST -------------------------

generate_all_scripts <- function() {
  c(
    base_scripts,
    clinical_scripts,
    proteomics_scripts,
    lipidomics_scripts,
    multiomics_scripts,
    sideproject_scripts
  )
}

# The canonical pipeline order. No manual maintenance required.
all_scripts <- generate_all_scripts()

# 3. SYNC FUNCTION -------------------------------

sync <- function(files = NULL, all = FALSE, publish = TRUE, preview = FALSE, message = NULL, verbose = TRUE) {

  # --- 1. Resolve Input Files ---
  if (all) {
    # If all = TRUE, we let workflowr determine the files later.
    target_files <- NULL
  } else if (is.null(files)) {
    # No argument: look for current_file in the global environment
    if (!exists("current_file", envir = .GlobalEnv)) {
      stop("Error: 'current_file' is not defined in the environment.", call. = FALSE)
    }
    target_files <- get("current_file", envir = .GlobalEnv)
  } else {
    target_files <- files
  }

  # --- 2. Validation and Formatting (if specific files are targeted) ---
  if (!is.null(target_files)) {
    # Check if files exist
    missing_files <- target_files[!file.exists(target_files)]
    if (length(missing_files) > 0) {
      stop(paste("Error: The following files do not exist:\n", paste(missing_files, collapse = "\n")), call. = FALSE)
    }

    # Deduplicate and reorder according to canonical all_scripts
    target_files <- unique(target_files)

    # Check for invalid files not in all_scripts
    invalid_files <- target_files[!target_files %in% all_scripts]
    if (length(invalid_files) > 0) {
      stop(paste("Error: Files not found in canonical all_scripts list:\n", paste(invalid_files, collapse = "\n")), call. = FALSE)
    }

    # Reorder based on the master list
    target_files <- all_scripts[all_scripts %in% target_files]

    # --- 3. Confirmation for all_scripts ---
    if (identical(target_files, all_scripts) && !preview) {
      cat(sprintf("You are about to publish %d analysis files.\n", length(target_files)))
      ans <- readline("Continue? [y/n]: ")
      if (tolower(ans) != "y") {
        stop("Sync cancelled by user.", call. = FALSE)
      }
    }
  }

  # --- 4. Handle Commit Message ---
  if (!preview) {
    if (is.null(message)) {
      message <- readline("Commit message: ")
      if (trimws(message) == "") {
        stop("Error: Commit message cannot be empty.", call. = FALSE)
      }
    } else if (isFALSE(message)) {
      message <- "Automated sync update" # Fallback if opt-out
    }
  }

  # --- 5. Preview Mode ---
  if (preview) {
    cat("\n=== PREVIEW MODE ===\n")
    if (publish) {
      cat("Files to be PUBLISHED:\n")
      if (all) {
        cat("- [workflowr automatically detected changes]\n")
      } else {
        cat(paste("-", target_files, collapse = "\n"), "\n")
      }
    } else {
      cat("PUBLISHING is skipped (publish = FALSE).\n")
    }
    cat("\nFiles to be COMMITTED:\n- All remaining uncommitted project changes (wflow_git_commit(all = TRUE))\n")
    cat("\nTarget to PUSH:\n- GitHub (Current branch)\n")
    cat("====================\n\n")
    return(invisible(TRUE))
  }

  # --- 6. Execution ---

  # PUBLISH
  if (publish) {
    if (verbose) cat("\nPublishing analyses...\n")

    tryCatch({
      if (all) {
        workflowr::wflow_publish(all = TRUE, message = message)
      } else {
        workflowr::wflow_publish(target_files, message = message)
      }
      if (verbose) cat("✓ Publishing completed\n")
    }, error = function(e) {
      stop(paste("Publishing failed with workflowr error:\n", e$message), call. = FALSE)
    })
  }

  # COMMIT ALL REMAINING
  if (verbose) cat("\nCommitting remaining project changes...\n")
  tryCatch({
    # Catch any supporting files (like 00_config.R) that publishing missed
    workflowr::wflow_git_commit(all = TRUE, message = paste(message, "(Catch-all commit)"))
    if (verbose) cat("✓ Commit completed\n")
  }, error = function(e) {
    stop(paste("Commit failed:\n", e$message), call. = FALSE)
  })

  # PUSH
  if (verbose) cat("\nPushing to GitHub...\n")
  tryCatch({
    workflowr::wflow_git_push()
    if (verbose) cat("✓ Push completed\n")
  }, error = function(e) {
    stop(paste("Push failed, but commits remain intact. Error:\n", e$message), call. = FALSE)
  })
}

# 4. SYNC STATUS FUNCTION --------------------


sync_status <- function() {
  cat("\n=========================================\n")
  cat("           PROJECT SYNC STATUS           \n")
  cat("=========================================\n\n")

  # --- 1. Workflowr Status ---
  cat("## WORKFLOWR STATUS ##\n")
  w_stat <- workflowr::wflow_status()
  status_df <- w_stat$status

  # Extract specific states
  scratch <- rownames(status_df)[status_df$scratch == TRUE]
  unpublished <- rownames(status_df)[status_df$published == FALSE & status_df$scratch == FALSE]
  modified <- rownames(status_df)[status_df$published == TRUE & status_df$mod_Rmd == TRUE]

  if (length(scratch) > 0) {
    cat("\n\u26A0\uFE0F Scratch (Untracked by Git, never published):\n")
    cat(paste("-", scratch, collapse = "\n"), "\n")
  }

  if (length(unpublished) > 0) {
    cat("\n\u26A0\uFE0F Unpublished (Tracked by Git, never published):\n")
    cat(paste("-", unpublished, collapse = "\n"), "\n")
  }

  if (length(modified) > 0) {
    cat("\n\u26A0\uFE0F Modified (Rmd file differs from existing HTML):\n")
    cat(paste("-", modified, collapse = "\n"), "\n")
  }

  if (length(scratch) == 0 && length(unpublished) == 0 && length(modified) == 0) {
    cat("\n\u2705 All analyses are up to date and published.\n")
  }

  # --- Workflowr Recommendations ---
  cat("\n\U0001F4A1 WORKFLOWR ACTIONS:\n")
  if (length(scratch) > 0 || length(modified) > 0 || length(unpublished) > 0) {
    cat("- Run `sync(all = TRUE)` to publish all detected changed Rmd files\n")
    cat("- Or target specific files: `sync(c(\"file1.Rmd\", \"file2.Rmd\"))`\n")
  } else {
    cat("- No action required.\n")
  }

  cat("\n-----------------------------------------\n")

  # --- 2. Git Status ---
  cat("## GIT STATUS ##\n")

  get_git <- function(cmd) {
    suppressWarnings(system(paste("git", cmd), intern = TRUE, ignore.stderr = TRUE))
  }

  branch <- get_git("rev-parse --abbrev-ref HEAD")
  cat(sprintf("\nBranch:\n%s\n", branch))

  get_git("fetch --quiet")
  ahead_behind <- get_git("rev-list --left-right --count HEAD...@{u}")
  is_ahead <- FALSE

  if (length(ahead_behind) > 0) {
    ab_split <- strsplit(ahead_behind, "\t")[[1]]
    is_ahead <- as.numeric(ab_split[1]) > 0
    is_behind <- as.numeric(ab_split[2]) > 0
    if (is_ahead | is_behind) { cat(sprintf("\n⚠️ ️GitHub status:\nahead by %s commits\nbehind by %s commits\n", ab_split[1], ab_split[2]))
    } else {cat(sprintf("\n\u2705 GitHub status:\nahead by %s commits & behind by %s commits\n", ab_split[1], ab_split[2]))}
  } else {
    cat("\nGitHub status: No upstream branch configured.\n")
  }

  status_raw <- get_git("status --porcelain")
  has_git_changes <- length(status_raw) > 0
  untracked_clean <- character(0)

  if (!has_git_changes) {
    cat("\nWorking tree:\nclean\n")
  } else {
    untracked <- status_raw[grepl("^\\?\\?", status_raw)]
    modified_git <- status_raw[!grepl("^\\?\\?", status_raw)]

    if (length(untracked) > 0) {
      untracked_clean <- sub("^\\?\\?\\s+", "", untracked)
      cat("\n\u26A0\uFE0F Untracked non-analysis files (Git has never seen these):\n")
      cat(paste("-", untracked_clean, collapse = "\n"), "\n")
    }
    if (length(modified_git) > 0) {
      cat("\n\u26A0\uFE0F Modified tracked non-analysis files (Tracked files with new changes):\n")
      cat(paste("-", sub("^\\s*[A-Z]+\\s+", "", modified_git), collapse = "\n"), "\n")
    }
  }

  # --- Git Recommendations ---
  cat("\n\U0001F4A1 GIT ACTIONS:\n")
  if (length(untracked_clean) > 0) {
    cat("- Run `untracked_files <- sync_status()` to assign the untracked files to a vector.\n")
    cat("- You can then print `untracked_files`, copy the ones you need, and commit them.\n")
  }
  if (has_git_changes) {
    cat("- Run `sync(publish = FALSE)` to commit non-analysis file changes (e.g. .R files).\n")
  }
  if (is_ahead && !has_git_changes && length(modified) == 0) {
    cat("- Run `workflowr::wflow_git_push()` to send your local commits to GitHub.\n")
  }
  if (!has_git_changes && !is_ahead) {
    cat("- No action required.\n")
  }

  # --- 3. Enriched Last Commit ---
  cat("\n-----------------------------------------\n")
  cat("## LAST COMMIT ##\n")
  # Format: Hash | YYYY-MM-DD HH:MM:SS | Author - Message
  last_commit <- get_git("log -1 --format=\"%h | %ci | %an - %s\"")
  cat(last_commit, "\n")

  cat("\n=========================================\n")

  # Invisibly return the vector of untracked files so it doesn't clutter the terminal output
  # but remains available for variable assignment.
  return(invisible(untracked_clean))
}

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
                             "Subset_diffuse",
                             "Puffy fingers",
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
                                         "Subset_diffuse",
                                         "Puffy fingers",
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

clinical_domains <- list(
  Demographics_and_Scores = c(
    "Age", "Sex", "Sample_Age", "cohort_group", "Total Score",
    "eustarAI", "ACTIVE_AI", "Active_our", "Inflamm_active", "fibrotic_active", "vascular_active"
  ),
  Cutaneous_and_Vascular = c(
    "Subset_diffuse",
    "Modified rodnan skin score", "mRSS_a", "mRSS Worsening (1yr)",
    "Raynaud phenomenon", "Score : Raynaud phenomenon", "Raynaud's present", "Raynaud VAS",
    "Digital Tip Ulcers", "Digital ulcer", "Fingertip pitting scars", "Pitting scars on finger tips",
    "Telangiectasia", "Telangiectasia (any)", "Abnormal Capillaroscopy", "Capillaroscopy Score",
    "Capillaroscopy Pattern", "Giant capillaries", "Hemorrhages", "Capillary loss", "Ramified bushy capillaries",
    "Puffy fingers", "Sclerodactyly", "Skin Thickening (Ext. MCP)", "Skin Thickening (Prox. MCP)",
    "Subcutaneous Calcinosis", "Gangrene", "log_Modified_rodnan_skin_score"
  ),
  Musculoskeletal = c(
    "Joint involvement", "Joint Contractures", "Tendon friction rubs",
    "Proximal muscle weakness", "Myalgia", "Muscle atrophy", "Activity of arthritis during the past week"
  ),
  Pulmonary = c(
    "PAH/ILD", "ILD diagnosed via HRCT", "FVC %pred", "Forced Vital Capacity (ml)",
    "DLCO (SB) (% predicted)", "DLCO/VA (% predicted)", "TLC %pred",
    "O2-saturation at rest (%)", "Min SpO2 (Exercise)",
    "Dyspnea NYHA Stage", "Max Borg Dyspnea", "Distance in m", "New PAH"
  ),
  Cardiac = c(
    "New Cardiac Manifest.",
    "Right bundle branch block", "Right axis deviation", "Right ventricular hypertrophy",
    "Ventricular arrhythmias", "Auricular arrhythmias", "Conduction blocks", "Arrhythmias requiring therapy",
    "Right atrium area (cm2)", "Right ventricular area (cm2)", "TAPSE (cm)",
    "Pericardial effusion on echo", "Diastolic Dysf. (Echo)", "LVEF %","log_Right_atrium_area__cm2", "log_TAPSE__cm"
  ),
  Gastrointestinal = c(
    "Stomach Symptoms", "GAVE", "Intestinal Symptoms", "Malabsorption syndrome", "Proximal Dysphagia"
  ),
  Labs_and_Biomarkers = c(
    "CRP", "Creatinine", "NT-proBNP (pg/ml)", "Uric Acid (mg/dl)",
    "Hemoglobin (g/dl)", "CK value in serum (U/L)", "Proteinuria (>300mg/d)",
    "SSc Antibodies", "SSc Antibody Score",  "log_CRP",
    "log_CK_value_in_serum__U_L", "log_Creatinine", "log_NT_proBNP__pg_ml", "log_Uric_Acid__mg_dl"
  ),
  Medications = c(
    "Med_rituximab",
    "Med_tocilizumab",                            "Med_mmf",                                    "Med_mtx",
    "Med_jaki",                                   "Med_abatacept",                              "Med_sildenafil",
    "Med_tadalafil",                              "Med_bosentan",                               "Med_iloprost",
    "Med_macitentan",                             "Med_ambrisentan",                            "Med_ppi",
    "Med_acei",                                   "Med_nintedanib",                             "Med_asa",
    "Med_ccb",                                    "Med_arb",                                    "Med_oac",
    "Immunosuppressants", "Major vascular", "Steroids", "Other_medication"
    )
)

# NEW: MASTER COLOR DICTIONARY
# A single, unified function that can map colors for any variable in the project
get_project_colors <- function(requested_levels, custom_map = NULL) {
  requested_levels <- as.character(requested_levels)

  base_map <- setNames(
    c(COLOR_ACTIVE, COLOR_NON_ACTIVE, COLOR_CONTROL, COLOR_MALE, COLOR_FEMALE),
    c(GRP_ACTIVE, GRP_NON_ACTIVE, GRP_CONTROL, SEX_MALE, SEX_FEMALE)
  )

  # Combine with custom overrides (custom takes precedence)
  if (!is.null(custom_map)) {
    full_map <- c(custom_map, base_map)
    full_map <- full_map[!duplicated(names(full_map))]
  } else {
    full_map <- base_map
  }

  out_colors <- full_map[requested_levels]

  # Smart Fallback for completely unknown levels (prevents dark grey blobs)
  unknown_idx <- is.na(out_colors)
  if (any(unknown_idx)) {
    num_unknown <- sum(unknown_idx)
    safe_pal <- if(exists("FALLBACK_CAT_PALETTE")) FALLBACK_CAT_PALETTE else c("#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")
    out_colors[unknown_idx] <- rep(safe_pal, length.out = num_unknown)
  }

  names(out_colors) <- requested_levels
  return(out_colors)
}

# GLOBAL GGPLOT THEME OVERRIDE
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

why_is_this_done <- function(...) {
  points <- list(...)
  li_items <- paste0("    <li>", unlist(points), "</li>", collapse = "\n")
  html <- sprintf('
<div class="why-box">
  <strong>Why is this done?</strong>
  <ul>
%s
  </ul>
</div>', li_items)
  return(html)
}
