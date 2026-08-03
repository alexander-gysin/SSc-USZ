# scratch directory
scratch_dir <- "~/scratch/SSc-USZ"

# GLOBAL CONFIGURATION FILE (00_config.R)
# Project: Systemic Sclerosis Disease Activity Multi-Omics



# 3. SYNC FUNCTION -------------------------------
sync <- function(files = NULL, all = FALSE, publish = FALSE, preview = FALSE) {

  # 1. Environment Health Checks
  is_git <- suppressWarnings(system("git rev-parse --is-inside-work-tree", intern = TRUE, ignore.stderr = TRUE)[1] == "true")
  is_wflow <- file.exists("_workflowr.yml") || file.exists("analysis/_site.yml")
  is_renv <- file.exists("renv.lock") && requireNamespace("renv", quietly = TRUE)

  if (any(!c(is_git, is_wflow, is_renv)) && !preview) {
    cat("\n\u26A0\uFE0F ENVIRONMENT WARNINGS:\n")
    if (!is_git) cat("- Git is missing. Run `git init` or check your repository path.\n")
    if (!is_wflow) cat("- Workflowr is missing. Run `workflowr::wflow_start()`.\n")
    if (!is_renv) cat("- renv is missing. Run `renv::init()`.\n")
    cat("-----------------------------------------\n")
  }

  # 2. Input Resolution
  if (all) publish <- TRUE
  target_files <- files

  # If publishing without passing specific files, grab current_file
  if (is.null(target_files) && publish && !all) {
    if (!exists("current_file", envir = .GlobalEnv)) {
      stop("Error: 'current_file' is undefined. Either pass a file, or use all = TRUE.", call. = FALSE)
    }
    target_files <- get("current_file", envir = .GlobalEnv)
  }

  # Validate specific files & enforce all_scripts order if it exists
  if (!is.null(target_files)) {
    missing <- target_files[!file.exists(target_files)]
    if (length(missing) > 0) stop(sprintf("Missing files:\n%s", paste(missing, collapse = "\n")), call. = FALSE)

    target_files <- unique(target_files)
    if (exists("all_scripts", envir = .GlobalEnv)) {
      not_in_list <- target_files[!target_files %in% all_scripts]

      if (length(not_in_list) > 0) {
        cat("\n\U0001F4A1 HINT: The following files are not in the canonical `all_scripts` list but will be processed:\n")
        cat(paste("-", not_in_list, collapse = "\n"), "\n")
      }

      # Reorder canonical files according to all_scripts, and append the non-canonical ones
      in_list <- all_scripts[all_scripts %in% target_files]
      target_files <- c(in_list, not_in_list)

      if (identical(target_files, all_scripts) && !preview) {
        if (tolower(readline("Publish all analyses? [y/n]: ")) != "y") stop("Sync cancelled.", call. = FALSE)
      }
    }
  }

  # 3. Message Prompting
  if (!preview) {
    msg_prompt <- readline("Commit message (Leave blank for default): ")
    if (trimws(msg_prompt) != "") {
      commit_msg <- trimws(msg_prompt)
    } else {
      # Smart defaults
      if (!is.null(target_files) && length(target_files) == 1) {
        commit_msg <- sprintf("Update %s", basename(target_files))
      } else if (!is.null(target_files) && length(target_files) > 1) {
        commit_msg <- sprintf("Update %d files", length(target_files))
      } else {
        commit_msg <- "routine sync"
      }
    }
  } else {
    commit_msg <- "[Preview Mode Message]"
  }

  # 4. Preview Mode
  if (preview) {
    cat("\n=== PREVIEW MODE ===\n")
    if (publish) {
      cat(sprintf("Action: PUBLISH, COMMIT & PUSH\nTarget: %s\n", if (all) "All detected changes" else paste(target_files, collapse = ", ")))
    } else {
      cat(sprintf("Action: COMMIT & PUSH (Hourly Backup)\nTarget: %s + all modified tracked files\n", if(is.null(target_files)) "None" else paste(target_files, collapse = ", ")))
    }
    cat(sprintf("Message: \"%s\"\n====================\n\n", commit_msg))
    return(invisible(TRUE))
  }

  # 5. Execution
  tryCatch({
    if (publish) {
      cat("\n\U0001F680 Publishing and committing...\n")
      if (all) {
        workflowr::wflow_publish(all = TRUE, message = commit_msg)
      } else {
        workflowr::wflow_publish(target_files, message = commit_msg)
      }
    } else {
      cat("\n\U0001F4BE Committing changes...\n")
      workflowr::wflow_git_commit(files = target_files, all = TRUE, message = commit_msg)
    }
  }, error = function(e) stop(sprintf("Commit/Publish failed:\n%s", e$message), call. = FALSE))

  # 6. Push with Error Handling
  tryCatch({
    cat("\u2601\uFE0F Pushing to GitHub...\n")
    workflowr::wflow_git_push()
    cat("\u2705 Sync complete!\n")
  }, error = function(e) {
    cat("\n\u274C PUSH FAILED:\n", e$message, "\n\n")
    cat("\U0001F4A1 INSTRUCTIONS: Your files were successfully committed, but the push failed.\n")
    cat("1. Check your internet connection.\n")
    cat("2. If this is a new repository/branch, you may need to set the upstream:\n")
    cat("   Run: `system(\"git push -u origin HEAD\")`\n")
    cat("3. If there is a remote conflict, pull changes first: `workflowr::wflow_git_pull()`\n")
    cat("4. If you are constantly asked for your username/password on the cluster,\n")
    cat("   run this in your terminal to save them: `git config --global credential.helper store`\n")
  })
}

# 4. SYNC STATUS FUNCTION --------------------
sync_status <- function() {
  cat("\n=========================================\n")
  cat("           PROJECT SYNC STATUS           \n")
  cat("=========================================\n")

  untracked_clean <- character(0)

  # Helper
  get_git <- function(cmd) suppressWarnings(system(paste("git", cmd), intern = TRUE, ignore.stderr = TRUE))

  # 1. Workflowr Status
  cat("\n## WORKFLOWR STATUS ##\n")
  is_wflow <- file.exists("_workflowr.yml") || file.exists("analysis/_site.yml")

  if (!is_wflow) {
    cat("\u26A0\uFE0F Workflowr not initialized.\n")
  } else {
    w_stat <- tryCatch(workflowr::wflow_status(), error = function(e) NULL)
    if (is.null(w_stat)) {
      cat("\u26A0\uFE0F Error reading Workflowr status.\n")
    } else {
      st <- w_stat$status
      scratch <- rownames(st)[st$scratch == TRUE]
      unpub <- rownames(st)[st$published == FALSE & st$scratch == FALSE]
      mod <- rownames(st)[st$published == TRUE & st$mod_Rmd == TRUE]

      if (length(scratch) > 0) cat(sprintf("\n\u26A0\uFE0F Scratch (Untracked, no HTML):\n- %s\n", paste(scratch, collapse = "\n- ")))
      if (length(unpub) > 0) cat(sprintf("\n\u26A0\uFE0F Unpublished (Tracked, no HTML):\n- %s\n", paste(unpub, collapse = "\n- ")))
      if (length(mod) > 0) cat(sprintf("\n\u26A0\uFE0F Modified (Rmd differs from HTML):\n- %s\n", paste(mod, collapse = "\n- ")))
      if (sum(length(scratch), length(unpub), length(mod)) == 0) cat("\n\u2705 All analyses published.\n")

      cat("\n\U0001F4A1 WORKFLOWR ACTIONS:\n")
      if (length(scratch) > 0) {
        cat("- Track a scratch file (no HTML): `sync(\"analysis/your_file.Rmd\")`\n")
        cat("- Publish a scratch file: `sync(\"analysis/your_file.Rmd\", publish = TRUE)`\n")
      }
      if (length(mod) > 0 || length(unpub) > 0) {
        cat("- Publish all tracked Rmd changes: `sync(all = TRUE)`\n")
      }
    }
  }

  cat("\n-----------------------------------------\n")

  # 2. Git Status
  cat("## GIT STATUS ##\n")
  is_git <- length(get_git("rev-parse --is-inside-work-tree")) > 0 && get_git("rev-parse --is-inside-work-tree")[1] == "true"

  if (!is_git) {
    cat("\u26A0\uFE0F Git not initialized.\n")
  } else {
    cat(sprintf("\nBranch: %s\n", get_git("rev-parse --abbrev-ref HEAD")))

    get_git("fetch --quiet")
    ab <- get_git("rev-list --left-right --count HEAD...@{u}")
    is_ahead <- FALSE

    if (length(ab) > 0) {
      ab_split <- strsplit(ab, "\t")[[1]]
      is_ahead <- as.numeric(ab_split[1]) > 0
      icon <- if (is_ahead || as.numeric(ab_split[2]) > 0) "\u26A0\uFE0F" else "\u2705"
      cat(sprintf("\n%s GitHub status: ahead by %s, behind by %s commits\n", icon, ab_split[1], ab_split[2]))
    } else {
      cat("\n\u26A0\uFE0F GitHub status: No upstream branch configured.\n")
    }

    stat_raw <- get_git("status --porcelain")
    if (length(stat_raw) == 0) {
      cat("\nWorking tree: clean\n")
    } else {
      untracked <- stat_raw[grepl("^\\?\\?", stat_raw)]
      modified_git <- stat_raw[!grepl("^\\?\\?", stat_raw)]

      if (length(untracked) > 0) {
        untracked_clean <- sub("^\\?\\?\\s+", "", untracked)
        cat(sprintf("\n\u26A0\uFE0F Untracked files:\n- %s\n", paste(untracked_clean, collapse = "\n- ")))
      }
      if (length(modified_git) > 0) {
        cat(sprintf("\n\u26A0\uFE0F Modified files:\n- %s\n", paste(sub("^\\s*[A-Z]+\\s+", "", modified_git), collapse = "\n- ")))
      }
    }

    # Check for Git credential helper config
    cred_helper <- get_git("config --global credential.helper")
    if (length(cred_helper) > 0 && grepl("store|cache|manager", cred_helper[1])) {
      cat("\n\u2705 Git credential helper is configured.\n")
    } else {
      cat("\n\u26A0\uFE0F Git credential helper not found. If pushing asks for a password,\n")
      cat("   run this in your terminal to save them the next time you are asked: \n`git config --global credential.helper store`\n
        KEEP IN MIND THAT THIS SAVES YOUR PAT IN PLAIN TEXT IN YOUR HOME DIRECTORY\n
        NOBODY ELSE SHOULD HAVE ACCESS TO YOUR HOME DIRECTORY TO KEEP YOUR PAT PRIVATE\n")
    }

    cat("\n\U0001F4A1 GIT ACTIONS:\n")
    if (length(untracked_clean) > 0) cat("- Untracked files returned invisibly. Assign with `untracked <- sync_status()`.\n")
    if (length(stat_raw) > 0) cat("- Hourly backup (commit all tracked changes): `sync()`\n")
    if (is_ahead && length(stat_raw) == 0) cat("- Push local commits: `workflowr::wflow_git_push()`\n")

    cat("\n## LAST COMMIT ##\n")
    cat(get_git("log -1 --format=\"%h | %ci | %an - %s\""), "\n")
  }

  cat("\n-----------------------------------------\n")

  # 3. Renv Status
  cat("## RENV STATUS ##\n")
  is_renv <- file.exists("renv.lock") && requireNamespace("renv", quietly = TRUE)

  if (!is_renv) {
    cat("\u26A0\uFE0F renv not initialized.\n")
  } else {
    r_out <- capture.output(suppressMessages(tryCatch(renv::status(), error = function(e) "Error")))
    if (length(r_out) == 0 || any(grepl("is in sync|is synchronized|no issues found", r_out, ignore.case = TRUE))) {
      cat("\n\u2705 renv environment is synchronized.\n")
    } else {
      cat(sprintf("\n\u26A0\uFE0F renv is out of sync:\n%s\n", paste(r_out, collapse = "\n")))
      cat("\n\U0001F4A1 RENV ACTIONS:\n")
      cat("- Update lockfile: `renv::snapshot()`\n")
      cat("- Revert environment: `renv::restore()`\n")
    }
  }

  cat("\n=========================================\n")
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
FALLBACK_CAT_PALETTE  <- c("#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#D55E00", "#0072B2", "#F0E442", "#999999")
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
