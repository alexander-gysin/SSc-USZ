# 05_sRNA_engines.R
# Engine functions for miRNA-Seq QC and Preprocessing via Command-Line Tools

# Preprocessing Wrapper ------------------------------------------

#' Master Orchestrator: Runs FastQC, UMI-Tools, Cutadapt, and MultiQC
#' @param fastq_files Character vector of absolute paths to raw R1 fastq files
#' @param config Nested list of configurations from RMarkdown
#' @param scratch_dir Directory on the scratch drive for heavy I/O files
#' @param docs_dir Directory in the workflowr project to copy lightweight MultiQC HTMLs
#' @return List containing diagnostics and absolute paths to processed scratch files
wrap_preprocessing <- function(fastq_files, config, scratch_dir, docs_dir) {

  # --- Absolute paths to isolated Miniconda tools ---
  conda_bin    <- path.expand("~/miniconda3/envs/mirna_env/bin")
  cmd_fastqc   <- file.path(conda_bin, "fastqc")
  cmd_multiqc  <- file.path(conda_bin, "multiqc")
  cmd_umi      <- file.path(conda_bin, "umi_tools")
  cmd_cutadapt <- file.path(conda_bin, "cutadapt")

  # 1. Directory Management
  dirs <- worker_setup_directories(scratch_dir, docs_dir)

  # 2. Compute Safe Parallel Cores
  n_cores <- worker_calculate_cores(fastq_files, config)

  message(sprintf("\nStarting pipeline: %d files. Allocating %d cores", length(fastq_files), n_cores))
  message(sprintf("(Limits - Hard: %d, Relative: %d, File-driven: %d)",
                  config$system$max_cores,
                  max(1, floor(parallel::detectCores() * config$system$core_ratio)),
                  max(1, floor(length(fastq_files) / 2))))

  # 3. Phase 1: Raw FastQC
  message("\n--- Phase 1: Raw FastQC ---")
  raw_qc_res <- pbmcapply::pbmclapply(fastq_files, function(fq) {
    worker_fastqc(fq, dirs$raw_qc, cmd_fastqc)
  }, mc.cores = n_cores)

  # 4. Phase 2: Raw MultiQC
  message("\n--- Phase 2: Raw MultiQC ---")
  worker_multiqc(dirs$raw_qc, "raw_multiqc_report.html", docs_dir, cmd_multiqc)

  # 5. Phase 3 & 4: UMI Extraction and Adapter Trimming (Per Sample)
  message("\n--- Phase 3 & 4: UMI Extract & Cutadapt ---")
  process_res <- pbmcapply::pbmclapply(fastq_files, function(fq) {
    worker_umi_cutadapt(fq, dirs$umi, dirs$trim, config, cmd_umi, cmd_cutadapt)
  }, mc.cores = n_cores, mc.preschedule = FALSE) # DYNAMIC TASK SCHEDULING APPLIED HERE

  # Extract successful trimmed files and capture errors
  processed_files <- sapply(process_res, function(x) x$processed)[sapply(process_res, function(x) x$status == "success")]
  error_logs <- sapply(process_res, function(x) x$error)[sapply(process_res, function(x) x$status != "success")]

  # 6. Phase 5: Processed FastQC
  message("\n--- Phase 5: Processed FastQC ---")
  pro_qc_res <- pbmcapply::pbmclapply(processed_files, function(fq) {
    worker_fastqc(fq, dirs$pro_qc, cmd_fastqc)
  }, mc.cores = n_cores)

  # 7. Phase 6: Processed MultiQC
  message("\n--- Phase 6: Processed MultiQC ---")
  worker_multiqc(dirs$pro_qc, "processed_multiqc_report.html", docs_dir, cmd_multiqc)

  # Return absolute paths pointing to the scratch drive and diagnostics
  return(list(
    status = "complete",
    diagnostics = list(
      total_samples = length(fastq_files),
      successful_processed = length(processed_files),
      errors = error_logs
    ),
    paths = list(
      processed_fastq_dir = dirs$trim
    )
  ))
}

# Legacy Alias Support
wrap_preprocessing_pipeline <- wrap_preprocessing

# Preprocessing Workers -----------------------------------------------

#' Worker: Directory Management
#' @param scratch_dir Directory on the scratch drive
#' @param docs_dir Directory for lightweight MultiQC HTMLs
#' @return Named list of directory paths
worker_setup_directories <- function(scratch_dir, docs_dir) {
  dirs <- list(
    raw_qc = file.path(scratch_dir, "01_raw_fastqc"),
    umi    = file.path(scratch_dir, "02_umi_extracted"),
    trim   = file.path(scratch_dir, "03_trimmed_reads"),
    pro_qc = file.path(scratch_dir, "04_processed_fastqc")
  )

  # Create scratch directories
  lapply(dirs, function(x) {
    if (!dir.exists(x)) dir.create(x, recursive = TRUE)
  })

  # Create docs directory for lightweight HTMLs
  if (!dir.exists(docs_dir)) dir.create(docs_dir, recursive = TRUE)

  return(dirs)
}

#' Worker: Compute Safe Parallel Cores
#' @param fastq_files Character vector of fastq files
#' @param config Nested list of configurations
#' @return Integer representing the number of cores to use
worker_calculate_cores <- function(fastq_files, config) {
  available_cores <- parallel::detectCores()

  hard_cap     <- config$system$max_cores
  relative_cap <- max(1, floor(available_cores * config$system$core_ratio))
  file_cap     <- max(1, floor(length(fastq_files) / 2)) # Half-file rule

  # The lowest limit wins
  return(min(hard_cap, relative_cap, file_cap))
}

#' Worker: Run FastQC on a single file
#' @param fq Path to fastq file
#' @param out_dir Output directory for FastQC results
#' @param cmd_fastqc Path to FastQC executable
#' @return List with file path and status
worker_fastqc <- function(fq, out_dir, cmd_fastqc) {
  args <- c(fq, "-o", out_dir, "-t", "1") # 1 thread per job, parallelized over samples
  res <- system2(cmd_fastqc, args, stdout = TRUE, stderr = TRUE)
  return(list(
    file = fq,
    status = ifelse(is.null(attr(res, "status")), "success", "failed")
  ))
}

#' Worker: Run MultiQC and copy report
#' @param target_dir Directory containing FastQC results and where MultiQC will run
#' @param report_name Name of the output HTML report
#' @param docs_dir Destination directory for the HTML report
#' @param cmd_multiqc Path to MultiQC executable
worker_multiqc <- function(target_dir, report_name, docs_dir, cmd_multiqc) {
  multiqc_args <- c(target_dir, "-f", "-n", report_name, "-o", target_dir)
  system2(cmd_multiqc, multiqc_args)

  # Copy ONLY the HTML file to the workflowr docs folder
  file.copy(
    from = file.path(target_dir, report_name),
    to = file.path(docs_dir, report_name),
    overwrite = TRUE
  )
}

#' Worker: UMI Extraction and Adapter Trimming for a single file
#' @param fq Path to raw fastq file
#' @param dir_umi Directory for UMI output
#' @param dir_trim Directory for trimmed output
#' @param config Configuration list
#' @param cmd_umi Path to umi_tools executable
#' @param cmd_cutadapt Path to cutadapt executable
#' @return List with input, processed path, status code, and error log
worker_umi_cutadapt <- function(fq, dir_umi, dir_trim, config, cmd_umi, cmd_cutadapt) {
  base_name <- tools::file_path_sans_ext(basename(fq))
  if(grepl("\\.fastq\\.gz$", basename(fq))) base_name <- tools::file_path_sans_ext(base_name)

  umi_out  <- file.path(dir_umi, paste0(base_name, "_umi.fastq.gz"))
  trim_out <- file.path(dir_trim, paste0(base_name, "_processed.fastq.gz"))

  err_log <- ""

  # 5a. UMI-Tools Extract (Single-threaded by default)
  umi_args <- c("extract",
                "--extract-method=regex",
                paste0("--bc-pattern=", shQuote(config$umi$pattern)),
                "--stdin", fq,
                "--stdout", umi_out)
  umi_status <- system2(cmd_umi, umi_args, stdout = TRUE, stderr = TRUE)

  # Aggressive check: Did UMI tools succeed and create a valid file?
  if (is.null(attr(umi_status, "status")) && file.exists(umi_out) && file.info(umi_out)$size > 0) {

    # 5b. Cutadapt with Strict Thread Pinning (-j 1), Quality Trimming, and Max/Min Length
    cut_args <- c("-a", config$cutadapt$adapter,
                  "-m", config$cutadapt$min_length,
                  "-M", config$cutadapt$max_length,
                  "-q", config$cutadapt$quality_cutoff,
                  "-j", "1", # Force 1 core per cutadapt instance
                  "-o", trim_out,
                  umi_out)
    cut_status <- system2(cmd_cutadapt, cut_args, stdout = TRUE, stderr = TRUE)

    # Aggressive check: Did Cutadapt succeed and create a valid file?
    if (is.null(attr(cut_status, "status")) && file.exists(trim_out) && file.info(trim_out)$size > 0) {
      status_code <- "success"
    } else {
      status_code <- "failed_cutadapt"
      err_log <- paste(cut_status, collapse = "\n")
    }

  } else {
    status_code <- "failed_umi"
    err_log <- paste(umi_status, collapse = "\n")
  }

  return(list(input = fq, processed = trim_out, status = status_code, error = err_log))
}


