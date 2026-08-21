# 05_sRNA_engines.R
# Engine functions for miRNA-Seq QC and Preprocessing via Command-Line Tools

# Packages -------------------------------------------------------
library(Rsamtools)
library(Rsubread)

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
                  max(1, ceiling(length(fastq_files) / 2))))

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
                  "-l", config$cutadapt$hard_crop,
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

# Alignment and Feature Counting Wrapper -------------------------

#' Master Orchestrator: Runs Alignment, BAM Sorting, and UMI Deduplication
#' @param fastq_files Character vector of absolute paths to processed fastq files
#' @param config Nested list of configurations from RMarkdown
#' @param scratch_dir Directory on the scratch drive for heavy I/O files
#' @return List containing diagnostics and absolute paths to deduplicated BAM files
wrap_alignment_dedup <- function(fastq_files, config, scratch_dir) {

  # Absolute path to isolated Miniconda tools
  conda_bin <- path.expand("~/miniconda3/envs/mirna_env/bin")
  cmd_umi   <- file.path(conda_bin, "umi_tools")

  # 1. Directory Management
  dirs <- worker_setup_align_directories(scratch_dir)

  # 2. Compute Safe Parallel Cores
  n_cores <- worker_calculate_cores(fastq_files, config)

  message(sprintf("\nStarting Alignment & Dedup: %d files. Allocating %d cores", length(fastq_files), n_cores))

  # 3. Phase 1-3: Align, Sort, and Dedup (Per Sample)
  message("\n--- Executing Align -> Sort -> Dedup ---")

  process_res <- pbmcapply::pbmclapply(fastq_files, function(fq) {

    # 3a. Align
    align_res <- worker_rsubread_align(fq, dirs$raw, config)
    if(align_res$status != "success") return(align_res)

    # 3b. Sort and Index
    sort_res <- worker_sort_index_bam(align_res$processed, dirs$sorted)
    if(sort_res$status != "success") return(sort_res)

    # 3c. UMI Dedup
    dedup_res <- worker_umi_dedup(sort_res$processed, dirs$dedup, cmd_umi)
    return(dedup_res)

  }, mc.cores = n_cores, mc.preschedule = FALSE)

  # Extract successful deduplicated BAM files and capture errors
  dedup_files <- sapply(process_res, function(x) x$processed)[sapply(process_res, function(x) x$status == "success")]
  error_logs  <- sapply(process_res, function(x) x$error)[sapply(process_res, function(x) x$status != "success")]

  return(list(
    status = "complete",
    diagnostics = list(
      total_samples = length(fastq_files),
      successful_processed = length(dedup_files),
      errors = error_logs
    ),
    paths = list(
      dedup_bam_files = dedup_files,
      dedup_bam_dir   = dirs$dedup
    )
  ))
}

# Alignment and Feature Counting Workers -------------------------------

#' Worker: Directory Management for Alignment
#' @param scratch_dir Directory on the scratch drive
#' @return Named list of directory paths
worker_setup_align_directories <- function(scratch_dir) {
  dirs <- list(
    raw    = file.path(scratch_dir, "05_bam_raw"),
    sorted = file.path(scratch_dir, "06_bam_sorted"),
    dedup  = file.path(scratch_dir, "07_bam_dedup")
  )

  lapply(dirs, function(x) {
    if (!dir.exists(x)) dir.create(x, recursive = TRUE)
  })

  return(dirs)
}

#' Worker: Align single FastQ using Rsubread
#' @param fq Path to processed fastq file
#' @param out_dir Directory for raw BAM output
#' @param config Configuration list containing reference and align parameters
#' @return List with input, processed path, status code, and error log
worker_rsubread_align <- function(fq, out_dir, config) {
  base_name <- tools::file_path_sans_ext(basename(fq))
  if(grepl("\\.fastq\\.gz$", basename(fq))) base_name <- tools::file_path_sans_ext(base_name)
  base_name <- sub("_processed$", "", base_name)

  bam_out <- file.path(out_dir, paste0(base_name, "_raw.bam"))
  err_log <- ""
  status_code <- "failed_align"

  tryCatch({
    # Strictly limit to 1 thread here to prevent core overloading in parallel loop
    Rsubread::align(
      index          = config$reference$index_path,
      readfile1      = fq,
      output_file    = bam_out,
      nthreads       = 1,
      type           = config$align$type,
      maxMismatches  = config$align$maxMismatches,
      unique         = config$align$unique,
      nBestLocations = config$align$nBestLocations
    )

    if (file.exists(bam_out) && file.info(bam_out)$size > 0) {
      status_code <- "success"
    } else {
      err_log <- "Align complete but BAM file missing or empty."
    }
  }, error = function(e) {
    err_log <<- paste("Rsubread::align failed:", e$message)
  })

  return(list(input = fq, processed = bam_out, status = status_code, error = err_log))
}

#' Worker: Sort and Index BAM file
#' @param bam_in Path to raw BAM file
#' @param out_dir Directory for sorted BAM output
#' @return List with input, processed path, status code, and error log
worker_sort_index_bam <- function(bam_in, out_dir) {
  base_name <- tools::file_path_sans_ext(basename(bam_in))
  base_name <- sub("_raw$", "", base_name)

  sort_prefix <- file.path(out_dir, paste0(base_name, "_sorted"))
  sorted_bam  <- paste0(sort_prefix, ".bam")
  err_log     <- ""
  status_code <- "failed_sort"

  tryCatch({
    Rsamtools::sortBam(file = bam_in, destination = sort_prefix)
    Rsamtools::indexBam(file = sorted_bam)

    if (file.exists(sorted_bam) && file.exists(paste0(sorted_bam, ".bai"))) {
      status_code <- "success"
    } else {
      err_log <- "Sort/Index complete but BAM or BAI missing."
    }
  }, error = function(e) {
    err_log <<- paste("Rsamtools sort/index failed:", e$message)
  })

  return(list(input = bam_in, processed = sorted_bam, status = status_code, error = err_log))
}

#' Worker: UMI Deduplication via UMI-tools
#' @param sorted_bam Path to sorted and indexed BAM file
#' @param out_dir Directory for deduplicated BAM output
#' @param cmd_umi Path to umi_tools executable
#' @return List with input, processed path, status code, and error log
worker_umi_dedup <- function(sorted_bam, out_dir, cmd_umi) {
  base_name <- tools::file_path_sans_ext(basename(sorted_bam))
  base_name <- sub("_sorted$", "", base_name)

  dedup_out <- file.path(out_dir, paste0(base_name, "_dedup.bam"))

  umi_args <- c("dedup",
                "--stdin", sorted_bam,
                "--stdout", dedup_out)

  umi_status <- system2(cmd_umi, umi_args, stdout = TRUE, stderr = TRUE)

  if (is.null(attr(umi_status, "status")) && file.exists(dedup_out) && file.info(dedup_out)$size > 0) {
    status_code <- "success"
    err_log <- ""
  } else {
    status_code <- "failed_umi_dedup"
    err_log <- paste(umi_status, collapse = "\n")
  }

  return(list(input = sorted_bam, processed = dedup_out, status = status_code, error = err_log))
}

#' Worker: Run FeatureCounts on all deduplicated BAMs
#' @param bam_files Character vector of absolute paths to deduplicated BAM files
#' @param config Configuration list containing feature_counts parameters
#' @return List containing the final Rsubread featureCounts object
worker_feature_counts <- function(bam_files, config) {
  n_cores <- worker_calculate_cores(bam_files, config)

  message(sprintf("Running featureCounts on %d BAM files using %d cores...", length(bam_files), n_cores))

  fc_res <- Rsubread::featureCounts(
    files                  = bam_files,
    annot.ext              = config$reference$gtf_path,
    isGTFAnnotationFile    = config$feature_counts$isGTFAnnotationFile,
    GTF.featureType        = config$feature_counts$GTF.featureType,
    GTF.attrType           = config$feature_counts$GTF.attrType,
    nthreads               = n_cores,
    strandSpecific         = config$feature_counts$strandSpecific,
    countMultiMappingReads = config$feature_counts$countMultiMappingReads,
    fraction               = config$feature_counts$fraction
  )

  return(fc_res)
}

#' Worker: Download External Reference Databases
#' @param url Character URL of the file
#' @param dest_file Character absolute path to the destination file
#' @return Logical TRUE if successful or already exists
worker_download_reference <- function(url, dest_file) {
  if (file.exists(dest_file)) {
    message(sprintf("File already exists: %s", dest_file))
    return(TRUE)
  }

  dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
  message(sprintf("Downloading %s ...", basename(dest_file)))

  old_timeout <- getOption("timeout")
  options(timeout = max(3600, old_timeout))

  tryCatch({
    download.file(url, destfile = dest_file, mode = "wb")
    options(timeout = old_timeout)
    return(file.exists(dest_file))
  }, error = function(e) {
    options(timeout = old_timeout)
    stop(paste("Download failed:", e$message))
  })
}
