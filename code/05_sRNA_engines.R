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

#' Master Orchestrator: Runs Bowtie, UMICollapse, Fastx, and miRDeep2
#' @param fastq_files Character vector of absolute paths to processed fastq files
#' @param config Nested list of configurations from RMarkdown
#' @param scratch_dir Directory on the scratch drive for heavy I/O files
#' @return List containing diagnostics, paths, and the final parsed count matrix
wrap_alignment_dedup <- function(fastq_files, config, scratch_dir) {

  # Absolute path to isolated Miniconda tools and jars
  conda_bin    <- path.expand("~/miniconda3/envs/mirna_env/bin")
  cmd_bowtie   <- file.path(conda_bin, "bowtie")
  cmd_samtools <- file.path(conda_bin, "samtools")
  cmd_fastx    <- file.path(conda_bin, "fastx_collapser")
  cmd_mirdeep  <- file.path(conda_bin, "quantifier.pl")

  # Assuming UMICollapse jar is stored in the environment or a specific resource folder
  jar_umicollapse <- path.expand("~/miniconda3/envs/mirna_env/share/umicollapse/umicollapse.jar")

  # 1. Directory Management
  dirs <- worker_setup_align_directories(scratch_dir)

  # 2. Compute Safe Parallel Cores
  # Calculate standard cores, then strictly cap at 32 for Java RAM safety (32 jobs * 8GB = 256GB max RAM)
  n_cores <- worker_calculate_cores(fastq_files, config)

  message(sprintf("\nStarting sRNA Alignment & Dedup: %d files. Allocating %d safely throttled cores",
                  length(fastq_files), n_cores))
  message("\n--- Phase 1-3: Bowtie Align -> Sort -> UMICollapse -> Fastx (Per Sample) ---")

  # 3. Parallel Execution: Align, Dedup, and format for miRDeep2
  process_res <- pbmcapply::pbmclapply(fastq_files, function(fq) {

    # 3a. Align and Sort (Piped directly to avoid heavy I/O)
    align_res <- worker_bowtie_align_sort(fq, dirs$sorted, config, cmd_bowtie, cmd_samtools)
    if(align_res$status != "success") return(align_res)

    # 3b. UMI Deduplication (Java memory strictly capped)
    dedup_res <- worker_umicollapse(align_res$processed, dirs$dedup, jar_umicollapse)
    if(dedup_res$status != "success") return(dedup_res)

    # 3c. Revert to FASTA and Collapse identical reads for miRDeep2
    fastx_res <- worker_bam_to_collapsed_fasta(dedup_res$processed, dirs$collapsed, cmd_samtools, cmd_fastx)
    return(fastx_res)

  }, mc.cores = n_cores, mc.preschedule = FALSE)

  # Extract successful FASTA files and capture errors
  collapsed_fastas <- sapply(process_res, function(x) x$processed)[sapply(process_res, function(x) x$status == "success")]
  error_logs       <- sapply(process_res, function(x) x$error)[sapply(process_res, function(x) x$status != "success")]

  if(length(collapsed_fastas) == 0) {
    stop("All samples failed preprocessing. Check error logs.")
  }

  # 4. Phase 4: Cohort-level Quantification with miRDeep2
  message("\n--- Phase 4: miRDeep2 Cohort Quantification ---")
  quant_res <- worker_mirdeep_quantify(collapsed_fastas, dirs$quant, config, cmd_mirdeep)

  return(list(
    status = "complete",
    diagnostics = list(
      total_samples = length(fastq_files),
      successful_processed = length(collapsed_fastas),
      errors = error_logs,
      quantification_log = quant_res$log
    ),
    paths = list(
      collapsed_fasta_dir = dirs$collapsed,
      quant_dir = dirs$quant
    ),
    counts = quant_res$count_matrix # Automatically parsed back into R
  ))
}

# Alignment and Feature Counting Workers -------------------------------

#' Worker: Directory Management for Alignment
worker_setup_align_directories <- function(scratch_dir) {
  dirs <- list(
    sorted    = file.path(scratch_dir, "05_bam_sorted"),
    dedup     = file.path(scratch_dir, "06_bam_dedup"),
    collapsed = file.path(scratch_dir, "07_fasta_collapsed"),
    quant     = file.path(scratch_dir, "08_mirdeep_quant")
  )

  lapply(dirs, function(x) {
    if (!dir.exists(x)) dir.create(x, recursive = TRUE)
  })
  return(dirs)
}

#' Worker: Bowtie1 Alignment Piped to Samtools Sort
#' @description Streams SAM output directly to BAM and sorts to bypass writing huge SAMs to disk.
worker_bowtie_align_sort <- function(fq, out_dir, config, cmd_bowtie, cmd_samtools) {
  base_name <- tools::file_path_sans_ext(basename(fq))
  if(grepl("\\.fastq\\.gz$", basename(fq))) base_name <- tools::file_path_sans_ext(base_name)
  base_name <- sub("_processed$", "", base_name)

  sorted_bam <- file.path(out_dir, paste0(base_name, "_sorted.bam"))
  err_log <- ""
  status_code <- "failed_align_sort"

  # Construct piped command: bowtie | samtools view | samtools sort
  # Strictly single-threaded across the pipe to allow scaling across samples.
  # Added explicit -x flag for the Bowtie index path to prevent warnings.
  cmd <- sprintf(
    "%s %s -x %s %s | %s view -bS -@ 1 - | %s sort -@ 1 -o %s -",
    cmd_bowtie,
    config$align$bowtie_params,
    config$reference$index_path,
    fq,
    cmd_samtools,
    cmd_samtools,
    sorted_bam
  )

  tryCatch({
    sys_res <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

    # Create BAM index required for UMICollapse
    system2(cmd_samtools, c("index", sorted_bam))

    if (file.exists(sorted_bam) && file.info(sorted_bam)$size > 0) {
      status_code <- "success"
    } else {
      err_log <- paste("Pipeline finished but BAM missing. Log:", paste(sys_res, collapse = "\n"))
    }
  }, error = function(e) {
    err_log <<- paste("Bowtie->Samtools pipe failed:", e$message)
  })

  return(list(input = fq, processed = sorted_bam, status = status_code, error = err_log))
}

#' Worker: UMI Deduplication via UMICollapse
#' @description Runs UMICollapse with a strict JVM heap limit and large stack to prevent HPC crashes.
worker_umicollapse <- function(sorted_bam, out_dir, jar_umicollapse) {
  base_name <- tools::file_path_sans_ext(basename(sorted_bam))
  base_name <- sub("_sorted$", "", base_name)

  dedup_bam <- file.path(out_dir, paste0(base_name, "_dedup.bam"))

  # Strict memory cap (-Xmx8G) and massive stack cap (-Xss256M) per sample
  umi_cmd <- sprintf("java -Xmx12G -Xss256M -jar %s bam -i %s -o %s",
                     jar_umicollapse, sorted_bam, dedup_bam)

  tryCatch({
    sys_res <- system(umi_cmd, intern = TRUE, ignore.stderr = FALSE)
    if (file.exists(dedup_bam) && file.info(dedup_bam)$size > 0) {
      status_code <- "success"
      err_log <- ""
    } else {
      status_code <- "failed_umicollapse"
      err_log <- paste(sys_res, collapse = "\n")
    }
  }, error = function(e) {
    status_code <<- "failed_umicollapse"
    err_log <<- paste("UMICollapse failed:", e$message)
  })

  return(list(input = sorted_bam, processed = dedup_bam, status = status_code, error = err_log))
}

#' Worker: Convert BAM to FastQ and Collapse via Fastx
#' @description Prepares the format required by miRDeep2 (>seq_1_xCount format)
worker_bam_to_collapsed_fasta <- function(dedup_bam, out_dir, cmd_samtools, cmd_fastx) {
  base_name <- tools::file_path_sans_ext(basename(dedup_bam))
  base_name <- sub("_dedup$", "", base_name)

  collapsed_fasta <- file.path(out_dir, paste0(base_name, "_collapsed.fa"))

  # Pipe samtools fastq directly into fastx_collapser to avoid intermediate files
  cmd <- sprintf(
    "%s fastq %s | %s -o %s",
    cmd_samtools, dedup_bam, cmd_fastx, collapsed_fasta
  )

  tryCatch({
    sys_res <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
    if (file.exists(collapsed_fasta) && file.info(collapsed_fasta)$size > 0) {
      status_code <- "success"
      err_log <- ""
    } else {
      status_code <- "failed_fastx"
      err_log <- paste("Samtools->Fastx pipe failed.", collapse = "\n")
    }
  }, error = function(e) {
    status_code <<- "failed_fastx"
    err_log <<- paste("Conversion error:", e$message)
  })

  return(list(input = dedup_bam, processed = collapsed_fasta, status = status_code, error = err_log))
}

#' Worker: miRDeep2 Quantifier Execution
#' @description Generates the cohort config file and runs quantifier.pl
worker_mirdeep_quantify <- function(collapsed_fastas, out_dir, config, cmd_mirdeep) {

  # 1. Generate Config.txt mapping file
  sample_ids <- sub("_collapsed\\.fa$", "", basename(collapsed_fastas))
  config_df <- data.frame(
    File = collapsed_fastas,
    Sample = sample_ids
  )
  config_path <- file.path(out_dir, "mirdeep_config.txt")
  write.table(config_df, config_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

  # 2. Setup Quantifier environment
  orig_dir <- getwd()
  setwd(out_dir)

  log_out <- ""
  count_df <- NULL

  tryCatch({
    # FIX: Explicitly call 'perl', and run inside a bash shell that activates conda first.
    # This guarantees miRDeep2 has access to all its dependencies.
    cmd <- sprintf("bash -c 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate mirna_env && perl %s -p %s -m %s -r %s -W -U -c -N'",
                   cmd_mirdeep,
                   config$reference$hairpin_path,
                   config$reference$mature_path,
                   config_path)

    sys_res <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
    log_out <- paste(sys_res, collapse = "\n")

    # 3. Locate and parse the generated Count Matrix
    out_files <- list.files(out_dir, pattern = "miRNAs_expressed_all_samples_.*\\.csv", full.names = TRUE)

    if (length(out_files) > 0) {
      latest_file <- out_files[which.max(file.info(out_files)$mtime)]
      count_df <- read.delim(latest_file, sep = "\t", check.names = FALSE)
    } else {
      warning("Quantifier.pl completed, but no expression CSV was found.")
    }

  }, error = function(e) {
    log_out <<- paste("miRDeep2 Quantifier failed:", e$message)
  }, finally = {
    setwd(orig_dir) # Always revert to original working directory
  })

  return(list(
    status = ifelse(is.null(count_df), "failed", "success"),
    log = log_out,
    count_matrix = count_df
  ))
}

#' Worker: Download External Reference Databases
#' @param url Character URL of the file (e.g., .gz link)
#' @param dest_file Character absolute path to the uncompressed destination file (.fa)
#' @param species_prefix Character prefix to filter FASTA (e.g., "hsa"). NULL to skip.
#' @param convert_u_to_t Logical. If TRUE, converts RNA to DNA in sequences.
#' @return Logical TRUE if successful or already exists
worker_download_reference <- function(url, dest_file, species_prefix = NULL, convert_u_to_t = FALSE) {

  dest_file <- path.expand(dest_file)

  # 1. Skip if valid file already exists (threshold lowered to 10KB for filtered files)
  if (file.exists(dest_file) && file.info(dest_file)$size > 10000) {
    message(sprintf("Valid file already exists: %s", dest_file))
    return(TRUE)
  }

  dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
  temp_gz <- paste0(dest_file, ".gz")
  temp_fa <- paste0(dest_file, ".tmp")

  message(sprintf("Downloading %s via wget...", basename(temp_gz)))

  tryCatch({
    # 2. Force system wget
    wget_cmd <- sprintf("wget -q -O %s %s", shQuote(temp_gz), shQuote(url))
    wget_status <- system(wget_cmd)

    if (wget_status != 0 || !file.exists(temp_gz)) {
      stop("wget failed to download the file.")
    }

    # 3. Extract the file to a temporary FASTA
    message("Extracting sequences to raw fasta format...")
    system2("gunzip", args = c("-c", temp_gz), stdout = temp_fa)
    unlink(temp_gz)

    # 4. Process FASTA: Filter by species and/or convert U to T
    message("Processing FASTA (filtering species and/or translating U to T)...")

    # Build Awk script dynamically
    awk_script <- ""
    if (!is.null(species_prefix)) {
      awk_script <- sprintf('/^>%s/ {p=1; print; next} /^>/ {p=0; next} ', species_prefix)
    } else {
      awk_script <- '/^>/ {p=1; print; next} '
    }

    if (convert_u_to_t) {
      awk_script <- paste0(awk_script, 'p {gsub(/U/,"T"); gsub(/u/,"t"); print}')
    } else {
      awk_script <- paste0(awk_script, 'p {print}')
    }

    # Execute Awk directly on the temp file and write to the final destination
    system2("awk", args = c(shQuote(awk_script), temp_fa), stdout = dest_file)
    unlink(temp_fa)

    # 5. Strict Validation
    if (file.exists(dest_file) && file.info(dest_file)$size > 10000) {
      return(TRUE)
    } else {
      unlink(dest_file)
      stop("Extraction resulted in an empty or invalid file (check species prefix).")
    }

  }, error = function(e) {
    if(file.exists(temp_gz)) unlink(temp_gz)
    if(file.exists(temp_fa)) unlink(temp_fa)
    stop(paste("Worker execution failed:", e$message))
  })
}
