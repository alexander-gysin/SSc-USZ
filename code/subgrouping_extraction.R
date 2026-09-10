library(dplyr)
library(tibble)
library(writexl)

omics_mat <- readRDS(here::here("output", "01a_data_import_and_qc", "clean_olink_matrix.rds"))
master_spine  <- readRDS(here::here("output", "01b_clinical_feature_engineering", "curated_clinical_spine.rds"))

# -----------------------------
# 1. Select the proteins of interest
# -----------------------------

# Get all CXCL proteins plus IDO1 and CASP4
proteins_keep <- rownames(omics_mat)[
  grepl("^CXCL", rownames(omics_mat))
]

proteins_keep <- unique(c("IDO1", proteins_keep, "CASP4", "IFNG"))

# Keep only proteins that actually exist in the matrix
proteins_keep <- proteins_keep[proteins_keep %in% rownames(omics_mat)]

# -----------------------------
# 2. Subset the Olink matrix
# -----------------------------

olink_subset <- omics_mat[proteins_keep, ]

# -----------------------------
# 3. Transpose so patients become rows
# -----------------------------

olink_df <- as.data.frame(t(olink_subset))

# Add Subject_ID as a column
olink_df <- rownames_to_column(olink_df, "Subject_ID")

# -----------------------------
# 4. Add patient information
# -----------------------------

final_df <- master_spine %>%
  select(Subject_ID, inCentraxx, Entnahmedatum, ACTIVE_AI, Age, Sex) %>%
  inner_join(olink_df, by = "Subject_ID")

clust_res <- readRDS(here::here("output", "03c_proteomics_clustering", "Clustering", "job_svc", "job_svc_Master_Assignments.rds"))

# ---> NEW: Read the updated list-based RDS from 03d, extract the raw matrix, and keep ONLY the 50/50 split
sub_res_list <- readRDS(here::here("output", "03d_proteomics_subgroups", "non_active_splits_raw_and_resid.rds"))
sub_res <- sub_res_list$raw %>% select(Subject_ID, Split_50_50)

# join clustering result to the final_df
joined_df <- inner_join(final_df, clust_res, by = "Subject_ID")
joined_df <- inner_join(joined_df, sub_res, by = "Subject_ID")

# compare clustering results (Updated to check the 50/50 split)
joined_df <- joined_df %>%
  mutate(
    Concordance = ifelse(
      (kmeans_Cluster == "Cluster 1" & Split_50_50 == "Control-Like") |
        (kmeans_Cluster == "Cluster 2" & Split_50_50 == "Active-Like"),
      1,
      0
    )
  )

# Subset to nonactive only
export_df <- joined_df %>%
  filter(grepl("Nonact", Subject_ID))

# -----------------------------
# 5. Export to Excel
# -----------------------------

# Output paths
current_file <- "subgrouping_extraction"
output_dir_data <- here::here("output", current_file)
if (!dir.exists(output_dir_data)) dir.create(output_dir_data, recursive = TRUE)

write_xlsx(
  export_df,
  file.path(output_dir_data, "Nonact_Protein_Expression_Subgrouping_Extraction_age_resid_50_50_split.xlsx")
)
