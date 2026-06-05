library(dplyr)
library(writexl)
library(here)
library(tidyr)

# PARAMETERS -------------------------------------------------------------------
input_path  <- "~/SSc-USZ/output/01b_clinical_feature_engineering/curated_clinical_spine.rds"
output_dir  <- "~/SSc-USZ/output/randomizer"
file_name   <- "randomized_groups_for_ligation.xlsx"
id_col_name <- "Subject_ID"
strata_col  <- "cohort_group"
group_size  <- 8

# DIRECTORY SETUP --------------------------------------------------------------
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
final_output_path <- file.path(output_dir, file_name)

# PROCESSING -------------------------------------------------------------------
df <- readRDS(input_path)

set.seed(321)

assign_blocks <- function(data, size) {
  data %>%
    sample_frac(1) %>%
    mutate(group_number = ceiling(row_number() / size))
}

act_df    <- df %>% filter(!!sym(strata_col) == "Active")    %>% assign_blocks(3)
nonact_df <- df %>% filter(!!sym(strata_col) == "Non-Active") %>% assign_blocks(3)
ctrl_df   <- df %>% filter(!!sym(strata_col) == "Control")   %>% assign_blocks(2)

result_df <- bind_rows(act_df, nonact_df, ctrl_df) %>%
  # 1. Shuffle all rows completely
  sample_frac(1) %>%
  # 2. Arrange ONLY by group_number. Since the data is shuffled,
  # the order within each group will now be entirely random.
  arrange(group_number) %>%
  select(ID = all_of(id_col_name), !!sym(strata_col), group_number)

# Create the distribution table for the second sheet
summary_df <- result_df %>%
  group_by(group_number, !!sym(strata_col)) %>%
  tally() %>%
  pivot_wider(names_from = !!sym(strata_col), values_from = n, values_fill = 0) %>%
  # Adding a total column for easy checking in Excel
  mutate(total_in_group = rowSums(across(where(is.numeric), ~ .x), na.rm = TRUE))

# EXPORT -----------------------------------------------------------------------

# Create a named list: "Name of Sheet" = Data_Frame
export_list <- list(
  "Randomized_Groups" = result_df,
  "Distribution_Check" = summary_df
)

write_xlsx(export_list, path = final_output_path)

cat("Script Complete! File saved with two sheets to:", final_output_path)
