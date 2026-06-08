library(tidyverse)
library(writexl)
library(here)

# 1. Load the old text-based mapping file (Using read_tsv for better encoding handling)
old_txt_path <- here::here("data", "clinical_data_type_dictionary.txt")
raw_map <- read_tsv(old_txt_path, col_names = TRUE, show_col_types = FALSE)
colnames(raw_map) <- c("Raw_Name", "Old_Type")

# BULLETPROOFING: Force valid UTF-8 encoding and clean out any corrupt hidden characters
raw_map$Raw_Name <- iconv(raw_map$Raw_Name, to = "UTF-8", sub = "")
raw_map$Old_Type <- replace_na(as.character(raw_map$Old_Type), "")

# 2. Recreate your hardcoded renaming dictionary (Old Name = Raw, New Name = Clean)
rename_vec <- c(
  "of the fingers as sclerodactyly (distal to the metacarpophalangeal joints but proximal to the proximal interphalangeal joints)" = "Sclerodactyly",
  "Score: Scleroderma related antibodies (any of anti-centromere, Anti-topoisomerase I [ anti-SCL 70], anti-RNA polymerase III" = "SSc Antibody Score",
  "Scleroderma related antibodies (any of anti-centromere, Anti-topoisomerase I [ anti-SCL 70], anti-RNA polymerase III" = "SSc Antibodies",
  "of the fingers of both hands extending proximally to the metacarpophalangeal joints" = "Skin Thickening (Ext. MCP)",
  "Skin thickening of the fingers of both hands extending proximal to the MCP joints" = "Skin Thickening (Prox. MCP)",
  "Does the patient have pulmonary arterial hypertension since the last visit?" = "New PAH",
  "Pulmonary arterial hypertension and/or Interstitial lung Disease" = "PAH/ILD",
  "Worst modified Borg dyspnea score during the test (from 0 to 10)" = "Max Borg Dyspnea",
  "Does the patient have cardiac manifestation since last visit?" = "New Cardiac Manifest.",
  "Diastolic function abnormal(on echo, E/a less than 10cm/sec)" = "Diastolic Dysf. (Echo)",
  "Physicians global assessment of disease activity VAS Unknown" = "MD Global Activity (NA)",
  "Patients global assessment of disease activity VAS Unknown" = "Pt Global Activity (NA)",
  "Was there a worsening of mRSS over the last 12 months?" = "mRSS Worsening (1yr)",
  "Intestinal symptoms (diarrhea, bloating, constipation)" = "Intestinal Symptoms",
  "Activity of arthritis during the past week Unknown" = "Arthritis Activity (NA)",
  "2013 EULAR/ACR classification criteria fulfilled" = "EULAR/ACR 2013 Criteria",
  "Tricuspid annular plane systolic excursion (cm)" = "TAPSE (cm)",
  "Stomach symptoms (early satiety, vomiting)" = "Stomach Symptoms",
  "Forced Vital Capacity (FVC- % predicted)" = "FVC %pred",
  "Score: Abnormal nailfold capillaroscopy" = "Capillaroscopy Score",
  "Total lung capacity (TLC- % predicted)" = "TLC %pred",
  "Gastric Antral Vascular Ectasia (GAVE)" = "GAVE",
  "Left ventricular ejection fraction %" = "LVEF %",
  "Worst O2-saturation at exercise (%)" = "Min SpO2 (Exercise)",
  "Capillaroscopy scleroderma pattern" = "Capillaroscopy Pattern",
  "Joint contractures (at any joint)" = "Joint Contractures",
  "Abnormal nailfold capillaroscopy" = "Abnormal Capillaroscopy",
  "Any non-Raynaud symptom present?" = "Non-Raynaud Symptom",
  "of the fingers as puffy fingers" = "Puffy Fingers",
  "OR serum creatinine (micromoles/L)" = "Creatinine",
  "CRP (mg/dl)" = "CRP",
  "Right atrium area (cm²)" = "Right atrium area (cm2)",
  "Right ventricular area (cm²)" = "Right ventricular area (cm2)"
)

# 3. Define the old exclusion lists
regex_patterns <- c("performed", "measured", "done", "velocity\\?", "ratio\\?", "area\\?", "excursion\\?", "fraction %\\?", "done\\?")

# 4. Build the new Master Mapping dataframe
master_mapping <- raw_map %>%
  mutate(
    # A. Apply Clean Names
    Clean_Name = ifelse(Raw_Name %in% names(rename_vec), rename_vec[Raw_Name], Raw_Name),

    # B. Initialize Data_Type and Exclusion_Reason
    Data_Type = ifelse(Old_Type == "", "UNKNOWN", Old_Type),
    Exclusion_Reason = ""
  ) %>%
  mutate(
    # C. Apply Exclusions (Using str_to_lower instead of tolower to prevent crashes)
    Exclusion_Reason = case_when(
      str_detect(str_to_lower(Raw_Name), paste(str_to_lower(regex_patterns), collapse = "|")) ~ "Regex Drop",
      Data_Type == "x" ~ "Mapped Drop",
      TRUE ~ ""
    ),
    Data_Type = ifelse(Exclusion_Reason != "", "x", Data_Type)
  ) %>%
  mutate(
    # D. Implement targeted Datatype fixes
    Data_Type = case_when(
      # 1. Date logic (Catches "Date", "datum", and the specific "onset_nonRP")
      Raw_Name == "Geburtsdatum" ~ "n", # Explicitly override Geburtsdatum (just a year) to numeric
      str_detect(str_to_lower(Raw_Name), "date|datum") ~ "d",
      Raw_Name == "onset_nonRP" ~ "d",

      # 2. Character-to-Numeric fixes
      str_detect(Clean_Name, "VAS") ~ "n",
      Clean_Name %in% c("Right atrium area (cm2)", "Right ventricular area (cm2)", "TLC %pred", "Arthritis Activity (NA)") ~ "n",

      # 3. Character-to-Categorical fixes
      Clean_Name %in% c("Diastolic Dysf. (Echo)", "PAH/ILD", "Right bundle branch block", "Geschlecht") ~ "c",

      # Keep whatever it already was
      TRUE ~ Data_Type
    )
  ) %>%
  # E. Clean up columns for the final Excel file
  select(Raw_Name, Clean_Name, Data_Type, Exclusion_Reason) %>%
  arrange(Data_Type, Clean_Name)

# 5. Export to Excel
output_path <- here::here("data", "clinical_master_mapping.xlsx")
write_xlsx(master_mapping, output_path)

cat("Successfully generated:", output_path, "\n")
