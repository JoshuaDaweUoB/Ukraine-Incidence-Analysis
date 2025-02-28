# load packages
pacman::p_load(arsenal, survival, readxl, dplyr)

## set wd
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Ukraine data/Data")

# load data
analysis_data_hiv_clean <- read_excel("HIV_data_clean.xlsx")

# baseline characteristics sex work
analysis_data_hiv_clean <- analysis_data_hiv_clean %>%
  group_by(ID) %>%
  mutate(id_seq = row_number())

analysis_data_hiv_bl <- analysis_data_hiv_clean %>%
  group_by(ID) %>%
  mutate(id_seq = row_number(),
         rec_sw_sell_90d = ifelse(all(is.na(rec_sw_sell_90d)), NA_real_, max(rec_sw_sell_90d, na.rm = TRUE)),
         rec_incarc_6m = ifelse(all(is.na(rec_incarc_6m)), NA_real_, max(rec_incarc_6m, na.rm = TRUE)),
         rec_oat_6m = ifelse(all(is.na(rec_oat_6m)), NA_real_, max(rec_oat_6m, na.rm = TRUE)),
         rec_homeless = ifelse(all(is.na(rec_homeless)), NA_real_, max(rec_homeless, na.rm = TRUE))) %>%
  ungroup() %>%
  subset(id_seq == 1)

# convert numeric to factor variables
analysis_data_hiv_bl$sex <- factor(analysis_data_hiv_bl$sex)
analysis_data_hiv_bl$rec_homeless <- factor(analysis_data_hiv_bl$rec_homeless)
analysis_data_hiv_bl$rec_incarc_6m <- factor(analysis_data_hiv_bl$rec_incarc_6m)
analysis_data_hiv_bl$rec_oat_6m <- factor(analysis_data_hiv_bl$rec_oat_6m)
analysis_data_hiv_bl$rec_sw_sell_90d <- factor(analysis_data_hiv_bl$rec_sw_sell_90d)

# convert numeric to factor variables
analysis_data_hiv_clean$sex <- factor(analysis_data_hiv_clean$sex)
analysis_data_hiv_clean$rec_homeless <- factor(analysis_data_hiv_clean$rec_homeless)
analysis_data_hiv_clean$rec_incarc_6m <- factor(analysis_data_hiv_clean$rec_incarc_6m)
analysis_data_hiv_clean$rec_oat_6m <- factor(analysis_data_hiv_clean$rec_oat_6m)

# table
tab_bl_sw_sex_hiv <- tableby(sex ~ age + rec_homeless + rec_incarc_6m + rec_oat_6m + rec_sw_sell_90d + inj_dur + rec_inj_freq_month, data=analysis_data_hiv_bl)
summary(tab_bl_sw_sex_hiv, text=TRUE)

## HIV analysis ##

## incidence rate calculations

# overall incidence rate
total_days_hiv <- sum(analysis_data_hiv_clean$days_risk)
total_cases <- sum(analysis_data_hiv_clean$hiv_rslt)
incidence_rate <- (total_cases / total_days_hiv) * 365.25 * 100

# Calculate 95% confidence intervals
incidence_rate_se <- sqrt(total_cases) / total_days_hiv * 365.25 * 100
ci_lower <- incidence_rate - 1.96 * incidence_rate_se
ci_upper <- incidence_rate + 1.96 * incidence_rate_se

cat("Incidence rate of HIV per 100 person years:", incidence_rate, "\n")
cat("95% CI:", ci_lower, "-", ci_upper, "\n")

### sex work HIV incidence rate calculations ###

# Filter data for males
analysis_data_hiv_men <- subset(analysis_data_hiv_clean, sex == 1)

### sex work HIV incidence rate calculations ###

# Function to manually calculate rate ratio for recent sex work with additional 0.5 cases
calculate_manual_rate_ratio <- function(data, group_label) {
  # Calculate incidence rate for recent sex work exposure with additional 0.5 cases
  total_days_hiv_sw_recent <- sum(data$days_risk[data$sw_time_bin_recent == 1])
  total_cases_sw_recent <- sum(data$hiv_rslt[data$sw_time_bin_recent == 1])
  incidence_rate_sw_recent <- ((total_cases_sw_recent + 0.5) / total_days_hiv_sw_recent) * 365.25 * 100

  # Calculate 95% confidence intervals for sex workers using exact Poisson method
  ci_sw <- poisson.test(total_cases_sw_recent, total_days_hiv_sw_recent / 365.25)
  ci_lower_sw_recent <- ci_sw$conf.int[1] * 100
  ci_upper_sw_recent <- ci_sw$conf.int[2] * 100

  cat("Incidence rate of HIV per 100 person years among sex workers (", group_label, "):", incidence_rate_sw_recent, "\n")
  cat("95% CI:", ci_lower_sw_recent, "-", ci_upper_sw_recent, "\n")
  cat("Number of cases among sex workers (", group_label, "):", total_cases_sw_recent, "\n")
  cat("Person years among sex workers (", group_label, "):", total_days_hiv_sw_recent / 365.25, "\n")

  # Calculate incidence rate for unexposed with additional 0.5 cases
  total_days_hiv_nosw_recent <- sum(data$days_risk[data$sw_time_bin_recent == 0])
  total_cases_nosw_recent <- sum(data$hiv_rslt[data$sw_time_bin_recent == 0])
  incidence_rate_nosw_recent <- ((total_cases_nosw_recent + 0.5) / total_days_hiv_nosw_recent) * 365.25 * 100

  # Calculate 95% confidence intervals for non-sex workers using exact Poisson method
  ci_nosw <- poisson.test(total_cases_nosw_recent, total_days_hiv_nosw_recent / 365.25)
  ci_lower_nosw_recent <- ci_nosw$conf.int[1] * 100
  ci_upper_nosw_recent <- ci_nosw$conf.int[2] * 100

  cat("Incidence rate of HIV per 100 person years among non-sex workers (", group_label, "):", incidence_rate_nosw_recent, "\n")
  cat("95% CI:", ci_lower_nosw_recent, "-", ci_upper_nosw_recent, "\n")
  cat("Number of cases among non-sex workers (", group_label, "):", total_cases_nosw_recent, "\n")
  cat("Person years among non-sex workers (", group_label, "):", total_days_hiv_nosw_recent / 365.25, "\n")

  # Calculate rate ratio and its 95% confidence interval
  rate_ratio_recent <- incidence_rate_sw_recent / incidence_rate_nosw_recent
  rate_ratio_recent_se <- sqrt((1 / (total_cases_sw_recent + 0.5)) + (1 / (total_cases_nosw_recent + 0.5)))
  ci_lower_rr_recent <- exp(log(rate_ratio_recent) - 1.96 * rate_ratio_recent_se)
  ci_upper_rr_recent <- exp(log(rate_ratio_recent) + 1.96 * rate_ratio_recent_se)

  cat("Rate ratio of HIV (sex workers vs non-sex workers) (", group_label, "):", rate_ratio_recent, "\n")
  cat("95% CI:", ci_lower_rr_recent, "-", ci_upper_rr_recent, "\n")
}

# Calculate for recent exposure
analysis_data_hiv_clean <- analysis_data_hiv_clean %>%
  mutate(sw_time_bin_recent = ifelse(is.na(rec_sw_sell_90d), 0, rec_sw_sell_90d))

# Filter data for females
analysis_data_hiv_women <- subset(analysis_data_hiv_clean, sex == 2)

# Calculate incidence rates and rate ratios for females
cat("### Analysis for Females ###\n")
calculate_manual_rate_ratio(analysis_data_hiv_women, "Recent Sex Work (Females)")

# Filter data for males
analysis_data_hiv_men <- subset(analysis_data_hiv_clean, sex == 1)

# Calculate incidence rates and rate ratios for males
cat("### Analysis for Males ###\n")
calculate_manual_rate_ratio(analysis_data_hiv_men, "Recent Sex Work (Males)")

# Calculate incidence rates and rate ratios for both
cat("### Analysis for Both ###\n")
calculate_manual_rate_ratio(analysis_data_hiv_clean, "Recent Sex Work (Both)")