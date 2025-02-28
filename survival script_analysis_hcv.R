# load packages
pacman::p_load(arsenal, survival, readxl, dplyr)

## set wd
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Ukraine data/Data")

# load data
analysis_data_hcv_clean <- read_excel("HCV_data_clean.xlsx")

# baseline characteristics sex work
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  group_by(ID) %>%
  mutate(id_seq = row_number())

analysis_data_hcv_bl <- analysis_data_hcv_clean %>%
  group_by(ID) %>%
  mutate(id_seq = row_number(),
         rec_sw_sell_90d = ifelse(all(is.na(rec_sw_sell_90d)), NA_real_, max(rec_sw_sell_90d, na.rm = TRUE)),
         rec_incarc_6m = ifelse(all(is.na(rec_incarc_6m)), NA_real_, max(rec_incarc_6m, na.rm = TRUE)),
         rec_oat_6m = ifelse(all(is.na(rec_oat_6m)), NA_real_, max(rec_oat_6m, na.rm = TRUE)),
         rec_homeless = ifelse(all(is.na(rec_homeless)), NA_real_, max(rec_homeless, na.rm = TRUE))) %>%
  ungroup() %>%
  subset(id_seq == 1)

# table
tab_bl_sw_sex_hcv <- tableby(sex ~ age + rec_homeless + rec_incarc_6m + rec_oat_6m + rec_sw_sell_90d + inj_dur + rec_inj_freq_month, data=analysis_data_hcv_bl)
summary(tab_bl_sw_sex_hcv, text=TRUE)

## HCV analysis ##

## incidence rate calculations

# overall incidence rate
total_days_hcv <- sum(analysis_data_hcv_clean$days_risk)
total_cases <- sum(analysis_data_hcv_clean$hcv_rslt)
incidence_rate <- (total_cases / total_days_hcv) * 365.25 * 100

# Calculate 95% confidence intervals
incidence_rate_se <- sqrt(total_cases) / total_days_hcv * 365.25 * 100
ci_lower <- incidence_rate - 1.96 * incidence_rate_se
ci_upper <- incidence_rate + 1.96 * incidence_rate_se

cat("Incidence rate of HCV per 100 person years:", incidence_rate, "\n")
cat("95% CI:", ci_lower, "-", ci_upper, "\n")

### sex work HCV incidence rate calculations ###

# Function to calculate incidence rates and rate ratios
calculate_incidence_and_rate_ratio <- function(data, time_bin, group_label) {
  # selling sex work incidence rate
  total_days_hcv_sw <- sum(data$days_risk[data[[time_bin]] == 1])
  total_cases_sw <- sum(data$hcv_rslt[data[[time_bin]] == 1])
  incidence_rate_sw <- (total_cases_sw / total_days_hcv_sw) * 365.25 * 100

  # Calculate 95% confidence intervals for sex workers
  incidence_rate_sw_se <- sqrt(total_cases_sw) / total_days_hcv_sw * 365.25 * 100
  ci_lower_sw <- incidence_rate_sw - 1.96 * incidence_rate_sw_se
  ci_upper_sw <- incidence_rate_sw + 1.96 * incidence_rate_sw_se

  cat("Incidence rate of HCV per 100 person years among sex workers (", group_label, "):", incidence_rate_sw, "\n")
  cat("95% CI:", ci_lower_sw, "-", ci_upper_sw, "\n")

  # no sex work incidence rate
  total_days_hcv_nosw <- sum(data$days_risk[data[[time_bin]] == 0])
  total_cases_nosw <- sum(data$hcv_rslt[data[[time_bin]] == 0])
  incidence_rate_nosw <- (total_cases_nosw / total_days_hcv_nosw) * 365.25 * 100

  # Calculate 95% confidence intervals for non-sex workers
  incidence_rate_nosw_se <- sqrt(total_cases_nosw) / total_days_hcv_nosw * 365.25 * 100
  ci_lower_nosw <- incidence_rate_nosw - 1.96 * incidence_rate_nosw_se
  ci_upper_nosw <- incidence_rate_nosw + 1.96 * incidence_rate_nosw_se

  cat("Incidence rate of HCV per 100 person years among non-sex workers (", group_label, "):", incidence_rate_nosw, "\n")
  cat("95% CI:", ci_lower_nosw, "-", ci_upper_nosw, "\n")

  # Calculate rate ratio and its 95% confidence interval
  rate_ratio <- incidence_rate_sw / incidence_rate_nosw
  rate_ratio_se <- sqrt((1 / total_cases_sw) + (1 / total_cases_nosw))
  ci_lower_rr <- exp(log(rate_ratio) - 1.96 * rate_ratio_se)
  ci_upper_rr <- exp(log(rate_ratio) + 1.96 * rate_ratio_se)

  cat("Rate ratio of HCV (sex workers vs non-sex workers) (", group_label, "):", rate_ratio, "\n")
  cat("95% CI:", ci_lower_rr, "-", ci_upper_rr, "\n")

  # Create a summary dataset for Poisson regression
  summary_data <- data %>%
    group_by(!!sym(time_bin), rec_oat_6m, rec_incarc_6m, rec_inj_freq_month, inj_dur) %>%
    summarise(
      total_cases = sum(hcv_rslt),
      total_days = sum(days_risk)
    ) %>%
    mutate(
      rate = total_cases / total_days * 365.25 * 100
    )

  # Fit Poisson regression model controlling for rec_oat_6m and rec_incarc_6m
  poisson_model1 <- glm(total_cases ~ get(time_bin) + rec_oat_6m + rec_incarc_6m + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio1 <- exp(coef(poisson_model1)[2])
  ci1 <- exp(confint(poisson_model1)[2, ])

  cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for rec_oat_6m and rec_incarc_6m (", group_label, "):", rate_ratio1, "\n")
  cat("95% CI:", ci1[1], "-", ci1[2], "\n")

  # Fit Poisson regression model controlling for rec_oat_6m, rec_incarc_6m, rec_inj_freq_month, and inj_dur
  poisson_model2 <- glm(total_cases ~ get(time_bin) + rec_oat_6m + rec_incarc_6m + rec_inj_freq_month + inj_dur + offset(log(total_days)), 
                        family = poisson(link = "log"), 
                        data = summary_data)

  # Extract rate ratio and confidence intervals
  rate_ratio2 <- exp(coef(poisson_model2)[2])
  ci2 <- exp(confint(poisson_model2)[2, ])

  cat("Rate ratio of HCV (sex workers vs non-sex workers) controlling for rec_oat_6m, rec_incarc_6m, rec_inj_freq_month, and inj_dur (", group_label, "):", rate_ratio2, "\n")
  cat("95% CI:", ci2[1], "-", ci2[2], "\n")
}

# Calculate for recent exposure
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(sw_time_bin_recent = ifelse(is.na(rec_sw_sell_90d), 0, rec_sw_sell_90d))

# Filter data for females
analysis_data_hcv_women <- subset(analysis_data_hcv_clean, sex == 2)

# Calculate incidence rates and rate ratios for females
cat("### Analysis for Females ###\n")
calculate_incidence_and_rate_ratio(analysis_data_hcv_women, "sw_time_bin_recent", "Recent Sex Work (Females)")

# Filter data for males
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, sex == 1)

# Manually calculate rate ratio for recent sex work with additional 0.5 cases
# Calculate incidence rate for recent sex work exposure with additional 0.5 cases
total_days_hcv_sw_recent <- sum(analysis_data_hcv_men$days_risk[analysis_data_hcv_men$sw_time_bin_recent == 1])
total_cases_sw_recent <- sum(analysis_data_hcv_men$hcv_rslt[analysis_data_hcv_men$sw_time_bin_recent == 1]) + 0.5
incidence_rate_sw_recent <- (total_cases_sw_recent / total_days_hcv_sw_recent) * 365.25 * 100

# Calculate 95% confidence intervals for sex workers
incidence_rate_sw_recent_se <- sqrt(total_cases_sw_recent) / total_days_hcv_sw_recent * 365.25 * 100
ci_lower_sw_recent <- incidence_rate_sw_recent - 1.96 * incidence_rate_sw_recent_se
ci_upper_sw_recent <- incidence_rate_sw_recent + 1.96 * incidence_rate_sw_recent_se

cat("Incidence rate of HCV per 100 person years among sex workers (recent exposure):", incidence_rate_sw_recent, "\n")
cat("95% CI:", ci_lower_sw_recent, "-", ci_upper_sw_recent, "\n")

# Calculate incidence rate for unexposed males with additional 0.5 cases
total_days_hcv_nosw_recent <- sum(analysis_data_hcv_men$days_risk[analysis_data_hcv_men$sw_time_bin_recent == 0])
total_cases_nosw_recent <- sum(analysis_data_hcv_men$hcv_rslt[analysis_data_hcv_men$sw_time_bin_recent == 0]) + 0.5
incidence_rate_nosw_recent <- (total_cases_nosw_recent / total_days_hcv_nosw_recent) * 365.25 * 100

# Calculate 95% confidence intervals for non-sex workers
incidence_rate_nosw_recent_se <- sqrt(total_cases_nosw_recent) / total_days_hcv_nosw_recent * 365.25 * 100
ci_lower_nosw_recent <- incidence_rate_nosw_recent - 1.96 * incidence_rate_nosw_recent_se
ci_upper_nosw_recent <- incidence_rate_nosw_recent + 1.96 * incidence_rate_nosw_recent_se

cat("Incidence rate of HCV per 100 person years among non-sex workers (recent exposure):", incidence_rate_nosw_recent, "\n")
cat("95% CI:", ci_lower_nosw_recent, "-", ci_upper_nosw_recent, "\n")

# Calculate rate ratio and its 95% confidence interval
rate_ratio_recent <- incidence_rate_sw_recent / incidence_rate_nosw_recent
rate_ratio_recent_se <- sqrt((1 / total_cases_sw_recent) + (1 / total_cases_nosw_recent))
ci_lower_rr_recent <- exp(log(rate_ratio_recent) - 1.96 * rate_ratio_recent_se)
ci_upper_rr_recent <- exp(log(rate_ratio_recent) + 1.96 * rate_ratio_recent_se)

cat("Rate ratio of HCV (sex workers vs non-sex workers) (recent exposure):", rate_ratio_recent, "\n")
cat("95% CI:", ci_lower_rr_recent, "-", ci_upper_rr_recent, "\n")