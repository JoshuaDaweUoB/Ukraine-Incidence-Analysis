# load packages
pacman::p_load(dplyr, arsenal, survival)

## HIV analysis ##

## incidence rate calculations

# overall incidence rate
total_months_hiv <- sum(analysis_data_hiv_clean$months_risk)
total_cases <- sum(analysis_data_hiv_clean$hiv_rslt)
incidence_rate <- (total_cases / total_months_hiv) * 12 *100

cat("Incidence rate of HCV per 100 person years:", incidence_rate)

# selling sex work incidence rate
analysis_data_hiv_clean$sw_time_bin <- analysis_data_hiv_clean$rec_sw_sell_90d
analysis_data_hiv_clean$sw_time_bin[is.na(analysis_data_hiv_clean$sw_time_bin)] <- 0
total_months_hiv_sw <- sum(analysis_data_hiv_clean$months_risk[analysis_data_hiv_clean$sw_time_bin == 1])
total_cases_sw <- sum(analysis_data_hiv_clean$hiv_rslt[analysis_data_hiv_clean$sw_time_bin == 1])
incidence_rate_sw <- (total_cases_sw / total_months_hiv_sw) * 12 *100

cat("Incidence rate of HIV per 100 person years among sex workers:", incidence_rate_sw)

# no sex work incidence rate
analysis_data_hiv_clean$sw_time_bin <- analysis_data_hiv_clean$rec_sw_sell_90d
analysis_data_hiv_clean$sw_time_bin <- ifelse(is.na(analysis_data_hiv_clean$rec_sw_sell_90d), 1, analysis_data_hiv_clean$rec_sw_sell_90d)
total_months_hiv_nosw <- sum(analysis_data_hiv_clean$months_risk[analysis_data_hiv_clean$sw_time_bin == 0])
total_cases_nosw <- sum(analysis_data_hiv_clean$hiv_rslt[analysis_data_hiv_clean$sw_time_bin == 0])
incidence_rate_nosw <- (total_cases_nosw / total_months_hiv_nosw) * 12 *100

cat("Incidence rate of HIV per 100 person years among non-sex workers:", incidence_rate_nosw)

## missing sex work rate
analysis_data_hiv_clean$sw_time_bin <- analysis_data_hiv_clean$rec_sw_sell_90d
analysis_data_hiv_clean$sw_time_bin <- ifelse(is.na(analysis_data_hiv_clean$rec_sw_sell_90d), NA, analysis_data_hiv_clean$rec_sw_sell_90d)
total_months_hiv_miss <- sum(analysis_data_hiv_clean$months_risk[is.na(analysis_data_hiv_clean$sw_time_bin)])
total_cases_miss <- sum(analysis_data_hiv_clean$hiv_rslt[is.na(analysis_data_hiv_clean$sw_time_bin)])
incidence_rate_miss <- (total_cases_miss / total_months_hiv_miss) * 12 *100

cat("Incidence rate of HIV per 100 person years among non-sex workers:", incidence_rate_miss)

## HR/IRR calculcations

# incidence rate ratio
irr_sw <- (incidence_rate_sw / incidence_rate_nosw)
print(irr_sw)

# buying sex work incidence rate
analysis_data_hcv_clean$sw_time_bin <- analysis_data_hcv_clean$rec_sw_buy_90d
analysis_data_hcv_clean$sw_time_bin[is.na(analysis_data_hcv_clean$sw_time_bin)] <- 0
total_months_hcv_sw <- sum(analysis_data_hcv_clean$months_risk[analysis_data_hcv_clean$sw_time_bin == 1])
total_cases_sw <- sum(analysis_data_hcv_clean$hcv_rslt[analysis_data_hcv_clean$sw_time_bin == 1])
incidence_rate_sw <- (total_cases_sw / total_months_hcv_sw) * 12 *100

cat("Incidence rate of HCV per 100 person years among sex workers:", incidence_rate_sw)

# no sex work incidence rate
analysis_data_hcv_clean$sw_time_bin <- analysis_data_hcv_clean$rec_sw_buy_90d
analysis_data_hcv_clean$sw_time_bin <- ifelse(is.na(analysis_data_hcv_clean$rec_sw_buy_90d), 1, analysis_data_hcv_clean$rec_sw_buy_90d)
total_months_hcv_nosw <- sum(analysis_data_hcv_clean$months_risk[analysis_data_hcv_clean$sw_time_bin == 0])
total_cases_nosw <- sum(analysis_data_hcv_clean$hcv_rslt[analysis_data_hcv_clean$sw_time_bin == 0])
incidence_rate_nosw <- (total_cases_nosw / total_months_hcv_nosw) * 12 *100

cat("Incidence rate of HCV per 100 person years among non-sex workers:", incidence_rate_nosw)

## HR/IRR calculcations

# incidence rate ratio
irr_sw <- (incidence_rate_sw / incidence_rate_nosw)
print(irr_sw)

#### hazard ratio selling sex work ####

# unadjusted hazard ratio selling sw HCV recent - men and women
sw_model_hcv_crude = coxph(
  Surv(time = months_start, time2 = months_end, event = hcv_rslt) ~ rec_sw_sell_90d + city + age + sex, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_crude)

# unadjusted hazard ratio selling sw HCV recent - men
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, sex == 1)

sw_model_hcv_crude_men = coxph(
  Surv(time = months_start, time2 = months_end, event = hcv_rslt) ~ rec_sw_sell_90d, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_crude_men)

# unadjusted hazard ratio selling sw HCV recent - women
analysis_data_hcv_women <- subset(analysis_data_hcv_clean, sex == 2)

sw_model_hcv_crude_women = coxph(
  Surv(time = months_start, time2 = months_end, event = hcv_rslt) ~ rec_sw_sell_90d, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_crude_women)

#### hazard ratio buying sex work ####


#### hazard ratio buying sex work ####

# unadjusted hazard ratio selling sw HCV recent - men and women
sw_model_hcv_crude = coxph(
  Surv(time = months_start, time2 = months_end, event = hcv_rslt) ~ rec_sw_buy_90d, 
  data = analysis_data_hcv_clean
)

summary(sw_model_hcv_crude)

# unadjusted hazard ratio selling sw HCV recent - men
analysis_data_hcv_men <- subset(analysis_data_hcv_clean, sex == 1)

sw_model_hcv_crude_men = coxph(
  Surv(time = months_start, time2 = months_end, event = hcv_rslt) ~ rec_sw_buy_90d, 
  data = analysis_data_hcv_men
)

summary(sw_model_hcv_crude_men)

# unadjusted hazard ratio selling sw HCV recent - women
analysis_data_hcv_women <- subset(analysis_data_hcv_clean, sex == 2)

sw_model_hcv_crude_women = coxph(
  Surv(time = months_start, time2 = months_end, event = hcv_rslt) ~ rec_sw_buy_90d, 
  data = analysis_data_hcv_women
)

summary(sw_model_hcv_crude_women)
