## load packages
pacman::p_load(dplyr, haven, tidyr, readxl, writexl, survival)

## set wd
# setwd("C:/Users/joshua.dawe/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Tijuana data/Data")
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Ukraine data/Data")

## open data for visits 1-4
visit_1 <- read_excel("Ukraine data_analysis_bl.xlsx")
hiv_rslt_summary_v1 <- table(visit_1$hiv_rslt)
print(hiv_rslt_summary_v1) ## no positive at baseline, 2156 negs

visit_2 <- read_excel("Ukraine data_analysis_6m.xlsx")
hiv_rslt_summary_v2 <- table(visit_2$hiv_rslt)
print(hiv_rslt_summary_v2) ## 26 positives 1808 negs

visit_3 <- read_excel("Ukraine data_analysis_12m.xlsx")
hiv_rslt_summary_v3 <- table(visit_3$hiv_rslt)
print(hiv_rslt_summary_v3) ## 36 positives 1879 negs

visit_4 <- read_excel("Ukraine data_analysis_18m.xlsx")
hiv_rslt_summary_v4 <- table(visit_4$hiv_rslt)
print(hiv_rslt_summary_v4) ## 40 positives 1864 negs

# Note: data stored separately for each visit

# Combine the dataframes into a single longitudinal dataframe
analysis_data_hiv_long <- bind_rows(visit_1, visit_2, visit_3, visit_4) %>%
  arrange(ID, visit)

hiv_rslt_summary_long <- table(analysis_data_hiv_long$hiv_rslt)
print(hiv_rslt_summary_long) ## 102 positives 7707 negs

write_xlsx(analysis_data_hiv_long,"Ukraine data_analysis_hiv_long.xlsx")

# remove rows where hiv test result is missing
analysis_data_hiv_long <- analysis_data_hiv_long[!is.na(analysis_data_hiv_long$hiv_rslt), ]

hiv_rslt_summary_long <- table(analysis_data_hiv_long$hiv_rslt)
print(hiv_rslt_summary_long) ## still 102 positives and 7707 negs

# check unique values of hiv_rslt
unique_values_before <- unique(analysis_data_hiv_long$hiv_rslt)
print(unique_values_before)

# recode hiv_rslt to sensible binary numbers
analysis_data_hiv_long$hiv_rslt <- ifelse(analysis_data_hiv_long$hiv_rslt == "2", "1", 
                                          ifelse(analysis_data_hiv_long$hiv_rslt == "1", "0", 
                                                 analysis_data_hiv_long$hiv_rslt))

hiv_rslt_summary_long <- table(analysis_data_hiv_long$hiv_rslt)
print(hiv_rslt_summary_long)

# find first negative hiv test per participant
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  filter(hiv_rslt == 0) %>%
  group_by(ID) %>%
  summarise(visit_first_neg = min(visit, na.rm = TRUE)) %>%
  left_join(analysis_data_hiv_long, ., by = "ID")

# drop participants who never have negative hiv test
analysis_data_hiv_long <- analysis_data_hiv_long[!is.na(analysis_data_hiv_long$visit_first_neg), ]

hiv_rslt_summary_long <- table(analysis_data_hiv_long$hiv_rslt)
print(hiv_rslt_summary_long) ## no drops

# remove participants positive at baseline
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  filter(!(hiv_rslt == 1 & visit == 1)) %>%
  group_by(ID)
  
# find first positive hiv test
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  filter(hiv_rslt == 1) %>%
  group_by(ID) %>%
  summarise(visit_first_pos = min(visit, na.rm = TRUE)) %>%
  left_join(analysis_data_hiv_long, ., by = "ID")

# remove tests that occurred after first positive test if visit_first_pos is not missing
analysis_data_hiv_long <- analysis_data_hiv_long %>%
  filter(is.na(visit_first_pos) | visit <= visit_first_pos)

hiv_rslt_summary_visit <- table(analysis_data_hiv_long$hiv_rslt)
print(hiv_rslt_summary_visit) ## now 55 positives

## testing data

# keep columns of interest
testing_df <- subset(analysis_data_hiv_long, select = c(ID, visit, hiv_rslt)) 

testing_df <- testing_df %>%
  mutate(months_start = case_when(
    visit == "1" ~ "0",
    visit == "2" ~ "6",
    visit == "3" ~ "12",
    visit == "4" ~ "18",
    TRUE ~ NA_character_
  ))

# create lag of test (using lead function)
testing_df <- testing_df %>%
  arrange(ID, visit) %>%  
  group_by(ID) %>%
  mutate(months_end = lead(months_start),
         hiv_rslt = lead(hiv_rslt)) 


# months at risk
testing_df <- testing_df %>%
  mutate(
    hiv_rslt = as.numeric(hiv_rslt),
    months_start = as.numeric(months_start),
    months_end = as.numeric(months_end)
  )

testing_df <- testing_df %>%
  mutate(months_risk = months_end-months_start)

testing_df <- testing_df[!is.na(testing_df$months_end), ]

write_xlsx(testing_df,"HIV_test_data.xlsx")

## exposure data

# keep columns of interest
exposure_df <- subset(analysis_data_hiv_long, select = c(ID, age, sex, city, rec_sw_sell_90d, rec_sw_buy_90d, rec_oat_6m, rec_incarc_6m, rec_homeless, lifetime_oat, lifetime_incarc, inj_dur, rec_inj_freq_month)) 

exposure_df <- exposure_df %>%
  mutate(
    city = as.character(city)
  )

exposure_df <- exposure_df %>%
  arrange(ID) %>%
  group_by(ID) %>%
  slice(-n()) 

# recode exposures to sensible binary numbers
sw_sell <- table(exposure_df$rec_sw_sell_90d)
print(sw_sell)

exposure_df$rec_sw_sell_90d[exposure_df$rec_sw_sell_90d == "0"] <- NA
exposure_df$rec_sw_sell_90d[exposure_df$rec_sw_sell_90d == "2"] <- "0"
exposure_df$rec_sw_sell_90d[!(exposure_df$rec_sw_sell_90d %in% c("0", "1"))] <- NA

sw_sell <- table(exposure_df$rec_sw_sell_90d)
print(sw_sell)

sw_buy <- table(exposure_df$rec_sw_buy_90d)
print(sw_buy)

exposure_df$rec_sw_buy_90d[exposure_df$rec_sw_buy_90d == "0"] <- NA
exposure_df$rec_sw_buy_90d[exposure_df$rec_sw_buy_90d == "2"] <- "0"
exposure_df$rec_sw_buy_90d[!(exposure_df$rec_sw_buy_90d %in% c("0", "1"))] <- NA

sw_buy <- table(exposure_df$rec_sw_buy_90d)
print(sw_buy)

# create analysis df

# sequence by id for merge
exposure_df <- exposure_df %>%
  arrange(ID) %>%
  mutate(id_seq = row_number())

testing_df <- testing_df %>%
  arrange(ID) %>%
  mutate(id_seq = row_number())

# merge
analysis_data_hiv_clean <- left_join(exposure_df, testing_df, by = c("ID", "id_seq"))

# Remove rows where rec_sw_sell_90d is NA
analysis_data_hiv_clean <- analysis_data_hiv_clean %>%
  filter(!is.na(rec_sw_sell_90d))

# Create days_risk variable
analysis_data_hiv_clean <- analysis_data_hiv_clean %>%
  mutate(days_risk = months_risk * 30.4)

hiv_rslt_summary <- table(analysis_data_hiv_long$hiv_rslt)
print(hiv_rslt_summary)

# save data
write_xlsx(analysis_data_hiv_clean,"HIV_data_clean.xlsx")




