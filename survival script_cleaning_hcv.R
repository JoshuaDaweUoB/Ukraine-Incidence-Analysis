## load packages
pacman::p_load(dplyr, haven, tidyr, readxl, writexl, survival)

## set wd
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Ukraine data/Data")

## open data for visits 1-4
visit_1 <- read_excel("Ukraine data_analysis_bl.xlsx")
hcv_rslt_summary_v1 <- table(visit_1$hcv_rslt)
print(hcv_rslt_summary_v1) ## 1149 positive at baseline

visit_2 <- read_excel("Ukraine data_analysis_6m.xlsx")
visit_3 <- read_excel("Ukraine data_analysis_12m.xlsx")
visit_4 <- read_excel("Ukraine data_analysis_18m.xlsx")
# Note: data stored separately for each visit

# Combine the dataframes into a single longitudinal dataframe
analysis_data_hcv_long <- bind_rows(visit_1, visit_2, visit_3, visit_4) %>%
  arrange(ID, visit)

write_xlsx(analysis_data_hcv_long,"Ukraine data_analysis_hcv_long.xlsx")

# recode hcv_rslt to sensible binary numbers
analysis_data_hcv_long$hcv_rslt <- ifelse(analysis_data_hcv_long$hcv_rslt == "2", "1", 
                                          ifelse(analysis_data_hcv_long$hcv_rslt == "1", "0", 
                                                 analysis_data_hcv_long$hcv_rslt))

hcv_rslt_summary_long <- table(analysis_data_hcv_long$hcv_rslt)
print(hcv_rslt_summary_long) ## 3770 positives across entire dataset

# remove rows where hcv test result is missing
analysis_data_hcv_long <- analysis_data_hcv_long[!is.na(analysis_data_hcv_long$hcv_rslt), ]

hcv_rslt_summary_long <- table(analysis_data_hcv_long$hcv_rslt)
print(hcv_rslt_summary_long) ## still 3770 positives 

# check unique values of hcv_rslt
unique_values_before <- unique(analysis_data_hcv_long$hcv_rslt)
print("Unique values before recoding:")
print(unique_values_before)

hcv_rslt_summary_long <- table(analysis_data_hcv_long$hcv_rslt)
print(hcv_rslt_summary_long) ## still 3770 positives 

# find first negative hcv test per participant
analysis_data_hcv_long <- analysis_data_hcv_long %>%
  filter(hcv_rslt == 0) %>%
  group_by(ID) %>%
  summarise(visit_first_neg = min(visit, na.rm = TRUE)) %>%
  left_join(analysis_data_hcv_long, ., by = "ID")

# look at visit data
visit_data <- subset(analysis_data_hcv_long, select = c(ID, visit, visit_first_neg, hcv_rslt))

# drop participants who never have negative hcv test
analysis_data_hcv_long <- analysis_data_hcv_long[!is.na(analysis_data_hcv_long$visit_first_neg), ]

# remove positive tests that occurred before first negative test
visit_data <- visit_data %>%
  filter(visit >= visit_first_neg)

hcv_rslt_summary_visit <- table(visit_data$hcv_rslt)
print(hcv_rslt_summary_visit) ## now 657 positives

analysis_data_hcv_long <- analysis_data_hcv_long %>%
  filter(visit >= visit_first_neg)

# find first positive hcv test
visit_data <- visit_data %>%
  filter(hcv_rslt == 1) %>%
  group_by(ID) %>%
  summarise(visit_first_pos = min(visit, na.rm = TRUE)) %>%
  left_join(visit_data, ., by = "ID")

analysis_data_hcv_long <- analysis_data_hcv_long %>%
  filter(hcv_rslt == 1) %>%
  group_by(ID) %>%
  summarise(visit_first_pos = min(visit, na.rm = TRUE)) %>%
  left_join(analysis_data_hcv_long, ., by = "ID")

# remove tests that occurred after first positive test if visit_first_pos is not missing
visit_data <- visit_data %>%
  filter(is.na(visit_first_pos) | visit <= visit_first_pos)

hcv_rslt_summary_visit <- table(visit_data$hcv_rslt)
print(hcv_rslt_summary_visit) ## now 325 positives

analysis_data_hcv_long <- analysis_data_hcv_long %>%
  filter(is.na(visit_first_pos) | visit <= visit_first_pos)

## testing data

# keep columns of interest
testing_df <- subset(analysis_data_hcv_long, select = c(ID, visit, hcv_rslt)) 

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
         hcv_rslt = lead(hcv_rslt))  

# months at risk
testing_df <- testing_df %>%
  mutate(
    hcv_rslt = as.numeric(hcv_rslt),
    months_start = as.numeric(months_start),
    months_end = as.numeric(months_end)
  )

testing_df <- testing_df %>%
  mutate(months_risk = months_end - months_start - ifelse(hcv_rslt == 1, 3, 0))

testing_df <- testing_df[!is.na(testing_df$months_end), ]

write_xlsx(testing_df,"HCV_test_data.xlsx")

## exposure data

# keep columns of interest
exposure_df <- subset(analysis_data_hcv_long, select = c(ID, age, sex, city, rec_sw_sell_90d, rec_sw_buy_90d, rec_oat_6m, rec_incarc_6m, rec_homeless, lifetime_oat, lifetime_incarc, inj_dur, rec_inj_freq_month)) 

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
analysis_data_hcv_clean <- left_join(exposure_df, testing_df, by = c("ID", "id_seq"))

# Remove rows where rec_sw_sell_90d is NA
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  filter(!is.na(rec_sw_sell_90d))

# Create days_risk variable
analysis_data_hcv_clean <- analysis_data_hcv_clean %>%
  mutate(days_risk = months_risk * 30.4)

# convert numeric to factor variables
analysis_data_hcv_bl$sex <- factor(analysis_data_hcv_bl$sex)
analysis_data_hcv_bl$rec_homeless <- factor(analysis_data_hcv_bl$rec_homeless)
analysis_data_hcv_bl$rec_incarc_6m <- factor(analysis_data_hcv_bl$rec_incarc_6m)
analysis_data_hcv_bl$rec_oat_6m <- factor(analysis_data_hcv_bl$rec_oat_6m)
analysis_data_hcv_bl$rec_sw_sell_90d <- factor(analysis_data_hcv_bl$rec_sw_sell_90d)

hcv_rslt_summary <- table(analysis_data_hcv_long$hcv_rslt)
print(hcv_rslt_summary)

# save data
write_xlsx(analysis_data_hcv_clean,"HCV_data_clean.xlsx")

