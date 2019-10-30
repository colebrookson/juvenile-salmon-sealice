library(tidyverse)
library(googlesheets)
library(tidyverse)
library(here)
library(hakaiApi)

survey_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/survey_data.csv")
seine_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/seine_data.csv")
fish = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/fish_field_data.csv")
mainlice = read.csv("Hakai_lice_data_CB_edits.csv")
newlice = read.csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sealice_field.csv')

`%notin%` = negate(`%in%`)


survey_data = survey_data %>% 
  filter(survey_date > as.Date("2019-01-01")) %>% 
  drop_na(ctd_bout) %>% 
  mutate(hakai_id = paste(site_id, survey_date, ctd_bout, sep="_")) %>% 
  select(survey_id, survey_date, site_id, hakai_id, secchi) %>% 
  # Finding most mismatches are due to the wrong sampling bout being entered into the portal, so create a more general linkage of date_site
  mutate(date_site = paste(survey_date, site_id, sep="_")) %>% 
  filter(site_id %in% mainlice$site.id)

seine_data = seine_data %>% 
  filter(so_taken > 4 & cu_taken > 4 & pi_taken > 4) %>% 
  filter(survey_id %in% survey_data$survey_id)

survey_data = survey_data %>% 
  filter(survey_id %in% seine_data$survey_id)

fish = fish %>% 
  filter(seine_id %in% seine_data$seine_id)


newlice_2019 = newlice %>% 
  filter(ufn %notin% mainlice$ufn)
