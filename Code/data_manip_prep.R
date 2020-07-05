##### TITLE: Differential Infection of Parasitic Sea Lice on Juvenile Pacific Salmon in British Columbia, Canada
##### CREATOR: Cole B. Brookson
##### INITIALIZATION DATE: 2019-11-21

library(tidyverse)
library(cowplot)
library(bbmle)
library(lubridate)
library(rsample)
library(tibble)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(MuMIn)
library(here)
library(mapdata)
library(maps)
library(maptools)
library(ggmap)
library(ggrepel)
library(raster)
library(ggthemes)
library(ggsn)
library(rgeos)
library(rgdal)
library(data.table)
library(stringr)
library(doParallel)
library(here)

`%notin%` = negate(`%in%`) #creates inverse operator of the %in% operator for use in filtering data

#NOTE: the most up-to-date versions of this data can be found at the following URL:
#https://github.com/HakaiInstitute/jsp-data/tree/master/data
survey_data = read_csv(here("./data/survey_data.csv"))
seine_data = read_csv(here("./data/seine_data.csv"))
fish = read_csv(here('./data/fish_field_data.csv'), guess_max = 300000)
field = read_csv(here('./data/sealice_field.csv'))
lab_fs = read_csv(here('./data/sealice_lab_fs.csv'))
lab_mot = read_csv(here('./data/sealice_lab_mot.csv'))
fish_lab = read_csv(here('./data/fish_lab_data.csv'), guess_max = 300000)
stock_data = read_csv(here('./data/stock_id.csv'), guess_max = 300000)
salmonsites = read_csv('./data/spatial_data/site_coordinates.csv')

#first, only keep the survey's where more than 4 of each focal species were kept
# proper_surveys = seine_data %>% 
#   filter(so_taken > 4) %>% 
#   filter(pi_taken > 4) %>% 
#   filter(cu_taken > 4)
# proper_fish = fish %>% 
#   filter(seine_id %in% proper_seines$seine_id)

field = as.data.table(field)
fish = as.data.table(fish)
seine_data = as.data.table(seine_data)
survey_data = as.data.table(survey_data)
lab_fs = as.data.table(lab_fs)
lab_mot = as.data.table(lab_mot)

#motile lab
setkey(lab_mot,ufn)
setkey(fish,ufn)
newlice_mot = merge(lab_mot,fish, all.x=TRUE)
newlice_mot = left_join(newlice_mot, fish_lab, by = 'ufn'); newlice_mot = as.data.table(newlice_mot)
setkey(newlice_mot, seine_id)
setkey(seine_data, seine_id)
newlice_mot = merge(newlice_mot, seine_data, all.x = TRUE)
setkey(newlice_mot, survey_id)
setkey(survey_data, survey_id)
newlice_mot = merge(newlice_mot, survey_data, all.x = TRUE)

#lab fine scale
setkey(lab_fs,ufn)
setkey(fish,ufn)
newlice_fs = merge(lab_fs,fish, all.x=TRUE)
newlice_fs = left_join(newlice_fs, fish_lab, by = 'ufn'); newlice_fs = as.data.table(newlice_fs)
setkey(newlice_fs, seine_id)
setkey(seine_data, seine_id)
newlice_fs = merge(newlice_fs, seine_data, all.x = TRUE)
setkey(newlice_fs, survey_id)
setkey(survey_data, survey_id)
newlice_fs = merge(newlice_fs, survey_data, all.x = TRUE)

#lab fine scale
setkey(field,ufn)
setkey(fish,ufn)
newlice_field = merge(field,fish, all.x=TRUE)
newlice_field = left_join(newlice_field, fish_lab, by = 'ufn'); newlice_field = as.data.table(newlice_field)
setkey(newlice_field, seine_id)
setkey(seine_data, seine_id)
newlice_field = merge(newlice_field, seine_data, all.x = TRUE)
setkey(newlice_field, survey_id)
setkey(survey_data, survey_id)
newlice_field = merge(newlice_field, survey_data, all.x = TRUE)

#default to keeping the fine scale over the motile, motile over field
#no fish in motile that are in fine scale
#881 fish in motile that are also in field - exclude those fish from 'field' 
mot_and_field_ufn = intersect(lab_mot$ufn, field$ufn)
newlice_field = newlice_field %>% 
  filter(ufn %notin% mot_and_field_ufn)

#now get counts of total motiles, motile leps and motile cal
newlice_field = newlice_field %>% 
  rowwise() %>% 
  mutate(all.cal = sum(cal_mot_field, cgf_field), 
         all.leps = sum(lpam_field, lpaf_field, laf_field, lam_field, lgf_field), 
         all.lice = sum(all.cal, all.leps))
newlice_mot = newlice_mot %>% 
  rowwise() %>% 
  mutate(all.cal = sum(cpaf_lab, caf_lab, cgf_lab, cm_lab), 
         all.leps = sum(lpaf_lab,lpam_lab,lam_lab, laf_lab, lgf_lab), 
         all.lice = sum(all.cal, all.leps))
newlice_fs = newlice_fs %>% 
  rowwise() %>% 
  mutate(all.cal = sum(cal_grav_f, cal_a_f, cal_a_m, cal_pa_f, cal_pa_m), 
         all.leps = sum(lep_pa_m_1, lep_pa_m_2, lep_pa_f_1, lep_pa_f_2, lep_pa_unid, lep_a_m, lep_a_f, lep_grav_f), 
         all.lice = sum(all.cal, all.leps))
#now keep only the relevant columns and get rid of the rest
#with fork length
newlice_field_fork = newlice_field %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% #unique place/time identifier for random effect
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice, fork_length, seine_id) %>% 
  rename(spp = species) %>% 
  filter(!is.na(fork_length))
newlice_mot_fork = newlice_mot %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% 
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice, fork_length, seine_id) %>% 
  rename(spp = species) %>% 
  filter(!is.na(fork_length))
newlice_fs_fork = newlice_fs %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% 
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice, fork_length, seine_id) %>% 
  rename(spp = species) %>% 
  filter(!is.na(fork_length))

#bind and remove duplicates
mainlice_fork = rbind(newlice_field_fork, newlice_fs_fork, newlice_mot_fork)
mainlice_fork = mainlice_fork[!duplicated(mainlice_fork$ufn), ]
mainlice_fork = mainlice_fork %>% 
  dplyr::select(-pi_taken, -so_taken, -cu_taken)
mainlice_fork = mainlice_fork %>% 
  filter(spp %in% c('SO', 'PI', 'CU'))

#now keep only collections (seines) that have min. 5 of all 3 focal species
collections_all_fork = c(unique(mainlice_fork$collection))
collections_proper_fork = c()
collections_bad_fork = c()
suppressWarnings(for(i in collections_all_fork) {
  temp = mainlice_fork %>% 
    filter(collection == i) %>% 
    group_by(spp) %>% 
    summarize(n = n())
  if(nrow(temp) == 3) {
    if(temp$n[1] > 4 && temp$n[2] > 4 && temp$n[3] > 4) {
      collections_proper_fork = c(collections_proper_fork, i)
    } else{
      collections_bad_fork = c(collections_bad_fork, i)
    }
    
  } else {
    collections_bad_fork = c(collections_bad_fork, i)
  }
}
)
mainlice_fork = mainlice_fork %>% 
  filter(collection %in% collections_proper_fork)

#without fork length
newlice_field = newlice_field %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% #unique place/time identifier for random effect
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice, fork_length, seine_id) %>% 
  rename(spp = species) 
newlice_mot = newlice_mot %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% 
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice, fork_length, seine_id) %>% 
  rename(spp = species) 
newlice_fs = newlice_fs %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% 
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice, fork_length, seine_id) %>% 
  rename(spp = species) 

#bind and remove duplicates
mainlice = rbind(newlice_field, newlice_fs, newlice_mot)
mainlice = mainlice[!duplicated(mainlice$ufn), ]
mainlice = mainlice %>% 
  dplyr::select(-pi_taken, -so_taken, -cu_taken)
mainlice = mainlice %>% 
  filter(spp %in% c('SO', 'PI', 'CU'))

#now keep only collections (seines) that have min. 5 of all 3 focal species
collections_all = c(unique(mainlice$collection))
collections_proper = c()
collections_bad = c()
suppressWarnings(for(i in collections_all) {
  temp = mainlice %>% 
    filter(collection == i) %>% 
    group_by(spp) %>% 
    summarize(n = n())
  if(nrow(temp) == 3) {
    if(temp$n[1] > 4 && temp$n[2] > 4 && temp$n[3] > 4) {
      collections_proper = c(collections_proper, i)
    } else{
      collections_bad = c(collections_bad, i)
    }
    
  } else {
    collections_bad = c(collections_bad, i)
  }
}
)
mainlice = mainlice %>% 
  filter(collection %in% collections_proper)

#mainlice_nofork = mainlice

#join up the stock data via UFN to get the stock ID into our dataset for sockeye
#stock_data = read.csv('stock_id.csv')
main_stock = left_join(mainlice, stock_data, by = 'ufn')  

main_stock =  main_stock %>% 
  filter(spp == 'SO') %>% 
  filter(prob_1 !=is.na(prob_1)) %>% 
  mutate(stock_1 = tolower(stock_1)) %>% 
  group_by(stock_1) %>% 
  summarize(count = length(stock_1))
main_stock$stock_1 = str_to_title(main_stock$stock_1)


sum_stock = main_stock %>% 
  filter(count > 10)
sum(sum_stock$count)
sum(main_stock$count)

#make vars into factors
mainlice$year = as.factor(mainlice$year);mainlice$collection = as.factor(mainlice$collection)

#####do some throw away calculations for paper/reviewers

#check number of observations in a dataset if collections were filtered to 10 fish
collections_check_10 = c(unique(mainlice$collection))
collections_check_10_proper = c()
collections_check_10_bad = c()
suppressWarnings(for(i in collections_check_10) {
  temp = mainlice %>% 
    filter(collection == i) %>% 
    group_by(spp) %>% 
    summarize(n = n())
  if(nrow(temp) == 3) {
    if(temp$n[1] > 9 && temp$n[2] > 9 && temp$n[3] > 9) {
      collections_check_10_proper = c(collections_check_10_proper, i)
    } else{
      collections_bad = c(collections_bad, i)
    }
    
  } else {
    collections_check_10_bad = c(collections_check_10_bad, i)
  }
}
)
mainlice_check_10 = mainlice %>% 
  filter(collection %in% collections_check_10_proper)

#get mean and sd of lengths for the fish we do have length measurement for 
#check sizes of fish we do have measurements for
size_test = mainlice_fork %>% 
  group_by(spp) %>% 
  summarize(mean(fork_length))
size_sd_test = mainlice_fork %>% 
  group_by(spp) %>% 
  summarize(sd(fork_length))

