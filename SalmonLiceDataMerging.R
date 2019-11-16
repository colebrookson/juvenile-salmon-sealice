library(tidyverse)
library(googlesheets)
library(data.table)


survey_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/survey_data.csv")
seine_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/seine_data.csv")
fish = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/fish_field_data.csv")
mainlice = read.csv("Hakai_lice_data_CB_edits.csv")
field = read.csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sealice_field.csv')

`%notin%` = negate(`%in%`)

#Process: get all the fish from 2019, only keep the ones from the same sites(? - didn't actually do this but might need to?), 
# keep only seines that are balanced collections

#keep only 2019 surveys
survey_data = survey_data %>% 
  filter(survey_date > as.Date("2019-01-01"))



#get the seines from the relevant surveys
seine_data = seine_data %>% 
  filter(survey_id %in% survey_data$survey_id)

seine_data = seine_data %>% 
  filter(so_taken > 4) %>% 
  filter(cu_taken > 4) %>% 
  filter(pi_taken > 4)

#only keep the fish with the right seine
fish = fish %>% 
  filter(seine_id %in% seine_data$seine_id)

#get the lice measurements for only those fish
field = field %>% 
  filter(ufn %in% fish$ufn)

#pull relevant info into a new table from the lice data and the relevant fish data
field = as.data.table(field)
fish = as.data.table(fish)
seine_data = as.data.table(seine_data)
survey_data = as.data.table(survey_data)

## pull fish data to the field data
setkey(field,ufn)
setkey(fish,ufn)
newlice = merge(field,fish, all.x=TRUE)
## pull the seiene data to the newlice data
setkey(newlice,seine_id)
setkey(seine_data,seine_id)
newlice = merge(newlice,seine_data, all.x = TRUE)
## pull the survey data to the newlice data
setkey(newlice,survey_id)
setkey(survey_data, survey_id)
newlice = merge(newlice, survey_data, all.x = TRUE)


#get rid of unnecessary columns
newlice_2019 = newlice[, c(1:3, 8:14, 28, 64, 65)]

#make new collection numbers for the 2019 data
newlice_2019 = newlice_2019 %>% 
  mutate(collection = paste(site_id, survey_date, sep="_")) 

collections = data.table(cbind(unique(newlice_2019$collection), c(53, 54, 55, 56, 57, 58, 59, 60))) %>% 
  rename(collection = V1)

newlice_2019$collection = as.factor(newlice_2019$collection)
newlice_2019 = as.data.table(newlice_2019)

setkey(newlice_2019, collection)
setkey(collections, collection)
newlice_2019 = merge(newlice_2019, collections, all.x = TRUE)

newlice_2019 = newlice_2019 %>% 
  select(-collection) %>% 
  rename(collection = V2)

#make all the columns I need
newlice_2019 = newlice_2019 %>% 
  rowwise() %>% 
  mutate(all.cal = sum(cal_mot_field, cgf_field), all.leps = sum(lpam_field, lpaf_field, laf_field, lam_field, lgf_field), all.lice = sum(all.cal, all.leps))
newlice_2019 = newlice_2019 %>% 
  mutate(site.region = paste(site_id[1]),
         year = survey_date)
newlice_2019$site.region<- substr(newlice_2019$site.region, 0, 1)
newlice_2019$year = format(as.Date(newlice_2019$survey_date, format="%m/%d/%Y"),"%Y")

newlice_2019 = newlice_2019 %>% 
  rename(spp = species)
mainlice = mainlice %>% 
  rename(site_id = site.id)

mainlice = mainlice %>% 
  select(year, spp, site.region, site_id, collection, ufn, all.cal, all.leps, all.lice)
newlice_2019_sub = newlice_2019 %>% 
  select(year, spp, site.region, site_id, collection, ufn, all.cal, all.leps, all.lice)

mainlice = rbind(mainlice, newlice_2019_sub)

#join up the stock data via UFN to get the stock ID into our dataset for sockeye

stock_data = read.csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/stock_id.csv')
library(stringr)
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

n_distinct(mainlice$)

