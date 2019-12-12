library(tidyverse)
library(googlesheets)
library(data.table)
library(stringr)

survey_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/survey_data.csv")
seine_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/seine_data.csv")
fish = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/fish_field_data.csv")
field = read_csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sealice_field.csv')
lab_fs = read_csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sealice_lab_fs.csv')
lab_mot = read_csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sealice_lab_mot.csv')
fish_lab = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/fish_lab_data.csv")

#mainlice = read_csv("Hakai_lice_data_CB_edits.csv")
#keep = mainlice
#write.csv(keep, 'Hakai_lice_CB_edits_analysis_from_first_draft_to_coauthors.csv') keeping this for posterity

`%notin%` = negate(`%in%`)

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
newlice_field = newlice_field %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% #unique place/time identifier for random effect
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice) %>% 
  rename(spp = species)
newlice_mot = newlice_mot %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% 
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice) %>% 
  rename(spp = species)
newlice_fs = newlice_fs %>% 
  mutate(year = format(as.Date(survey_date, format="%m/%d/%Y"),"%Y"), 
         site.region = substr(site_id, 0, 1),
         collection = paste(site_id, survey_date, sep="_")) %>% 
  dplyr::select(year, species, site.region, site_id, collection, ufn, 
                all.cal, all.leps, so_taken, cu_taken, pi_taken, all.lice) %>% 
  rename(spp = species)

#now keep only collections (seines) that have min. 5 of all 3 focal species
newlice_field = newlice_field %>% 
  filter(so_taken > 4) %>% 
  filter(cu_taken > 4) %>% 
  filter(pi_taken > 4)
newlice_fs = newlice_fs %>% 
  filter(so_taken > 4) %>% 
  filter(cu_taken > 4) %>% 
  filter(pi_taken > 4)
newlice_mot = newlice_mot %>% 
  filter(so_taken > 4) %>% 
  filter(cu_taken > 4) %>% 
  filter(pi_taken > 4)

#bind and remove duplicates
mainlice = rbind(newlice_field, newlice_fs, newlice_mot)
mainlice = mainlice[!duplicated(mainlice$ufn), ]
mainlice = mainlice %>% 
  dplyr::select(-pi_taken, -so_taken, -cu_taken)



##### comented code below will show how my initial dataset had some lice from all of the different licing protocols
##### but not all and I don't know why - therefore, when pulling directly from Hakai data, I have 3092 fish instead of 
##### ~2100 fish which is obviously a big difference. Hopefully results will stay consistent. Looking at the collections
##### previously I had around 60 collections, now I have 129, and there are more than double the sites than previously
##### so that's where the discrepancy is coming from. 

# t = rbind(newlice_field, newlice_fs, newlice_mot)
# n_distinct(t$ufn)
# x = t %>% 
#   filter(year != 2019)
# int = intersect(keep$ufn, x$ufn)
# int_mot = intersect(newlice_mot$ufn, int)
# int_field = intersect(newlice_field$ufn, int)
# int_fs = intersect(newlice_fs$ufn, int)


#join up the stock data via UFN to get the stock ID into our dataset for sockeye
stock_data = read.csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/stock_id.csv')
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

t = fish %>% 
  filter(ufn %in% mainlice$ufn)
summary(t$fork_length_field)
