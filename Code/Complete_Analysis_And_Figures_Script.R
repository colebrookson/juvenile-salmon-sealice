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
library(googlesheets)
library(stringr)
library(doParallel)

`%notin%` = negate(`%in%`)

#NOTE: the most up-to-date versions of this data can be found at the following URL:
#https://github.com/HakaiInstitute/jsp-data/tree/master/data

survey_data = survey_data
seine_data = seine_data
fish = fish_field_data
field = fish_field_data
lab_fs = sealice_lab_fs
lab_mot = sealice_lab_mot
fish_lab = fish_lab_data

survey_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/survey_data.csv")
seine_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/seine_data.csv")
fish = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/fish_field_data.csv")
field = read_csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sealice_field.csv')
lab_fs = read_csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sealice_lab_fs.csv')
lab_mot = read_csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sealice_lab_mot.csv')
fish_lab = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/fish_lab_data.csv")

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
stock_data = read.csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/sample_results/stock_id.csv')
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
mainlice$year <- as.factor(mainlice$year);mainlice$collection <- as.factor(mainlice$collection)

## Figure 1: Map
salmonsites <- read_csv('C:/Users/brookson/Documents/GitHub/jsp-data/data/site_coordinates.csv')

unique(mainlice$site_id)
salmonsites = salmonsites %>% 
  filter(site_id %in% unique(mainlice$site_id)) %>% 
  dplyr::select(site_id, site_lat, site_long)
salmonsites = salmonsites[!duplicated(salmonsites$site_id), ]

BC_shp = readOGR('C:/Users/brookson/Documents/GitHub/SalmonWork/SpatialData/COAST_TEST2.shp')
coords <- data.frame(cbind(salmonsites$site_lat,salmonsites$site_long))
colnames(coords) <- c('lat', 'long')

sort(unique(ggplot2::map_data("world")$region))
westcoast <- map_data('canada')

#make two themes, one for the inset, one for the big map
fte_theme_map_small <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_rect(colour = 'black')) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 12, color = color.axis.text, angle = 90)) + 
    theme(axis.text.y = element_text(size = 12, color = color.axis.text)) + 
    theme(axis.title.x = element_text(size = 14, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 14, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.ticks = element_line(colour = 'black')) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) 
}
fte_theme_map_big <- function(){
  color.background = 'grey75'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill = 'white', color = 'white')) +
    theme(plot.background = element_rect(fill=color.background,color = color.background)) +
    theme(panel.border = element_rect(colour = 'black')) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 15, vjust = 1.25)) +
    theme(axis.text.x = element_blank()) + 
    theme(axis.text.y = element_blank()) + 
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_blank()) +
    theme(plot.title = element_blank()) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) 
}
#create the two maps
states    <- c('Washington')
provinces <- c("British Columbia")

us <- getData("GADM",country="USA",level=1)
canada <- getData("GADM",country="CAN",level=1)

us.states <- us[us$NAME_1 %in% states,]
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

biggermap = ggplot()+
  geom_polygon(data = us.states,aes(x=long,y=lat,group=group), colour = 'black', size = 0.01, fill = 'grey75')+
  geom_polygon(data=ca.provinces, aes(x=long,y=lat,group=group), colour = 'black', size = 0.01, fill = 'grey75')+
  coord_cartesian(xlim = c(-128.5,-119.5), ylim = c(48.1,51.0)) +
  fte_theme_map_big()+
  annotate("rect", xmin = -127.3, xmax = -125, ymin = 49.8, ymax = 50.9, alpha = .7)+
  annotate('text', x = -121, y = 50.7, label = 'British Columbia', size = 4)+
  annotate('text', x = -120.7, y = 48.5, label = 'Washington', size = 4)+
  annotate('text', x = -126.7, y = 48.5, label = 'Vancouver Island', size = 4)+
  annotate('segment',x=-126.7, y=48.6, xend=-125.5, yend=49.5, arrow=arrow(length = unit(0.16, "npc")), 
           alpha = 0.8, size=1.2, color="black")
biggermap

smallermap = ggplot()+
  geom_polygon(data=ca.provinces,aes(x=long,y=lat,group=group), colour = 'black', size = 0.01, fill = 'grey75')+
  coord_cartesian(xlim = c(-127.0,-125), ylim = c(50,50.75))+
  fte_theme_map_small()+
  labs(x = 'Longitude', y = 'Latitude')+
  annotate("rect", xmin = -125.51, xmax = -125.05, ymin = 50.1, ymax = 50.5, alpha = .65)+
  annotate("rect", xmin = -126.9, xmax = -126.55, ymin = 50.45, ymax = 50.7, alpha = .65)+
  annotate('text', x= -125.8, y = 50.33, label = 'Discovery Islands', size = 4)+
  annotate('text', x = -126.6, y = 50.38, label = 'Johnstone Strait', size = 4)+
  geom_point(data = coords, aes(long,lat), color = 'black', size = 4, shape = 21, fill = 'black')+
  annotate('segment',x=-125.8, y=50.302, xend=-125.565, yend=50.27, arrow=arrow(length = unit(0.04, "npc")), 
           alpha = 0.8, size=1.1, color="black")+
  annotate('segment',x=-126.6, y=50.4, xend=-126.65, yend=50.44, arrow=arrow(length = unit(0.04, "npc")), 
           alpha = 0.8, size=1.1, color="black")
smallermap

#make the one plot inset with the other
insetmap = ggdraw()+
  draw_plot(smallermap) + 
  draw_plot(biggermap, x=0.105, y=0.161, width=0.5, height=0.3) 
insetmap
ggsave('study_map.png', plot = insetmap,
       width = 8, height = 7.5,
       dpi = 300)

## Figure 2: 3x2 of sal and lice species
stde <- function(x) sd(x)/sqrt(length(x))
fig_2_data = mainlice %>% 
  group_by(year, spp) %>% 
  summarize(mean_cal = mean(all.cal), mean_lep = mean(all.leps), 
            low_cal = (mean(all.cal) - stde(all.cal)),
            up_cal = (mean(all.cal) + stde(all.cal)),
            low_lep = (mean(all.leps) - stde(all.leps)),
            up_lep = (mean(all.leps) + stde(all.leps)))
fig_2_allyear_data = mainlice %>% 
  group_by(spp) %>% 
  summarize(mean_cal = mean(all.cal), mean_lep = mean(all.leps),
            low_cal = (mean(all.cal) - stde(all.cal)),
            up_cal = (mean(all.cal) + stde(all.cal)),
            low_lep = (mean(all.leps) - stde(all.leps)),
            up_lep = (mean(all.leps) + stde(all.leps))) %>% 
  mutate(year = 'All')
fig_2_allyear_data = fig_2_allyear_data[,c(8,1,2,3,4,5,6,7)]
fig_2_data = rbind(data.frame(fig_2_data), data.frame(fig_2_allyear_data))
lep_y_title = expression(paste("Mean Motile ", italic("L. salmonis"), " Abundance"))
cal_y_title = expression(paste("Mean Motile ", italic("C. clemensi"), " Abundance"))

#make each plot then stitch them together
fig_2_1 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'PI')), aes(x = year, y = mean_lep), 
           colour = 'black', fill = '#ff9999', width = 0.7) +
  labs(x = 'Year', y = lep_y_title, title = 'Pink') +
  geom_errorbar(data = (fig_2_data %>% filter(spp == 'PI')),
                aes(x = year, ymin=low_lep, ymax = up_lep,width = 0), size = 0.78 ,colour = 'Black')+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 0.5))
fig_2_2 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'CU')), aes(x = year, y = mean_lep), 
           colour = 'black', fill = '#59AE7F', width = 0.7) +
  geom_errorbar(data = (fig_2_data %>% filter(spp == 'CU')),
                aes(x = year, ymin=low_lep, ymax = up_lep,width = 0), size = 0.78 ,colour = 'Black')+
  labs(x = 'Year', y = lep_y_title, title = 'Chum') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0.0, 0.5))
fig_2_3 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'SO')), aes(x = year, y = mean_lep), 
           colour = 'black', fill = '#23359d', width = 0.7) +
  geom_errorbar(data = (fig_2_data %>% filter(spp == 'SO')),
                aes(x = year, ymin=low_lep, ymax = up_lep,width = 0), size = 0.78 ,colour = 'Black')+
  labs(x = 'Year', y = lep_y_title, title = 'Sockeye') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0.0, 0.5))
fig_2_4 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'PI')), aes(x = year, y = mean_cal), 
           colour = 'black', fill = '#ff9999', width = 0.7) +
  geom_errorbar(data = (fig_2_data %>% filter(spp == 'PI')),
                aes(x = year, ymin=low_cal, ymax = up_cal,width = 0), size = 0.78 ,colour = 'Black')+
  labs(x = 'Year', y = cal_y_title, title = 'Pink') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1))
fig_2_5 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'CU')), aes(x = year, y = mean_cal), 
           colour = 'black', fill = '#59AE7F', width = 0.7) +
  geom_errorbar(data = (fig_2_data %>% filter(spp == 'CU')),
                aes(x = year, ymin=low_cal, ymax = up_cal,width = 0), size = 0.78 ,colour = 'Black')+
  labs(x = 'Year', y = cal_y_title, title = 'Chum') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1)) 
fig_2_6 = ggplot() +
  geom_col(data = (fig_2_data %>% filter(spp == 'SO')), aes(x = year, y = mean_cal), 
           colour = 'black', fill = '#23359d', width = 0.7) +
  geom_errorbar(data = (fig_2_data %>% filter(spp == 'SO')),
                aes(x = year, ymin=low_cal, ymax = up_cal,width = 0), size = 0.78 ,colour = 'Black')+
  labs(x = 'Year', y = cal_y_title, title = 'Sockeye') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1)) 

fig_2 = plot_grid(fig_2_1, fig_2_2, fig_2_3,
                  fig_2_4, fig_2_5, fig_2_6, nrow = 2, rel_heights = c(0.85, 1), rel_widths = c(1, 0.8, 0.8))
ggsave('lice_per_fish_sp.png', plot = fig_2,
       width = 8, height = 8.0,
       dpi = 300)

## Figure 3: average number of lice per fish (fish sp and by year)
fig_3_data = mainlice %>% 
  group_by(year, spp) %>% 
  summarize(mean_lice = mean(all.lice))

fig_3_1 = ggplot() +
  geom_col(data = (fig_3_data %>% filter(spp == 'PI')), aes(x = year, y = mean_lice), 
           colour = 'black', fill = '#ff9999', width = 0.7) +
  labs(x = 'Year', y = 'Mean Number of Motile Lice Per Fish', title = 'Pink') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1.15))
fig_3_2 = ggplot() +
  geom_col(data = (fig_3_data %>% filter(spp == 'CU')), aes(x = year, y = mean_lice), 
           colour = 'black', fill = '#59AE7F', width = 0.7) +
  labs(x = 'Year', y = cal_y_title, title = 'Chum') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1.15)) 
fig_3_3 = ggplot() +
  geom_col(data = (fig_3_data %>% filter(spp == 'SO')), aes(x = year, y = mean_lice), 
           colour = 'black', fill = '#23359d', width = 0.7) +
  labs(x = 'Year', y = cal_y_title, title = 'Sockeye') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.0, 1.15))

fig_3 = plot_grid(fig_3_1, fig_3_2, fig_3_3, nrow = 1,
                  rel_widths = c(1, 0.8, 0.8))
ggsave('lice_per_fish.png', plot = fig_3,
       scale = 1, width = 8, height = 6.7,
       dpi = 300)

#get the necessary model selection things:
lepmod.crossedp <- glmmTMB(all.leps ~ spp * site.region + spp * year + 
                            site.region * year + (1 | collection),
                            ziformula = ~1,
                          data = mainlice, family=poisson)
lepmod.crossednb <- glmmTMB(all.leps ~ spp * site.region + spp * year + 
                            site.region * year + (1 | collection),
                          ziformula = ~1,
                          data = mainlice, family=nbinom2)
AICtab(lepmod.crossedp, lepmod.crossednb)
calmod.crossednb <- glmmTMB(all.cal ~ spp * site.region + spp * year + 
                            site.region * year + (1 | collection), 
                          data = mainlice, family=nbinom2)
calmod.crossedp <- glmmTMB(all.cal ~ spp * site.region + spp * year + 
                            site.region * year + (1 | collection), 
                          data = mainlice, family=poisson)
AICtab(calmod.crossedp, calmod.crossednb)

#check models
lep_p_sr <- simulateResiduals(lepmod.crossedp)
cal_nb_sr <- simulateResiduals(calmod.crossednb)
plot(lep_p_sr)
plot(cal_nb_sr)
#NOTE: despite the fact that the nb one fit better (dAIC = 2.2), not all of the submodels converged with the nb, 
# so using the poisson one instead so they actually fit 
lepmod.crossed_dredge = MuMIn::dredge(lepmod.crossedp, subset = (`cond(site.region)` && `cond(year)`))
calmod.crossed_dredge = MuMIn::dredge(calmod.crossednb, subset = (`cond(site.region)` && `cond(year)`))
 
#try models with forklength included -- leps model won't converge and cal model isn't any better 
lepmod_forkp = glmmTMB(all.leps ~ spp * site.region + spp * year + 
                        site.region * year + fork_length*year + fork_length*site.region +
                        (1 | collection), 
                      data = mainlice_fork, family=poisson)
lepmod_forknb = glmmTMB(all.leps ~ spp * site.region + spp * year + 
                         site.region * year + fork_length*year + fork_length*site.region +
                         (1 | collection), 
                       data = mainlice_fork, family=nbinom2)
#neither distribution fits
calmod_fork = glmmTMB(all.cal ~ spp * site.region + spp * year + 
                        site.region * year + fork_length*spp + fork_length*site.region +
                        (1 | collection), 
                      data = mainlice_fork, family=nbinom2)
calmod_nofork = glmmTMB(all.cal ~ spp * site.region + spp * year + 
                        site.region * year + (1 | collection), 
                      data = mainlice_fork, family=nbinom2)
AICtab(calmod_nofork, calmod_fork)

#so the goal here is to bootstrap the data (parametric) by resampling the hierarchical levels, then 
#run the model averaging process with the new data and use that to get our model-averaged CI's. Essentially, 
#we're wrapping our esitmator of mu in a function and bootstrapping it

bootlice =  mainlice %>% 
  dplyr::select(all.cal, all.leps, spp, site.region, collection, year, ufn)

sock2015D <- bootlice %>% filter(spp == 'SO' & year == '2015' & site.region == 'D')
sock2015J <- bootlice %>% filter(spp == 'SO' & year == '2015' & site.region == 'J')
sock2016D <- bootlice %>% filter(spp == 'SO' & year == '2016' & site.region == 'D')
sock2016J <- bootlice %>% filter(spp == 'SO' & year == '2016' & site.region == 'J')
sock2017D <- bootlice %>% filter(spp == 'SO' & year == '2017' & site.region == 'D')
sock2017J <- bootlice %>% filter(spp == 'SO' & year == '2017' & site.region == 'J')
sock2018D <- bootlice %>% filter(spp == 'SO' & year == '2018' & site.region == 'D')
sock2018J <- bootlice %>% filter(spp == 'SO' & year == '2018' & site.region == 'J')
sock2019D <- bootlice %>% filter(spp == 'SO' & year == '2019' & site.region == 'D')
sock2019J <- bootlice %>% filter(spp == 'SO' & year == '2019' & site.region == 'J')

chum2015D <- bootlice %>% filter(spp == 'CU' & year == '2015' & site.region == 'D')
chum2015J <- bootlice %>% filter(spp == 'CU' & year == '2015' & site.region == 'J')
chum2016D <- bootlice %>% filter(spp == 'CU' & year == '2016' & site.region == 'D')
chum2016J <- bootlice %>% filter(spp == 'CU' & year == '2016' & site.region == 'J')
chum2017D <- bootlice %>% filter(spp == 'CU' & year == '2017' & site.region == 'D')
chum2017J <- bootlice %>% filter(spp == 'CU' & year == '2017' & site.region == 'J')
chum2018D <- bootlice %>% filter(spp == 'CU' & year == '2018' & site.region == 'D')
chum2018J <- bootlice %>% filter(spp == 'CU' & year == '2018' & site.region == 'J')
chum2019D <- bootlice %>% filter(spp == 'CU' & year == '2019' & site.region == 'D')
chum2019J <- bootlice %>% filter(spp == 'CU' & year == '2019' & site.region == 'J')

pink2015D <- bootlice %>% filter(spp == 'PI' & year == '2015' & site.region == 'D')
pink2015J <- bootlice %>% filter(spp == 'PI' & year == '2015' & site.region == 'J')
pink2016D <- bootlice %>% filter(spp == 'PI' & year == '2016' & site.region == 'D')
pink2016J <- bootlice %>% filter(spp == 'PI' & year == '2016' & site.region == 'J')
pink2017D <- bootlice %>% filter(spp == 'PI' & year == '2017' & site.region == 'D')
pink2017J <- bootlice %>% filter(spp == 'PI' & year == '2017' & site.region == 'J')
pink2018D <- bootlice %>% filter(spp == 'PI' & year == '2018' & site.region == 'D')
pink2018J <- bootlice %>% filter(spp == 'PI' & year == '2018' & site.region == 'J')
pink2019D <- bootlice %>% filter(spp == 'PI' & year == '2019' & site.region == 'D')
pink2019J <- bootlice %>% filter(spp == 'PI' & year == '2019' & site.region == 'J')

bootintervalcal = matrix(nrow = 30, ncol = 10000)
bootintervallep = matrix(nrow = 30, ncol = 10000)

#Run the Loop to populate the matrix - CAUTION:: THIS TAKES OVER 2 FULL DAYS TO RUN IN PARALLEL ON A POWERFUL DESKTOP COMPUTER
cores = detectCores()
cl = makeCluster(cores[1]-1)
registerDoParallel(cl)

start = Sys.time()

y = foreach(i = 1:10000, .packages = c('tidyverse', 'rsample', 'tibble', 'glmmTMB', 'MuMIn', 'ggeffects')) %dopar% {
  
  sock2015Dboot = matrix(nrow = nrow(sock2015D), ncol = 7)
  n = 1
  for(k in unique(sock2015D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2015D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2015Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2015Jboot = matrix(nrow = nrow(sock2015J), ncol = 7)
  n = 1
  for(k in unique(sock2015J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2015J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2015Jboot[rows,] = res
    n = n + nrow(res)
  }
  sock2016Dboot = matrix(nrow = nrow(sock2016D), ncol = 7)
  n = 1
  for(k in unique(sock2016D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2016D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2016Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2016Jboot = matrix(nrow = nrow(sock2016J), ncol = 7)
  n = 1
  for(k in unique(sock2016J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2016J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2016Jboot[rows,] = res
    n = n + nrow(res)
  }
  sock2017Dboot = matrix(nrow = nrow(sock2017D), ncol = 7)
  n = 1
  for(k in unique(sock2017D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2017D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2017Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2017Jboot = matrix(nrow = nrow(sock2017J), ncol = 7)
  n = 1
  for(k in unique(sock2017J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2017J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2017Jboot[rows,] = res
    n = n + nrow(res)
  }
  sock2018Dboot = matrix(nrow = nrow(sock2018D), ncol = 7)
  n = 1
  for(k in unique(sock2018D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2018D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2018Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2018Jboot = matrix(nrow = nrow(sock2018J), ncol = 7)
  n = 1
  for(k in unique(sock2018J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2018J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2018Jboot[rows,] = res
    n = n + nrow(res)
  }  
  sock2019Dboot = matrix(nrow = nrow(sock2019D), ncol = 7)
  n = 1
  for(k in unique(sock2019D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2019D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2019Dboot[rows,] = res
    n = n + nrow(res)
  }
  sock2019Jboot = matrix(nrow = nrow(sock2019J), ncol = 7)
  n = 1
  for(k in unique(sock2019J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(sock2019J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    sock2019Jboot[rows,] = res
    n = n + nrow(res)
  } 
  
  #pink
  
  chum2015Dboot = matrix(nrow = nrow(chum2015D), ncol = 7)
  n = 1
  for(k in unique(chum2015D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2015D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2015Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2015Jboot = matrix(nrow = nrow(chum2015J), ncol = 7)
  n = 1
  for(k in unique(chum2015J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2015J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2015Jboot[rows,] = res
    n = n + nrow(res)
  }
  chum2016Dboot = matrix(nrow = nrow(chum2016D), ncol = 7)
  n = 1
  for(k in unique(chum2016D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2016D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2016Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2016Jboot = matrix(nrow = nrow(chum2016J), ncol = 7)
  n = 1
  for(k in unique(chum2016J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2016J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2016Jboot[rows,] = res
    n = n + nrow(res)
  }
  chum2017Dboot = matrix(nrow = nrow(chum2017D), ncol = 7)
  n = 1
  for(k in unique(chum2017D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2017D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2017Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2017Jboot = matrix(nrow = nrow(chum2017J), ncol = 7)
  n = 1
  for(k in unique(chum2017J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2017J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2017Jboot[rows,] = res
    n = n + nrow(res)
  }
  chum2018Dboot = matrix(nrow = nrow(chum2018D), ncol = 7)
  n = 1
  for(k in unique(chum2018D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2018D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2018Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2018Jboot = matrix(nrow = nrow(chum2018J), ncol = 7)
  n = 1
  for(k in unique(chum2018J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2018J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2018Jboot[rows,] = res
    n = n + nrow(res)
  }  
  chum2019Dboot = matrix(nrow = nrow(chum2019D), ncol = 7)
  n = 1
  for(k in unique(chum2019D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2019D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2019Dboot[rows,] = res
    n = n + nrow(res)
  }
  chum2019Jboot = matrix(nrow = nrow(chum2019J), ncol = 7)
  n = 1
  for(k in unique(chum2019J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(chum2019J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    chum2019Jboot[rows,] = res
    n = n + nrow(res)
  } 
  
  #pink
  
  pink2015Dboot = matrix(nrow = nrow(pink2015D), ncol = 7)
  n = 1
  for(k in unique(pink2015D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2015D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2015Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2015Jboot = matrix(nrow = nrow(pink2015J), ncol = 7)
  n = 1
  for(k in unique(pink2015J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2015J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2015Jboot[rows,] = res
    n = n + nrow(res)
  }
  pink2016Dboot = matrix(nrow = nrow(pink2016D), ncol = 7)
  n = 1
  for(k in unique(pink2016D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2016D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2016Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2016Jboot = matrix(nrow = nrow(pink2016J), ncol = 7)
  n = 1
  for(k in unique(pink2016J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2016J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2016Jboot[rows,] = res
    n = n + nrow(res)
  }
  pink2017Dboot = matrix(nrow = nrow(pink2017D), ncol = 7)
  n = 1
  for(k in unique(pink2017D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2017D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2017Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2017Jboot = matrix(nrow = nrow(pink2017J), ncol = 7)
  n = 1
  for(k in unique(pink2017J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2017J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2017Jboot[rows,] = res
    n = n + nrow(res)
  }
  pink2018Dboot = matrix(nrow = nrow(pink2018D), ncol = 7)
  n = 1
  for(k in unique(pink2018D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2018D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2018Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2018Jboot = matrix(nrow = nrow(pink2018J), ncol = 7)
  n = 1
  for(k in unique(pink2018J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2018J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2018Jboot[rows,] = res
    n = n + nrow(res)
  }  
  pink2019Dboot = matrix(nrow = nrow(pink2019D), ncol = 7)
  n = 1
  for(k in unique(pink2019D$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2019D %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2019Dboot[rows,] = res
    n = n + nrow(res)
  }
  pink2019Jboot = matrix(nrow = nrow(pink2019J), ncol = 7)
  n = 1
  for(k in unique(pink2019J$collection)) { #for each collection
    #resample the observations in the collection with replacement for the number of rows in that given collection
    res <- as.matrix(pink2019J %>%
                       group_by(collection) %>% 
                       filter(collection == k) %>% 
                       sample_n(., n(), replace = TRUE))
    rows = c(n:(n+nrow(res)-1))
    pink2019Jboot[rows,] = res
    n = n + nrow(res)
  }  
  
  #bind the matrices so we can have our resampled dataframe
  bootdata = data.frame(rbind(sock2015Dboot,sock2015Jboot,sock2016Dboot,sock2016Jboot,sock2017Dboot,sock2017Jboot,sock2018Dboot,sock2018Jboot,sock2019Dboot,sock2019Jboot,
                              chum2015Dboot,chum2015Jboot,chum2016Dboot,chum2016Jboot,chum2017Dboot,chum2017Jboot,chum2018Dboot,chum2018Jboot,chum2019Dboot,chum2019Jboot,
                              pink2015Dboot,pink2015Jboot,pink2016Dboot,pink2016Jboot,pink2017Dboot,pink2017Jboot,pink2018Dboot,pink2018Jboot,pink2019Dboot,pink2019Jboot)) %>% 
    rename(all.cal = X1, all.leps = X2, spp = X3, site.region = X4, collection = X5, year = X6, ufn = X7)
  bootdata$all.cal = as.integer(as.character(bootdata$all.cal))
  bootdata$all.leps = as.integer(as.character(bootdata$all.leps))
  
  #now run our set of models       
  #omited the last two models in each set since their weights were 0. 
  lep1 = glmmTMB(all.leps ~ site.region + year + spp + 
                   site.region * year + 
                   (1 | collection), data = bootdata, family = poisson)
  lep2 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * site.region + site.region * year + 
                   (1 | collection), data = bootdata, family = poisson) 
  lep3 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * site.region + spp * year + site.region * year + 
                   (1 | collection), data = bootdata, family = poisson) 
  lep4 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * year + site.region * year + 
                   (1 | collection), data = bootdata, family = poisson) 
  lep5 = glmmTMB(all.leps ~ site.region + year + spp + 
                   (1 | collection), data = bootdata, family = poisson) 
  lep6 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * site.region + 
                   (1 | collection), data = bootdata, family = poisson) 
  lep7 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * year +  
                   (1 | collection), data = bootdata, family = poisson) 
  lep8 = glmmTMB(all.leps ~ site.region + year + spp + 
                   spp * site.region + spp * year +  
                   (1 | collection), data = bootdata, family = poisson) 

  
  cal1 = glmmTMB(all.cal ~ site.region + year + spp + 
                   spp * site.region + site.region * year + spp * year + 
                   (1 | collection), data = bootdata, family = nbinom2)
  cal2 = glmmTMB(all.cal ~ site.region + year + spp + 
                   spp * year + site.region * year + 
                   (1 | collection), data = bootdata, family = nbinom2)
  cal3 = glmmTMB(all.cal ~ site.region + year + spp + 
                   site.region * year + site.region * spp +
                   (1 | collection), data = bootdata, family = nbinom2) 
  cal4 = glmmTMB(all.cal ~ site.region + year + spp + 
                   site.region * year + 
                   (1 | collection), data = bootdata, family = nbinom2) 
  cal5 = glmmTMB(all.cal ~ site.region + year + spp + 
                   spp * year + 
                   (1 | collection), data = bootdata, family = nbinom2)
  cal6 = glmmTMB(all.cal ~ site.region + year + spp + 
                   spp * site.region + spp * year +  
                   (1 | collection), data = bootdata, family = nbinom2) 
  cal7 = glmmTMB(all.cal ~ site.region + year + spp + 
                   site.region * spp +  
                   (1 | collection), data = bootdata, family = nbinom2)
  cal8 = glmmTMB(all.cal ~ site.region + year + spp + 
                   (1 | collection), data = bootdata, family = nbinom2)
 

  
  #get the predictions of the estiamtes
  lep1pred <- ggpredict(lep1, terms = c('spp', 'year', 'site.region'))
  lep2pred <- ggpredict(lep2, terms = c('spp', 'year', 'site.region')) 
  lep3pred <- ggpredict(lep3, terms = c('spp', 'year', 'site.region')) 
  lep4pred <- ggpredict(lep4, terms = c('spp', 'year', 'site.region')) 
  lep5pred <- ggpredict(lep5, terms = c('spp', 'year', 'site.region')) 
  lep6pred <- ggpredict(lep6, terms = c('spp', 'year', 'site.region')) 
  lep7pred <- ggpredict(lep7, terms = c('spp', 'year', 'site.region')) 
  lep8pred <- ggpredict(lep8, terms = c('spp', 'year', 'site.region')) 
  
  cal1pred <- ggpredict(cal1, terms = c('spp', 'year', 'site.region'))
  cal2pred <- ggpredict(cal2, terms = c('spp', 'year', 'site.region')) 
  cal3pred <- ggpredict(cal3, terms = c('spp', 'year', 'site.region')) 
  cal4pred <- ggpredict(cal4, terms = c('spp', 'year', 'site.region')) 
  cal5pred <- ggpredict(cal5, terms = c('spp', 'year', 'site.region')) 
  cal6pred <- ggpredict(cal6, terms = c('spp', 'year', 'site.region')) 
  cal7pred <- ggpredict(cal7, terms = c('spp', 'year', 'site.region')) 
  cal8pred <- ggpredict(cal8, terms = c('spp', 'year', 'site.region')) 
  
  ###start by getting them all in one dataframe with the weights
  
  #pull the predicted values from each one
  lepallpred = data.frame(cbind(lep1pred$predicted, lep2pred$predicted, lep3pred$predicted, lep4pred$predicted,
                                lep5pred$predicted, lep6pred$predicted, lep7pred$predicted, lep8pred$predicted)) %>% 
    rename(lep1 = X1, lep2 = X2, lep3 = X3, lep4 = X4, lep5 = X5, lep6 = X6, lep7 = X7, lep8 = X8)
  calallpred = data.frame(cbind(cal1pred$predicted, cal2pred$predicted, cal3pred$predicted, cal4pred$predicted,
                                cal5pred$predicted, cal6pred$predicted, cal7pred$predicted, cal8pred$predicted)) %>% 
    rename(cal1 = X1, cal2 = X2, cal3 = X3, cal4 = X4, cal5 = X5, cal6 = X6, cal7 = X7, cal8 = X8)
  
  #add the weights from the model selection object
  lepallpred = lepallpred %>% 
    mutate(w1 = rep(lepmod.crossed_dredge$weight[1], nrow(lepallpred)), 
           w2 = rep(lepmod.crossed_dredge$weight[2], nrow(lepallpred)),
           w3 = rep(lepmod.crossed_dredge$weight[3], nrow(lepallpred)),
           w4 = rep(lepmod.crossed_dredge$weight[4], nrow(lepallpred)),
           w5 = rep(lepmod.crossed_dredge$weight[5], nrow(lepallpred)),
           w6 = rep(lepmod.crossed_dredge$weight[6], nrow(lepallpred)),
           w7 = rep(lepmod.crossed_dredge$weight[7], nrow(lepallpred)),
           w8 = rep(lepmod.crossed_dredge$weight[8], nrow(lepallpred)))
  
  calallpred = calallpred %>% 
    mutate(w1 = rep(calmod.crossed_dredge$weight[1], nrow(calallpred)), 
           w2 = rep(calmod.crossed_dredge$weight[2], nrow(calallpred)),
           w3 = rep(calmod.crossed_dredge$weight[3], nrow(calallpred)),
           w4 = rep(calmod.crossed_dredge$weight[4], nrow(calallpred)),
           w5 = rep(calmod.crossed_dredge$weight[5], nrow(calallpred)),
           w6 = rep(calmod.crossed_dredge$weight[6], nrow(calallpred)),
           w7 = rep(calmod.crossed_dredge$weight[7], nrow(calallpred)),
           w8 = rep(calmod.crossed_dredge$weight[8], nrow(calallpred)))
  
  #now make averaged predictions!
  lepallpred = lepallpred %>% 
    mutate(lep1w = lep1*w1, lep2w = lep2*w2, lep3w = lep3*w3, lep4w = lep4*w4, lep5w = lep5*w5, lep6w = lep6*w6, lep7w = lep7*w7, lep8w = lep8*w8) %>% 
    mutate(avg = lep1w + lep2w + lep3w + lep4w + lep5w + lep6w + lep7w + lep8w)
  
  calallpred = calallpred %>% 
    mutate(cal1w = cal1*w1, cal2w = cal2*w2, cal3w = cal3*w3, cal4w = cal4*w4, cal5w = cal5*w5, cal6w = cal6*w6, cal7w = cal7*w7, cal8w = cal8*w8) %>% 
    mutate(avg = cal1w + cal2w + cal3w + cal4w + cal5w + cal6w + cal7w + cal8w)
  
  #keep just the averaged predictions and the relevant grouping info 
  lepavgpred = lepallpred %>% 
    dplyr::select(avg) %>% 
    mutate(sal = lep1pred$x, reg = lep1pred$facet, yr = lep1pred$group)
  #lepavgpred$sal = factor(lepavgpred$sal, levels = c(1, 2, 3), labels = c('CU', 'PI', 'SO'))
  
  calavgpred = calallpred %>% 
    dplyr::select(avg) %>% 
    mutate(sal = cal1pred$x, reg = cal1pred$facet, yr = cal1pred$group)
  #calavgpred$sal = factor(calavgpred$sal, levels = c(1, 2, 3), labels = c('CU', 'PI', 'SO'))
  
 # bootintervalcal[,i] = calavgpred$avg
 # bootintervallep[,i] = lepavgpred$avg
  c(calavgpred$avg, lepavgpred$avg)

}
y_copy = y 

for(i in 1:10000) {
  bootintervalcal[,i] = y[[i]][1:30]
  bootintervallep[,i] = y[[i]][31:60]
}


stopCluster(cl)
end = Sys.time()
run_time = end - start


boot_int_cal = data.frame(bootintervalcal)
boot_int_lep = data.frame(bootintervallep)
write_csv(boot_int_cal, 'boot_int_cal.csv')
write_csv(boot_int_lep, 'boot_int_lep.csv')


interval_cal_long = as.data.frame(t(bootintervalcal))
interval_lep_long = as.data.frame(t(bootintervallep))


#name the columns, sort them, and transpose them
names_lep = lepavgpred %>% 
  unite(., col = 'names', sal:yr, sep = '_')
names_cal = calavgpred %>% 
  unite(., col = 'names', sal:yr, sep = '_')

names_lep = as.vector(names_lep$names)
names_cal = as.vector(names_cal$names)

colnames(interval_cal_long) = names_cal
colnames(interval_lep_long) = names_lep

interval_cal_long_sorted <- apply(interval_cal_long,2,sort,decreasing=F)
interval_lep_long_sorted <- apply(interval_lep_long,2,sort,decreasing=F)

upci_cal = interval_cal_long_sorted[9750, ]
loci_cal = interval_cal_long_sorted[250, ]
upci_lep = interval_lep_long_sorted[9750, ]
loci_lep = interval_lep_long_sorted[250, ]

#fit the models and do the model averaging process again to get the actual estimates themselves
lep1 = glmmTMB(all.leps ~ site.region + year + spp + 
                 site.region * year + 
                 (1 | collection), data = mainlice, family = poisson)
lep2 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * site.region + site.region * year + 
                 (1 | collection), data = mainlice, family = poisson) 
lep3 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * site.region + spp * year + site.region * year + 
                 (1 | collection), data = mainlice, family = poisson) 
lep4 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * year + site.region * year + 
                 (1 | collection), data = mainlice, family = poisson) 
lep5 = glmmTMB(all.leps ~ site.region + year + spp + 
                 (1 | collection), data = mainlice, family = poisson) 
lep6 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * site.region + 
                 (1 | collection), data = mainlice, family = poisson) 
lep7 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * year +  
                 (1 | collection), data = mainlice, family = poisson) 
lep8 = glmmTMB(all.leps ~ site.region + year + spp + 
                 spp * site.region + spp * year +  
                 (1 | collection), data = mainlice, family = poisson) 


cal1 = glmmTMB(all.cal ~ site.region + year + spp + 
                 spp * site.region + site.region * year + spp * year + 
                 (1 | collection), data = mainlice, family = nbinom2)
cal2 = glmmTMB(all.cal ~ site.region + year + spp + 
                 spp * year + site.region * year + 
                 (1 | collection), data = mainlice, family = nbinom2)
cal3 = glmmTMB(all.cal ~ site.region + year + spp + 
                 site.region * year + site.region * spp +
                 (1 | collection), data = mainlice, family = nbinom2) 
cal4 = glmmTMB(all.cal ~ site.region + year + spp + 
                 site.region * year + 
                 (1 | collection), data = mainlice, family = nbinom2) 
cal5 = glmmTMB(all.cal ~ site.region + year + spp + 
                 spp * year + 
                 (1 | collection), data = mainlice, family = nbinom2)
cal6 = glmmTMB(all.cal ~ site.region + year + spp + 
                 spp * site.region + spp * year +  
                 (1 | collection), data = mainlice, family = nbinom2) 
cal7 = glmmTMB(all.cal ~ site.region + year + spp + 
                 site.region * spp +  
                 (1 | collection), data = mainlice, family = nbinom2)
cal8 = glmmTMB(all.cal ~ site.region + year + spp + 
                 (1 | collection), data = mainlice, family = nbinom2)

#get the predictions of the estiamtes
lep1pred <- ggpredict(lep1, terms = c('spp', 'year', 'site.region'))
lep2pred <- ggpredict(lep2, terms = c('spp', 'year', 'site.region')) 
lep3pred <- ggpredict(lep3, terms = c('spp', 'year', 'site.region')) 
lep4pred <- ggpredict(lep4, terms = c('spp', 'year', 'site.region')) 
lep5pred <- ggpredict(lep5, terms = c('spp', 'year', 'site.region')) 
lep6pred <- ggpredict(lep6, terms = c('spp', 'year', 'site.region')) 
lep7pred <- ggpredict(lep7, terms = c('spp', 'year', 'site.region')) 
lep8pred <- ggpredict(lep8, terms = c('spp', 'year', 'site.region')) 

cal1pred <- ggpredict(cal1, terms = c('spp', 'year', 'site.region'))
cal2pred <- ggpredict(cal2, terms = c('spp', 'year', 'site.region')) 
cal3pred <- ggpredict(cal3, terms = c('spp', 'year', 'site.region')) 
cal4pred <- ggpredict(cal4, terms = c('spp', 'year', 'site.region')) 
cal5pred <- ggpredict(cal5, terms = c('spp', 'year', 'site.region')) 
cal6pred <- ggpredict(cal6, terms = c('spp', 'year', 'site.region')) 
cal7pred <- ggpredict(cal7, terms = c('spp', 'year', 'site.region')) 
cal8pred <- ggpredict(cal8, terms = c('spp', 'year', 'site.region')) 

###start by getting them all in one dataframe with the weights

#pull the predicted values from each one
lepallpred = data.frame(cbind(lep1pred$predicted, lep2pred$predicted, lep3pred$predicted, lep4pred$predicted,
                              lep5pred$predicted, lep6pred$predicted, lep7pred$predicted, lep8pred$predicted)) %>% 
  rename(lep1 = X1, lep2 = X2, lep3 = X3, lep4 = X4, lep5 = X5, lep6 = X6, lep7 = X7, lep8 = X8)
calallpred = data.frame(cbind(cal1pred$predicted, cal2pred$predicted, cal3pred$predicted, cal4pred$predicted,
                              cal5pred$predicted, cal6pred$predicted, cal7pred$predicted, cal8pred$predicted)) %>% 
  rename(cal1 = X1, cal2 = X2, cal3 = X3, cal4 = X4, cal5 = X5, cal6 = X6, cal7 = X7, cal8 = X8)

#add the weights from the model selection object
lepallpred = lepallpred %>% 
  mutate(w1 = rep(lepmod.crossed_dredge$weight[1], nrow(lepallpred)), 
         w2 = rep(lepmod.crossed_dredge$weight[2], nrow(lepallpred)),
         w3 = rep(lepmod.crossed_dredge$weight[3], nrow(lepallpred)),
         w4 = rep(lepmod.crossed_dredge$weight[4], nrow(lepallpred)),
         w5 = rep(lepmod.crossed_dredge$weight[5], nrow(lepallpred)),
         w6 = rep(lepmod.crossed_dredge$weight[6], nrow(lepallpred)),
         w7 = rep(lepmod.crossed_dredge$weight[7], nrow(lepallpred)),
         w8 = rep(lepmod.crossed_dredge$weight[8], nrow(lepallpred)))

calallpred = calallpred %>% 
  mutate(w1 = rep(calmod.crossed_dredge$weight[1], nrow(calallpred)), 
         w2 = rep(calmod.crossed_dredge$weight[2], nrow(calallpred)),
         w3 = rep(calmod.crossed_dredge$weight[3], nrow(calallpred)),
         w4 = rep(calmod.crossed_dredge$weight[4], nrow(calallpred)),
         w5 = rep(calmod.crossed_dredge$weight[5], nrow(calallpred)),
         w6 = rep(calmod.crossed_dredge$weight[6], nrow(calallpred)),
         w7 = rep(calmod.crossed_dredge$weight[7], nrow(calallpred)),
         w8 = rep(calmod.crossed_dredge$weight[8], nrow(calallpred)))

#now make averaged predictions!
lepallpred = lepallpred %>% 
  mutate(lep1w = lep1*w1, lep2w = lep2*w2, lep3w = lep3*w3, lep4w = lep4*w4, lep5w = lep5*w5, lep6w = lep6*w6, lep7w = lep7*w7, lep8w = lep8*w8) %>% 
  mutate(avg = lep1w + lep2w + lep3w + lep4w + lep5w + lep6w + lep7w + lep8w)

calallpred = calallpred %>% 
  mutate(cal1w = cal1*w1, cal2w = cal2*w2, cal3w = cal3*w3, cal4w = cal4*w4, cal5w = cal5*w5, cal6w = cal6*w6, cal7w = cal7*w7, cal8w = cal8*w8) %>% 
  mutate(avg = cal1w + cal2w + cal3w + cal4w + cal5w + cal6w + cal7w + cal8w)

#keep just the averaged predictions and the relevant grouping info 
lepavgpred = lepallpred %>% 
  dplyr::select(avg) %>% 
  mutate(sal = lep1pred$x, reg = lep1pred$facet, yr = lep1pred$group)

calavgpred = calallpred %>% 
  dplyr::select(avg) %>% 
  mutate(sal = cal1pred$x, reg = cal1pred$facet, yr = cal1pred$group)

#put the up and lo CI's into the df's 
sal_order = c("CU","PI", "SO")
yr_order = c('2015', '2016', '2017', '2018')
calavgpred = calavgpred[order(factor(calavgpred$sal, levels = sal_order), 
                              factor(calavgpred$yr, levels = yr_order)),]

lepavgpred = lepavgpred[order(factor(lepavgpred$sal, levels = sal_order),
                              factor(lepavgpred$yr, levels = yr_order)),]

calavgpred$conf.high = upci_cal
calavgpred$conf.low = loci_cal
lepavgpred$conf.high = upci_lep
lepavgpred$conf.low = loci_lep

#Make the plots
fte_theme1 <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 16, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 13, color = color.axis.text)) + 
    theme(axis.text.y = element_text(size = 13, color = color.axis.text)) + 
    theme(axis.title.x = element_text(size = 15, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 15, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.ticks.x = element_line(colour = 'black')) +
    theme(axis.line.x = element_line(color="black", size = 0.15),
          axis.line.y = element_line(color="black", size = 0.15)) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 13))+
    theme(legend.position = c(0.6,0.82),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13))
} 

leg_title <- 'Salmon Species'
lep_y_effects = expression(paste("Mean Motile  ", italic("L. salmonis"), " Abundance"))
cal_y_effects = expression(paste("Mean Motile  ", italic("C. clemensi"), " Abundance"))
lepavgpred$sal <- factor(lepavgpred$sal,levels = c('PI', 'CU', 'SO'))
calavgpred$sal <- factor(calavgpred$sal,levels = c('PI', 'CU', 'SO'))
lepavgpred$yr = factor(lepavgpred$yr, levels = c('2015', '2016', '2017', '2018', '2019'))
calavgpred$yr = factor(calavgpred$yr, levels = c('2015', '2016', '2017', '2018', '2019'))


lepsmodplot_avg <- lepavgpred %>% 
  group_by(., yr,sal,reg) %>% 
  ggplot(aes(x = sal, y = avg, colour = sal, shape = reg)) +
  scale_shape_manual(values = c(15,17), labels = c('Discovery Islands', 'Johnstone Strait')) +
  geom_errorbar(aes(ymin=conf.low, ymax = conf.high,width = 0), size = 0.78, position = position_dodge(width = 0.835),colour = 'Black')+
  geom_point(size = 4.7,position = position_dodge(width = 0.8)) +
  facet_wrap(~yr,nrow=1,strip.position = "bottom")+
  theme(strip.background = element_blank(), strip.placement = "outside") + 
  scale_color_manual(leg_title,values=c('#ff9999','#59AE7F','#23359d'), labels = c('Pink', 'Chum', 'Sockeye'))+
  labs(x = 'Salmon Species/Year', y = lep_y_effects) +
  guides(shape = guide_legend(title = 'Region', override.aes = list(shape = c(0,2)), type = 'b')) +
  scale_y_continuous(limits = c(0,1.5), breaks = c(0.5,1.0,1.5))+
  fte_theme1()
lepsmodplot_avg
ggsave('model_ests_lep.png', plot = lepsmodplot_avg,
       width = 8, height = 7.5,
       dpi = 300)
calmodplot_avg <- calavgpred %>% 
  group_by(., yr,sal,reg) %>% 
  ggplot(aes(x = sal, y = avg, colour = sal, shape = reg)) +
  scale_shape_manual(values = c(15,17), labels = c('Discovery Islands', 'Johnstone Strait')) +
  geom_errorbar(aes(ymin=conf.low, ymax = conf.high,width = 0), size = 0.78, position = position_dodge(width = 0.835),colour = 'Black')+
  geom_point(size = 4.7,position = position_dodge(width = 0.8)) +
  facet_wrap(~yr,nrow=1,strip.position = "bottom")+
  theme(strip.background = element_blank(), strip.placement = "outside") + 
  scale_colour_manual(leg_title,values=c('#ff9999','#59AE7F','#23359d'), labels = c('Pink', 'Chum', 'Sockeye'))+
  labs(x = 'Salmon Species/Year', y = cal_y_effects) +
  guides(shape = guide_legend(title = 'Region', override.aes = list(shape = c(0,2)), type = 'b')) +
  fte_theme1()
calmodplot_avg

ggsave('model_ests_cal.png', plot = calmodplot_avg,
       width = 8, height = 7.5,
       dpi = 300)


# Relative Importance Extraction
sw(lepmod.crossed_dredge)
sw(calmod.crossed_dredge)


# Get average estimated lice 

cal_allyears = glmmTMB(all.cal ~ spp +
                 (1 | collection), data = mainlice, family = nbinom2)
lep_allyears = glmmTMB(all.leps ~ spp +
                         (1 | collection), data = mainlice, family = nbinom2)
cal_allyears_predict = ggpredict(cal_allyears, terms = c('spp'), ci.lvl = 0.95)
lep_allyears_predict = ggpredict(lep_allyears, terms = c('spp'), ci.lvl = 0.95)

# Make histograms of fish caught to show how many of our focal species there are
fish = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/fish_field_data.csv")
seine_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/seine_data.csv")
survey_data = read_csv("C:/Users/brookson/Documents/GitHub/jsp-data/data/survey_data.csv")


fish = left_join(fish, seine_data, by = "seine_id")
fish = left_join(fish, survey_data, by = 'survey_id')
fish$year = substr(fish$survey_date, 0, 4)
fish$year = as.factor(fish$year)

#Make the plots
fte_theme3 <- function(){
  color.background = 'white'
  color.grid.major = 'black'
  color.axis.text = 'black'
  color.axis.title = 'black'
  color.title = 'black'
  theme_bw(base_size = 9) + 
    theme(panel.background = element_rect(fill=color.background,color = color.background)) +
    theme(plot.background = element_rect(fill = color.background, color = color.background)) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank()) + 
    theme(panel.grid.minor = element_blank()) + 
    theme(axis.ticks = element_blank()) +
    theme(plot.title = element_text(color = color.title, size = 16, vjust = 1.25)) +
    theme(axis.text.x = element_text(size = 13, color = color.axis.text)) + 
    theme(axis.text.y = element_text(size = 13, color = color.axis.text)) + 
    theme(axis.title.x = element_text(size = 15, color = color.axis.title, vjust = 0)) +
    theme(axis.title.y = element_text(size = 15, color = color.axis.title, vjust = 1.25)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.ticks.x = element_line(colour = 'black')) +
    theme(axis.line.x.bottom = element_line(color="black", size = 0.15),
          axis.line.y.left = element_line(color="black", size = 0.15),
          axis.line.x.top = element_line(color="black", size = 0.15),
          axis.line.y.right = element_line(color="black", size = 0.15)) +
    theme(strip.background = element_blank(),
          strip.placement = 'outside',
          strip.text = element_text(size = 13))+
    theme(legend.position = c(0.8,0.82),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13))
} 

leg_title = 'Fish Speices'
fish_caught = ggplot(data = fish, aes(x = year)) +
  geom_bar(aes(x = year, fill = species), colour = 'black', position = 'dodge', stat = 'count')+
  scale_fill_manual(leg_title,values=c('skyblue1', 'goldenrod3', '#59AE7F', 'red4', '#ff9999','#23359d'), labels = c('Chinook', 'Coho', 'Chum', 'Herring', 'Pink', 'Sockeye'))+
  fte_theme3() +
  labs(x = 'Year', y = 'Number of Fish Caught')


ggsave('all_fish_caught.png', plot = fish_caught,
       width = 8, height = 7.5,
       dpi = 300)

sites_year_spp = mainlice %>% 
  group_by(year, site_id, spp) %>% 
  summarize(obs = n())
write.table(sites_year_spp, file = "sites_year_spp.txt", sep = ",", quote = FALSE, row.names = F)





