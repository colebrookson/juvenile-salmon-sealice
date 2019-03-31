regionlice <- read.csv('Hakai_lice_data_all_fish.csv')

library(lubridate)
library(tidyverse)
library(glmmADMB)

#add in the summed lice columns
regionlice <- regionlice %>%
    rowwise() %>%
    mutate(all.cal = sum(c(cal.mot.notgf,cal.gf), na.rm = TRUE))
regionlice <- regionlice %>%
    rowwise() %>%
    mutate(all.leps = sum(c(lep.pam,lep.paf,lep.am,lep.af,lep.gf,lep.mot.unid.notgf), na.rm = TRUE))

#add in the site region column
regionlice <- regionlice %>%
  mutate(site.region = site.id)
regionlice$site.region<- substr(regionlice$site.region, 0, 1)

#add a column for week in the year
regionlice <- regionlice %>%
  mutate(week = date)
regionlice$week <- strftime(regionlice$week, format = "%V")

#now add the year in
regionlice <- regionlice %>%
   mutate(year = date)
regionlice$year <- year(regionlice$year)


## now we can write the models - we'll write one for each fish species
#make some sub-dfs first
regionlice$week <- as.factor(regionlice$week); regionlice$year <- as.factor(regionlice$year)
chum.region <- regionlice %>% 
  filter(spp == 'CU')
pink.region <- regionlice %>% 
  filter(spp == 'PI')
sock.region <- regionlice %>% 
  filter(spp == 'SO')
sock.region2 <- sock.region[sample(nrow(sock.region),3100),]
sock.region3 <- sock.region[sample(nrow(sock.region),2500),]

#models with year as a random effect
chumrmod.cal1 <- glmmadmb(all.cal ~ site.region - 1 + (1|week) + (1|year), 
                          data = chum.region, family = 'nbinom', zeroInflation = FALSE) 
chumrmod.leps1 <- glmmadmb(all.leps ~ site.region - 1 + (1|week) + (1|year), 
                           data = chum.region, family = 'nbinom', zeroInflation = FALSE)#this one doesn't run
pinkrmod.cal1 <- glmmadmb(all.cal ~ site.region - 1  + (1|year) + (1|week), 
                          data = pink.region, family = 'nbinom', zeroInflation = FALSE) 
pinkrmod.leps1 <- glmmadmb(all.leps ~ site.region - 1 + (1|week) + (1|year), 
                           data = pink.region, family = 'nbinom', zeroInflation = FALSE)
sockrmod.cal1 <- glmmadmb(all.cal ~ site.region - 1 + (1|week) + (1|year), 
                          data = sock.region, family = 'nbinom', zeroInflation = FALSE)
sockrmod.leps1 <- glmmadmb(all.leps ~ site.region - 1 + (1|week) + (1|year), 
                           data = sock.region2, family = 'nbinom', zeroInflation = FALSE)#runs with reduced observations sock.region2, not with sock.region

summary(chumrmod.cal1)
summary(chumrmod.leps1)
summary(pinkrmod.cal1)
summary(pinkrmod.leps1)
summary(sockrmod.cal1)
summary(sockrmod.leps1)

#models with year as fixed effect

chumrmod.cal <- glmmadmb(all.cal ~ site.region + year -1 + (1|week), 
                         data = chum.region, family = 'nbinom', zeroInflation = FALSE)
chumrmod.leps <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                          data = chum.region, family = 'nbinom', zeroInflation = FALSE)
pinkrmod.cal <- glmmadmb(all.cal ~ site.region + year - 1  + (1|year), 
                         data = pink.region, family = 'nbinom', zeroInflation = FALSE)
pinkrmod.leps <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                          data = pink.region, family = 'nbinom', zeroInflation = FALSE)
sockrmod.cal <- glmmadmb(all.cal ~ site.region + year - 1 + (1|week), 
                         data = sock.region, family = 'nbinom', zeroInflation = FALSE)
sockrmod.leps <- glmmadmb(all.leps ~ site.region + year - 1 + (1|week), 
                          data = sock.region, family = 'nbinom', zeroInflation = FALSE) #only runs with reduced number of observation (such as in sock.region3)

summary(chumrmod.cal)
summary(chumrmod.leps)
summary(pinkrmod.cal)
summary(pinkrmod.leps)
summary(sockrmod.cal)
summary(sockrmod.leps)

## Things I've tried to fix it

# 1. running with admbControl(shess=FALSE,noinit=FALSE) set
# 2. running with the extra command: '...zeroInflation = FALSE, extra.args= '-ndi 100000'')
# 3. running with zeroInflation = TRUE


