##### TITLE: Differential Infection of Parasitic Sea Lice on Juvenile Pacific Salmon in British Columbia, Canada
##### CREATOR: Cole B. Brookson
##### INITIALIZATION DATE: 2019-11-21

library(here)
source(here('./code/figures_pre_model_runs.R'))


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
y = foreach(i = 1:10, .packages = c('tidyverse', 'rsample', 'tibble', 'glmmTMB', 'MuMIn', 'ggeffects'),
            .export = c(lepallpred, calallpred, lepavgpred, calavgpred)) %dopar% {
  
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

for(i in 1:10) {
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

#pull shape parameters and fiex-effects estimates
fixed_cal1 = fixef(cal1)
unlist(fixed_cal1)
print.default(fixed_cal1)

fixed_lep1 = fixef(lep1)
unlist(fixed_lep1)
print.default(fixed_lep1)

cal_shape = sigma(cal1)
lep_shape = sigma(lep1)

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

