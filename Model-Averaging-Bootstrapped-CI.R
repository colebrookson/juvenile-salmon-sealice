library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(MuMIn)
library(cowplot)
library(AICcmodavg)
library(latexpdf)

mainlice <- read.csv("Hakai_lice_data_CB_edits.csv")

#make vars into factors
mainlice$year <- as.factor(mainlice$year);mainlice$collection <- as.factor(mainlice$collection)