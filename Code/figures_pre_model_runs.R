##### TITLE: Differential Infection of Parasitic Sea Lice on Juvenile Pacific Salmon in British Columbia, Canada
##### CREATOR: Cole B. Brookson
##### INITIALIZATION DATE: 2019-11-21

library(here)
source(here('./code/data_manip_prep.R'))

coords <- data.frame(cbind(salmonsites$site_lat,salmonsites$site_long))
colnames(coords) <- c('lat', 'long')

#sort(unique(ggplot2::map_data("world")$region))
#westcoast <- ggplot2::map_data('Canada') ---- DEPRECATED

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
ggsave(here('./figures/study_map.png'), plot = insetmap,
       width = 8, height = 7.5,
       dpi = 200)
ggsave(here('./figures/study_map_fullres.png'), plot = insetmap,
       width = 8, height = 7.5,
       dpi = 1200)

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
ggsave(here('./figures/lice_per_fish_sp.png'), plot = fig_2,
       width = 8, height = 8.0,
       dpi = 300)
ggsave(here('./figures/lice_per_fish_sp_fullres.png'), plot = fig_2,
       width = 8, height = 8.0,
       dpi = 1200)

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
ggsave(here('./figures/lice_per_fish.png'), plot = fig_3,
       scale = 1, width = 8, height = 6.7,
       dpi = 300)
ggsave(here('./figures/lice_per_fish_fullres.png'), plot = fig_3,
       scale = 1, width = 8, height = 6.7,
       dpi = 1200)
