#### Distance Decay & Mean Latitudinal Range ####

#### Load dependencies ####

rm(list = ls(all=T))

library(tidyverse)

devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datalistFunctions.R")
devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datatableFunctions.R")

##### Read processed tables #####

datalist_Pacific <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok") %>%
  mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
  filter_station_datalist(Counts_Total > 8000) %>%
  mutate_meta_datalist(Depth_Grp = c(20, 40, 60, 100, 200, 300, 500)[findInterval(Depth, c(20, 40, 60, 100, 200, 300, 500))]) %>%
  mutate_meta_datalist(Depth_Grp_2 = ifelse(.$Meta_Data$Depth_Grp < 100, "Epi", "Meso")) %>%
  filter_abundance()

##### Core Analysis Functions #####

distance_decay <- function(datalist) {
  
  geographic_distance <- function(long1, lat1, long2, lat2) {
    
    dlong <- long1 * pi/180 - long2 * pi/180
    dlat <- lat1 * pi/180 - lat2 * pi/180
    
    R <- 6371
    a <- sin(dlat/2)^2 + cos(lat1 * pi/180) * cos(lat2 * pi/180) * sin(dlong/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    
    distance <- R * c
    
    return(distance)
    
  }
  
  DD_datatable <- datalist$Count_Data %>%
    select_if(is.numeric) %>%
    t(.) %>% 
    vegan::vegdist(.) %>% 
    t(.) %>% 
    as.matrix(.) %>%
    reshape2::melt() %>%
    dplyr::rename("From" = "Var1", "Sample_ID" = "Var2", !!dist_measure := "value") %>%
    left_join(., datalist$Meta_Data %>%
                select(Sample_ID, Latitude, Longitude, Province, Depth_Grp), by = "Sample_ID") %>%
    dplyr::rename("To" = "Sample_ID", "To_Latitude" = "Latitude", "To_Province" = "Province", 
           "To_Longitude" = "Longitude", "Sample_ID" = "From", "To_Depth" = "Depth_Grp") %>%
    left_join(., datalist$Meta_Data %>%
                select(Sample_ID, Latitude, Longitude, Province, Depth_Grp), by = "Sample_ID") %>%
    dplyr::rename("From" = "Sample_ID", "From_Latitude" = "Latitude", "From_Longitude" = "Longitude", 
           "From_Province" = "Province", "From_Depth" = "Depth_Grp") %>%
    mutate(Distance = geographic_distance(To_Longitude, To_Latitude, From_Longitude, From_Latitude)) %>%
    as_tibble()
  
  return(DD_datatable)
  
}

##### Data-wrangling Pipeline #####

datalist_Pacific_range <- datalist_Pacific

datalist_Pacific_range$Count_Data <- datalist_Pacific_range %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  with(., apply(., 1, function(x) ifelse(x > 0, filter(datalist_Pacific_range$Meta_Data, x > 0) %>% with(., range(.$Latitude)) %>% diff, 0))) %>%
  t() %>%
  as_tibble() %>%
  with(., bind_cols(select_if(datalist_Pacific_range$Count_Data, is.character), .)) 

datatable_Pacific_range <- datalist_Pacific_range %>%
  create_datatable(grpBy = OTU_ID, otherThreshold = 0) %>% 
  filter(Abundance > 0) %>%
  mutate(Size_Fraction = ordered(ifelse(Size_Fraction == 0.22, "FL",
                                 ifelse(Size_Fraction == 3, "SPA", "LPA")), levels = c("FL", "SPA", "LPA"))) %>%
  group_by(Latitude, Depth_Grp_2, Size_Fraction) %>%
  summarize(Abundance = mean(Abundance))

DD_Pacific_table <- rbind(datalist_Pacific %>%
                            filter_station_datalist(Size_Fraction == 0.22) %>%
                            distance_decay() %>%
                            mutate(Size_Fraction = "FL"),
                          datalist_Pacific %>%
                            filter_station_datalist(Size_Fraction == 3) %>%
                            distance_decay() %>%
                            mutate(Size_Fraction = "SPA"),
                          datalist_Pacific %>%
                            filter_station_datalist(Size_Fraction == 8) %>%
                            distance_decay() %>%
                            mutate(Size_Fraction = "LPA")) %>%
  mutate(Size_Fraction = ordered(Size_Fraction, levels = c("FL", "SPA", "LPA"))) %>%
  filter(From_Depth == To_Depth) %>%
  mutate(To_Province = ordered(To_Province, levels = unique(arrange(., To_Latitude)$To_Province))) %>%
  mutate(From_Province = ordered(From_Province, levels = unique(arrange(., From_Latitude)$From_Province))) %>%
  mutate(Province_Difference = abs(as.numeric(From_Province)-as.numeric(To_Province))) %>%
  filter(Distance > 0) %>%
  mutate(From_Depth = ordered(From_Depth, level = sort(unique(From_Depth)), label = paste(sort(unique(From_Depth)), "m", sep = ""))) %>%
  mutate(From_Depth_Grp = ifelse(From_Depth < "100m", "Epi", "Meso")) %>%
  
#### Distance-Decay Analysis ####
  
ggplot(DD_Pacific_table, aes(x = Distance, y = Bray_Curtis)) +
  geom_point(aes(col = as.factor(Province_Difference)), alpha = 0.3) +
  geom_smooth(se = F, col = "black") +
  facet_grid(From_Depth_Grp~Size_Fraction) +
  labs(y = "Bray-Curtis Dissimilarity", x = "Distance [km]", col = "No. of Provinces\napart") +
  guides(col = guide_legend(nrow = 6, title.position = "top", override.aes = list(alpha = 1))) +
  theme(panel.spacing = unit(1.2, "lines"))

#### Mean Latitudinal Range of ASVs ####
  
ggplot(datatable_Pacific_range, aes(x = Latitude, y = Abundance, col = Size_Fraction)) +
  geom_point() +
  geom_smooth(se = F) +
  scale_color_discrete() +
  labs(y = "Mean Latitudinal Range of ASVs [° Lat]", x = "Latitude", col = "Size Fraction") +
  facet_wrap(~Depth_Grp_2, ncol = 1, strip.position = "right")
