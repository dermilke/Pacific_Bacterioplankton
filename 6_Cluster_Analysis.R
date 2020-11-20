### Global Cluster Visualization

##### Load dependencies #####

rm(list = ls(all=T))

library(tidyverse)
library(ggsci)

devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datalistFunctions.R")
devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datatableFunctions.R")

##### Core-Analysis Function #####

cluster_analysis <- function(datalist, Depth_Filter = "Epi", SF_Filter = 0.22, tax_lvl = Family, k_cluster = NULL, 
                             label_order = NULL, title = NULL) {
  
  tax_lvl_en <- enquo(tax_lvl)
  
  datatable_summarized <- datalist %>%
    filter_station_datalist(Size_Fraction == !!SF_Filter) %>%
    filter_station_datalist(Depth_Grp == !!Depth_Filter) %>%
    summarize_by_taxa(!!tax_lvl_en) %>%
    create_datatable(., grpBy = !!tax_lvl_en, otherThreshold = 0) %>%
    group_by(Group) %>%
    filter(Group != "Others") %>%
    mutate(Abundance = (Abundance - mean(Abundance))/sd(Abundance)) %>%
    mutate(Abundance = loess(Abundance ~ Latitude, span = .7)$fitted) 
  
  Count_Matrix <- datatable_summarized %>%
    select(Group, Latitude, Abundance) %>%
    reshape2::dcast(., Group ~ Latitude)
  
  if (is.null(k_cluster)) {
    k_cluster <- factoextra::fviz_nbclust(select_if(Count_Matrix, is.numeric), factoextra::hcut, method = "silhouette", k.max = min(nrow(Count_Matrix)-1, 9))$data$y %>%
      which.max()
  }
  
  datatable_clustered <- Count_Matrix %>%
    select_if(is.numeric) %>%
    vegan::vegdist(method = "euclid") %>%
    hclust(method = "ward.D2") %>%
    cutree(., k = k_cluster) %>%
    tibble(Cluster = .) %>%
    mutate(Group = Count_Matrix$Group) %>%
    left_join(., datatable_summarized, by = "Group") %>%
    filter(!is.na(Cluster)) %>%
    mutate(Cluster = if (!is.null(label_order)) as.numeric(factor(Cluster, levels = label_order)) else Cluster)
  
  datatable_mean <- datatable_clustered %>%
    group_by(Latitude, Depth_Grp, Cluster) %>%
    summarize(Abundance = mean(Abundance)) %>%
    filter(!is.na(Cluster))
  
  ggplot(datatable_clustered, aes(x = Latitude, y = Abundance, col = Group)) +
    geom_line() +
    geom_smooth(aes(x = Latitude, y = Abundance), col = "black", size = 2, se = F, span = .7, data = datatable_mean) +
    geom_smooth(aes(x = Latitude, y = Abundance), col = "red", se = F, span = .7, data = datatable_mean) +
    theme(legend.position = "none") +
    facet_wrap(~Cluster, nrow = 1) +
    labs(title = if (!is.null(title)) title else "Clustered abundance patterns", y = "Normalized Abundance")
}

##### Read processed tables #####

datalist_Pacific <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok") %>%
  mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
  mutate_meta_datalist(Latitude = round(Latitude),
                       Depth_Grp = ifelse(Depth < 100, "Epi", "Meso")) %>%
  filter_station_datalist(Counts_Total > 2000) %>%
  filter_abundance() %>%
  correct_ambiguous(.)

##### Data-wrangling Pipeline #####

datalist_Pacific_summarized <- datalist_Pacific %>%
  summarize_by_param(Latitude, Depth_Grp, Size_Fraction)

##### Cluster-Analysis #####

p1 <- cluster_analysis(datalist_Pacific_summarized, Depth_Filter = "Epi", SF_Filter = 0.22, OTU_ID, label_order = c(5,3,1,2,4), 
                       title = "Epipelagic layer\nFL Community")
p2 <- cluster_analysis(datalist_Pacific_summarized, Depth_Filter = "Epi", SF_Filter = 3, OTU_ID, label_order = c(4,5,1,2,3),
                       title = "\nSPA Community")
p3 <- cluster_analysis(datalist_Pacific_summarized, Depth_Filter = "Epi", SF_Filter = 8, OTU_ID, label_order = c(4,3,1,2,5),
                       title = "\nLPA Community")

p4 <- cluster_analysis(datalist_Pacific_summarized, Depth_Filter = "Meso", SF_Filter = 0.22, OTU_ID, k_cluster = 4, label_order = c(4,2,3,1),
                       title = "Mesopelagic layer\nFL Community")
p5 <- cluster_analysis(datalist_Pacific_summarized, Depth_Filter = "Meso", SF_Filter = 3, OTU_ID, label_order = c(2,1,4,3),
                       title = "\nSPA Community")
p6 <- cluster_analysis(datalist_Pacific_summarized, Depth_Filter = "Meso", SF_Filter = 8, OTU_ID, label_order = c(4,1,2,3),
                       title = "\nLPA Community")

cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2)
