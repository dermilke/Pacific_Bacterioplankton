### Diversity analysis

##### Load dependencies #####

rm(list = ls(all=T))

library(tidyverse)
library(VennDiagram)

devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datalistFunctions.R")
devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datatableFunctions.R")

##### Read processed tables #####

datalist_Pacific <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok") %>%
  mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
  filter_station_datalist(Counts_Total > 8000) %>%
  mutate_meta_datalist(Depth_Grp = c(20, 40, 60, 100, 200, 300, 500)[findInterval(Depth, c(20, 40, 60, 100, 200, 300, 500))]) %>%
  mutate_meta_datalist(Depth_Grp_2 = ifelse(Depth_Grp < 100, "Epi", "Meso")) %>%
  filter_abundance()

datalist_Pacific_unfiltered <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok") %>%
  mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
  filter_station_datalist(Counts_Total > 8000)

datalist_Pacific_unassigned <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete_NA/", "Prok") %>%
  mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
  filter_station_datalist(Counts_Total > 8000) %>%
  filter_abundance()

##### Core Analysis Functions #####

unassigned_taxonomy <- function(datalist, tax_grp) {
  
  tax_grp <- enquo(tax_grp)
  
  datalist %>%
    mutate_count_datalist(function(x) ifelse(x > 0, 1, 0))  %>%
    .$Count_Data %>%
    mutate(!!tax_grp := ifelse(is.na(select(., !!tax_grp)), "Unassigned", "Assigned")) %>%
    group_by(!!tax_grp) %>%
    summarize_if(is.numeric, sum) %>%
    t() %>%
    .[-1,] %>%
    as_tibble(rownames = "Sample_ID") %>%
    dplyr::rename("Assigned" = "V1", "Unassigned" = "V2") %>%
    mutate(Assigned = as.numeric(Assigned)) %>%
    mutate(Unassigned = as.numeric(Unassigned)) %>%
    right_join(., datalist$Meta_Data) %>%
    mutate(Unassigned_Rel = Unassigned/(Unassigned + Assigned)) 
  
}

diversity_datatable <- function(datalist) {
  
  diversity <- tibble(Richness = apply(select_if(datalist$Count_Data, is.numeric) > 0, 2, sum),
                      Shannon = apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::diversity),
                      Evenness = Shannon/log(apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::specnumber)),
                      Simpson = apply(select_if(datalist$Count_Data, is.numeric), 2, vegan::diversity, index = "simpson")) %>%
    bind_cols(datalist$Meta_Data, .)
  
  return(diversity)
  
} 

##### Data-wrangling Pipeline #####

Pacific_Diversity <- datalist_Pacific %>%
  rarefy_datalist(., rare_lim = 8000) %>%
  diversity_datatable(.) %>%
  mutate(Size_Fraction = ordered(ifelse(Size_Fraction == 0.22, "FL",
                                 ifelse(Size_Fraction == 3, "SPA", "LPA")), levels = c("FL", "SPA", "LPA")))

Pacific_Diversity_unassigned <- datalist_Pacific_unassigned %>%
  rarefy_datalist(., rare_lim = 8000) %>%
  diversity_datatable(.) 

unassigned_Pacific <- datalist_Pacific_unassigned %>%
  unassigned_taxonomy(tax_grp = Genus) %>%
  mutate(Size_Fraction = ordered(ifelse(Size_Fraction == 0.22, "FL", 
                                 ifelse(Size_Fraction == 3, "SPA", "LPA")), levels = c("FL", "SPA", "LPA"))) %>%
  mutate(Richness = Pacific_Diversity_unassigned$Richness) %>%
  mutate(Shannon = Pacific_Diversity_unassigned$Shannon) %>%
  mutate(Evenness = Pacific_Diversity_unassigned$Evenness)

VD_list_filtered <- datalist_Pacific %>%
  summarize_by_param(Size_Fraction) %>%
  with(., list(FL  = .$Count_Data$OTU_ID[rowSums(select_if(select_if(.$Count_Data, is.numeric), .$Meta_Data$Sample_ID == 0.22)) > 0],
               SPA = .$Count_Data$OTU_ID[rowSums(select_if(select_if(.$Count_Data, is.numeric), .$Meta_Data$Sample_ID == 3)) > 0],
               LPA = .$Count_Data$OTU_ID[rowSums(select_if(select_if(.$Count_Data, is.numeric), .$Meta_Data$Sample_ID == 8)) > 0])
  )

VD_list_unfiltered <- datalist_Pacific_unfiltered %>%
  summarize_by_param(Size_Fraction) %>%
  with(., list(FL  = .$Count_Data$OTU_ID[rowSums(select_if(select_if(.$Count_Data, is.numeric), .$Meta_Data$Sample_ID == 0.22)) > 0],
               SPA = .$Count_Data$OTU_ID[rowSums(select_if(select_if(.$Count_Data, is.numeric), .$Meta_Data$Sample_ID == 3)) > 0],
               LPA = .$Count_Data$OTU_ID[rowSums(select_if(select_if(.$Count_Data, is.numeric), .$Meta_Data$Sample_ID == 8)) > 0])
  )

##### Significant Correlation - Environment with Richness ####

summary(lm(log(Pacific_Diversity$Richness)~poly(Pacific_Diversity$Pot_Temperature, 2)))
summary(lm(log(Pacific_Diversity$Richness)~Pacific_Diversity$Salinity))

ggplot(Pacific_Diversity, aes(x = Salinity, y = log(Richness), col = as.factor(Size_Fraction))) +
  geom_point() +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
  geom_smooth(method = "lm", formula = y~x, se = F)

#### Richness Pattern along Transect ####

ggplot(Pacific_Diversity, aes(x = Latitude, y = log(Richness), col = as.factor(Size_Fraction))) +
  geom_point(size = 2.2) +
  geom_smooth(span = .6, size = 2) +
  facet_wrap(~Depth_Grp_2) +
  labs(col = "Size Fraction") +
  scale_color_discrete(labels = c("FL", "SPA", "LPA")) +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", colour = "black"),
        strip.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.text = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#### Richness Pattern along Depth ####

ggplot(Pacific_Diversity, aes(x = Depth_Grp, y = log(Richness), group = Depth_Grp)) +
  geom_smooth(aes(x = Depth_Grp, y = log(Richness), group = NULL), se = F, method = "lm", formula = y~poly(x, 2)) +
  geom_boxplot(aes(fill = as.factor(Size_Fraction))) +
  coord_flip() +
  ylim(c(4,6)) +
  scale_x_continuous(trans = "reverse", name = "Depth [m]") +
  scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
  facet_wrap(~Size_Fraction) +
  theme(legend.position = "none",
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", colour = "black"),
        strip.text = element_text(size = 12, face = "bold", colour = "black"))

#### Richness correlation with proportion of unassigned taxa ####

unassigned_Pacific %>%
  filter(Size_Fraction == "FL") %>%
  with(., cor.test(.$Unassigned_Rel, log(.$Richness)))

unassigned_Pacific %>%
  filter(Size_Fraction == "SPA") %>%
  with(., cor.test(.$Unassigned_Rel, log(.$Richness)))

unassigned_Pacific %>%
  filter(Size_Fraction == "LPA") %>%
  with(., cor.test(.$Unassigned_Rel, log(.$Richness)))

ggplot(unassigned_Pacific, aes(x = log(Richness), y = Unassigned_Rel*100, col = Size_Fraction)) +
  geom_smooth(method = "lm", se = F, col = "black", size = 1.3) +
  geom_point() +
  facet_wrap(~ Size_Fraction) +
  xlim(c(4,6)) +
  scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
  labs(y = "Proportion of unassigned genera [%]", x = "log(Richness)",col = "Size Fraction")

#### Effect of Abundance Filter - Venn-Diagrams ####

venn.diagram(x = VD_list_filtered,  category.names = c("FL", "SPA", "LPA"),
             filename = "~/PhD/Statistics/Plots/2020_11_Paper_Pacific/VennDiagram_SizeFraction_filtered.png",
             output = F, imagetype="png", height = 800, width = 800, resolution = 200, compression = "lzw",
             lwd = 2, lty = 'blank', fill = c("#F8766D", "#00BA38", "#619CFF"),
             print.mode = "percent", fontfamily = "sans", fontface = "bold", 
             cat.fontfamily = "sans", cat.fontface = "bold")

venn.diagram(x = VD_list_unfiltered,  category.names = c("FL", "SPA", "LPA"),
             filename = "~/PhD/Statistics/Plots/2020_11_Paper_Pacific/VennDiagram_SizeFraction_Unfiltered.png",
             output = F, imagetype="png", height = 800, width = 800, resolution = 200, compression = "lzw",
             lwd = 2, lty = 'blank', fill = c("#F8766D", "#00BA38", "#619CFF"),
             print.mode = "percent", fontfamily = "sans", fontface = "bold", 
             cat.fontfamily = "sans", cat.fontface = "bold")
