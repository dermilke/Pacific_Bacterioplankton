#### Barplots & NMDS ####

#### Load dependencies ####

rm(list = ls(all=T))

library(tidyverse)

devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datalistFunctions.R")
devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datatableFunctions.R")

#### Read processed tables ####

datalist_Pacific <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok") %>%
  mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
  filter_station_datalist(Counts_Total > 2000) %>%
  mutate_meta_datalist(Depth_Grp = c(20, 40, 60, 100, 200, 300, 500)[findInterval(Depth, c(20, 40, 60, 100, 200, 300, 500))]) %>%
  mutate_meta_datalist(Depth_Grp_2 = ifelse(Depth < 100, "Epi", "Meso")) %>%
  correct_ambiguous(.) %>%
  filter_abundance(.)

Province_Colour <- read_csv("https://raw.githubusercontent.com/dermilke/Pacific_Bacterioplankton/master/Colors/Province_Colour.csv")

##### Data-wrangling Pipeline #####

datatable_Pacific_summarized <- datalist_Pacific %>%
  summarize_by_param(Province, Depth_Grp_2, Size_Fraction) %>%
  create_datatable(., grpBy = Family, otherThreshold = 0.008, addColorScheme = T) 

colorScheme <- datatable_Pacific_summarized$color

datatable_Pacific <- datatable_Pacific_summarized$table %>%
  mutate(Group = plyr::mapvalues(Group, levels(Group), str_replace(levels(Group), pattern = ";", replacement = " - "))) %>%
  mutate(Province = ordered(Province, levels = unique(arrange(datalist_Pacific$Meta_Data, Latitude)$Province)))

NMDS_list <- datalist_Pacific %>%
  summarize_by_param(Latitude, Depth_Grp_2, Size_Fraction) %>%
  mutate_meta_datalist(Province = {left_join(.$Meta_Data, distinct(select(datalist_Pacific$Meta_Data, Province, Latitude)), by = "Latitude")$Province}) %>%
  mutate_meta_datalist(Province = ordered(Province, levels = unique(arrange(.$Meta_Data, Latitude)$Province))) %>%
  mutate_count_datalist(make_proportion) %>%
  mutate_count_datalist(function(x) sqrt(x)) %>%
  sizefraction_communities() %>%
  purrr::map(NMDS_ordination_datalist)

NMDS_data <- NMDS_list %>%
  purrr::map("table") %>%
  bind_rows() %>%
  mutate(Size_Fraction = ordered(ifelse(Size_Fraction == 0.22, paste("FL - Stress:", round(NMDS_list$Fraction_022$stress, digits = 2)), 
                                        ifelse(Size_Fraction == 3, paste("SPA - Stress:", round(NMDS_list$Fraction_3$stress, digits = 2)),
                                               paste("LPA - Stress:", round(NMDS_list$Fraction_8$stress, digits = 2)))), 
                                 levels = c(paste("FL - Stress:", round(NMDS_list$Fraction_022$stress, digits = 2)), 
                                            paste("SPA - Stress:", round(NMDS_list$Fraction_3$stress, digits = 2)), 
                                            paste("LPA - Stress:", round(NMDS_list$Fraction_8$stress, digits = 2)))))

#### Barplot of Provinces ####

ggplot(datatable_Pacific, aes(x = Province, y = Abundance*100, fill = Group)) +
  geom_col() +
  scale_fill_manual(values = colorScheme) +
  facet_grid(Depth_Grp_2~Size_Fraction) +
  labs(fill = "Class - Family", x = "Latitude", y = "Abundance [%]") +
  theme(legend.position = "right",
       legend.text = element_text(size = 8),
        legend.key.size = unit(.8, "line"),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  guides(fill = guide_legend(nrow = length(unique(datatable_Pacific$Group)), 
                             title.position = "top"))

#### NMDS Biplots ####

ggplot(NMDS_data, aes(x = MDS1, y = MDS2, fill = Province, shape = Depth_Grp_2, size = Depth_Grp_2)) +
  geom_point() +
  scale_fill_manual(values = Province_Colour$Colour) +
  scale_shape_manual(values = c(21, 24), name = "Depth") +
  scale_size_manual(values = c(4,3)) +
  facet_wrap(~Size_Fraction) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
         shape = guide_legend(override.aes = list(fill = "black", size = 4)),
         size = FALSE) +
  theme(panel.spacing.x = unit(1.5, "lines"))
