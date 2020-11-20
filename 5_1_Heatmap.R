##### Load dependencies #####

rm(list = ls(all=T))

library(tidyverse)
library(pheatmap)
library(gridExtra)

devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datalistFunctions.R")
devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datatableFunctions.R")

##### Read processed tables #####

datalist_Pacific <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok") %>%
  mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
  filter_station_datalist(Counts_Total > 2000) %>%
  mutate_meta_datalist(Depth_Grp = c(20, 40, 60, 100, 200, 300, 500)[findInterval(Depth, c(20, 40, 60, 100, 200, 300, 500))]) %>%
  mutate_meta_datalist(Depth_Grp_2 = ifelse(.$Meta_Data$Depth_Grp < 100, "Epi", "Meso")) %>%
  correct_ambiguous(.) %>%
  filter_abundance()

##### Core Analysis Functions #####

plot_heatmap <- function(datalist, taxa_filter = NULL, envir_filter = NULL, cluster_cols = F, row_label = NULL,
                         scale_method = make_proportion, col_label = NULL, title = NULL, cluster_order = NULL,
                         gaps_col = F) {
  
  taxa_filter <- enquo(taxa_filter) 
  envir_filter <- enquo(envir_filter)
  col_label <- enquo(col_label)
  row_label <- enquo(row_label)
  
  datalist_filtered <- datalist %>%
    with(., if (!rlang::quo_is_null(taxa_filter)) filter_taxa_datalist(., !!taxa_filter) else .) %>%
    with(., if(!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .)
  
  gaps_col = if (gaps_col == T & length(unique(datalist_filtered$Meta_Data$Size_Fraction)) == 1) NULL else {
    cumsum(table(datalist_filtered$Meta_Data$Size_Fraction))
  }
  
  Count_Matrix <- datalist_filtered %>%
    .$Count_Data %>%
    select_if(is.numeric) %>%
    select(datalist_filtered$Meta_Data %>% 
             with(., 
                  if (!rlang::quo_is_null(col_label)) {
                    if ("Size_Fraction" %in% names(datalist$Meta_Data)) {
                      arrange(., Size_Fraction, !!col_label)
                    } else {
                      arrange(., !!col_label)
                    }
                  } else {
                    if ("Size_Fraction" %in% names(datalist$Meta_Data)) {
                      arrange(., Size_Fraction, Latitude)
                    } else {
                      arrange(., Latitude)
                    }
                  }) %>% .$Sample_ID) %>%
    as.matrix() %>%
    apply(., 1, function(x) scale_method(x)) %>% 
    t()
  
  rownames(Count_Matrix) <- datalist_filtered$Count_Data %>%
    select(OTU_ID) %>%
    deframe()
  
  row_label <- if(!rlang::quo_is_null(row_label)) {
    datalist_filtered$Count_Data %>%
      select(!!row_label) %>%
      deframe()
  } else NULL
  
  col_label <- if (!rlang::quo_is_null(col_label)) {
    
    if ("Size_Fraction" %in% names(datalist$Meta_Data)) {
      datalist_filtered$Meta_Data %>%
        arrange(Size_Fraction, !!col_label) %>%
        select(!!col_label) %>%
        mutate(!!col_label := round(!!col_label, digits = 0)) %>%
        deframe() %>%
        as.character()
    } else {
      datalist_filtered$Meta_Data %>%
        arrange(!!col_label) %>%
        mutate(!!col_label := if (is.numeric(!!col_label)) round(!!col_label, digits = 0) else !!col_label) %>%
        select(!!col_label) %>%
        deframe() %>%
        as.character()
    }
    
  } else NULL
  
  clust_num <- factoextra::fviz_nbclust(Count_Matrix, cluster::pam, method = "silhouette", k.max = min(nrow(Count_Matrix)-1, 8))$data$y %>%
    which.max()
  
  cluster_obj <- Count_Matrix %>%
    vegan::vegdist(method = "euclidian") %>%
    hclust(., method = "ward.D2") %>%
    with(., if ("Latitude" %in% names(datalist_filtered$Meta_Data)) {
      reorder(., wts = abs(datalist_filtered$Meta_Data$Latitude), agglo.FUN = "uwmean")
    } else .)
  
  row_annotation <- data.frame(Cluster = paste("Cluster", cutree(cluster_obj, k = clust_num)), row.names = rownames(Count_Matrix))
  row_annotation_color <- list(Cluster = structure(ggsci::pal_rickandmorty(palette = c("schwifty"))(clust_num), names = unique(row_annotation$Cluster)))
  
  if (is.null(title)) {
    title <- paste(str_split_fixed(as_label(taxa_filter), pattern = "==", 2)[1] %>% str_replace(., " ", ""), ": ",
                   str_split_fixed(as_label(taxa_filter), pattern = "==", 2)[2] %>% str_replace_all(., "\\\"", "") %>% str_replace(., " ", ""), sep = "")
  }
  
  pheatmap(Count_Matrix, cluster_cols = cluster_cols, cluster_rows = cluster_obj, cutree_rows = clust_num, 
           main = title, border_color = "grey60",
           legend_breaks = c(min(Count_Matrix), max(Count_Matrix)),legend_labels = c("Min  ", "Max  "),
           angle_col = "90", annotation_legend = F, labels_col = col_label,
           labels_row = row_label, show_rownames = ifelse(!is.null(row_label),T,F),
           annotation_row = row_annotation, fontsize_row = 6, fontsize_col = 7,
           annotation_colors = row_annotation_color, fontsize = 7,
           gaps_col = gaps_col, 
           clustering_callback = callback)
  
}

##### Data-wrangling Pipeline #####

Latitude_Resolution = 12

datalist_Pacific_Summarized <- datalist_Pacific %>%
  mutate_meta_datalist(Latitude_Interval = round(seq(min(Latitude),max(Latitude), length.out = Latitude_Resolution)[findInterval(Latitude, seq(min(datalist_Pacific$Meta_Data$Latitude),max(datalist_Pacific$Meta_Data$Latitude), length.out = Latitude_Resolution))])) %>%
  mutate_meta_datalist(Latitude = round(Latitude, digits = 1)) %>%
  summarize_by_param(Latitude, Depth_Grp_2, Size_Fraction)

#### Create heatmaps ####

plot_heatmap(datalist_Pacific_Summarized, taxa_filter = Family == "Clade I",
             envir_filter = Depth_Grp_2 == "Epi", col_label = Latitude, row_label = Genus,
             scale_method = function(x) (x-mean(x))/sd(x), gaps_col = T,
             title = "Family: SAR11 Clade I - Depth: Epipelagic\nFL-Community                     SPA-Community                     LPA-Community")
