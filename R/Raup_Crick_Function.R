#### Raup-Crick-Metric by Stegen et al. 2013 (ISME) ####

library(picante)
library(parallel)
library(doParallel)
library(tidyverse)
library(ape)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

raup_crick_abundance = function(datalist, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, 
                                as.distance.matrix=TRUE, report_similarity=FALSE, use.cores = 4){
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  
  # By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 
  # instead of 0 to 1). 
  
  # classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. 
  
  # split_ties (defaults to TRUE) adds half of the number of null observations that are equal 
  # to the observed number of shared species to the calculation- this is highly recommended.  
  
  # report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate 
  # as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick 
  # originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity 
  # (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges 
  # from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). 
  # If ties are not split (and there are ties between the observed and expected shared number of species) 
  # this conversion will not work. 
  
  # reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  
  
  # set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their 
  # occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species 
  # is drawn from the null model	
  
  registerDoParallel(cores = use.cores)
  
  spXsite <- datalist %>%
    .$Count_Data %>%
    select_if(is.numeric) %>%
    as.data.frame() %>%
    magrittr::set_rownames(datalist$Count_Data$OTU_ID) %>%
    t()
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites <- nrow(spXsite)
  gamma <- ncol(spXsite)
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur <- apply(spXsite > 0, MARGIN = 2, FUN = sum)
  
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance <- apply(spXsite, MARGIN = 2, FUN = sum)
  
  ##make_null
  com_shuffle <- matrix(0, nrow = n_sites, ncol = gamma)
  seq_gamma <- 1:gamma
  
  richness_per_sample <- apply(spXsite > 0, MARGIN = 1, FUN = sum)
  counts_per_sample <- apply(spXsite, MARGIN = 1, FUN = sum)
  
  start_time <- Sys.time()
  
  null_bray_curtis <- foreach(i = 1:reps) %dopar% {
    
    # seq_gamma = vector mit ASV indizes
    # spXsite.inc = CountTable mit absence/presence (samples in rows)
    # occur = count-sums für jeden ASV
        
    # step1: zieh zufaellige ASVs mit anzahl = tatsaechliche Anzahl an ASVs in der Probe. Jeder ASV hat Prob = occur
    # dann setze diese ASVs in unserer kuenstlichen Count-Table = 1
        
    shuffled_com <- map(1:n_sites, function(x) {
      com_shuffle_tmp <- com_shuffle
      ASV_inds <- sample(seq_gamma, richness_per_sample[x], replace = F, prob = occur)        
      counts <- sample(ASV_inds, (counts_per_sample[x] - richness_per_sample[x]), replace = T, prob = abundance[ASV_inds]) %>%
        table() %>%
        as.data.frame(stringsAsFactors = F) %>%
        dplyr::rename("ASVs" = ".") %>%
        mutate(ASVs = as.numeric(ASVs))
      com_shuffle_tmp[x,counts$ASVs] <- counts$Freq
      com_shuffle_tmp[x,]
    }) %>%
    cbind.data.frame() %>%
    magrittr::set_colnames(1:n_sites) %>%
    t()
        
    vegan::vegdist(shuffled_com) %>%
      as.matrix()
  }
  
  end_time <- Sys.time()
  cat((end_time - start_time)/reps)
  
  obs_bray_curtis <- vegan::vegdist(spXsite) %>%
    as.matrix()
      
  ##how many null observations is the observed value tied with?
  num_exact_matching_in_null = lapply(null_bray_curtis, function(x) x == obs_bray_curtis) %>%
      reduce(`+`)
  
  ##how many null values are smaller than the observed *dissimilarity*?
  num_less_than_in_null =  lapply(null_bray_curtis, function(x) x < obs_bray_curtis) %>%
    reduce(`+`)
      
  if (split_ties) {
    rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
  } else {
    rc = (num_less_than_in_null )/reps
  }  
      
  if(!classic_metric){
    ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
    rc = (rc-.5)*2
  }
  
  rc <- rc %>%
    magrittr::set_rownames(rownames(obs_bray_curtis)) %>%
    magrittr::set_colnames(colnames(obs_bray_curtis))
  
  if (as.distance.matrix){ 
    results <- as.dist(round(rc,digits=2))
  }	else {
    results <- round(rc,digits=2)
  }
  
  return(results)
  
}

bootstrap_reps = 9999
parallel_cores = 6

datalist <- import_data("data/", kingdom = "Euk", min_counts = 500, abundance_filter = T) #%>%
  #filter_taxa_datalist(Kingdom == "Archaea")
  
#rc_FL <- datalist %>% 
#  filter_station_datalist(colSums(select_if(datalist$Count_Data, is.numeric) > 0) >= 2) %>%
#  filter_station_datalist(Size_Fraction == 0.22) %>%
#  
#  raup_crick_abundance(., classic_metric = F, split_ties = T, reps = bootstrap_reps, set_all_species_equal = F,
#                       as.distance.matrix = F, report_similarity = F, use.cores = parallel_cores)
#
#write.csv(rc_FL, "output/Raup_Crick_Archaea_FL.csv")

rc_SPA <- datalist %>% 
  filter_station_datalist(colSums(select_if(datalist$Count_Data, is.numeric) > 0) >= 2) %>%
  filter_station_datalist(Size_Fraction == 3) %>%
  raup_crick_abundance(., classic_metric = F, split_ties = T, reps = bootstrap_reps, set_all_species_equal = F,
                       as.distance.matrix = F, report_similarity = F, use.cores = parallel_cores)

write.csv(rc_SPA, "output/Raup_Crick_Eukaryotes_SPA.csv")

rc_LPA <- datalist %>% 
  filter_station_datalist(colSums(select_if(datalist$Count_Data, is.numeric) > 0) >= 2) %>%
  filter_station_datalist(Size_Fraction == 8) %>%
  raup_crick_abundance(., classic_metric = F, split_ties = T, reps = bootstrap_reps, set_all_species_equal = F,
                       as.distance.matrix = F, report_similarity = F, use.cores = parallel_cores)

write.csv(rc_LPA, "output/Raup_Crick_Eukaryotes_LPA.csv")
