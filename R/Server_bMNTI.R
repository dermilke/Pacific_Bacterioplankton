library(picante)
library(parallel)
library(doParallel)
library(tidyverse)
library(ape)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

get_weighted_bNTI <- function(datalist, tree, boot_num = 999, core_num = 4) {
  
  require(picante)
  require(doParallel)
  
  registerDoParallel(cores = core_num)
  
  otu <- datalist$Count_Data %>%
    select_if(is.numeric) %>%
    as.data.frame() %>%
    magrittr::set_rownames(datalist$Count_Data$OTU_ID)
  
  match.phylo.otu <- suppressMessages(picante::match.phylo.data(tree, otu))
  
  cophenetic_phylo <- cophenetic(match.phylo.otu$phy)
  trans_otus <- t(match.phylo.otu$data)
  
  beta.mntd.weighted <- picante::comdistnt(trans_otus, cophenetic_phylo, abundance.weighted = T) %>%
    as.matrix()
  
  rand.weighted.bMNTD.comp <- parallel::mclapply(seq(1, boot_num), function(x) {
    tmp <- picante::comdistnt(trans_otus, taxaShuffle(cophenetic_phylo), abundance.weighted = T) %>%
      as.matrix()
    
    cat(c(date(), x, "\n"))
    
    return(tmp)
  }, mc.cores = core_num) %>%
    unlist() %>%
    array(., dim = c(ncol(match.phylo.otu$data), ncol(match.phylo.otu$data), boot_num))
  
  weighted.bNTI <- (beta.mntd.weighted - apply(rand.weighted.bMNTD.comp, c(1,2), mean)) / apply(rand.weighted.bMNTD.comp, c(1,2), sd)
  
  diag(weighted.bNTI) <- 0
  
  rownames(weighted.bNTI) = colnames(match.phylo.otu$data)
  colnames(weighted.bNTI) = colnames(match.phylo.otu$data)
  
  return(weighted.bNTI)
  
}

use_cores <- 6
beta.reps <- 999

my_tree_18S <- ape::read.tree("data/Eukaryotes_SSUalign_FastTree.tree")
my_tree <- ape::read.tree("data/Prokaryotes_SINA_FastTree.tree")

datalist <- import_data("data/", kingdom = "Prok", min_counts = 2000, abundance_filter = T) %>%
  filter_taxa_datalist(Family != "Mitochondria")
datalist_18S <- import_data("data/", kingdom = "Euk", min_counts = 500, abundance_filter = T) %>%
  filter_taxa_datalist(Phylum != "Metazoa")

weighted_bNTI_FL <- datalist %>%
  filter_station_datalist(Size_Fraction == 0.22) %>%
  get_weighted_bNTI(., my_tree, boot_num = beta.reps, core_num = use_cores)

write.csv(weighted_bNTI_FL,"output/Prokaryotes_weighted_bNTI_FL.csv")

weighted_bNTI_SPA <- datalist %>%
  filter_station_datalist(Size_Fraction == 3) %>%
  get_weighted_bNTI(., my_tree, boot_num = beta.reps, core_num = use_cores)

write.csv(weighted_bNTI_SPA,"output/Prokaryotes_weighted_bNTI_SPA.csv")

weighted_bNTI_LPA <- datalist %>%
  filter_station_datalist(Size_Fraction == 8) %>%
  get_weighted_bNTI(., my_tree, boot_num = beta.reps, core_num = use_cores)

write.csv(weighted_bNTI_LPA,"output/Prokaryotes_weighted_bNTI_LPA.csv")

