library(picante)
library(parallel)
library(doParallel)
library(tidyverse)
library(ape)

source("R/Datalist_Wrangling_Functions.R")
source("R/Import_Data.R")

get_unweighted_bNTI <- function(datalist, tree, boot_num = 999, core_num = 4) {
  
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
  
  beta.mntd.weighted <- picante::comdistnt(trans_otus, cophenetic_phylo, abundance.weighted = F) %>%
    as.matrix()
  
  rand.weighted.bMNTD.comp <- parallel::mclapply(seq(1, boot_num), function(x) {
    tmp <- picante::comdistnt(trans_otus, taxaShuffle(cophenetic_phylo), abundance.weighted = F) %>%
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

use_cores <- 9
beta.reps <- 999

my_tree <- ape::read.tree("output/FastTree_Tree/Prok_Combined.tree")

cluster <- read_csv("output/SparCC/SparCC_Cluster.csv") %>%
  select(OTU_ID, Cluster) %>%
  mutate(Abundance = 1) %>%
  reshape2::dcast(OTU_ID ~ Cluster, fill = 0) %>%
  magrittr::set_colnames(c("OTU_ID", as.character(1:10)))

datalist_cluster <- list(Count_Data = cluster, Meta_Data = NULL)

weighted_bNTI_SparCC_Cluster <- datalist_cluster %>%
  get_unweighted_bNTI(., my_tree, boot_num = beta.reps, core_num = use_cores)

write.csv(weighted_bNTI_SparCC_Cluster,"output/SparCC_Cluster_bMNTI.csv")

bNTI_mod <- weighted_bNTI_SparCC_Cluster %>%
  magrittr::set_colnames(rownames(.)) %>%
  as.matrix() 

bNTI_mod[lower.tri(bNTI_mod)] <- NA
bNTI_mod[diag(bNTI_mod)] <- NA

bNTI_mod[abs(bNTI_mod) < 2] <- NA
bNTI_mod[abs(bNTI_mod) >= 2] <- 1

pheatmap::pheatmap(bNTI_mod, cluster_cols = F, cluster_rows = F, kmeans_k = NA)

detailled <- bNTI_mod %>%
  reshape2::melt() %>%
  magrittr::set_colnames(c("Cluster", "To_Cluster", "bNTI")) %>%
  filter(!is.na(bNTI)) %>%
  mutate(Mechanism = ifelse(bNTI > 2, "Heterogeneous Selection", 
                            ifelse(bNTI < -2, "Homogeneous Selection", "Not significant"))) %>%
  filter(Cluster != To_Cluster) %>%
  mutate(bNTI_mod = ifelse(abs(bNTI) < 2, 0, bNTI))
