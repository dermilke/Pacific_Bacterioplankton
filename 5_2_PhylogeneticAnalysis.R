##### Load dependencies #####

rm(list = ls(all=T))

library(tidyverse)
library(msa)
library(phyloseq)

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

phylogenetic_distance <- function(datalist, fasta_location = NULL) {
  
  faster_readLines = function(fname) {
    s <- file.info(fname)$size 
    buf <- readChar(fname, s, useBytes=T)
    strsplit(buf, "\n",fixed=T, useBytes=T)[[1]]
  }
  
  subset <- datalist$Count_Data$OTU_ID
  Sequence <- NULL
  
  for (i in 1:length(subset)) {
    tmp_index <- grep(subset[i], faster_readLines(fasta_location))
    Sequence <- c(Sequence, readLines(fasta_location)[tmp_index + 1])
  }
  
  tibble(OTU_ID = subset, Sequence = Sequence) %>%
    mutate(OTU_ID = paste(">", OTU_ID, sep = "")) %>%
    with(., do.call(rbind, lapply(seq(nrow(.)), function(i) t(.[i, ])))) %>%
    write.table(., "tmp.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  DNA_set <- Biostrings::readDNAStringSet("tmp.fasta", format = "fasta")
  DNA_align <- msa(DNA_set, method = "ClustalW") %>%
    msaConvert(type="seqinr::alignment")
  
  Phylo_Distances <- seqinr::dist.alignment(DNA_align)
  
  Count_Matrix <- datalist$Count_Data %>%
    dplyr::slice(match(attr(Phylo_Distances, "Labels"), datalist$Count_Data$OTU_ID)) %>%
    select_if(is.numeric) %>%
    data.frame(row.names = attr(Phylo_Distances, "Labels"), check.names = F) %>%
    as.matrix()
  
  Taxo_Distances <- Count_Matrix %>%
    apply(., 1, function(x) (x-mean(x))/sd(x)) %>%
    t() %>%
    vegan::vegdist(method = "euclidean")
  
  Tax_Matrix <-  datalist$Count_Data %>%
    dplyr::slice(match(attr(Phylo_Distances, "Labels"), datalist$Count_Data$OTU_ID)) %>%
    dplyr::rename("Domain" = "Kingdom") %>%
    select_if(is.character) %>%
    select(-OTU_ID) %>%
    data.frame(row.names = attr(Phylo_Distances, "Labels")) %>%
    as.matrix()
  
  Meta_Matrix <- datalist$Meta_Data %>%
    data.frame(row.names = .$Sample_ID)
  
  Unifrac_Distances <- phyloseq(otu_table(Count_Matrix, taxa_are_rows = T), 
                                tax_table(Tax_Matrix), 
                                sample_data(Meta_Matrix), 
                                refseq(DNA_set), 
                                ape::nj(Phylo_Distances)) %>%
    UniFrac(., weighted=F, normalized=TRUE, parallel=TRUE, fast=TRUE)
  
  return(list(Phylogenetic_Distance = Phylo_Distances,
              Unifrac_Distance = Unifrac_Distances))
  
}

##### Data-wrangling Pipeline #####

phylo_dist <- datalist_Pacific %>%
  filter_taxa_datalist(Family == "Clade I") %>%
  phylogenetic_distance(., fasta_location = "~/PhD/Data_Storage/Fasta/Pacific/V4V5_Primerset/Prok/Sequences_Pacific_Complete.fasta")

phy_DD_table <- reshape2::melt(as.matrix(phylo_dist$Unifrac_Distance)) %>%
  filter(!is.nan(value)) %>%
  dplyr::rename( "From" = "Var1", "Sample_ID" = "Var2", "Unifrac" = "value") %>%
  left_join(., datalist_Pacific$Meta_Data %>%
              select(Sample_ID, Latitude, Longitude, Province, Depth_Grp, Size_Fraction), by = "Sample_ID") %>%
  dplyr::rename("To" = "Sample_ID", "To_Latitude" = "Latitude", "To_Province" = "Province" , 
                "To_Longitude"  = "Longitude", "Sample_ID" = "From", "To_Depth" = "Depth_Grp", "To_Size_Fraction" = "Size_Fraction") %>%
  left_join(., datalist_Pacific$Meta_Data %>%
              select(Sample_ID, Latitude, Longitude, Province, Depth_Grp, Size_Fraction), by = "Sample_ID") %>%
  dplyr::rename("From" = "Sample_ID", "From_Latitude" = "Latitude", "From_Longitude" = "Longitude", 
                "From_Province" = "Province", "From_Depth" = "Depth_Grp", "From_Size_Fraction" = "Size_Fraction") %>%
  mutate(Distance = geographic_distance(To_Longitude, To_Latitude, From_Longitude, From_Latitude)) %>%
  as_tibble() %>%
  mutate(From_Depth = ifelse(From_Depth < 100, "Epi", "Meso")) %>%
  mutate(To_Depth = ifelse(To_Depth < 100, "Epi", "Meso")) %>%
  filter(abs(From_Latitude) == min(abs(From_Latitude))) %>%
  filter(From_Depth == To_Depth) %>%
  filter(From_Size_Fraction == To_Size_Fraction) %>%
  filter(Unifrac > 0 & Unifrac < 1)
  
#### Phylogenetic divergence from equatorial communities ####  
  
phy_DD_table %>%
  filter(From_Size_Fraction == 0.22) %>%
  
  ggplot(., aes(x = To_Latitude, y = Unifrac)) +
    geom_point() +
    lims(y = c(0,1)) +
    geom_smooth(formula = y~poly(x, 2), method = lm) +
    labs(title = "Phylogenetic distances between SAR11 Clade I Communities", y = "unweighted UniFrac-Distance", x = "Latitude",
         subtitle = "Distances calculated from equatorial communities") +
    facet_wrap(~From_Depth)
