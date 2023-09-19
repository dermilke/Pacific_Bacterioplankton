#### Import Function ####
# Start of every Analysis Pipeline
# 
# Define:
# file_ASV = Location of folder structure (top level of substructure -> see Github ReadMe)
# kingdom = Prokaryotes, Chloroplasts, Eukaryotes (Prok, Chloroplast, Euk)
# rare_lim = Integer defining the rarefying level. If NULL no rarefying will be done
# drop_rare = Logical defining if samples below rare_lim will be dropped from analysis
# abundance_filter = Logical defining if abundance filter after Milici et al. 2016 should be applied
# min_count = Integer defining the minimum count number a sample should have. NULL for no filtering

import_data <- function(file_ASV, kingdom = "Prok", rare_lim = NULL, drop_rare = T, 
                        abundance_filter = F, min_counts = NULL) {
  
  data_select <- function(file_ASV, kingdom = "Prok") {
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Processed/", kingdom, "/Full_", kingdom,"_Count.tsv", sep = ""), 
                                              del = "\t")) %>%
      select(-grep("[Mm]ock", names(.))) %>%
      select(-grep("NC", names(.))) 
    
    Meta_Data <- suppressMessages(read_delim(paste(file_ASV, "Meta_Data/", kingdom, "/Meta_Data.tsv", sep = ""), 
                                             del = "\t"))
    
    return(list(Meta_Data = Meta_Data, 
                Count_Data = Count_Data))
    
  }
  
  filter_abundance <- function(datalist) {
    
    counts <- datalist$Count_Data
    prop <- datalist$Count_Data %>%
      mutate_if(is.numeric,
                function(x) x/sum(x))
    
    counts_filtered <- counts %>%
      filter(((rowSums(select_if(counts, is.numeric))/sum(rowSums(select_if(counts, is.numeric)))) > 0.00001) &
               (apply(select_if(prop, is.numeric),1,max) > 0.01) |
               ((apply(select_if(prop, is.numeric) > 0.001, 1, sum) / ncol(select_if(prop, is.numeric))) > 0.02) |
               ((apply(select_if(prop, is.numeric) > 0, 1, sum) / ncol(select_if(prop, is.numeric))) > 0.05))
    
    datalist$Count_Data <- counts_filtered
    
    return(datalist)
    
  }
  
  rarefy_datalist <- function(datalist, rare_lim, drop = F) {
    
    count_rared <- datalist$Count_Data %>%
      select_if(is.numeric) %>%
      select_if(colSums(.) >= ifelse(drop, rare_lim, 0)) %>%
      t() %>% 
      vegan::rrarefy(., rare_lim) %>%
      t() %>%
      as_tibble() %>%
      bind_cols(select_if(datalist$Count_Data, is.character), .) %>%
      filter(rowSums(select_if(., is.numeric)) > 0)
    
    meta_subset <- datalist$Meta_Data %>%
      dplyr::slice(match(names(select_if(count_rared, is.numeric)), datalist$Meta_Data$Sample_ID))
    
    datalist$Count_Data <- count_rared
    datalist$Meta_Data <- meta_subset
    
    return(datalist)
    
  }
  
  correct_ambiguous <- function(datalist, fromTaxLvl = 8) {
    
    replacer <- function(Count_Data, taxLvl, replaceLvl, pattern) {
      Count_Data[grep(pattern, x = as_vector(Count_Data[, taxLvl])), taxLvl] <- paste0("Unknown ", as_vector(Count_Data[grep(pattern, x = as_vector(Count_Data[, taxLvl])), replaceLvl]))
      return(Count_Data)
    }
    
    for (taxLvl in fromTaxLvl:1) {
      
      for (i in (taxLvl-1):1) {
        
        datalist$Count_Data <- replacer(datalist$Count_Data, taxLvl, replaceLvl = i, pattern = 'bacterium$') %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "uncultured") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "metagenome") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "unidentified") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "Ambiguous_taxa") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "unknown") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "unidentified marine bacterioplankton") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "^Unknown.*bacterium$") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "^Uncultured.*bacterium$") 
        
      }
      
    }  
    
    return(datalist)
    
  }
  
  data_import <- data_select(file_ASV, kingdom = kingdom) %>%
    mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
    with(., if (!is.null(min_counts)) filter_station_datalist(., Counts_Total > !!min_counts) else .) %>%
    with(., if (!is.null(rare_lim)) rarefy_datalist(., rare_lim, drop_rare) else .) %>%
    with(., if (abundance_filter) filter_abundance(.) else .) %>%
    correct_ambiguous()
  
  return(data_import)
  
}
