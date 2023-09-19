#### Datalist Wrangling Functions ####

filter_station_datalist <- function(.datalist, ...,  removeEmpty = T) {
  
  exprFilter <- enexprs(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    filter(!!! exprFilter)
  
  .datalist$Count_Data <- .datalist$Count_Data %>%
    select_if(is.numeric) %>%
    select_if(names(.) %in% .datalist$Meta_Data$Sample_ID) %>%
    bind_cols(select_if(.datalist$Count_Data, is.character), .) %>%
    with(., if(removeEmpty) filter(., {select_if(., is.numeric) %>% rowSums(.)} > 0) else .)
  
  return(.datalist)
  
}

dissimilarity_datalist <- function(datalist, dist_mat = NULL, method = "Bray_Curtis") {
  
  geographic_distance <- function(long1, lat1, long2, lat2) {
    
    dlong <- long1 * pi/180 - long2 * pi/180
    dlat <- lat1 * pi/180 - lat2 * pi/180
    
    R <- 6371
    a <- sin(dlat/2)^2 + cos(lat1 * pi/180) * cos(lat2 * pi/180) * sin(dlong/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    
    distance <- R * c
    
    return(distance)
    
  }
  
  if (is_null(dist_mat)) {
    
    source("R/Functions_Similarity_Indices.R")
    
    Sample_Data <- datalist$Meta_Data %>%
      as.data.frame() %>%
      magrittr::set_rownames(datalist$Meta_Data$Sample_ID)
    
    Count_Table <- datalist$Count_Data %>%
      select_if(is.numeric) %>%
      as.matrix() %>%
      Matrix::Matrix() %>%
      magrittr::set_rownames(datalist$Count_Data$OTU_ID)
    
    my.ps <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(Count_Table), taxa_are_rows=T), 
                                phyloseq::sample_data(Sample_Data))
    
    dist_mat <- distance_wrapper(my.ps, method = method, size.thresh = 1, pseudocount = 10^-6, nblocks = 50, use.cores = 4, 
                                 cor.use = "na.or.complete") %>%
      magrittr::set_rownames(datalist$Count_Data$OTU_ID) %>%
      magrittr::set_colnames(datalist$Count_Data$OTU_ID)
  }
  
  result <- dist_mat %>%
    as.matrix() %>%
    reshape2::melt() %>%
    magrittr::set_colnames(c("From", "Sample_ID", method)) %>%
    left_join(., datalist$Meta_Data, by = "Sample_ID") %>%
    magrittr::set_names(c("Sample_ID", "To", "Similarity", paste0("To_", names(.)[-c(1:3)]))) %>%
    left_join(., datalist$Meta_Data, by = "Sample_ID") %>%
    magrittr::set_names(c("From", "To", "Similarity", 
                          names(.)[grepl(pattern = "^To_", x = names(.))],
                          paste0("From_", names(.)[!grepl(pattern = "^To_", x = names(.))])[-c(1:3)])) %>%
    mutate(Distance = geographic_distance(To_Longitude, To_Latitude, From_Longitude, From_Latitude)) %>%
    as_tibble()
  
  return(result)
  
}

filter_taxa_datalist <- function(.datalist, ...) {
  
  exprFilter <- enquos(...)
  
  .datalist$Count_Data <- filter(.datalist$Count_Data, !!!(exprFilter))
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    mutate(Counts_Total = colSums(select_if(.datalist$Count_Data, is.numeric))) 
  
  .datalist$Count_Data <- .datalist$Count_Data %>%
    select_if(is.numeric) %>%
    select_if(.datalist$Meta_Data$Counts_Total > 0) %>%
    with(., bind_cols(select_if(.datalist$Count_Data, is.character), .))
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    filter(Sample_ID %in% names(.datalist$Count_Data))
  
  return(.datalist)
  
}

singleton_filter <- function(datalist, min_count = 5, min_station = 2) {
  
  Count_Data <- datalist$Count_Data %>%
    filter(select_if(., is.numeric) %>% rowSums(.>0) >= min_station) %>%
    filter(select_if(., is.numeric) %>% rowSums() >= min_count)
  
  datalist$Count_Data <- Count_Data
  
  return(datalist)
  
}

slice_station_datalist <- function(.datalist, ...) {
  
  exprSlice <- enexprs(...)
  
  .datalist$Count_Data <- .datalist$Count_Data %>%
    select_if(is.numeric) %>%
    select(!!! exprSlice) %>%
    bind_cols(select_if(.datalist$Count_Data, is.character), .) %>%
    filter(select_if(., is.numeric) %>% rowSums(.) > 0)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    dplyr::slice(!!! exprSlice)
  
  return(.datalist)
  
}

select_meta_datalist <- function(.datalist, ...) {
  
  exprFunc <- enquos(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    select(!!!exprFunc)
  
  return(.datalist)
  
}

mutate_count_datalist <- function(.datalist, func, ...) {
  
  exprFunc <- rlang::enexpr(func)
  arguments <- rlang::list2(...)
  
  .datalist$Count_Data <- .datalist$Count_Data %>%
    select_if(is.numeric) %>%
    mutate_all(eval(exprFunc), !!! arguments) %>%
    cbind(select_if(.datalist$Count_Data, is.character), .) %>%
    as_tibble()
  
  return(.datalist)
  
}

mutate_meta_datalist <- function(.datalist, ...) {
  
  exprFunc <- enquos(...)
  
  .datalist$Meta_Data <- .datalist$Meta_Data %>%
    mutate(!!!exprFunc)
  
  return(.datalist)
  
}

summarize_by_taxa <- function(datalist, tax_lvl = Species) {
  
  tax_lvl <- enquo(tax_lvl)
  
  if (is.name(rlang::quo_get_expr(tax_lvl))) {
    
    datalist$Count_Data <- datalist$Count_Data %>%
      group_by(!! tax_lvl) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup()
    
  } else {
    
    datalist$Count_Data <- datalist$Count_Data %>%
      group_by_at(2:(tax_lvl+1)) %>%
      summarize_if(is.numeric, sum) %>%
      ungroup()
    
  }
  
  return(datalist)
  
}

arrange_datalist <- function(datalist, ...) {
  
  param <- enquos(...)
  
  datalist$Meta_Data <- datalist$Meta_Data %>%
    dplyr::arrange(!!!param)
  
  datalist$Count_Data <- bind_cols(select_if(datalist$Count_Data, is.character),
                                   select(datalist$Count_Data, datalist$Meta_Data$Sample_ID))
  
  return(datalist)
  
}

summarize_by_param <- function(datalist, ...) {
  
  param <- enquos(...)
  
  datalist <- datalist %>%
    mutate_count_datalist(function(x) x/sum(x))
  
  tmp_final <- select_if(datalist$Count_Data, is.character)
  
  tmp_obj <- unique(select(datalist$Meta_Data, !!!param)) %>%
    mutate_if(is.ordered, as.character)
  
  for (j in 1:nrow(tmp_obj)) {
    
    tmp <- datalist 
    
    for (i in 1:ncol(tmp_obj)) {
      
      tmp <- tmp %>%
        filter_station_datalist(!!param[[i]] == !!deframe(tmp_obj[j,i]), removeEmpty = F)
      
    }
    
    tmp_final <- tmp_final %>%
      mutate(!! paste(tmp_obj[j,], collapse = "_") := 
               rowSums(select_if(tmp$Count_Data, is.numeric)))
  } 
  
  Meta_Data_new <- datalist$Meta_Data %>%
    group_by(!!!param) %>%
    summarize_if(is.numeric, mean) %>%
    ungroup() %>%
    mutate(Sample_ID = deframe(unite(select(., !!! param), ""))) %>%
    dplyr::slice(match(names(select_if(tmp_final, is.numeric)),Sample_ID))
  
  
  datalist_return <- list(Count_Data = tmp_final, Meta_Data = Meta_Data_new) %>%
    mutate_count_datalist(function(x) x/sum(x))
  
  return(datalist_return)
  
}

create_datatable <- function(datalist, grpBy = Class, upperGroup = NULL, otherThreshold = 0.01, selectGrps = NULL, 
                             addColorScheme = F, addOthers = T) {
  
  make_proportion <- function(x) {
    result <- x/sum(x, na.rm = T)
    return(ifelse(is.finite(result), result, 0))
  }
  
  grpBy <- enquo(grpBy)
  upperGroup <- enquo(upperGroup)
  
  if (addColorScheme) {
    
    datalist$Count_Data <- mutate(datalist$Count_Data, !!quo_name(grpBy) := paste(Class, !!grpBy, sep = " - "))
    
  }
  
  if (is.null(selectGrps)) {
    
    if (addOthers) {
      
    tmp_Class <- datalist$Count_Data %>%
      group_by(!!upperGroup) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T)
    
    tmp_Family <- datalist$Count_Data %>%
      group_by(!!grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(apply(select(., -1), 2, make_proportion) %>% apply(., 1, mean, na.rm = T) > otherThreshold) 
    
    tmp_both <- tmp_Family %>%
      mutate(Class = left_join(tmp_Family, distinct(select(datalist$Count_Data, Class, !!grpBy)))$Class, .before = 1) 
    
    tmp_Class_red <- tmp_both %>%
      group_by(Class) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T)
    
    tmp <- (select_if(filter(tmp_Class, Class %in% tmp_Class_red$Class), is.numeric) - 
                        select_if(tmp_Class_red, is.numeric)) %>%
      as_tibble() %>%
      mutate(Family = paste0(tmp_Class_red$Class, " - Others"), .before = 1) %>%
      rbind(., tmp_Family) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0) %>%
      add_case(summarise(., across(where(is.numeric), sum)) %>% 
                 mutate_if(is.numeric, function(x) 1-x) %>%
                 mutate(Family = "Other Classes")) %>%
      filter(rowSums(select_if(., is.numeric)) > 0)
    
    } else {
      
      tmp <- datalist$Count_Data %>%
        group_by(!!grpBy) %>%
        select_if(is.numeric) %>%
        summarize_all(sum, na.rm = T) %>%
        filter(apply(select(., -1), 2, make_proportion) %>% apply(., 1, mean, na.rm = T) > otherThreshold) %>%
        rbind(., c("Others", apply(select_if(datalist$Count_Data, is.numeric), 2, sum) - apply(select(., -1), 2, sum))) %>%
        mutate_at(1, ~replace(., is.na(.), "")) %>%
        mutate_at(2:ncol(.), as.numeric) %>%
        select_if(~sum(is.na(.)) == 0)
      
    }
      
    
  } else {
    
    tmp <- datalist$Count_Data %>%
      group_by(!! grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(!! grpBy %in% selectGrps) %>%
      rbind(., c(paste0("Others (<", otherThreshold * 100,"%)"), apply(select_if(datalist$Count_Data, is.numeric), 2, sum) - apply(select(., -1), 2, sum))) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0)
    
  }
  
  if (addOthers) {
    
    group_label <- NA
  
    for (i in 1:length(tmp_Class_red$Class)) {
      
      group_label <- c(group_label, filter(tmp_both, Class == tmp_Class_red$Class[i])$Family, 
                       paste0(tmp_Class_red$Class[i], " - Others"))
      
    }
    
    group_label <- c(group_label[-1], "Other Classes")
    
    datatable <- tmp %>%
      reshape2::melt(.) %>%
      dplyr::rename("Sample_ID" = "variable", "Abundance" = "value") %>%
      left_join(., distinct(datalist$Meta_Data), by = "Sample_ID") %>%
      as_tibble() %>%
      dplyr::rename("Group" = as_label(grpBy)) %>%
      mutate(Group = ordered(Group, levels = group_label)) %>%
      arrange(Group)
    
  } else {
    
    tmp_richness <- datalist$Count_Data %>%
      mutate_if(is.numeric, function(x) as.numeric(x > 0)) %>%
      group_by(!!grpBy) %>%
      select_if(is.numeric) %>%
      summarize_all(sum, na.rm = T) %>%
      filter(!!grpBy %in% deframe(select(tmp, as_label(grpBy)))) %>%
      rbind(., c("Others", apply(select_if(datalist$Count_Data, is.numeric) > 0, 2, sum) - apply(select(., -1), 2, sum))) %>%
      mutate_at(1, ~replace(., is.na(.), "")) %>%
      mutate_at(2:ncol(.), as.numeric) %>%
      select_if(~sum(is.na(.)) == 0) 
    
    datatable <- tmp %>%
      reshape2::melt(.) %>%
      dplyr::rename("Sample_ID" = "variable", "Abundance" = "value") %>%
      left_join(., {tmp_richness %>% reshape2::melt(.) %>% dplyr::rename("Sample_ID" = "variable", "Richness" = "value")}, by = c("Sample_ID", as_label(grpBy))) %>%
      left_join(., distinct(datalist$Meta_Data), by = "Sample_ID") %>%
      as_tibble() %>%
      dplyr::rename("Group" = as_label(grpBy)) %>%
      mutate(Group = ordered(Group, levels = c(unique(Group)[-which(Group == "Others")], "Others"))) %>%
      arrange(Group)
    
  }
  
  if (addColorScheme) {
    
    #colorTaxonomy <- read_csv("https://raw.githubusercontent.com/dermilke/ExCom/master/Data/Colors/Taxonomy_Colour.csv")
    colorTaxonomy <- read_csv("~/PhD/SoftwareBuilds/ExCom/Data/Colors/Taxonomy_Colour.csv")
    
    derivative_color = str_split_fixed(datatable$Group, " - ", 2) %>%
      unique() %>%
      .[,1] %>%
      table() %>%
      as_tibble() %>%
      dplyr::rename("Group" = ".") %>%
      mutate(Group = ordered(Group, levels = c(unique(Group[Group != "Other Classes"]), "Other Classes"))) %>%
      arrange(Group) %>%
      mutate(main_color = colorTaxonomy$Colour[match(Group, colorTaxonomy$Group)] %>%
               ifelse(is.na(.), "grey", .)) %>%
      mutate(rowNumber = 1:n()) %>%
      with(., {
        derivative_color <- NULL
        for (i in rowNumber) {
          
          derivative_color <- c(derivative_color,
                                colorspace::lighten(main_color[i], 
                                                    amount = {if (n[i] == 1) 0
                                                      else if (n[i] > 6) (seq_len(n[i])/n[i])*1.6-0.8
                                                      else (seq_len(n[i])/n[i])-0.5}
                                )
          )
          
        }
        derivative_color
      })
    
    result <- list(table = datatable, color = derivative_color)
    
  } else {
    
    result <- datatable
    
  }
  
  return(result)
  
}
