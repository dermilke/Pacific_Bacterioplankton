distance_decay_core <- function(datalist) {
  
  geographic_distance <- function(long1, lat1, long2, lat2) {
    
    dlong <- long1 * pi/180 - long2 * pi/180
    dlat <- lat1 * pi/180 - lat2 * pi/180
    
    R <- 6371
    a <- sin(dlat/2)^2 + cos(lat1 * pi/180) * cos(lat2 * pi/180) * sin(dlong/2)^2
    c <- 2 * atan2(sqrt(a), sqrt(1-a))
    
    distance <- R * c
    
    return(distance)
    
  }
  
  DD_datatable <- datalist$Count_Data %>%
    select_if(is.numeric) %>%
    t(.) %>% 
    vegan::vegdist(.) %>% 
    t(.) %>% 
    as.matrix(.) %>%
    reshape2::melt() %>%
    dplyr::rename("From" = "Var1", "Sample_ID" = "Var2", "Bray_Curtis" = "value") %>%
    left_join(., datalist$Meta_Data, by = "Sample_ID") %>%
    magrittr::set_names(c("Sample_ID", "To", "Bray_Curtis", paste0("To_", names(.)[-c(1:3)]))) %>%
    left_join(., datalist$Meta_Data, by = "Sample_ID") %>%
    magrittr::set_names(c("From", "To", "Bray_Curtis", 
                          names(.)[grepl(pattern = "^To_", x = names(.))],
                          paste0("From_", names(.)[!grepl(pattern = "^To_", x = names(.))])[-c(1:3)])) %>%
    mutate(Distance = geographic_distance(To_Longitude, To_Latitude, From_Longitude, From_Latitude)) %>%
    as_tibble() %>%
    filter(From_Size_Fraction == To_Size_Fraction) %>%
    filter(From_Depth_Grp == To_Depth_Grp)
  
  return(DD_datatable)
  
}