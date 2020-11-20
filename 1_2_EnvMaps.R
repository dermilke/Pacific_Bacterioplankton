### Environment Maps ###

rm(list = ls(all=T))

library(tidyverse)
library(gridExtra)

devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datalistFunctions.R")
devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datatableFunctions.R")

##### Read processed tables #####

Meta_Data <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok") %>%
  .$Meta_Data %>%
  dplyr::rename("Temperature" = "Pot_Temperature")

datatable_CTD <- read_delim("~/PhD/Data_Storage/Meta_Data/Pacific/Raw/SO248_Data/CTD/SO248_phys_oce.tab", del = "\t", skip = 43) %>%
  dplyr::rename("Temperature" = "Tpot [°C]", "Depth" = "Depth water [m]", "Salinity" = "Sal (CTD, SEA-BIRD SBE 911plus, SN...)",
                "Fluorescence" = "Fluores [µg/l]", "Turbidity" = "Turbidity [NTU]", "Oxygen" = "OXYGEN [µmol/kg]") %>%
  mutate(Depth = seq(0, 600, 5)[findInterval(Depth, seq(0, 600, 5))]) %>%
  group_by(Latitude, Depth) %>%
  summarize_if(is.numeric, mean) %>%
  select(Latitude, Depth, Temperature, Salinity, Fluorescence) %>%
  full_join(., distinct(select(Meta_Data, Latitude, Temperature, Depth, Salinity, Fluorescence)), by = c("Latitude", "Depth")) %>%
  mutate(Temperature.x = ifelse(is.na(Temperature.y), Temperature.x, Temperature.y)) %>%
  mutate(Salinity.x = ifelse(is.na(Salinity.y), Salinity.x, Salinity.y)) %>%
  mutate(Fluorescence.x = ifelse(is.na(Fluorescence.y), Fluorescence.x, Fluorescence.y)) %>%
  select(-Temperature.y, -Salinity.y, -Fluorescence.y) %>%
  dplyr::rename("Temperature" = "Temperature.x", "Salinity" = "Salinity.x", "Fluorescence" = "Fluorescence.x")

##### Data-wrangling Pipeline #####

model_temp <- loess(Temperature ~ Latitude * Depth, data = datatable_CTD, span = 0.1)
model_sal <- loess(Salinity ~ Latitude * Depth, data = datatable_CTD, span = 0.1)
model_fluo <- loess(Fluorescence ~ Latitude * Depth, data = datatable_CTD, span = 0.1)

grid_resolution <- 1

data.fit <- expand.grid(Latitude = seq(min(datatable_CTD$Latitude), max(datatable_CTD$Latitude), grid_resolution),
                        Depth    = seq(min(datatable_CTD$Depth), max(datatable_CTD$Depth), grid_resolution*5))

temp_colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(300))
sal_colors <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(300))
fluo_colors <- colorRampPalette(c("white", "darkgreen"))(300)

datatable_predicted <-  predict(model_temp, newdata = data.fit) %>%
  reshape2::melt(., id.vars = c("Latitude", "Depth"), measure.vars = "Temperature", value.name = "Temperature") %>%
  as_tibble() %>%
  left_join(., {predict(model_sal, newdata = data.fit) %>%
      reshape2::melt(., id.vars = c("Latitude", "Depth"), measure.vars = "Salinity", value.name = "Salinity") %>%
      as_tibble()}, by = c("Latitude", "Depth")) %>%
  left_join(., {predict(model_fluo, newdata = data.fit) %>%
      reshape2::melt(., id.vars = c("Latitude", "Depth"), measure.vars = "Fluorescence", value.name = "Fluorescence") %>%
      as_tibble()}, by = c("Latitude", "Depth")) %>%
  mutate(Latitude = as.numeric(str_sub(Latitude, str_locate(Latitude, "=")[1,1] + 1))) %>%
  mutate(Depth =  as.numeric(str_sub(Depth, str_locate(Depth, "=")[1,1] + 1)))

#### Create Environment-Maps ####

p1 <- ggplot() +
  stat_contour(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Temperature, fill = ..level..), geom = "polygon") +
  geom_tile(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Temperature, fill = Temperature)) +
  stat_contour(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Temperature, fill = Temperature), 
               col = "black", bins = 20) +
  scale_fill_gradientn(colours = temp_colors, name = "Temperature [°C]") +
  geom_point(aes(y = Depth, x = Latitude), col = "black", alpha = .25, data = datalist_Pacific$Meta_Data, size = 2) +
  scale_x_continuous(position = "bottom", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  theme(legend.position = "right", legend.key.height = unit(1, "cm"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size = 1.5),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 15)) +
  guides(fill = guide_colorbar(title.position = "right", title.theme = element_text(angle = 90, hjust = .5)))

p2 <- ggplot() +
  stat_contour(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Salinity, fill = ..level..), geom = "polygon") +
  geom_tile(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Salinity, fill = Salinity)) +
  stat_contour(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Salinity, fill = Salinity), 
               col = "black", bins = 20) +
  scale_fill_gradientn(colours = sal_colors, name = "Salinity [PSU]") +
  geom_point(aes(y = Depth, x = Latitude), col = "black", alpha = .25, data = datalist_Pacific$Meta_Data, size = 2) +
  scale_x_continuous(position = "bottom", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  theme(legend.position = "right", legend.key.height = unit(1, "cm"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size = 1.5),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 15)) +
  guides(fill = guide_colorbar(title.position = "right", title.theme = element_text(angle = 90, hjust = .5)))

p3 <- ggplot() +
  stat_contour(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Fluorescence, fill = ..level..), geom = "polygon") +
  geom_tile(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Fluorescence, fill = Fluorescence)) +
  stat_contour(data = datatable_predicted, aes(y = Depth, x = Latitude, z = Fluorescence, fill = Fluorescence), 
               col = "black", bins = 20) +
  scale_fill_gradientn(colours = fluo_colors, name = "Fluorescence [V]") +
  geom_point(aes(y = Depth, x = Latitude), col = "black", alpha = .25, data = datalist_Pacific$Meta_Data, size = 2) +
  scale_x_continuous(position = "bottom", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  theme(legend.position = "right", legend.key.height = unit(1, "cm"),
        panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size = 1.5),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12), axis.title = element_text(size = 15)) +
  guides(fill = guide_colorbar(title.position = "right", title.theme = element_text(angle = 90, hjust = .5)))

cowplot::plot_grid(p1, p2, p3, ncol = 1)
