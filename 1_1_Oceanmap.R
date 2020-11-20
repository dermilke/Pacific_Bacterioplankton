#### Chlorophyll Oceanmap ####

#### Load dependencies ####

rm(list = ls(all=T))

library(tidyverse)
library(raster)
library(oceanmap)

devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datalistFunctions.R")
devtools::source_url("https://raw.githubusercontent.com/dermilke/ExCom/master/datatableFunctions.R")

##### Analysis Functions #####

pacific_centered <- function(raster) {
  
  x1 <- raster::crop(raster, extent(-180, 0, -90, 90))
  x2 <- raster::crop(raster, extent(0, 180, -90, 90))   
  raster::extent(x1) <- c(180, 360, -90, 90)
  raster.projected <- raster::merge(x1, x2)
  
  return(raster.projected)
  
}

#### Read station and satellite data and formatting ####

Stations <- data_select("~/PhD/Data_Storage/ASV_Tabs/Pacific/V4V5_Primerset/Complete/", "Prok")$Meta_Data %>%
  select(Latitude, Longitude, Province) %>%
  distinct() %>%
  arrange(Latitude) %>%
  left_join(., read_csv("https://raw.githubusercontent.com/dermilke/Pacific_Bacterioplankton/master/Colors/Province_Colour.csv"), by = "Province") %>%
  mutate(Longitude = ifelse(Longitude < 0, Longitude + 360, Longitude))

crop.vals <- c(lat = c(140,220), 
               lon = c(-75,75))

Chl <- ncdf4::nc_open("~/PhD/Data_Storage/Meta_Data/OceanMaps/Chl_Annual_4degr_mapped.nc") %>%
  oceanmap::nc2raster(., "chlor_a", lonname="lon", latname="lat", date=T) %>%
  raster::flip(., "y") %>%
  pacific_centered(.) %>%
  raster::crop(., extent(crop.vals))

##### Plotting Oceanmap with Chl a backgroung ####

par(mfrow=c(1,1))
dev.off()

oceanmap::v(Chl, cbpos = "b", pal = colorRampPalette(c("white", "darkgreen"))(300), zlim = c(0,1.5), 
            cb.xlab = expression("Annual chlorophyll-a (mg m"^-3*")"),
            bwd = 0.01, grid = F, replace.na = F, col.bg = "white", border = "#504f4f",
            cex.ticks = 1.5, axeslabels = F, figdim = c(4,5), show.colorbar = T)

points(y = Stations$Latitude, x = Stations$Longitude, pch = 21, col = "black", bg = Stations$Colour, cex = 2.2)

#### Plot Province Legend ####

par(mfrow = c(1,1))
par(xpd=T)

plot(0,0, type = "n", axes = F, ylab = "", xlab = "")
legend("center", pch = 21, pt.bg = unique(Meta %>% arrange(desc(Latitude)) %>% .$Colour), 
       legend = unique(Meta %>% arrange(desc(Latitude)) %>% .$Province),  bty = "n", pt.cex = 1.5,
       cex = .9, horiz = F)

par(xpd=F)
