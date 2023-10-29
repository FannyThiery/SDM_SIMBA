library(stars)
library(raster)


############
# Cleaning Extracted data

#### SST, PP, PAR, CHL

Sighting_DFO_SST_L4_2017 <- read.csv("~/data/GC/For_Dan_and_Fanny/Extracted_data/SST_L4_DFO/Sighting_SST_L4_2017.csv")

Sighting_DFO_SST_L4_2017$optional = NULL

Sighting_DFO_SST_L4_2017$SST_mean_8  <- rowMeans(Sighting_DFO_SST_L4_2017[, 78:85])
Sighting_DFO_SST_L4_2017$SST_mean_16 <- rowMeans(Sighting_DFO_SST_L4_2017[, 70:77])  
Sighting_DFO_SST_L4_2017$SST_mean_24 <- rowMeans(Sighting_DFO_SST_L4_2017[, 62:69])
Sighting_DFO_SST_L4_2017$SST_mean_31 <- rowMeans(Sighting_DFO_SST_L4_2017[, 54:61])

SST_DFO_2017_8days <- Sighting_DFO_SST_L4_2017[!(Sighting_DFO_SST_L4_2017$Effort_l == "off" | Sighting_DFO_SST_L4_2017$Effort_r == "off"), ]

SST_DFO_2017_8days$Sp_code = ifelse(SST_DFO_2017_8days$Sp_code == "Eg", 1, 0)
SST_DFO_2017_8days$Sp_code[is.na(SST_DFO_2017_8days$Sp_code)] = 0

SST_DFO_2017_8days = SST_DFO_2017_8days[SST_DFO_2017_8days$Lat > 45.5, ]

SST_DFO_2017_8days <- SST_DFO_2017_8days[, 86:89]

saveRDS(SST_DFO_2017_8days, '/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Extracted_Data_DFO_Clean/SST/SST_DFO_2017_8days.rds')



######## BATHY

Bathy_DFO_2020_8days = readRDS('/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Extracted_Data_DFO_Clean/Bathy/Bathy_DFO_2020_8days.rds')

XY = as.data.frame(st_coordinates(Bathy_DFO_2020_8days$geometry))
Bathy_DFO_2020_8days = cbind(Bathy_DFO_2020_8days$elevation, XY)
Bathy_DFO_2020_8days = Bathy_DFO_2020_8days[Bathy_DFO_2020_8days$Y > 45.5,]

saveRDS(Bathy_DFO_2020_8days,'/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Extracted_Data_DFO_Clean/Bathy/Bathy_DFO_2020_8days.rds')


######Sightings

Sightings_DFO_2017_8days = readRDS('/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Extracted_Data_DFO_Clean/Sightings/Sightings_DFO_2017_8days.rds')


Sightings_DFO_2017_8days = Sightings_DFO_2017_8days[Sightings_DFO_2017_8days$Lat > 45.5,]

saveRDS(Bathy_DFO_2017_8days,'/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Extracted_Data_DFO_Clean/Sightings/Sightings_DFO_2017_8days.rds')








