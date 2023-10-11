library(raster)
library(dplyr)
library(readr)
library(tidyr)
library(terra)
library(lubridate)
library(ncdf4)
library(data.table)
library(tibble)
library(stringr)
library(spatstat)
library(sp)
library(stars)
library(ncmeta)
library(readxl)



############### ---------------
#############################################
##################           VARIABLES ENV
############################################

#Ouverture des rasters avec le package stars
SST8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/SST/SST_L4_2020_8days.nc",
                         var = "sst", make_units = TRUE, make_time = TRUE, proxy = F)


PP8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/PP/PP_L3_2020_8days.nc",
                        var = c("pp"), make_units = TRUE, make_time = TRUE, proxy = F)


CDM8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/CDM/CDM_L3_2020_8days.nc",
                         var = c("CDM_mean"), make_units = TRUE, make_time = TRUE, proxy = FALSE)


PAR8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/_PAR_/_PAR__PAR__L3_2020_8days.nc",
                         var = "PAR_mean", make_units = TRUE, make_time = TRUE, proxy = F)


BBP8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/BBP/BBP_L3_2020_8days.nc",
                         var = "BBP_mean", make_units = TRUE, make_time = TRUE, proxy = F)


CHL8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/CHL-PCA/CHL-PCA_L3_2020_8days.nc",
                         var = "CHL-PCA_mean", make_units = TRUE, make_time = TRUE, proxy = F)


Bathy <- read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Bathymetry/REGRIDDED1km_gebco_2023_n54.0_s37.0_w-74.0_e-45.0.nc", var = 'elevation', make_units = TRUE, make_time = TRUE, proxy = FALSE)

#plot(Bathy)

BathyCRS=st_crs(SST8days2020)
st_crs(Bathy)=BathyCRS

#Bathy2 <- read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Bathymetry/NWA36_Bathymetry_for_GLORYS12V1_forcing2.nc", var = 'Bathymetry', make_units = TRUE, make_time = TRUE, proxy = FALSE)



############### ---------------
#############################################
##################           SIGHTINGS
############################################


############## CONSORTIUM
#########################
RIWHs_effort_occ_aerial_2015_2020 <- read.csv("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Insitu_data/NARWC/Dan & Fanny - RIWHs & effort - aerial 2015-2020.csv", sep=",")

#Extraction des données d'occurrences pour l'année 2020
RIWH_occ2020 = subset(RIWHs_effort_occ_aerial_2015_2020, YEAR == "2020") 

#Filtre des sightings a retenir suivant les codes LEGTYPE LEGSTAGE
RIWH_occ2020_filt <- subset(RIWH_occ2020, (RIWH_occ2020$LEGTYPE == 2 & RIWH_occ2020$LEGSTAGE %in% c("1", "2", "3", "4", "5", "7")))
RIWH_occ2020_filt2 <- subset(RIWH_occ2020, (RIWH_occ2020$LEGTYPE == 7 & RIWH_occ2020$LEGSTAGE %in% c("1", "2", "5")))
RIWH_occ2020_filt3 <- subset(RIWH_occ2020, (RIWH_occ2020$LEGTYPE == 9 & RIWH_occ2020$LEGSTAGE %in% c("1", "2")))
RIWH_occ2020_filt4 <- subset(RIWH_occ2020, (RIWH_occ2020$LEGTYPE == 5 & RIWH_occ2020$LEGSTAGE %in% c("1", "2")))
RIWH_occ2020_filt5 <- subset(RIWH_occ2020, (RIWH_occ2020$LEGTYPE == 6 & RIWH_occ2020$LEGSTAGE %in% c("1", "2")))


RIWH_occ2020_filt <- rbind(RIWH_occ2020_filt, RIWH_occ2020_filt2, RIWH_occ2020_filt3, RIWH_occ2020_filt4, RIWH_occ2020_filt5)
RIWH_occ2020_filt$SPECCODE <- ifelse(RIWH_occ2020_filt$SPECCODE == "", 0, 1)


#
# Création de la colonne Date au meme format que les rasters YYYY-MM-DD
RIWH_occ2020_filt$date <- make_date(year = RIWH_occ2020_filt$YEAR, RIWH_occ2020_filt$MONTH, RIWH_occ2020_filt$DAY)

####Créer une colonne des dates des occurences qui correspondent aux 8days des variables ENV

# Colonne supplémentaire qui va matcher avec les 8days ENV
Date_8days_vect = (st_get_dimension_values(SST8days2020, "time"))

#Fonction pour faire matcher les dates des occurrences aux dates des 8days prédédents
RIWH_occ2020_filt$Date8days = floor_date(
  RIWH_occ2020_filt$date,
  unit = Date_8days_vect
  #change_on_boundary = NULL,
  #week_start = getOption("lubridate.week.start", 8)
)

# Transformer occurrences df en objet spatial pour l'extraction
RIWH_occ2020_sf = st_as_sf(RIWH_occ2020_filt, coords = c("LONG_DD", "LAT_DD"), crs = 4326)
###############
########
ggplot(data = RIWH_occ2020_filt) + 
  geom_point(mapping = aes (x = LONG_DD, y = LAT_DD, color = SPECCODE))



############## DFO - SL
#########################
RIWH_effort_occ_aerial_2020_DFO <- read.csv("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Insitu_data/DFO_sighting/Formatted_sighting_data/20200410_155125_20201125_201024_UTC_Platform_337 JOD_337 YOB_Partenavia_Twin Otter_Survey_NARW_2020.csv", sep=",")

RIWH_effort_occ_aerial_2020_DFO_2 <- read.csv("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Insitu_data/DFO_sighting/Formatted_sighting_data/20200829_092730_20201115_161846_UTC_Platform_Twin Otter_Survey_NARW_2020.csv", sep=",")

RIWH_effort_occ_aerial_2020_DFO_2$optional = NULL

RIWH_effort_occ_aerial_2020_DFO = rbind(RIWH_effort_occ_aerial_2020_DFO_1, RIWH_effort_occ_aerial_2020_DFO_2)

RIWH_2020_DFO_filt <- RIWH_effort_occ_aerial_2020_DFO[!(RIWH_effort_occ_aerial_2020_DFO$Effort_l == "off" | RIWH_effort_occ_aerial_2020_DFO$Effort_r == "off"), ]

RIWH_2020_DFO_filt$Sp_code = ifelse(RIWH_2020_DFO_filt$Sp_code == "Eg", 1, 0)
RIWH_2020_DFO_filt$Sp_code[is.na(RIWH_2020_DFO_filt$Sp_code)] = 0

# Domaine du SL
RIWH_2020_DFO_filt = RIWH_2020_DFO_filt[RIWH_2020_DFO_filt$Lat > 45.5 & RIWH_2020_DFO_filt$Long < -60, ]


# Création de la colonne Date au meme format que les rasters YYYY-MM-DD
RIWH_2020_DFO_filt$date <- as.Date(RIWH_2020_DFO_filt$DT_utc)

####Créer une colonne des dates des occurences qui correspondent aux 8days des variables ENV

# Colonne supplémentaire qui va matcher avec les 8days ENV
Date_8days_vect = (st_get_dimension_values(SST8days2020, "time"))

#Fonction pour faire matcher les dates des occurrences aux dates des 8days prédédents
RIWH_2020_DFO_filt$Date8days = floor_date(
  RIWH_2020_DFO_filt$date,
  unit = Date_8days_vect
  #change_on_boundary = NULL,
  #week_start = getOption("lubridate.week.start", 8)
)


# Renommer les colonnes qui vont servir a la liaison avec les donnees du CONSORTIUM
colnames(RIWH_2020_DFO_filt)[31] = "SPECCODE"

RIWH_2020_DFO_filt_sf = st_as_sf(RIWH_2020_DFO_filt, coords= c('Long', "Lat"), crs = 4326)


###############
########
ggplot(data = RIWH_2020_DFO_filt) + 
  geom_point(mapping = aes (x = Long, y = Lat, color = SPECCODE))


#############
################## OPPORTUNISTIC
#############

NARW_MAR_Oppotunistic_Sightings <- read_excel("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Insitu_data/DFO_sighting/NARW_MAR_Oppotunistic Sightings.xlsx")


# Domaine du SL
NARW_MAR_Oppotunistic_Sightings = NARW_MAR_Oppotunistic_Sightings[NARW_MAR_Oppotunistic_Sightings$LATITUDE > 45.5, ] 
                                                                  #& NARW_MAR_Oppotunistic_Sightings$LONGITUDE < -60, ]


# Remplacer les mois par des chiffres
NARW_MAR_Oppotunistic_Sightings$date = NARW_MAR_Oppotunistic_Sightings$DATE
NARW_MAR_Oppotunistic_Sightings$date <- gsub("JAN", "01", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("FEB", "02", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("MAR", "03", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("AVR", "04", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("MAY", "05", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("JUN", "06", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("JUL", "07", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("AUG", "08", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("SEP", "09", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("OCT", "10", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("NOV", "11", NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$date <- gsub("DEC", "12", NARW_MAR_Oppotunistic_Sightings$date)


# Création de la colonne Date au meme format que les rasters YYYY-MM-DD
NARW_MAR_Oppotunistic_Sightings$date = dmy(NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$Year <- year(NARW_MAR_Oppotunistic_Sightings$date)
NARW_MAR_Oppotunistic_Sightings$Month <- month(NARW_MAR_Oppotunistic_Sightings$date)


# Garder uniquement l'année d'etude
NARW_Oppotunistic_2020 = NARW_MAR_Oppotunistic_Sightings[NARW_MAR_Oppotunistic_Sightings$Year == '2020', ]


# Colonne supplémentaire qui va matcher avec les 8days ENV
Date_8days_vect = (st_get_dimension_values(SST8days2020, "time"))

#Fonction pour faire matcher les dates des occurrences aux dates des 8days prédédents
NARW_Oppotunistic_2020$Date8days = floor_date(
  NARW_Oppotunistic_2020$date,
  unit = Date_8days_vect
  #change_on_boundary = NULL,
  #week_start = getOption("lubridate.week.start", 8)
)

# Créer la colonne d'occurences = 1 (opportunistic donc sighting =1)
NARW_Oppotunistic_2020$SPECCODE = 1

# Transformer occurrences df en objet spatial pour l'extraction
NARW_Oppotunistic_2020_sf = st_as_sf(NARW_Oppotunistic_2020, coords= c('LONGITUDE', "LATITUDE"), crs = 4326)



#####################
#### Ajouter les sightings de DFO au CONSORTIUM

RIWH_occ2020_df = as.data.frame(RIWH_occ2020_filt[, c('SPECCODE', 'Date8days', 'LAT_DD', 'LONG_DD')])
RIWH_occ2020_DFO_df = as.data.frame(RIWH_2020_DFO_filt[, c('SPECCODE', 'Date8days', 'Lat', 'Long')])
colnames(RIWH_occ2020_df)[3] = "Lat"
colnames(RIWH_occ2020_df)[4] = "Long"

Sightings_2020 = rbind(RIWH_occ2020_df, RIWH_occ2020_DFO_df)

# Transformer occurrences df en objet spatial pour l'extraction
Sightings_2020_sf = st_as_sf(Sightings_2020, coords= c('Long', "Lat"), crs = 4326)




############### ---------------
#############################################
##################   SIGHTINGS + VARIABLES
############################################


###### Merger les stars des variables ENV
ENV8days2020 = c(SST8days2020, PP8days2020, PAR8days2020) #CDM8days2020, BBP8days2020,



# Extraction des var ENV suivant les donnes opp

extract_CHL_8days2020_opp = st_extract(CHL8days2020, at = NARW_Oppotunistic_2020_sf, time_column = "Date8days")
extract_Bathy_8days2020_opp = (st_extract(Bathy, NARW_Oppotunistic_2020_sf))
extract_ENV_8days2020_opp = (st_extract(ENV8days2020, at = NARW_Oppotunistic_2020_sf, time_column = "Date8days"))

# Binding des extractions 
extract_ENV_8days2020_opp_bind = cbind(extract_ENV_8days2020_opp, extract_Bathy_8days2020_opp$elevation, extract_CHL_8days2020_opp$`CHL-PCA_mean`)

OPP_ENV_8days2020 <- rename(extract_ENV_8days2020_opp_bind, bathy = 'extract_Bathy_8days2020_opp.elevation', CHL = 'extract_CHL_8days2020_opp..CHL.PCA_mean.')

OPP_ENV_8days2020 = OPP_ENV_8days2020[, c('sst', 'CHL', 'pp', 'CDM_mean', 'BBP_mean', 'PAR_mean', 'bathy', 'time', 'Date8days', 'geometry')]

OPP_ENV_8days2020_df = as.data.frame(OPP_ENV_8days2020, xy = TRUE)
OPP_ENV_8days2020_df$bathy = as.numeric(OPP_ENV_8days2020_df$bathy)



# Extraction des var ENV suivant les donnes occ

extract_CHL_8days2020 = st_extract(CHL8days2020, at = RIWH_2020_DFO_filt_sf, time_column = "Date8days")
extract_Bathy_8days2020 = (st_extract(Bathy, RIWH_2020_DFO_filt_sf))
extract_ENV_8days2020 = (st_extract(ENV8days2020, at = RIWH_2020_DFO_filt_sf, time_column = "Date8days"))

# Binding des extractions 
extract_ENV_8days2020_bind = cbind(extract_ENV_8days2020, extract_Bathy_8days2020$elevation, extract_CHL_8days2020$`CHL-PCA_mean`)

OCC_ENV_8days2020 <- rename(extract_ENV_8days2020_bind, bathy = 'extract_Bathy_8days2020.elevation', CHL = 'extract_CHL_8days2020..CHL.PCA_mean.')

OCC_ENV_8days2020 = OCC_ENV_8days2020[, c('sst', 'CHL', 'pp', 'CDM_mean', 'BBP_mean', 'PAR_mean', 'bathy', 'time', 'Date8days', 'geometry')]     

OCC_ENV_8days2020_df = as.data.frame(OCC_ENV_8days2020, xy = TRUE)
OCC_ENV_8days2020_df$bathy = as.numeric(OCC_ENV_8days2020_df$bathy)

saveRDS(OCC_ENV_8days2020, '/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Environmental_Var/8days/Surveys/St_Laurent/2020/Env_8days_Surveys_St_Laurent_2020.RDS')

saveRDS(RIWH_2020_DFO_filt, '/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Sightings/8days/Surveys/St_Laurent/2020/Sightings_8days_Surveys_St_Laurent_2020.RDS')


# Jeu de données pour les projections                   

nb_weeks = 24

#Juste le St lolo pour test
ENV8days2020_proj_8days = (ENV8days2020[ ,400:1700 , 350:1000, 17:40])
ENV8days2020_proj_8days_df = as.data.frame(ENV8days2020_proj_8days)
Bathy_proj = Bathy[, 400:1700 , 350:1000]
Bathy_proj_df = as.data.frame(Bathy_proj)
Bathy_proj_8days = replicate(nb_weeks, Bathy_proj_df, simplify = F)
Bathy_proj_8days_df = as.data.frame(do.call(rbind, Bathy_proj_8days)) 
CHL8days2020_proj_8days = (CHL8days2020[ ,400:1700 , 350:1000, 17:40])
CHL8days2020_proj_8days_df = as.data.frame(CHL8days2020_proj_8days)

ENV8days2020_proj_all_8days= cbind(ENV8days2020_proj_8days_df, Bathy_proj_8days_df$elevation, CHL8days2020_proj_8days_df$CHL.PCA_mean)
ENV8days2020_proj_all_8days_df = rename(ENV8days2020_proj_all_8days, bathy = 'Bathy_proj_8days_df$elevation', CHL = 'CHL8days2020_proj_8days_df$CHL.PCA_mean')
ENV8days2020_proj_all_8days_df$bathy = as.numeric(ENV8days2020_proj_all_8days_df$bathy)

# semaine 33
ENV8days2020_test2 = (ENV8days2020[ ,400:1300 , 400:1000, 33])
ENV8days2020_test2_df = as.data.frame(ENV8days2020_test2)
Bathy_test = Bathy[, 400:1300, 400:1000]
Bathy_test_df = as.data.frame(Bathy_test)

ENV8days2020_Bathy_test2_df = cbind(ENV8days2020_test2_df, Bathy_test_df$elevation)
ENV8days2020_Bathy_test2_df = rename(ENV8days2020_Bathy_test2_df, bathy = 'Bathy_test_df$elevation')
ENV8days2020_Bathy_test2_df$bathy = as.numeric(ENV8days2020_Bathy_test2_df$bathy)

#Semaine 30
ENV8days2020_test3 = (ENV8days2020[ ,400:1300 , 400:1000, 30])
ENV8days2020_test3_df = as.data.frame(ENV8days2020_test3)
Bathy_test = Bathy[, 400:1300, 400:1000]
Bathy_test_df = as.data.frame(Bathy_test)

ENV8days2020_Bathy_test3_df = cbind(ENV8days2020_test3_df, Bathy_test_df$elevation)
ENV8days2020_Bathy_test3_df = rename(ENV8days2020_Bathy_test3_df, bathy = 'Bathy_test_df$elevation')
ENV8days2020_Bathy_test3_df$bathy = as.numeric(ENV8days2020_Bathy_test3_df$bathy)

load("/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Environmental_Var/stars/Proj/8days/St_Laurent/2020/proj_8days_St_Laurent_2020.RData")

sst_df = as.data.frame(sst)
pp_df = as.data.frame(pp)
CHL_df = as.data.frame(CHL)
CDM_df = as.data.frame(CDM_mean)
BBP_df = as.data.frame(BBP_mean)
PAR_df = as.data.frame(PAR_Mean)
bathy_df = as.data.frame(bathy)

proj_env = cbind(sst_df, pp_df, CDM_df, BBP_df, PAR_df, CHL_df, bathy)

rm(SST8days2020, PP8days2020, BBP8days2020, CDM8days2020, PAR8days2020, CHL8days2020, Bathy)
save.image("/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Environmental_Var/stars/Proj/8days/St_Laurent/2020/proj_8days_St_Laurent_2020.RData")
