## Package selection ----

library(yaml)
library(dplyr)
library(marmap)
library(akima)
library(sf)
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(doParallel)
library(stars)


sst_path = paste0("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/SST/SST_L4_2020_8days.nc")
pp_path = paste0("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/PP/PP_L3_2020_8days.nc")
CHL_path = paste0("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/CHL-PCA_L4/CHL-PCA_L4_2020_8days.nc")
PAR_path = paste0("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/_PAR_/_PAR__PAR__L3_2020_8days.nc")
Bathy_path = paste0("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Bathymetry/REGRIDDED1km_gebco_2023_n54.0_s37.0_w-74.0_e-45.0.nc")



############### ---------------
##            LOAD VARIABLES ENV SAVED AS NETCDF
############### ---------------

#Ouverture des rasters avec le package stars
sst      <- stars::read_ncdf(sst_path,
                         var = "sst", 
                         make_units = TRUE, 
                         make_time = TRUE, 
                         proxy = F)


pp       <- stars::read_ncdf(pp_path,
                        var = c("pp"), 
                        make_units = TRUE, 
                        make_time = TRUE, 
                        proxy = F)


CHL      <- stars::read_ncdf(CHL_path,
                         var = c("CHL-PCA_mean"), 
                         make_units = TRUE, 
                         make_time = TRUE, 
                         proxy = FALSE)


PAR_mean <- stars::read_ncdf(PAR_path,
                         var = "PAR_mean", 
                         make_units = TRUE, 
                         make_time = TRUE, 
                         proxy = F)


bathy    <- stars::read_ncdf(Bathy_path, 
                             var = 'elevation', 
                             make_units = TRUE, 
                             make_time = TRUE, 
                             proxy = FALSE)


BathyCRS=st_crs(sst)
st_crs(bathy)=BathyCRS





##############---------
##        Load Model
#############----------

cfg = read_yaml('~/Documents/SDM_SIMBA/Code/SDM_Biomod2_Optim.yml')


Model_path = paste0("/home/thieryf/Documents/SDM_SIMBA/Models/8days/St_Laurent/2017to2020.2020/v1.001/NARW/NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out")
myBiomodEM = load(Model_path)



###############---------- 
##        Spatial Projections
##############---------

# Load data
Proj_path = paste0("/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Environmental_Var/stars/Proj/8days/St_Laurent/", cfg$Proj_Year, '/', "proj_8days_St_Laurent_", cfg$Proj_Year, ".RData")
proj_Env = load(Proj_path)


###########################
# One week in particular

week_nb = 25

# semaine 25
bathy_week = as.data.frame( bathy )
rm(bathy)

list = ls()
env <- bathy_week
new_columns = NULL
for (i in (list)) {
  objet <- get(i)  # Obtient l'objet en utilisant son nom
  if (class(objet)[1] == "stars") {  # Vérifie si l'objet est de la classe "stars"
    env_dataframe <- as.data.frame(objet[ , 400:1300 , 400:1000 , week_nb ])  # Convertit l'objet stars en dataframe
    colnames(env_dataframe)[4] = i
    env = cbind(env, env_dataframe[,4])
    new_columns = c(new_columns, i)
    colnames(env[length(colnames(env))]) = eval(i)
  }
}

colnames(env)[-c(1:4)] = new_columns
env = rename(env, bathy = "elevation")
env$bathy = as.numeric(env$bathy)



# #######################
# # All weeks combined
# 
# nb_weeks = 30
# 
# bathy_week = as.data.frame( bathy )
# rm(bathy)
# 
# list = ls()
# env <- do.call(rbind, replicate(nb_weeks, bathy_week, simplify = F))
# new_columns = NULL
# for (i in (list)) {
#   objet <- get(i)  # Obtient l'objet en utilisant son nom
#   if (class(objet)[1] == "stars") {  # Vérifie si l'objet est de la classe "stars"
#     env_dataframe <- as.data.frame(objet[ , , , 17:46])  # Convertit l'objet stars en dataframe
#     colnames(env_dataframe)[4] = i
#     env = cbind(env, env_dataframe[,4])
#     new_columns = c(new_columns, i)
#     colnames(env[length(colnames(env))]) = eval(i)
#   }
# }
# 
# colnames(env)[-c(1:4)] = new_columns
# env = rename(env, bathy = "elevation")
# env$bathy = as.numeric(env$bathy)



######################
## set data for projection ----

EMproj.model  <- NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out#_v1.0
EMproj.name   <- paste0('Proj', '.', cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run)
EMproj.env    <- env
EMproj.xy     <- env[, 1:2]
EMproj.chosen <- get_built_models( NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out,
                                   full.name = "NARW_EMcaByROC_mergedData_mergedRun_mergedAlgo" ) 


# Project ensemble models (building single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = EMproj.model,
                                             proj.name = EMproj.name,
                                             new.env = EMproj.env,
                                             models.chosen = EMproj.chosen,
                                             #metric.binary = 'all',
                                             #metric.filter = 'all',
                                             compress = TRUE,
                                             nb.cpu = 15,
                                             na.rm = TRUE)



# get information from the plot
proj.gg <- plot(
  x           = myBiomodEMProj, #a BIOMOD.projection.out object
  coord       = EMproj.xy,                  #a 2-columns data.frame containing the corresponding X and Y
  plot.output = "list",                #(optional, default facet) a character determining the type of output: with plot.output = 'list' the function will return a list of plots (one plot per model) ; with 'facet' ; with plot.output = 'facet' the function will return a single plot with all asked projections as facet.
  do.plot     = FALSE,                 #(optional, default TRUE) a boolean determining whether the plot should be displayed or just returned.
  #std         = TRUE,                  #(optional, default TRUE) a boolean controlling the limits of the color scales. With std = TRUE color scales are displayed between 0 and 1 (or 1000). With std = FALSE color scales are displayed between 0 and the maximum value observed.
  #scales      = "fixed",               #(optional, default fixed) a character determining whether x and y scales are shared among facet. Argument passed to facet_wrap. Possible values: 'fixed', 'free_x', 'free_y', 'free'.
  #size        = 0.3,                   #(optional, default 0.75) a numeric determing the size of points on the plots and passed to geom_point.
  #maxcell     = 5e+05,                 #maximum number of cells to plot. Argument transmitted to plot.
)


## --- save proj

#png( paste0( EMproj.name, '.', myBiomodEMProj@models.projected ), width = 1000, height = 800 )
proj.gg
#dev.off()


# get data & get rid of NA
# !!! WARNING !!! there is only one element in the output list
#                 because I specified with the option "models.chosen"
#                 that I want only 1 model; if you request more, you 
#                 will simply have a list with several elements


proj.na   <- is.na( proj.gg$data$pred )
proj.data <- proj.gg$data[!proj.na,]


Date8days = paste0("2020-07-11")

# Add the presences
pres.i      <- which(as.character(Eval_Sightings$Date8days) == Date8days & Eval_Sightings$SPECCODE == 1)
Eval_Sightings_sf = st_as_sf(Eval_Sightings, coords= c('Long', "Lat"), crs = 4326)
pres.xy     <- st_coordinates( Eval_Sightings_sf$geometry[pres.i] ) %>% 
  as.data.frame()

opp = readRDS("/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Sightings/8days/Opportunistic/St_Laurent/2020/Sightings_8days_Opportunistic_St_Laurent_2020.RDS")
opp.i =  which(as.character(opp$Date8days) == Date8days)
opp.xy     <- st_coordinates( opp$geometry[opp.i] ) %>% 
  as.data.frame()

#pres.xy_SL = pres.xy[pres.xy$X < -60 & pres.xy$Y > 45.5, ] 
#opp.xy_SL = opp.xy[opp.xy$X < -60 & opp.xy$Y > 45.5, ] 


# Add the absences
abs.i      <- which(as.character(Eval_Sightings$Date8days) == Date8days & Eval_Sightings$SPECCODE == 0)
Eval_Sightings_sf = st_as_sf(Eval_Sightings, coords= c('Long', "Lat"), crs = 4326)
abs.xy     <- st_coordinates( Eval_Sightings_sf$geometry[abs.i] ) %>% 
  as.data.frame()



# make our own ggplot of the projection
g.proj <- ggplot()                            +
  geom_point( data = proj.data,
              aes( x     = x,
                   y     = y,
                   color = pred*1e-3 ),
              size = 0.8 )                    +
  scale_color_viridis_c( option = "viridis",
                         limits = c(0,1),
                         na.value = "white" ) +
  geom_point( data = pres.xy,
              aes( x     = X,
                   y     = Y ),
              color = "red",
              shape = 3,
              size  = 3 )                   +
  geom_point( data = opp.xy,
              aes( x     = X,
                   y     = Y ),
              color = "cyan",
              shape = 3,
              size  = 3 )                   +
  labs( x        = "Longitude",
        y        = "Latitude",
        subtitle = "Predicted probability of presence",
        color    = "Prob.")

png( paste0( "ggplot", EMproj.name, '.', NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.outProj@models.projected[[n]], ".", week_nb, ".png"), width = 1500, height = 1000 )
g.proj
dev.off()




######Frequency of presences
proj_df = as.data.frame(cbind(EMproj.xy$x, EMproj.xy$y, NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.outProj@proj.out@val$pred[NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.outProj@proj.out@val$full.name == "NARW_EMcaByROC_mergedData_mergedRun_mergedAlgo"]))     #paste0(cfg$Algo)[[n]]]))

colnames(proj_df) = c('Long', 'Lat', 'Pred')

####### conversion en sf objet
proj_sf = st_as_sf(proj_df, coords= c('Long', "Lat"), crs = 4326)   # all projections/predictions
pres.xy_sf = st_as_sf(pres.xy, coords= c('X', "Y"), crs = 4326)     # actual sightings

###### conversion en star objet
proj_stars <- st_rasterize(proj_sf %>% dplyr::select(Pred, geometry))  # predictions as star object
extract_proj = st_extract(proj_stars, pres.xy_sf)  # extraction of the prob of prediction at sightings localizations

hist(extract_proj$Pred, xlim = c(0, 1000), xlab = 'Prediction value at Sightings spots', main = 'Frequency of probability of presence at Sightings spots')


######Frequency of absences
proj_df = as.data.frame(cbind(EMproj.xy$x, EMproj.xy$y, NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.outProj@proj.out@val$pred[NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.outProj@proj.out@val$full.name == "NARW_EMcaByROC_mergedData_mergedRun_mergedAlgo"]))     #paste0(cfg$Algo)[[n]]]))

colnames(proj_df) = c('Long', 'Lat', 'Pred')

####### conversion en sf objet
proj_sf = st_as_sf( proj_df,        # all projections/predictions
                    coords= c('Long', "Lat"), 
                    crs = 4326)   

abs.xy_sf = st_as_sf(abs.xy,        # actual sightings
                     coords= c('X', "Y"), 
                     crs = 4326)    

###### conversion en star objet
proj_stars <- st_rasterize(proj_sf %>% dplyr::select(Pred, geometry))  # predictions as star object

extract_proj_abs = st_extract(proj_stars, abs.xy_sf)  # extraction of the prob of prediction at sightings localizations

hist(extract_proj_abs$Pred, xlim = c(0, 1000), xlab = 'Prediction value at NO Sightings spots', main = 'Frequency of probability of absence at NO Sightings spots')


######Frequency of opportunistic
proj_df = as.data.frame(cbind(EMproj.xy$x, EMproj.xy$y, NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.outProj@proj.out@val$pred[NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.outProj@proj.out@val$full.name == "NARW_EMcaByROC_mergedData_mergedRun_mergedAlgo"]))     #paste0(cfg$Algo)[[n]]]))

colnames(proj_df) = c('Long', 'Lat', 'Pred')

####### conversion en sf objet
proj_sf = st_as_sf( proj_df,        # all projections/predictions
                    coords= c('Long', "Lat"), 
                    crs = 4326)   

opp.xy_sf = st_as_sf(opp.xy,        # actual sightings
                     coords= c('X', "Y"), 
                     crs = 4326)    

###### conversion en star objet
proj_stars <- st_rasterize(proj_sf %>% dplyr::select(Pred, geometry))  # predictions as star object

extract_proj_opp = st_extract(proj_stars, opp.xy_sf)  # extraction of the prob of prediction at sightings localizations


donnees_combines <- bind_rows(list(Presence = extract_proj, Absence = extract_proj_abs, Opp = extract_proj_opp), .id = 'Occ')

# Créez l'histogramme en spécifiant différentes esthétiques pour chaque série de données
histogramme = ggplot(
  donnees_combines, 
  aes(Pred, fill = Occ)) +
  
  geom_histogram(bins = 20, 
                 position = "identity", 
                 alpha = 0.5,
                 mapping = aes(y = stat(ncount)))

histogramme

#################
### Compare hist with random points 

# Créez une liste vide pour stocker les échantillons
random_pts <- list()

# Nombre de répétitions
n_repetitions <- 100

# Taille de chaque échantillon
sample_lenght <- 30

# Boucle for pour effectuer l'échantillonnage 100 fois
for (i in 1:n_repetitions) {
  # Sélectionnez aléatoirement 30 indices de lignes sans remplacement
  index <- sample(1:nrow(proj_df), sample_lenght, replace = FALSE)
  
  # Extrayez les lignes correspondantes de votre dataframe
  random_pts <- donnees[index, ]
  
  # Ajoutez l'échantillon à la liste
  random_pts[[i]] <- random_pts

}



################
### Calcul evaluation of prediction by months by year

## filter EMmean by ROC or by TSS
NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out_EMmean_evalpredictROC = NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out@models.prediction.eval@val$pred.eval[NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out@models.prediction.eval@val$algo == 'EMmean' & NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out@models.prediction.eval@val$filtered.by == 'ROC']
NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out_EMmean_evalpredictTSS = NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out@models.prediction.eval@val$pred.eval[NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out@models.prediction.eval@val$algo == 'EMmean' & NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out@models.prediction.eval@val$filtered.by == 'TSS']

## binding with the date and coordinates
predict_eval_SL_2017to2020 = cbind(OPP_ENV_2017to2020$LONGITUDE, OPP_ENV_2017to2020$LATITUDE, as.data.frame(NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out_EMmean_evalpredictROC), as.data.frame(NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out_EMmean_evalpredictTSS), OPP_ENV_2017to2020$Date8days)
predict_eval_SL_2017to2020$year = year(predict_eval_SL_2017to2020$`OPP_ENV_2017to2020$Date8days`)
predict_eval_SL_2017to2020$month = month(predict_eval_SL_2017to2020$`OPP_ENV_2017to2020$Date8days`)
predict_eval_SL_2017to2020 = rename(predict_eval_SL_2017to2020, Long = 'OPP_ENV_2017to2020$LONGITUDE', Lat = 'OPP_ENV_2017to2020$LATITUDE', EMmean_evalpredictROC = 'NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out_EMmean_evalpredictROC', EMmean_evalpredictTSS = 'NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out_EMmean_evalpredictTSS', Date8days = 'OPP_ENV_2017to2020$Date8days' )

## mean of the evaluation of the prediction by months by year
mean_predict_eval_SL_2017to2020 = predict_eval_SL_2017to2020  %>%
  group_by(month, year) %>%
  summarize(Mean_ROC = mean(EMmean_evalpredictROC)) 

meanTSS_predict_eval_SL_2017to2020= predict_eval_SL_2017to2020 %>%
  group_by(month, year) %>%
  summarize(Mean_TSS = mean(EMmean_evalpredictTSS))


## plot
ggplot(mean_predict_eval_SL_2017to2020, 
       aes(x = month, y = Mean_ROC, color = as.factor(year)))  +
  geom_point(size=4) +
  geom_smooth(method=lm, se=F)

##
capture.output(predict_eval_SL_2020,
               file=paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run, 'EM_Var_Imp.csv'))




png(paste0(directory, "/EM_Eval_mean.png"), width = 1024, height = 768)
bm_PlotEvalMean(bm.out = NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out, group.by = 'algo')
dev.off()

png(paste0(directory, "/EM_Var_Imp.png"))
bm_PlotVarImpBoxplot(bm.out = NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out, group.by = c('expl.var', 'algo', 'algo'))
dev.off()

png(paste0(directory, "/EMResp_curv.png"), width = 1024, height = 768)
bm_PlotResponseCurves(bm.out = NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out, 
                      models.chosen = get_built_models(NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out, full.name = c('RIWH_EMwmeanByROC_mergedData_mergedRun_mergedAlgo','RIWH_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo')),
                      fixed.var = 'median', 
                      do.bivariate = F)         
dev.off()

