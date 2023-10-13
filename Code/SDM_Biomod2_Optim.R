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
library(ggOceanMaps)
library(ggspatial)




#Load the yml file that corresponds
cfg = read_yaml('~/Documents/SDM_SIMBA/Code/SDM_Biomod2_Optim.yml')

my_WD = paste0( cfg$output_wd,
                cfg$Tempstemp,  '/',
                cfg$Domain          )

setwd(my_WD)

n = length(cfg$Year)
model_path = paste0( my_WD,         '/', 
                     cfg$Year[1],   'to', 
                     cfg$Year[n],   '.',
                     cfg$Eval_Year, '/' )
dir.create(model_path)

setwd(model_path)

out_config_file =  file.path(dir.create(cfg$run), paste0(cfg$Tempstemp, '/' , cfg$Domain, '/' , cfg$Year[1], 'to', cfg[n], '/', cfg$run))

file.copy("/home/thieryf/Documents/SDM_SIMBA/Code/SDM_Biomod2_Optim.yml", to = paste0(model_path, cfg$run))

setwd(paste0(model_path, cfg$run))

# Functions that will open and Row bind years for training dataset
  #Sightings
Read_Sightings = function( cfg ) {
  
  S_filename <- paste0( cfg$S_data,    '_',
                        cfg$Tempstemp, '_', 
                        cfg$Type,      '_', 
                        cfg$Domain,    '_',
                        cfg$Year,      '.', 
                        'RDS'              )
  
  S_filepath <- paste0( cfg$initial_wd, 
                        cfg$S_data,    '/', 
                        cfg$Tempstemp, '/', 
                        cfg$Type,      '/', 
                        cfg$Domain,    '/', 
                        cfg$Year,      '/', 
                        S_filename         )
  
    x = lapply( S_filepath, readRDS ) |> 
          dplyr::bind_rows()
    
    return( x )
}

  #Variables
Read_Env = function( cfg ) {
  
  E_filename <- paste0( 'Env',         '_',
                        cfg$Tempstemp, '_', 
                        cfg$Type,      '_', 
                        cfg$Domain,    '_', 
                        cfg$Year,      '.', 
                        'RDS'              )
  
  E_filepath <- paste0( cfg$initial_wd, 
                        cfg$E_data,    '/', 
                        cfg$Tempstemp, '/', 
                        cfg$Type,      '/', 
                        cfg$Domain,    '/', 
                        cfg$Year,      '/', 
                        E_filename         ) 
  
  x = lapply( E_filepath, readRDS ) |> 
    dplyr::bind_rows()
  
  return( x )
}
 

# Row binding data for Evaluation
  #Sightings
Read_EvalSightings = function( cfg ) {
  
  EvalS_filename <- paste0( cfg$S_data,    '_', 
                            cfg$Eval_Tempstemp, '_', 
                            cfg$Eval_Type,      '_', 
                            cfg$Eval_Domain,    '_', 
                            cfg$Eval_Year, '.', 
                            'RDS'              )
  
  EvalS_filepath <- paste0( cfg$initial_wd, 
                            cfg$Eval_S_data,    '/', 
                            cfg$Eval_Tempstemp, '/', 
                            cfg$Eval_Type,      '/', 
                            cfg$Eval_Domain,    '/', 
                            cfg$Eval_Year, '/', 
                            EvalS_filename     )
  
  x = lapply( EvalS_filepath, readRDS ) |> 
    dplyr::bind_rows()
  
  return( x )
}

#Variables
Read_EvalEnv = function( cfg ) {
  
  EvalE_filename <- paste0( 'Env',         '_',
                            cfg$Eval_Tempstemp, '_', 
                            cfg$Eval_Type,      '_', 
                            cfg$Eval_Domain,    '_', 
                            cfg$Eval_Year, '.', 
                            'RDS'              )
  
  EvalE_filepath <- paste0( cfg$initial_wd, 
                            cfg$Eval_E_data,    '/', 
                            cfg$Eval_Tempstemp, '/', 
                            cfg$Eval_Type,      '/', 
                            cfg$Eval_Domain,    '/', 
                            cfg$Eval_Year, '/', 
                            EvalE_filename    ) 
  
  x = lapply( EvalE_filepath, readRDS ) |> 
    dplyr::bind_rows()
  return( x )
}




# Put the file into the Sightings object
Sightings = Read_Sightings( cfg )
# Put the file into the Env_var object
Env_var <- Read_Env( cfg )


#Load Evaluation dataset
Eval_Sightings <- Read_EvalSightings( cfg )
Eval_Env <- Read_EvalEnv( cfg )





######################. -------------
###### --- Setting data for Biomod2
######################. -------------

# species name
myRespName <- cfg$MyRespName

# presence/absence data
myResp <- Sightings$SPECCODE

# XY coordinates
myXY <- st_coordinates( Env_var$geometry )

# environmental variables
env.var <- cfg$`Environmental Variables`

myExpl  <- Env_var[ , env.var ]
myExpl$geometry <- NULL

#Evaluation Dataset
myEvalResp  <- Eval_Sightings$SPECCODE
myEvalExpl  <- Eval_Env[ , env.var ]
myEvalExpl$geometry <- NULL
myEvalXY    <- st_coordinates(Eval_Env$geometry)


## Format training & evaluation data ----

## DATA 1.0) use ALL data for training & 2020 data for evaluation ----
myData.all <- BIOMOD_FormatingData( resp.var      = myResp,
                                    expl.var      = myExpl,
                                    resp.xy       = myXY,
                                    resp.name     = myRespName, 
                                    eval.resp.var = myEvalResp,
                                    eval.expl.var = myEvalExpl,
                                    eval.resp.xy  = myEvalXY   )

myData.all

## Plot datasets ----

plot(myData.all)


## Define models options ----

myBiomodOption <- bm_DefaultModelingOptions()


## XXXXXXXXXXXXXXXXXXXX ----
## --------------------   SINGLE MODELLING   -------------------- 


myBiomodModel <- BIOMOD_Modeling(
  bm.format         = myData.all,                                 #a BIOMOD.formated.data or BIOMOD.formated.data.PA object returned by the BIOMOD_FormatingData function
  bm.options        = myBiomodOption,
  modeling.id       = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run), #a character corresponding to the name (ID) of the simulation set (a random number by default)
  models            = cfg$Single_Models, #, "CTA", "GAM", "RF"),           #a vector containing model names to be computed, must be among GLM, GBM, GAM, CTA, ANN, SRE, FDA, MARS, RF, MAXENT.Phillips, MAXENT.Phillips.2
  models.pa         = NULL,                                       #a list containing for each model a vector defining which pseudo-absence datasets are to be used, must be among colnames(bm.format@PA.table)
  CV.strategy       = cfg$CV_strategy,                                   #a character corresponding to the cross-validation selection strategy, must be among random, kfold, block, strat, env or user.defined
  CV.nb.rep         = 10,                                          #if strategy = 'random' or strategy = 'kfold', an integer corresponding to the number of repetitions to be done for calibration/validation splitting
  CV.perc           = cfg$CV_perc,                                        #if strategy = 'random', a numeric between 0 and 1 corresponding to the percentage of data used to calibrate the models (calibration/validation splitting)
  CV.do.full.models = TRUE,                                       #a logical value defining whether models calibrated and evaluated over the whole dataset should be computed or not
  #CV.k              = 3,                                          #if strategy = 'kfold' or strategy = 'strat' or strategy = 'env', an integer corresponding to the number of partitions
  #CV.balance        = 'presences',                                #if strategy = 'strat' or strategy = 'env', a character corresponding to how data will be balanced between partitions, must be either presences or absences
  #CV.env.var        = 'all',                                      #if strategy = 'env', a character corresponding to the environmental variables used to build the partition. k partitions will be built for each environmental variables. By default the function uses all environmental variables available.
  #CV.strat          = 'both',                                     #if strategy = 'env', a character corresponding to how data will partitioned along gradient, must be among x, y, both
  #CV.user.table     = NULL,                                       #a matrix or data.frame defining for each repetition (in columns) which observation lines should be used for models calibration (TRUE) and validation (FALSE) (see BIOMOD_CrossValidation)
  #weights           = 0.6,                                       #a vector of numeric values corresponding to observation weights (one per observation, see Details)
  #prevalence        = 0.6,                                        #a numeric between 0 and 1 corresponding to the species prevalence to build 'weighted response weights' (see Details)
  metric.eval       = cfg$Metrics_Eval,                            #"KAPPA", "TSS","ACCURACY", BIAS, POD, FAR, POFD, SR, CSI, ETS, HK, HSS, OR, ORSS
  var.import        = 5,                                          #an integer corresponding to the number of permutations to be done for each variable to estimate variable importance
  scale.models      = FALSE,                                      #a logical value defining whether all models predictions should be scaled with a binomial GLM or not
  nb.cpu            = 15,                                         #a integer value corresponding to the number of computing resources to be used to parallelize the single models computation
  seed.val          = 42,                                         #an integer value corresponding to the new seed value to be set
  do.progress       = TRUE                                        #a logical value defining whether the progress bar is to be rendered or not
)








## Evaluation scores ----

# create ggplot2 object
# extract evaluations metrics for the "calibration" dataset
g.1 <- bm_PlotEvalMean( bm.out      = myBiomodModel, #a BIOMOD.models.out or BIOMOD.ensemble.models.out object that can be obtained with the BIOMOD_Modeling or BIOMOD_EnsembleModeling functions 
                        metric.eval = c("ROC", "TSS"),  #a vector containing evaluation metric names to be used, must be among ROC, TSS, KAPPA, ACCURACY, BIAS, POD, FAR, POFD, SR, CSI, ETS, HK, HSS, OR, ORSS
                        dataset     = "calibration",    #a character corresponding to the dataset upon which evaluation metrics have been calculated and that is to be represented, must be among calibration, validation, evaluation////
                        group.by    = "algo",           #a character corresponding to the way kept models will be combined to compute mean and sd evaluation scores, must be among full.name, PA, run, algo (if bm.out is a BIOMOD.models.out object), or full.name, merged.by.PA, merged.by.run, merged.by.algo (if bm.out is a BIOMOD.ensemble.models.out object)
                        do.plot     = FALSE )

# extract evaluations metrics for the "validation" dataset
g.2 <- bm_PlotEvalMean( bm.out      = myBiomodModel,
                        metric.eval = c("ROC", "TSS"),
                        dataset     = "validation",
                        group.by    = "algo",
                        do.plot     = FALSE )

# extract evaluations metrics for the "evaluation" dataset
g.3 <- bm_PlotEvalMean( bm.out      = myBiomodModel,
                        metric.eval = c("ROC", "TSS"),
                        dataset     = "evaluation",
                        group.by    = "algo",
                        do.plot     = FALSE )


# create data.frame for plotting
g.data_Eval          <- rbind( g.1$tab[1,], g.2$tab[1,], g.3$tab[1,] )

g.data_Eval[,"name"] <- c( "Calib.",
                          "Valid.",
                          "Eval."  )
g.data_Eval[,"low1"] <- g.data_Eval[,"mean1"] - g.data_Eval[,"sd1"]
g.data_Eval[,"up1"]  <- g.data_Eval[,"mean1"] + g.data_Eval[,"sd1"]
g.data_Eval[,"low2"] <- g.data_Eval[,"mean2"] - g.data_Eval[,"sd2"]
g.data_Eval[,"up2"]  <- g.data_Eval[,"mean2"] + g.data_Eval[,"sd2"]



# ggplot2 graphic
g.eval <- ggplot( data = g.data_Eval,
                  aes( x      = mean1,
                       y      = mean2,
                       colour = name) )   +
  geom_pointrange( aes( xmin = low1,
                        xmax = up1 ) )    +
  geom_linerange(  aes( ymin = low2,
                        ymax = up2 ) )    +
  labs( x = "ROC",
        y = "TSS",
        subtitle = "Evaluation metrics" ) +
  xlim( 0, 1 )                            +
  ylim( 0, 1 )                            +
  theme(axis.text    = element_text ( size = 14 ),
        axis.title   = element_text ( size = 16, face = "bold" ),
        legend.title = element_text ( size = 14, face = "bold" ),
        legend.text  = element_text ( size = 12) )


png( paste0("Eval_mean", "_", cfg$Single_Models ,".png"), width = 1500, height = 1000 )
g.eval
dev.off()

# save evaluation scores on hard drive
capture.output( myBiomodModel@models.evaluation@val,
                file = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run, '.', 'Eval_mean', '_', cfg$Single_Models, '.csv'))







## Importance of variables ----

# extract values for each variable's importance
g.data_Var <- bm_PlotVarImpBoxplot( bm.out   = myBiomodModel,
                                group.by = c('expl.var', 'algo', 'algo'), #a 3-length vector containing the way kept models will be represented, must be among full.name, PA, run, algo, expl.var (if bm.out is a BIOMOD.models.out object), or full.name, merged.by.PA, merged.by.run, merged.by.algo, expl.var (if bm.out is a BIOMOD.ensemble.models.out object)
                                do.plot  = FALSE                          #a logical value defining whether the plot is to be rendered or not
)

# ggplot2 graphic
g.VarImp <- ggplot( data = g.data_Var$tab,
                    aes( x = reorder( expl.var, var.imp ),
                         y = var.imp * 100 ) ) +
  geom_boxplot( colour = "firebrick" )         +
  labs( x        = "Explanatory variables",
        y        = "Importance (%)",
        subtitle = "Variable importance" )     +
  theme(axis.text    = element_text ( size = 14 ),
        axis.title   = element_text ( size = 16, face = "bold" ),
        legend.title = element_text ( size = 14, face = "bold" ),
        legend.text  = element_text ( size = 12) )

png( paste0("Var_Imp", "_", cfg$Single_Models ,".png"), width = 1500, height = 1000 )
g.VarImp
dev.off()

# save variable importance scores on hard drive
capture.output( g.data_Var$tab,
                file = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run, '.', 'Var_Imp', '_', cfg$Single_Models, '.csv'))







## --- Represent response curves 

png( "Resp_curv.png", width = 1000, height = 800 )
bm_PlotResponseCurves( bm.out        = myBiomodModel, 
                       models.chosen = get_built_models( myBiomodModel,
                                                         full.name = "NARW_allData_allRun_GBM" ), #a vector containing model names to be kept, must be either all or a sub-selection of model names that can be obtained with the get_built_models function
                       fixed.var     = 'median',                                                  #a character corresponding to the statistic to be used to fix as constant the remaining variables other than the one used to predict response, must be either mean, median, min, max
                       do.bivariate  = FALSE,                                                     #a logical value defining whether the response curves are to be represented in 3 dimensions (meaning 2 explanatory variables at a time) or not (meaning only 1)
                       do.plot       = TRUE )
dev.off()









###############---------- 
##        Spatial Projections
##############---------

# Load data
#Ouverture des rasters avec le package stars
SST8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/SST/SST_L4_2020_8days.nc",
                         var = "sst", make_units = TRUE, make_time = TRUE, proxy = F)


PP8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/PP/PP_L3_2020_8days.nc",
                        var = c("pp"), make_units = TRUE, make_time = TRUE, proxy = F)


PAR8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/_PAR_/_PAR__PAR__L3_2020_8days.nc",
                         var = "PAR_mean", make_units = TRUE, make_time = TRUE, proxy = F)


CHL8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/CHL-PCA/CHL-PCA_L3_2020_8days.nc",
                         var = "CHL-PCA_mean", make_units = TRUE, make_time = TRUE, proxy = F)


Bathy <- read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Bathymetry/REGRIDDED1km_gebco_2023_n54.0_s37.0_w-74.0_e-45.0.nc", var = 'elevation', make_units = TRUE, make_time = TRUE, proxy = FALSE)

BathyCRS=st_crs(SST8days2020)
st_crs(Bathy)=BathyCRS

#######################
#Set dataframe for the proj

nb_weeks = 24  #the number of week(s) you want to project on

#for the Gulf of SL

ENV8days2020_proj_8days = (ENV8days2020[ ,400:1700 , 350:1000, 17:40])     #crops the SL with the rights week(s)
ENV8days2020_proj_8days_df = as.data.frame(ENV8days2020_proj_8days)        #makes a dataframe 

Bathy_proj = Bathy[, 400:1700 , 350:1000]
Bathy_proj_df = as.data.frame(Bathy_proj)
Bathy_proj_8days = replicate(nb_weeks, Bathy_proj_df, simplify = F)        #crops the SL with the rights week(s)
Bathy_proj_8days_df = as.data.frame(do.call(rbind, Bathy_proj_8days))      #makes a dataframe

CHL8days2020_proj_8days = (CHL8days2020[ ,400:1700 , 350:1000, 17:40])     #crops the SL with the rights week(s)
CHL8days2020_proj_8days_df = as.data.frame(CHL8days2020_proj_8days)        #makes a dataframe

ENV8days2020_proj_all_8days= cbind(ENV8days2020_proj_8days_df, Bathy_proj_8days_df$elevation, CHL8days2020_proj_8days_df$CHL.PCA_mean)                                    #makes a dataframe
ENV8days2020_proj_all_8days_df = rename(ENV8days2020_proj_all_8days, bathy = 'Bathy_proj_8days_df$elevation', CHL = 'CHL8days2020_proj_8days_df$CHL.PCA_mean')
ENV8days2020_proj_all_8days_df$bathy = as.numeric(ENV8days2020_proj_all_8days_df$bathy)



######################
## set data for projection ----

proj.model  <- myBiomodModel #_v1.0
proj.name   <- paste0('Proj', '.', cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run, 'GAM', '.' ,'20to31')              # don't forget to change the name according to the model
proj.env    <- ENV8days2020_proj_all_8days_df[ , env.var ]
proj.xy     <- ENV8days2020_proj_all_8days_df[, 1:2]
proj.chosen <- get_built_models( myBiomodModel,
                                 full.name = "NARW_allData_allRun_GAM" )   # don't forget to change the name according to the model



# run the projection  
myBiomodProj_GAM_proj <- BIOMOD_Projection(
  bm.mod              = proj.model,  #a BIOMOD.models.out object returned by the BIOMOD_Modeling function
  proj.name           = proj.name,   #a character corresponding to the name (ID) of the projection set (a new folder will be created within the simulation folder with this name)
  new.env             = proj.env,    #a matrix, data.frame or SpatRaster object containing the new explanatory variables (in columns or layers, with names matching the variables names given to the BIOMOD_FormatingData function to build bm.mod) that will be used to project the species distribution model(s)
  new.env.xy          = proj.xy,     #if new.env is a matrix or a data.frame, a 2-columns matrix or data.frame containing the corresponding X and Y coordinates that will be used to project the species distribution model(s)
  models.chosen       = proj.chosen, #a vector containing model names to be kept, must be either all or a sub-selection of model names that can be obtained with the get_built_models function
  metric.binary       = NULL,        #a vector containing evaluation metric names to be used to transform prediction values into binary values based on models evaluation scores obtained with the BIOMOD_Modeling function. Must be among all (same evaluation metrics than those of bm.mod) or ROC, TSS, KAPPA, ACCURACY, BIAS, POD, FAR, POFD, SR, CSI, ETS, HK, HSS, OR, ORSS
  metric.filter       = NULL,        #a vector containing evaluation metric names to be used to transform prediction values into filtered values based on models evaluation scores obtained with the BIOMOD_Modeling function. Must be among all (same evaluation metrics than those of bm.mod) or ROC, TSS, KAPPA, ACCURACY, BIAS, POD, FAR, POFD, SR, CSI, ETS, HK, HSS, OR, ORSS
  compress            = TRUE,        #a logical or a character value defining whether and how objects should be compressed when saved on hard drive. Must be either TRUE, FALSE, xz or gzip (see Details)
  build.clamping.mask = F,        #a logical value defining whether a clamping mask should be built and saved on hard drive or not (see Details)
  nb.cpu              = 20,           #an integer value corresponding to the number of computing resources to be used to parallelize the single models computation
  seed.val            = 42           #an integer value corresponding to the new seed value to be set
)



###############################
# get information from the plot
proj.gg <- plot(
  x           = myBiomodProj_GAM_proj, #a BIOMOD.projection.out object
  #coord       = NULL,                  #a 2-columns data.frame containing the corresponding X and Y
  plot.output = "list",                #(optional, default facet) a character determining the type of output: with plot.output = 'list' the function will return a list of plots (one plot per model) ; with 'facet' ; with plot.output = 'facet' the function will return a single plot with all asked projections as facet.
  do.plot     = FALSE,                 #(optional, default TRUE) a boolean determining whether the plot should be displayed or just returned.
  #std         = TRUE,                  #(optional, default TRUE) a boolean controlling the limits of the color scales. With std = TRUE color scales are displayed between 0 and 1 (or 1000). With std = FALSE color scales are displayed between 0 and the maximum value observed.
  #scales      = "fixed",               #(optional, default fixed) a character determining whether x and y scales are shared among facet. Argument passed to facet_wrap. Possible values: 'fixed', 'free_x', 'free_y', 'free'.
  #size        = 0.3,                   #(optional, default 0.75) a numeric determing the size of points on the plots and passed to geom_point.
  #maxcell     = 5e+05,                 #maximum number of cells to plot. Argument transmitted to plot.
)


n = 1  # gives the model you want to project with

proj.na   <- is.na( proj.gg[[n]]$data$pred )
proj.data <- proj.gg[[n]]$data[!proj.na,]


#give all 8days periods of the according weeks of projection to select the right Sightings to add to the plot
date8days = paste0(c("2020-05-08", "2020-05-16", "2020-05-24" ,"2020-06-01", "2020-06-09", "2020-06-17", "2020-06-25", "2002-07-03", "2020-07-11", "2020-07-19", "2020-07-27", "2020-08-04", "2020-08-12", "2020-08-20", "2020-08-28", "2020-09-05", "2020-09-13", "2020-09-21", "2020-09-29", "2020-10-07", "2020-10-15", "2020-10-23", "2020-10-31", "2020-11-08"))


# Add the presences
condition1 <- Eval_Sightings$Date8days %in% date8days
condition2 <- Eval_Sightings$SPECCODE == 1

# Combinez les deux conditions avec l'opérateur logique '&'
combined_condition <- condition1 & condition2

# Utilisez 'which' pour obtenir les emplacements (indices) correspondants
pres.i <- which(combined_condition)
Eval_Sightings_sf = st_as_sf(Eval_Sightings, coords= c('Long', "Lat"), crs = 4326)
pres.xy     <- st_coordinates( Eval_Sightings_sf$geometry[pres.i] ) %>% 
  as.data.frame()




# Add the absences
condition1 <- Eval_Sightings$Date8days %in% date8days
condition2.1 <- Eval_Sightings$SPECCODE == 0

# Combinez les deux conditions avec l'opérateur logique '&'
combined_condition2 <- condition1 & condition2.1

# Utilisez 'which' pour obtenir les emplacements (indices) correspondants
abs.i <- which(combined_condition2)

Eval_Sightings_sf = st_as_sf(Eval_Sightings, coords= c('Long', "Lat"), crs = 4326)
abs.xy     <- st_coordinates( Eval_Sightings_sf$geometry[abs.i] ) %>% 
  as.data.frame()




# Add opportunistic sightings
opp = readRDS("/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Sightings/8days/Opportunistic/St_Laurent/2020/Sightings_8days_Opportunistic_St_Laurent_2020.RDS")
#opp.i =  which(as.character(opp$Date8days) == date8days)
condition1.1 <- opp$Date8days %in% date8days

opp.i <- which(condition1.1)

opp.xy     <- st_coordinates( opp$geometry[opp.i] ) %>% 
  as.data.frame()



# make our own ggplot of the projection
g.proj <-  basemap(
  limits = c(-70, -50, 45, 52), # sets the plot dimensions
  bathymetry = F,             # adds the bathymetry
  rotate = TRUE,                 # rotates the plot
  land.border.col = "black",
  land.col = "white",
  lon.interval = 10, 
  lat.interval = 5,
  bathy.border.col = NA) + 
  
  #ggplot()                            +
  geom_spatial_point( data = proj.data,      #geo_SPATIAL_point when using the basemap function
                      aes( x     = x,
                           y     = y,
                           color = pred*1e-3 ),
                      size = 0.8 )                    +
  scale_color_viridis_c( option = "plasma",
                         limits = c(0,1),
                         na.value = "white" ) +
  geom_spatial_point( data = pres.xy,
                      aes( x     = X,
                           y     = Y ),
                      color = "deepskyblue",
                      shape = 8,
                      size  = 5      )                   +
  geom_spatial_point( data = opp.xy,
                      aes( x     = X,
                           y     = Y ),
                      color = "green",
                      shape = 8,
                      size  = 5 )                   +
  theme(axis.text    = element_text ( size = 15 ),
        axis.title   = element_text ( size = 20, face = "bold" ),
        legend.title = element_blank()           ,
        legend.text  = element_text ( size = 15) ,
        legend.position = c(0.1, 0.5)            ,
        title        = element_text ( size = 20 ) ) +
  labs( x        = "Longitude",
        y        = "Latitude",
        #subtitle = "Predicted probability of presence",
        color    = "Prob.")                 


png( paste0( "ggplot", EMproj.name, '.', myBiomodEMProj@models.projected[[n]], ".", week_nb, ".png"), width = 1500, height = 1000 )
g.proj
dev.off()



###############-------
######Frequency of presences
proj_df = as.data.frame(cbind(EMproj.xy$x, EMproj.xy$y, myBiomodEMProj@proj.out@val$pred[myBiomodEMProj@proj.out@val$full.name == "NARW_EMcaByROC_mergedData_mergedRun_mergedAlgo"]))     #paste0(cfg$Algo)[[n]]])) subsets values of predictions of the selected model to a dataframe

colnames(proj_df) = c('Long', 'Lat', 'Pred')

### conversion en sf objet
proj_sf = st_as_sf(proj_df, coords= c('Long', "Lat"), crs = 4326)   # all projections/predictions
pres.xy_sf = st_as_sf(pres.xy, coords= c('X', "Y"), crs = 4326)     # actual sightings

###conversion en star objet
proj_stars <- st_rasterize(proj_sf %>% dplyr::select(Pred, geometry))  # predictions as star object
extract_proj = st_extract(proj_stars, pres.xy_sf)  # extraction of the prob of prediction at sightings localizations

#hist(extract_proj$Pred, xlim = c(0, 1000), xlab = 'Prediction value at Sightings spots', main = 'Frequency of probability of presence at Sightings spots')



###############-------
######Frequency of absences
abs.xy_sf = st_as_sf(abs.xy,        # actual sightings
                     coords= c('X', "Y"), 
                     crs = 4326)    

###conversion en star objet
extract_proj_abs = st_extract(proj_stars, abs.xy_sf)  # extraction of the prob of prediction at sightings localizations

#hist(extract_proj_abs$Pred, xlim = c(0, 1000), xlab = 'Prediction value at NO Sightings spots', main = 'Frequency of probability of absence at NO Sightings spots')



###############-------
######Frequency of opportunistic
opp.xy_sf = st_as_sf(opp.xy,        # actual sightings
                     coords= c('X', "Y"), 
                     crs = 4326)    

###conversion en star objet
extract_proj_opp = st_extract(proj_stars, opp.xy_sf)  # extraction of the prob of prediction at sightings localizations




###############-------
###### Créez l'histogramme en spécifiant différentes esthétiques pour chaque série de données
donnees_combines <- bind_rows(list(Presence = extract_proj, Absence = extract_proj_abs, Opp = extract_proj_opp), .id = 'Occ')


histogramme = ggplot(
  donnees_combines, 
  aes(Pred*1e-3 , fill = Occ)) +
  geom_histogram(bins = 20, 
                 position = "identity", 
                 alpha = 0.5,
                 mapping = aes(y = stat(ndensity))) +
  labs( x        = "Probability of presence",
        y        = "Density")        +
  theme(axis.text    = element_text ( size = 30 ),
        axis.title   = element_text ( size = 20, face = "bold" ),
        legend.title = element_blank()           ,
        legend.text  = element_text ( size = 15)  ,
        title        = element_text ( size = 20 ) )

histogramme





## XXXXXXXXXXXXXXXXXXXX ----
## --------------------   ENSEMBLE MODELS    --------------------


myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModel,  #a BIOMOD.models.out object returned by the BIOMOD_Modeling function
  models.chosen = "all",  #a vector containing model names to be kept, must be either all or a sub-selection of model names that can be obtained with the get_built_models function
  em.by = "all", #a character corresponding to the way kept models will be combined to build the ensemble models, must be among PA_dataset+repet, PA_dataset+algo, PA_dataset, algo, all
  metric.select = cfg$Metrics_Eval, # with metric.select.thresh to exclude single models "KAPPA", "TSS", "ROC"....
  metric.select.thresh = c(0.7, 0.3), #A vector of numeric values corresponding to the minimum scores (one for each metric.select) below which single models will be excluded from the ensemble model building
  metric.select.table = NULL, #If metric.select = 'user.defined', a data.frame containing evaluation scores calculated for each single models and that will be compared to metric.select.thresh values 
  #to exclude some of them from the ensemble model building, with metric.select rownames, and models.chosen colnames
  metric.eval = cfg$Metrics_Eval, #"KAPPA", "TSS", 
  var.import = 5, #An integer corresponding to the number of permutations to be done for each variable to estimate variable importance
  EMwmean.decay = 0.6, #A value defining the relative importance of the weights (if prob.mean.weight = TRUE). 
  #A high value will strongly discriminate good models from the bad ones (see Details), 
  #while proportional will attribute weights proportionally to the models evaluation scores
  em.algo = cfg$Algo,#, 'EMci'),
  nb.cpu = 15, #An integer value corresponding to the number of computing resources to be used to parallelize the single models computation
  seed.val = 42, #An integer value corresponding to the new seed value to be set
  do.progress = TRUE #A logical value defining whether the progress bar is to be rendered or not
)


## Evaluation scores ----

# create ggplot2 object
# extract evaluations metrics for the "calibration" dataset
g.1 <- bm_PlotEvalMean( bm.out      = myBiomodEM, #a BIOMOD.models.out or BIOMOD.ensemble.models.out object that can be obtained with the BIOMOD_Modeling or BIOMOD_EnsembleModeling functions 
                        metric.eval = c("ROC", "TSS"),  #a vector containing evaluation metric names to be used, must be among ROC, TSS, KAPPA, ACCURACY, BIAS, POD, FAR, POFD, SR, CSI, ETS, HK, HSS, OR, ORSS
                        dataset     = "calibration",    #a character corresponding to the dataset upon which evaluation metrics have been calculated and that is to be represented, must be among calibration, validation, evaluation////
                        group.by    = "algo",           #a character corresponding to the way kept models will be combined to compute mean and sd evaluation scores, must be among full.name, PA, run, algo (if bm.out is a BIOMOD.models.out object), or full.name, merged.by.PA, merged.by.run, merged.by.algo (if bm.out is a BIOMOD.ensemble.models.out object)
                        do.plot     = FALSE )

# extract evaluations metrics for the "validation" dataset
g.2 <- bm_PlotEvalMean( bm.out      = myBiomodEM,
                        metric.eval = c("ROC", "TSS"),
                        dataset     = "validation",
                        group.by    = "algo",
                        do.plot     = FALSE )

# extract evaluations metrics for the "evaluation" dataset
g.3 <- bm_PlotEvalMean( bm.out      = myBiomodEM,
                        metric.eval = c("ROC", "TSS"),
                        dataset     = "evaluation",
                        group.by    = "algo",
                        do.plot     = FALSE )


# create data.frame for plotting
g.data_Eval          <- rbind( g.1$tab[1,] , g.3$tab[1,] )

g.data_Eval[,"name"] <- c( "Calib.",
                           #"Valid.",
                           "Eval."  )
g.data_Eval[,"low1"] <- g.data_Eval[,"mean1"] - g.data_Eval[,"sd1"]
g.data_Eval[,"up1"]  <- g.data_Eval[,"mean1"] + g.data_Eval[,"sd1"]
g.data_Eval[,"low2"] <- g.data_Eval[,"mean2"] - g.data_Eval[,"sd2"]
g.data_Eval[,"up2"]  <- g.data_Eval[,"mean2"] + g.data_Eval[,"sd2"]



# ggplot2 graphic
g.eval <- ggplot( data = g.data_Eval,
                  aes( x      = mean1,
                       y      = mean2,
                       colour = name   )) +
         geom_point(size = 8)             +
  geom_pointrange( aes( xmin = low1,
                        xmax = up1  ) )   +
  geom_linerange(  aes( ymin = low2,
                        ymax = up2  ) )   +
  labs( x = "ROC",
        y = "TSS")                        +
        #subtitle = "Evaluation metrics")  +
  xlim( 0, 1 )                            +
  ylim( 0, 1 )                            +
  theme(axis.text    = element_text ( size = 30 ),
        axis.title   = element_text ( size = 35, face = "bold" ),
        legend.title = element_blank()           ,
        legend.text  = element_text ( size = 25) ,
        legend.position = c(0.1, 0.9)            ,
        title        = element_text ( size = 20 ) )


png( paste0("Eval_mean", "_", "Ensemble" , ".png"), width = 1500, height = 1000 )
g.eval
dev.off()


# save evaluation scores on hard drive
capture.output( myBiomodEM@models.evaluation@val,
                file = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run, '.', 'Eval_mean', '_', '.csv'))



## Importance of variables ----

# extract values for each variable's importance
g.data_Var <- bm_PlotVarImpBoxplot( bm.out   = myBiomodEM,
                                    group.by = c('expl.var', 'algo', 'algo'), #a 3-length vector containing the way kept models will be represented, must be among full.name, PA, run, algo, expl.var (if bm.out is a BIOMOD.models.out object), or full.name, merged.by.PA, merged.by.run, merged.by.algo, expl.var (if bm.out is a BIOMOD.ensemble.models.out object)
                                    do.plot  = FALSE                          #a logical value defining whether the plot is to be rendered or not
)

# ggplot2 graphic
g.VarImp <- ggplot( data = g.data_Var$tab,
                    aes( x = reorder( expl.var, var.imp ),
                         y = var.imp * 100 ) ) +
  geom_boxplot( colour = "firebrick", size = 1.5 )         +
  labs( x        = "Explanatory variables",
        y        = "Importance (%)")            +
        #subtitle = "Variable importance" )     +
  theme(axis.text    = element_text ( size = 30 ),
        axis.title   = element_text ( size = 35, face = "bold" ),
        legend.title = element_blank(),
        legend.text  = element_text ( size = 2),
        title        = element_text ( size = 20 ) )

png( paste0("Var_Imp",  "Ensemble" , "_" ,".png"), width = 1500, height = 1000 )
g.VarImp
dev.off()

# save variable importance scores on hard drive
capture.output( g.data_Var$tab,
                file = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run, '.', 'Var_Imp', '_', '.csv'))



## --- Represent response curves 

png( "Resp_curv.png", width = 1000, height = 800 )
bm_PlotResponseCurves( bm.out        = myBiomodEM, 
                       models.chosen = get_built_models( myBiomodEM,
                                                         full.name = c("NARW_EMmeanByROC_mergedData_mergedRun_mergedAlgo", 
                                                                       "NARW_EMcaByROC_mergedData_mergedRun_mergedAlgo" ) ), #a vector containing model names to be kept, must be either all or a sub-selection of model names that can be obtained with the get_built_models function
                       fixed.var     = 'median',                                                  #a character corresponding to the statistic to be used to fix as constant the remaining variables other than the one used to predict response, must be either mean, median, min, max
                       do.bivariate  = FALSE,                                                     #a logical value defining whether the response curves are to be represented in 3 dimensions (meaning 2 explanatory variables at a time) or not (meaning only 1)
                       do.plot       = TRUE )
dev.off()





###############---------- 
##        Spatial Projections
##############---------

# Load data
#Ouverture des rasters avec le package stars
SST8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/SST/SST_L4_2020_8days.nc",
                         var = "sst", make_units = TRUE, make_time = TRUE, proxy = F)


PP8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/PP/PP_L3_2020_8days.nc",
                        var = c("pp"), make_units = TRUE, make_time = TRUE, proxy = F)


PAR8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/_PAR_/_PAR__PAR__L3_2020_8days.nc",
                         var = "PAR_mean", make_units = TRUE, make_time = TRUE, proxy = F)


CHL8days2020 = read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/CHL-PCA/CHL-PCA_L3_2020_8days.nc",
                         var = "CHL-PCA_mean", make_units = TRUE, make_time = TRUE, proxy = F)


Bathy <- read_ncdf("/Ext_16T_andromede/0_ARCTUS_Projects/15_SmartWhales/data/GC/ATLNW/DATA/Bathymetry/REGRIDDED1km_gebco_2023_n54.0_s37.0_w-74.0_e-45.0.nc", var = 'elevation', make_units = TRUE, make_time = TRUE, proxy = FALSE)

BathyCRS=st_crs(SST8days2020)
st_crs(Bathy)=BathyCRS

#######################
#Set dataframe for the proj

nb_weeks = 24  #the number of week(s) you want to project on

#for the Gulf of SL

ENV8days2020_proj_8days = (ENV8days2020[ ,400:1700 , 350:1000, 17:40])     #crops the SL with the rights week(s)
ENV8days2020_proj_8days_df = as.data.frame(ENV8days2020_proj_8days)        #makes a dataframe 

Bathy_proj = Bathy[, 400:1700 , 350:1000]
Bathy_proj_df = as.data.frame(Bathy_proj)
Bathy_proj_8days = replicate(nb_weeks, Bathy_proj_df, simplify = F)        #crops the SL with the rights week(s)
Bathy_proj_8days_df = as.data.frame(do.call(rbind, Bathy_proj_8days))      #makes a dataframe

CHL8days2020_proj_8days = (CHL8days2020[ ,400:1700 , 350:1000, 17:40])     #crops the SL with the rights week(s)
CHL8days2020_proj_8days_df = as.data.frame(CHL8days2020_proj_8days)        #makes a dataframe

ENV8days2020_proj_all_8days= cbind(ENV8days2020_proj_8days_df, Bathy_proj_8days_df$elevation, CHL8days2020_proj_8days_df$CHL.PCA_mean)                                    #makes a dataframe
ENV8days2020_proj_all_8days_df = rename(ENV8days2020_proj_all_8days, bathy = 'Bathy_proj_8days_df$elevation', CHL = 'CHL8days2020_proj_8days_df$CHL.PCA_mean')
ENV8days2020_proj_all_8days_df$bathy = as.numeric(ENV8days2020_proj_all_8days_df$bathy)




######################
## set data for projection ----

EMproj.model  <- myBiomodEM #_v1.0
EMproj.name   <- paste0('EMProj', '.', cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run, '.', '20to31')             # don't forget to change the name according to the period of proj
EMproj.env    <- ENV8days2020_proj_all_8days_df[ , env.var ]
EMproj.xy     <- ENV8days2020_proj_all_8days_df[, 1:2]
EMproj.chosen <- get_built_models( NARW.8days.St_Laurent.2017to2020.2020.v1.001.ensemble.models.out,
                                full.name = "NARW_EMcaByROC_mergedData_mergedRun_mergedAlgo" )  # don't forget to change the name according to the model




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


myBiomodEMProj


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



# number of the chosen model to be projected
n = 1

## --- save proj

proj.na   <- is.na( proj.gg[[n]]$data$pred )
proj.data <- proj.gg[[n]]$data[!proj.na,]


# give all 8days periods of the according weeks of projection to select the right Sightings to add to the plot
date8days = paste0(c("2020-05-08", "2020-05-16", "2020-05-24" ,"2020-06-01", "2020-06-09", "2020-06-17", "2020-06-25", "2002-07-03", "2020-07-11", "2020-07-19", "2020-07-27", "2020-08-04", "2020-08-12", "2020-08-20", "2020-08-28", "2020-09-05", "2020-09-13", "2020-09-21", "2020-09-29", "2020-10-07", "2020-10-15", "2020-10-23", "2020-10-31", "2020-11-08"))



# Add the presences
condition1 <- Eval_Sightings$Date8days %in% date8days
condition2 <- Eval_Sightings$SPECCODE == 1

# Combinez les deux conditions avec l'opérateur logique '&'
combined_condition <- condition1 & condition2

# Utilisez 'which' pour obtenir les emplacements (indices) correspondants
pres.i <- which(combined_condition)
Eval_Sightings_sf = st_as_sf(Eval_Sightings, coords= c('Long', "Lat"), crs = 4326)
pres.xy     <- st_coordinates( Eval_Sightings_sf$geometry[pres.i] ) %>% 
  as.data.frame()



# Add the absences
condition1 <- Eval_Sightings$Date8days %in% date8days
condition2.1 <- Eval_Sightings$SPECCODE == 0

# Combinez les deux conditions avec l'opérateur logique '&'
combined_condition2 <- condition1 & condition2.1

# Utilisez 'which' pour obtenir les emplacements (indices) correspondants
abs.i <- which(combined_condition2)

Eval_Sightings_sf = st_as_sf(Eval_Sightings, coords= c('Long', "Lat"), crs = 4326)
abs.xy     <- st_coordinates( Eval_Sightings_sf$geometry[abs.i] ) %>% 
  as.data.frame()




# Add opportunistic sightings
opp = readRDS("/Ext_16T_andromede/Extern_workdir/thieryf_WD/DataFrame/Sightings/8days/Opportunistic/St_Laurent/2020/Sightings_8days_Opportunistic_St_Laurent_2020.RDS")
condition1.1 <- opp$Date8days %in% date8days

opp.i <- which(condition1.1)

opp.xy     <- st_coordinates( opp$geometry[opp.i] ) %>% 
  as.data.frame()




# make our own ggplot of the projection
g.proj <-  basemap(
  limits = c(-70, -50, 45, 52), # sets the plot dimensions
  bathymetry = F,             # adds the bathymetry
  rotate = TRUE,                 # rotates the plot
  land.border.col = "black",
  land.col = "white",
  lon.interval = 10, 
  lat.interval = 5,
  bathy.border.col = NA) + 
  
  #ggplot()                            +
  geom_spatial_point( data = proj.data,
              aes( x     = x,
                   y     = y,
                   color = pred*1e-3 ),
              size = 0.8 )                    +
  scale_color_viridis_c( option = "plasma",
                         limits = c(0,1),
                         na.value = "white" ) +
  geom_spatial_point( data = pres.xy,
              aes( x     = X,
                   y     = Y ),
              color = "deepskyblue",
              shape = 8,
              size  = 5      )                   +
  geom_spatial_point( data = opp.xy,
              aes( x     = X,
                   y     = Y ),
              color = "green",
              shape = 8,
              size  = 5 )                   +
  theme(axis.text    = element_text ( size = 15 ),
        axis.title   = element_text ( size = 20, face = "bold" ),
        legend.title = element_blank()           ,
        legend.text  = element_text ( size = 15) ,
        legend.position = c(0.1, 0.5)            ,
        title        = element_text ( size = 20 ) ) +
  labs( x        = "Longitude",
        y        = "Latitude",
        #subtitle = "Predicted probability of presence",
        color    = "Prob.")                 
  

png( paste0( "ggplot", EMproj.name, '.', myBiomodEMProj@models.projected[[n]], ".", week_nb, ".png"), width = 1500, height = 1000 )
g.proj
dev.off()



###############-------
######Frequency of presences
proj_df = as.data.frame(cbind(EMproj.xy$x, EMproj.xy$y, myBiomodEMProj@proj.out@val$pred[myBiomodEMProj@proj.out@val$full.name == "NARW_EMcaByROC_mergedData_mergedRun_mergedAlgo"]))     #paste0(cfg$Algo)[[n]]])) subsets values of predictions of the selected model to a dataframe

colnames(proj_df) = c('Long', 'Lat', 'Pred')

### conversion en sf objet
proj_sf = st_as_sf(proj_df, coords= c('Long', "Lat"), crs = 4326)   # all projections/predictions
pres.xy_sf = st_as_sf(pres.xy, coords= c('X', "Y"), crs = 4326)     # actual sightings

###conversion en star objet
proj_stars <- st_rasterize(proj_sf %>% dplyr::select(Pred, geometry))  # predictions as star object
extract_proj = st_extract(proj_stars, pres.xy_sf)  # extraction of the prob of prediction at sightings localizations

#hist(extract_proj$Pred, xlim = c(0, 1000), xlab = 'Prediction value at Sightings spots', main = 'Frequency of probability of presence at Sightings spots')



###############-------
######Frequency of absences
abs.xy_sf = st_as_sf(abs.xy,        # actual sightings
                     coords= c('X', "Y"), 
                     crs = 4326)    

###conversion en star objet
extract_proj_abs = st_extract(proj_stars, abs.xy_sf)  # extraction of the prob of prediction at sightings localizations

#hist(extract_proj_abs$Pred, xlim = c(0, 1000), xlab = 'Prediction value at NO Sightings spots', main = 'Frequency of probability of absence at NO Sightings spots')



###############-------
######Frequency of opportunistic
opp.xy_sf = st_as_sf(opp.xy,        # actual sightings
                     coords= c('X', "Y"), 
                     crs = 4326)    

###conversion en star objet
extract_proj_opp = st_extract(proj_stars, opp.xy_sf)  # extraction of the prob of prediction at sightings localizations




###############-------
###### Créez l'histogramme en spécifiant différentes esthétiques pour chaque série de données
donnees_combines <- bind_rows(list(Presence = extract_proj, Absence = extract_proj_abs, Opp = extract_proj_opp), .id = 'Occ')


histogramme = ggplot(
  donnees_combines, 
  aes(Pred*1e-3 , fill = Occ)) +
  geom_histogram(bins = 20, 
                 position = "identity", 
                 alpha = 0.5,
                 mapping = aes(y = stat(ndensity))) +
  labs( x        = "Probability of presence",
        y        = "Density")        +
  theme(axis.text    = element_text ( size = 30 ),
        axis.title   = element_text ( size = 20, face = "bold" ),
        legend.title = element_blank()           ,
        legend.text  = element_text ( size = 15)  ,
        title        = element_text ( size = 20 ) )

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
myBiomodEM_EMmean_evalpredictROC = myBiomodEM@models.prediction.eval@val$pred.eval[myBiomodEM@models.prediction.eval@val$algo == 'EMmean' & myBiomodEM@models.prediction.eval@val$filtered.by == 'ROC']
myBiomodEM_EMmean_evalpredictTSS = myBiomodEM@models.prediction.eval@val$pred.eval[myBiomodEM@models.prediction.eval@val$algo == 'EMmean' & myBiomodEM@models.prediction.eval@val$filtered.by == 'TSS']

## binding with the date and coordinates
predict_eval_SL_2017to2020 = cbind(OPP_ENV_2017to2020$LONGITUDE, OPP_ENV_2017to2020$LATITUDE, as.data.frame(myBiomodEM_EMmean_evalpredictROC), as.data.frame(myBiomodEM_EMmean_evalpredictTSS), OPP_ENV_2017to2020$Date8days)
predict_eval_SL_2017to2020$year = year(predict_eval_SL_2017to2020$`OPP_ENV_2017to2020$Date8days`)
predict_eval_SL_2017to2020$month = month(predict_eval_SL_2017to2020$`OPP_ENV_2017to2020$Date8days`)
predict_eval_SL_2017to2020 = rename(predict_eval_SL_2017to2020, Long = 'OPP_ENV_2017to2020$LONGITUDE', Lat = 'OPP_ENV_2017to2020$LATITUDE', EMmean_evalpredictROC = 'myBiomodEM_EMmean_evalpredictROC', EMmean_evalpredictTSS = 'myBiomodEM_EMmean_evalpredictTSS', Date8days = 'OPP_ENV_2017to2020$Date8days' )

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
bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'algo')
dev.off()

png(paste0(directory, "/EM_Var_Imp.png"))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'algo'))
dev.off()

png(paste0(directory, "/EMResp_curv.png"), width = 1024, height = 768)
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM, full.name = c('RIWH_EMwmeanByROC_mergedData_mergedRun_mergedAlgo','RIWH_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo')),
                      fixed.var = 'median', 
                      do.bivariate = F)         
dev.off()

