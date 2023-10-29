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



#Load the yml file that corresponds
cfg = read_yaml('~/Documents/SDM_SIMBA/Code/SDM_Biomod2_Optim_lag.yml')

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

# Functions that will open and Row bind years for training dataset
#Sightings
Read_Sightings = function( cfg ) {
  
  S_filename <- paste0( 'Sightings',    '_',
                        'DFO',          '_', 
                        cfg$Year,       '_', 
                        '8days',        '.',
                        'rds'              )
  
  S_filepath <- paste0( cfg$initial_wd, '/',
                        'Sightings',    '/', 
                        S_filename         )
  
  x = lapply( S_filepath, readRDS ) |> 
    dplyr::bind_rows()
  
  return( x )
}

#Variables

 Read_PP= function( cfg ) {
  
   E_filename <- paste0( 'PP', '_',
                        'DFO',                       '_', 
                        cfg$Year,                    '_', 
                        '8days',                     '.',
                        'rds'              )
  
   E_filepath <- paste0( cfg$initial_wd,              '/',
                         'PP', '/' ,
                         E_filename         ) 
  
   x = lapply( E_filepath, readRDS ) |> 
     dplyr::bind_rows()
  
  return( x )
}

 
 Read_CHL= function( cfg ) {
   
   E_filename <- paste0( 'CHL', '_',
                         'DFO',                       '_', 
                         cfg$Year,                    '_', 
                         '8days',                     '.',
                         'rds'              )
   
   E_filepath <- paste0( cfg$initial_wd,              '/',
                         'CHL', '/' ,
                         E_filename         ) 
   
   x = lapply( E_filepath, readRDS ) |> 
     dplyr::bind_rows()
   
   return( x )
 }
 

 Read_PAR= function( cfg ) {
   
   E_filename <- paste0( 'PAR', '_',
                         'DFO',                       '_', 
                         cfg$Year,                    '_', 
                         '8days',                     '.',
                         'rds'              )
   
   E_filepath <- paste0( cfg$initial_wd,              '/',
                         'PAR', '/' ,
                         E_filename         ) 
   
   x = lapply( E_filepath, readRDS ) |> 
     dplyr::bind_rows()
   
   return( x )
 }
 
 
 Read_SST= function( cfg ) {
   
   E_filename <- paste0( 'SST', '_',
                         'DFO',                       '_', 
                         cfg$Year,                    '_', 
                         '8days',                     '.',
                         'rds'              )
   
   E_filepath <- paste0( cfg$initial_wd,              '/',
                         'SST', '/' ,
                         E_filename         ) 
   
   x = lapply( E_filepath, readRDS ) |> 
     dplyr::bind_rows()
   
   return( x )
 }
 
 
 Read_Bathy= function( cfg ) {
   
   E_filename <- paste0( 'Bathy', '_',
                         'DFO',                       '_', 
                         cfg$Year,                    '_', 
                         '8days',                     '.',
                         'rds'              )
   
   E_filepath <- paste0( cfg$initial_wd,              '/',
                         'Bathy', '/' ,
                         E_filename         ) 
   
   x = lapply( E_filepath, readRDS ) |> 
     dplyr::bind_rows()
   
   return( x )
 }
 

# Row binding data for Evaluation
#Sightings
Read_EvalSightings = function( cfg ) {
  
  EvalS_filename <- paste0( 'Sightings',    '_',
                            'DFO',          '_', 
                            cfg$Eval_Year,  '_', 
                            '8days',        '.',
                            'rds'              )
  
  EvalS_filepath <- paste0( cfg$initial_wd, '/',
                            'Sightings',    '/', 
                            EvalS_filename         )
  
  x = lapply( EvalS_filepath, readRDS ) |> 
    dplyr::bind_rows()
  
  return( x )
}

#Variables
Read_EvalPP = function( cfg ) {
  
  EvalE_filename <- paste0( 'PP', '_',
                            'DFO',                       '_', 
                            cfg$Eval_Year,               '_', 
                            '8days',                     '.',
                            'rds'              )
  
  EvalE_filepath <- paste0( cfg$initial_wd,              '/',
                            'PP', '/' ,
                            EvalE_filename         ) 
  
  x = lapply( EvalE_filepath, readRDS ) |> 
    dplyr::bind_rows()
  return( x )
}


Read_EvalCHL = function( cfg ) {
  
  EvalE_filename <- paste0( 'CHL', '_',
                            'DFO',                       '_', 
                            cfg$Eval_Year,               '_', 
                            '8days',                     '.',
                            'rds'              )
  
  EvalE_filepath <- paste0( cfg$initial_wd,              '/',
                            'CHL', '/' ,
                            EvalE_filename         ) 
  
  x = lapply( EvalE_filepath, readRDS ) |> 
    dplyr::bind_rows()
  return( x )
}


Read_EvalPAR = function( cfg ) {
  
  EvalE_filename <- paste0( 'PAR', '_',
                            'DFO',                       '_', 
                            cfg$Eval_Year,               '_', 
                            '8days',                     '.',
                            'rds'              )
  
  EvalE_filepath <- paste0( cfg$initial_wd,              '/',
                            'PAR', '/' ,
                            EvalE_filename         ) 
  
  x = lapply( EvalE_filepath, readRDS ) |> 
    dplyr::bind_rows()
  return( x )
}


Read_EvalSST = function( cfg ) {
  
  EvalE_filename <- paste0( 'SST', '_',
                            'DFO',                       '_', 
                            cfg$Eval_Year,               '_', 
                            '8days',                     '.',
                            'rds'              )
  
  EvalE_filepath <- paste0( cfg$initial_wd,              '/',
                            'SST', '/' ,
                            EvalE_filename         ) 
  
  x = lapply( EvalE_filepath, readRDS ) |> 
    dplyr::bind_rows()
  return( x )
}


Read_EvalBathy = function( cfg ) {
  
  EvalE_filename <- paste0( 'Bathy', '_',
                            'DFO',                       '_', 
                            cfg$Eval_Year,               '_', 
                            '8days',                     '.',
                            'rds'              )
  
  EvalE_filepath <- paste0( cfg$initial_wd,              '/',
                            'Bathy', '/' ,
                            EvalE_filename         ) 
  
  x = lapply( EvalE_filepath, readRDS ) |> 
    dplyr::bind_rows()
  return( x )
}



# Put the file into the Sightings object
Sightings = Read_Sightings( cfg )
# Put the file into the Env object
pp    <- Read_PP  ( cfg )
par   <- Read_PAR ( cfg )
chl   <- Read_CHL ( cfg )
sst   <- Read_SST ( cfg )
bathy <- Read_Bathy ( cfg )

#Load Evaluation dataset
Eval_Sightings <- Read_EvalSightings ( cfg )
Eval_pp        <- Read_EvalPP  ( cfg )
Eval_par       <- Read_EvalPAR ( cfg )
Eval_chl       <- Read_EvalCHL ( cfg )
Eval_sst       <- Read_EvalSST  ( cfg )
Eval_bathy     <- Read_EvalBathy ( cfg )



######################. -------------
###### --- Setting data for Biomod2
######################. -------------
# Creating dataframe for saving TSS and ROC Ensemble modeling scores for each model
EvalEM_TSS = data.frame()
EvalEM_ROC = data.frame()

  for (i_sst in 1:ncol(sst)) {
    for (i_pp in 1:ncol(pp))   {
      for (i_chl in 1:ncol(chl)) {
        for (i_par in 1:ncol(par)) {
      
      
      # Créer le label du model ainsi que le modèle. ex: sst-1-pp-1
      label <- paste("sst", i_sst, "pp", i_pp, "chl", i_chl, "par", i_par, sep = "_")
      
# species name
myRespName <- cfg$MyRespName

# presence/absence data
myResp <- Sightings$Sp_code

# XY coordinates
myXY <- Sightings[, 1:2]

# environmental variables
# env.var <- c( "sst",
#               "CHL",
#               "pp",
#               # "CDM_mean",
#               # "BBP_mean",
#               "PAR_mean",
#               "bathy" )

 # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl <- cbind(bathy$elevation, sst[, i_sst], pp[ ,i_pp], chl[ ,i_chl], par[, i_par])

#myExpl  <- Env_var[ , env.var ]

#Evaluation Dataset
myEvalResp  <- Eval_Sightings$Sp_code
myEvalExpl  <- cbind(Eval_bathy$elevation, Eval_sst[, i_sst], Eval_pp[ ,i_pp], Eval_chl[,i_chl], Eval_par[, i_par])
myEvalXY    <- Eval_Sightings[, 1:2]


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


## RUN 1.0) ALL Train. & 2020 Eval. ----

myBiomodModel <- BIOMOD_Modeling(
  bm.format         = myData.all,                                 #a BIOMOD.formated.data or BIOMOD.formated.data.PA object returned by the BIOMOD_FormatingData function
  bm.options        = myBiomodOption,
  modeling.id       = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', label), #a character corresponding to the name (ID) of the simulation set (a random number by default)
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
  weights           = NULL,                                       #a vector of numeric values corresponding to observation weights (one per observation, see Details)
  prevalence        = 0.7,                                        #a numeric between 0 and 1 corresponding to the species prevalence to build 'weighted response weights' (see Details)
  metric.eval       = cfg$Metrics_Eval,                            #"KAPPA", "TSS","ACCURACY", BIAS, POD, FAR, POFD, SR, CSI, ETS, HK, HSS, OR, ORSS
  var.import        = 5,                                          #an integer corresponding to the number of permutations to be done for each variable to estimate variable importance
  scale.models      = FALSE,                                      #a logical value defining whether all models predictions should be scaled with a binomial GLM or not
  nb.cpu            = 15,                                          #a integer value corresponding to the number of computing resources to be used to parallelize the single models computation
  seed.val          = 42,                                         #an integer value corresponding to the new seed value to be set
  do.progress       = TRUE                                        #a logical value defining whether the progress bar is to be rendered or not
  #modeling.id       = paste('RIWH',"FirstModeling",sep="")
)



## Evaluation scores ----

# choose which model output to analyze
model.postproc <- myBiomodModel

# create ggplot2 object
# extract evaluations metrics for the "calibration" dataset
g.1 <- bm_PlotEvalMean( bm.out      = myBiomodModel, #a BIOMOD.models.out or BIOMOD.ensemble.models.out object that can be obtained with the BIOMOD_Modeling or BIOMOD_EnsembleModeling functions 
                        metric.eval = c("ROC", "TSS"),  #a vector containing evaluation metric names to be used, must be among ROC, TSS, KAPPA, ACCURACY, BIAS, POD, FAR, POFD, SR, CSI, ETS, HK, HSS, OR, ORSS
                        dataset     = "calibration",    #a character corresponding to the dataset upon which evaluation metrics have been calculated and that is to be represented, must be among calibration, validation, evaluation
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
g.eval <- ggplot( data = g.data_RF,
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
  ylim( 0, 1 )                         

png( paste0("Eval_mean", "_", cfg$Single_Models , label, ".png"), width = 1500, height = 1000 )
g.eval
dev.off()

# save evaluation scores
capture.output( myBiomodModel@models.evaluation@val,
                file = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', label, '.', 'Eval_mean', '_', cfg$Single_Models, '.csv'))




## Importance of variables ----

# extract values for each variable's importance
g.data <- bm_PlotVarImpBoxplot( bm.out   = myBiomodModel,
                                group.by = c('expl.var', 'algo', 'algo'), #a 3-length vector containing the way kept models will be represented, must be among full.name, PA, run, algo, expl.var (if bm.out is a BIOMOD.models.out object), or full.name, merged.by.PA, merged.by.run, merged.by.algo, expl.var (if bm.out is a BIOMOD.ensemble.models.out object)
                                do.plot  = FALSE                          #a logical value defining whether the plot is to be rendered or not
)

# ggplot2 graphic
g.VarImp <- ggplot( data = g.data$tab,
                    aes( x = reorder( expl.var, var.imp ),
                         y = var.imp * 100 ) ) +
  geom_boxplot( colour = "firebrick" )         +
  labs( x        = "Explanatory variables",
        y        = "Importance (%)",
        subtitle = "Variable importance" )

png( paste0("Var_Imp", "_", cfg$Single_Models, label,".png"), width = 1500, height = 1000 )
g.VarImp
dev.off()

# save variable importance scores on hard drive
capture.output( g.data_Var$tab,
                file = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', cfg$run, '.', 'Var_Imp', '_', cfg$Single_Models, '.csv'))





## --- Represent response curves 

png( "Resp_curv.png", label, width = 1000, height = 800 )
bm_PlotResponseCurves( bm.out        = myBiomodModel, 
                       models.chosen = get_built_models( model.postproc,
                                                         full.name = "RIWH_allData_allRun_GAM" ), #a vector containing model names to be kept, must be either all or a sub-selection of model names that can be obtained with the get_built_models function
                       fixed.var     = 'median',                                                  #a character corresponding to the statistic to be used to fix as constant the remaining variables other than the one used to predict response, must be either mean, median, min, max
                       do.bivariate  = FALSE,                                                     #a logical value defining whether the response curves are to be represented in 3 dimensions (meaning 2 explanatory variables at a time) or not (meaning only 1)
                       do.plot       = TRUE )
dev.off()






## XXXXXXXXXXXXXXXXXXXX ----
## --------------------   ENSEMBLE MODELS    --------------------


myBiomodEM <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModel.all_v1.0,  #a BIOMOD.models.out object returned by the BIOMOD_Modeling function
  models.chosen = "all",  #a vector containing model names to be kept, must be either all or a sub-selection of model names that can be obtained with the get_built_models function
  em.by = "all", #a character corresponding to the way kept models will be combined to build the ensemble models, must be among PA_dataset+repet, PA_dataset+algo, PA_dataset, algo, all
  metric.select = c("ROC", "TSS"), # with metric.select.thresh to exclude single models "KAPPA", "TSS", "ROC"....
  metric.select.thresh = c(0.7, 0.3), #A vector of numeric values corresponding to the minimum scores (one for each metric.select) below which single models will be excluded from the ensemble model building
  metric.select.table = NULL, #If metric.select = 'user.defined', a data.frame containing evaluation scores calculated for each single models and that will be compared to metric.select.thresh values 
  #to exclude some of them from the ensemble model building, with metric.select rownames, and models.chosen colnames
  metric.eval = c("ROC", "TSS"), #"KAPPA", "TSS", 
  var.import = 5, #An integer corresponding to the number of permutations to be done for each variable to estimate variable importance
  #EMmmean = TRUE, #A logical value defining whether to compute the mean probabilities across predictions or not
  #EMmedian = F, #A logical value defining whether to compute the median probabilities across predictions or not
  #EMcv = T, #A logical value defining whether to compute the coeff of variation across predictions or not
  #EMci = T, #A logical value defining whether to compute the confidence interval around the prob.mean ensemble model or not
  #EMci.alpha = 0.05, #A numeric value corresponding to the significance level to estimate confidence interval
  #EMca = TRUE, #A logical value defining whether to compute the committee averaging across predictions or not
  #EMwmean = TRUE, #A logical value defining whether to compute the weighted sum of probabilities across predictions or not
  EMwmean.decay = "proportional", #A value defining the relative importance of the weights (if prob.mean.weight = TRUE). 
  #A high value will strongly discriminate good models from the bad ones (see Details), 
  #while proportional will attribute weights proportionally to the models evaluation scores
  em.algo = cfg$Algo,
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


png( paste0("EMEval_mean", "_", label , ".png"), width = 1500, height = 1000 )
g.eval
dev.off()


# save evaluation scores on hard drive
capture.output( myBiomodEM@models.evaluation@val,
                file = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', label, '.', 'Eval_mean', '_', '.csv'))



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

png( paste0("EMVar_Imp",  label , "_" ,".png"), width = 1500, height = 1000 )
g.VarImp
dev.off()

# save variable importance scores on hard drive
capture.output( g.data_Var$tab,
                file = paste0(cfg$Tempstemp, '.', cfg$Domain, '.', cfg$Year[1], 'to', cfg$Year[n], '.', cfg$Eval_Year, '.', label, '.', 'Var_Imp', '_', '.csv'))



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


      }
    }
   }
  }


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
               file=file.path(myRespName, 
                              paste(myRespName,"_Ensmodels_prediction_evaluation.csv", sep="")))


png(paste0(directory, "/EMEval_mean.png"), width = 1024, height = 768)
bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'algo')
dev.off()

png(paste0(directory, "/EMImp_var.png"))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'algo'))
dev.off()

png(paste0(directory, "/EMResp_curv.png"), width = 1024, height = 768)
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM, full.name = c('RIWH_EMwmeanByROC_mergedData_mergedRun_mergedAlgo','RIWH_EMwmeanByTSS_mergedData_mergedRun_mergedAlgo')),
                      fixed.var = 'median', 
                      do.bivariate = F)         
dev.off()



library(readr)
library(tidyr)
library(terra)
library(stars)
library(ncmeta)
library(biomod2)
library(plotrix)
library(ggtext)
library(purrr)




# Modeling
## Preparing data

#All data must have same number of lines and environmental data must have the same number of colonnes as well.



sst_test = SST_2015to2020[, 1:7]
pp_test = PP_2015to2020[, 1:7]
cdm_test =  CDM_2015to2020[, 1:7]

Occ_test = Occ_2015to2020    #must be a vector  
XY_test = XY_2015to2020


## Single modeling function  - BIOMOD_Modeling



# Creating dataframe for saving TSS and ROC Ensemble modeling scores for each model
EvalEM_TSS = data.frame()
EvalEM_ROC = data.frame()

for (i_sst in 1:ncol(sst_test)) {
  for (i_pp in 1:ncol(pp_test)) {
    for (i_cdm in 1:ncol(cdm_test)) {
      
      
      # Créer le label du model ainsi que le modèle. ex: sst-1-pp-1
      label <- paste("sst", i_sst, "pp", i_pp, "cdm", i_cdm, sep = "_")
      
      
      
      
      # Select the name of the studied species
      myRespName <- 'NARW'
      
      # Get corresponding presence/absence data
      myResp <- Occ_test
      
      # Get corresponding XY coordinates
      myRespXY <- XY_test
      
      # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
      myExpl <- cbind(Bathy_2015to2020$Bathy, sst_test[, i_sst], pp_test[ ,i_pp], cdm_test[ ,i_cdm])
      
      #myExpl <- as.data.frame(cbind(sst_test, pp_test, cdm_test))
      
      
      # Format Data with true absences
      myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = myExpl,
                                           resp.xy = myRespXY,
                                           resp.name = myRespName)
      
      
      cat(myRespName,file="")  
      
      
      ### Options definition
      myBiomodOption <- BIOMOD_ModelingOptions()
      
      
      # function from 'purr' to skip error and go to the next iteration
      possibly_BIOMOD_Model = possibly(BIOMOD_Modeling, otherwise = NA)
      
      
      
      ###################################################################
      ## ----------------     -------------- Single Modelling 
      ###################################################################
      myBiomodModelOut <- possibly_BIOMOD_Model(
        bm.format = myBiomodData, #a BIOMOD.formated.data or BIOMOD.formated.data.PA object returned by the BIOMOD_FormatingData function
        bm.options = myBiomodOption,
        modeling.id = label, #a character corresponding to the name (ID) of the simulation set (a random number by default)
        models = c("GBM", "GAM"), #a vector containing model names to be computed, must be among GLM, GBM, GAM, CTA, ANN, SRE, FDA, MARS, RF, MAXENT.Phillips, MAXENT.Phillips.2
        models.pa = NULL,   #(optional, default NULL) A list containing for each model a vector defining which pseudo-absence datasets are to be used, must be among colnames(bm.format@PA.table)
        CV.strategy = "random",  #a character corresponding to the cross-validation selection strategy, must be among random, kfold, block, strat, env or user.defined
        CV.nb.rep = 2,  #an integer corresponding to the number of repetitions to be done for calibration/validation splitting
        CV.perc = 0.75, #a numeric between 0 and 1 corresponding to the percentage of data used to calibrate the models (calibration/validation splitting)
        #### CV.k = 3,   #(optional, default 0) If strategy = 'kfold' or strategy = 'strat' or strategy = 'env', an integer corresponding to the number of partitions
        #### CV.balance = 'presences',  #(optional, default 'presences') If strategy = 'strat' or strategy = 'env', a character corresponding to how data will be balanced between partitions, must be either presences or absences
        ####CV.env.var = 'all',  #(optional) If strategy = 'env', a character corresponding to the environmental variables used to build the partition. k partitions will be built for each environmental variables. By default the function uses all environmental variables available.
        #### CV.strat = 'both',  #(optional, default 'both') If strategy = 'env', a character corresponding to how data will partitioned along gradient, must be among x, y, both
        # CV.user.table = NULL, #. A matrix or data.frame defining for each repetition (in columns) which observation lines should be used for models calibration (TRUE) and validation (FALSE) (see BIOMOD_CrossValidation)
        #(if specified, nb.rep, data.split.perc and do.full.models will be ignored)
        CV.do.full.models = TRUE, #A logical value defining whether models calibrated and evaluated over the whole dataset should be computed or not
        weights = NULL, #A vector of numeric values corresponding to observation weights (one per observation, see Details)
        prevalence = 0.6, #A numeric between 0 and 1 corresponding to the species prevalence to build 'weighted response weights' (see Details)
        metric.eval = c("ROC", "TSS"),#"KAPPA", "TSS","ACCURACY", BIAS, POD, FAR, POFD, SR, CSI, ETS, HK, HSS, OR, ORSS
        var.import = 2, # An integer corresponding to the number of permutations to be done for each variable to estimate variable importance
        #save.output = TRUE, #A logical value defining whether all outputs should be saved on hard drive or not (! strongly recommended !)
        scale.models = FALSE, #A logical value defining whether all models predictions should be scaled with a binomial GLM or not
        nb.cpu = 15, #An integer value corresponding to the number of computing resources to be used to parallelize the single models computation
        seed.val = NULL,#An integer value corresponding to the new seed value to be set
        do.progress = TRUE, # A logical value defining whether the progress bar is to be rendered or not
        #modeling.id = paste('RIWH',"FirstModeling",sep="")
      )
      
      
      
      label <- paste("sst", i_sst, "pp", i_pp, "cdm", i_cdm, sep = "_")
      
      
      # Directory repository
      directory <- paste0("/Ext_16T_andromede/Extern_workdir/thieryf_WD/NARW/", label)
      dir.create(directory)
      
      # save evaluation
      capture.output(ev, file = paste0(directory, "/formal_eval_model.csv"))
      
      #save plots
      png(paste0(directory, "/Eval_mean.png"))
      bm_PlotEvalMean(bm.out = myBiomodModelOut, group.by = 'algo')
      dev.off()
      
      png(paste0(directory, "/Imp_var.png"))
      bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
      dev.off()
      
      png(paste0(directory, "/Resp_curv.png"))
      bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                            models.chosen = 'all',
                            fixed.var = 'median', 
                            do.bivariate = F)
      dev.off()

      
      
      ## Ensemble modeling function - BIOMOD_EnsembleModeling

      # function to skip error and go to the next iteration
      possibly_BIOMOD_EMModel = possibly(BIOMOD_EnsembleModeling, otherwise = NA)
      
      myBiomodEM <- possibly_BIOMOD_EMModel(
        bm.mod = myBiomodModelOut,  #a BIOMOD.models.out object returned by the BIOMOD_Modeling function
        models.chosen = "all",  #a vector containing model names to be kept, must be either all or a sub-selection of model names that can be obtained with the get_built_models function
        em.by = "all", #a character corresponding to the way kept models will be combined to build the ensemble models, must be among PA_dataset+repet, PA_dataset+algo, PA_dataset, algo, all
        metric.select = c("ROC", "TSS"), # with metric.select.thresh to exclude single models "KAPPA", "TSS", "ROC"....
        #metric.select.thresh = c(0,75, 0.5), #A vector of numeric values corresponding to the minimum scores (one for each metric.select) below which single models will be excluded from the ensemble model building
        #metric.select.table = NULL, #If metric.select = 'user.defined', a data.frame containing evaluation scores calculated for each single models and that will be compared to metric.select.thresh values 
        #to exclude some of them from the ensemble model building, with metric.select rownames, and models.chosen colnames
        metric.eval = c("ROC", "TSS"), #"KAPPA", "TSS", 
        var.import = 5, #An integer corresponding to the number of permutations to be done for each variable to estimate variable importance
        #EMmmean = TRUE, #A logical value defining whether to compute the mean probabilities across predictions or not
        #EMmedian = F, #A logical value defining whether to compute the median probabilities across predictions or not
        #EMcv = T, #A logical value defining whether to compute the coeff of variation across predictions or not
        #EMci = T, #A logical value defining whether to compute the confidence interval around the prob.mean ensemble model or not
        #EMci.alpha = 0.05, #A numeric value corresponding to the significance level to estimate confidence interval
        #EMca = TRUE, #A logical value defining whether to compute the committee averaging across predictions or not
        #EMwmean = TRUE, #A logical value defining whether to compute the weighted sum of probabilities across predictions or not
        EMwmean.decay = "proportional", #A value defining the relative importance of the weights (if prob.mean.weight = TRUE). 
        #A high value will strongly discriminate good models from the bad ones (see Details), 
        #while proportional will attribute weights proportionally to the models evaluation scores
        em.algo = c('EMmean', 'EMwmean', 'EMcv', 'EMca'),#, 'EMci'),
        nb.cpu = 20, #An integer value corresponding to the number of computing resources to be used to parallelize the single models computation
        seed.val = NULL, #An integer value corresponding to the new seed value to be set
        do.progress = TRUE #A logical value defining whether the progress bar is to be rendered or not
      )
      
      
      # Ensemble modeling mean committee averaging TSS and ROC scores 
      Eval_EMca_TSS = mean(myBiomodEM@models.evaluation@val[myBiomodEM@models.evaluation@val$algo== 'EMca' & myBiomodEM@models.evaluation@val$metric.eval == 'TSS', 'calibration'])
      Eval_EMca_ROC = mean(myBiomodEM@models.evaluation@val[myBiomodEM@models.evaluation@val$algo== 'EMca' & myBiomodEM@models.evaluation@val$metric.eval == 'ROC', 'calibration'])
      
      # Ensemble modeling mean weigthed mean TSS and ROC scores 
      Eval_EMwmean_TSS = mean(myBiomodEM@models.evaluation@val[myBiomodEM@models.evaluation@val$algo== 'EMwmean' & myBiomodEM@models.evaluation@val$metric.eval == 'TSS', 'calibration'])
      Eval_EMwmean_ROC = mean(myBiomodEM@models.evaluation@val[myBiomodEM@models.evaluation@val$algo== 'EMwmean' & myBiomodEM@models.evaluation@val$metric.eval == 'ROC', 'calibration'])  
      
      # Fill the dataframes
      EvalEM_TSS = rbind(EvalEM_TSS, (data.frame(Eval_EMca_TSS, Eval_EMwmean_TSS, label)))
      EvalEM_ROC = rbind(EvalEM_ROC, (data.frame(Eval_EMca_ROC, Eval_EMwmean_ROC, label)))   
      
      
      # Save outputs and plots
      capture.output(EvalEM_ROC, file = paste0(directory, "/EvalEM_ROC.csv"))
      capture.output(EvalEM_TSS, file = paste0(directory, "/EvalEM_TSS.csv"))
      
      capture.output(ev_EM, file = paste0(directory, "/formal_evalEM_model.csv"))
      
      
      png(paste0(directory, "/Eval_mean_EM.png"))
      bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'algo')
      dev.off()
      
      png(paste0(directory, "/Imp_var_EM.png"))
      bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'algo'))
      dev.off()
      
      png(paste0(directory, "/Resp_curv_EM.png"))
      bm_PlotResponseCurves(bm.out = myBiomodEM, 
                            models.chosen = 'all',
                            fixed.var = 'mean',
                            do.bivariate = F)
      dev.off()
      
      
    }
  }
}




