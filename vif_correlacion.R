## -*- coding: utf-8 -*-
#"""
#Created on Abr 2025
#
#@author: Carlos Munoz, Nerieth Leuro
#@edit: 
#"""

##Apply VIF analysis to the predictor variables prior to the model run. In addition, two correlation matrices are obtained  using layerCor y cor
##The corresponding VIF calculation code  was obtained from biomodelos-sdm scripts
#****************************************************************************************************************************
##Load librarys
#****************************************************************************************************************************
library(dplyr)
library(maps)
library(ggplot2)
library(sf)
library(terra)
library(raster)
library(ENMeval)
library(fmsb)
library(corrplot)

#****************************************************************************************************************************
##Define the working directory where the variables are located
#****************************************************************************************************************************
setwd("D:/biomodelos-sdm-master/modelling/Xenatros/Data/other/current")

#****************************************************************************************************************************
##Load data set
#****************************************************************************************************************************
vifsample <- read.csv("D:/biomodelos-sdm-master/modelling/Tapirus_Tayassus/occurrence_Tapirus_pinchaque.csv")

##Select only the coordinate columns
vifsample <-select(vifsample,decimalLatitude,decimalLongitude)

##Convert data into geographic points
vifsample <- vifsample |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs("EPSG:4326"))

##Generate the list of variables
lista.variables <- list.files(pattern='*.tif', full.names=TRUE)

##Convert the list to raster format
envars <- rast(lista.variables) 
print(names(envars))

## Assign unique names to variables
names(envars) <- make.unique(names(envars))
print(names(envars))

#****************************************************************************************************************************
#Function for applying VIF and correlations
#****************************************************************************************************************************

##Specify the VIF threshold value
vifdetails <- 10

##Execute function
vif_apply <- function(vifsample, envars, vifdetails){
  
  d <- terra::extract(envars, vifsample) %>% 
    na.omit()
  names(d) <- make.names(names(d))

  ## Extract a random sample of n non-zero pixels
  sample_vals <- spatSample(envars, size = 5000, method = "random", na.rm = TRUE)  
    
  ###Correlation with layerCor
  variables.cor <- layerCor(envars, na.rm=TRUE,fun = "pearson",maxcell=5000)

  ##Correlation with cor
  cor_matrix <- cor(sample_vals, method = "spearman", use = "pairwise.complete.obs")
 
  
  vifd <- vif_func(d, thresh = vifdetails)
  return(list(variables_filtradas = vifd,
              matriz_correlacion_alet = cor_matrix,
              matriz_correlacion_tot = variables.cor$pearson)) 
}

#' Calculates VIF values for a set of explanatory variables using linear regression models and iteratively
#' removes variables with high VIF values until all VIF values are below a specified threshold. 
#' 
#' @description This function computes the VIF values for a set of explanatory variables using a data frame 
#' (in_frame). The VIF is calculated by running a linear regression on each variable while using all other 
#' variables as predictors. The function returns the names of the variables that have a VIF value below a 
#' specified threshold (thresh).
#' 
#' @param in_frame data frame containing the data used to calculate the VIF.
#' @param thresh numeric value specifying the VIF threshold.
#' @param trace logical value indicating whether to print the output of each iteration.
#' @param ... additional arguments to be passed to the lm() function.

vif_func<-function(in_frame,thresh=10,trace=F,...){
  
  require(fmsb) 
  require(dismo)
  require(rgdal)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

#****************************************************************************************************************************
##Perform function call
#****************************************************************************************************************************
resultado <- vif_apply (vifsample, envars, vifdetails)

#****************************************************************************************************************************
##Show result
#****************************************************************************************************************************
##Print Variables selected by the VIF
print(paste('Variable seleccionada no colineal:',resultado$variables_filtradas))

##Export csv with Variables selected by the VIF
write.csv(resultado$variables_filtradas,"resultado$variables_filtradas.csv")

# Show layerCor correlation matrix
View(resultado$matriz_correlacion_tot)
write.csv(resultado$matriz_correlacion_tot,"matriz_correlacion_variables.cor.csv")

# Show cor correlation matrix
View(resultado$matriz_correlacion_alet )
write.csv(resultado$matriz_correlacion_alet,"matriz_correlacion_cor_matrix.csv")

#Adjust the names displayed on the graphs
colnames(resultado$matriz_correlacion_tot)<-c("mean","aetf","30average10","30average30","30average40","alpha","aspectcosine","aspectsine","Bare","class1","class10","class11","class12","class2"
                                              ,"class3","class4","class5","class6","class7","class8","class9","eviqa1","eviqa2","eviqa3","fpar8qa1","fpar8qa2","fpar8qa3","gppqa1"
                                              ,"gppqa2","gppqa3","lai8qa1","lai8qa2","lai8qa3","ndviqa1","ndviqa2","ndviqa3","dx","dxx","dy","dyy","eastness","geomflat"
                                              ,"geomfootslope","geomhollow","geompeak","geompit","geomridge","geomshoulder","geomslope","geomspur","geomvalley","northness","pcurv","roughness"
                                              ,"Short_vegetation","slope","tcurv","tpi","Tree","tri","vrm","bio1","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18",
                                              "bio2","bio3","bio4","bio5","bio6","bio7","elev","mean.1","mean.2","mean.3","mean.4","mean.5","mean.6","mean.7" )

rownames(resultado$matriz_correlacion_tot)<-c("mean","aetf","30average10","30average30","30average40","alpha","aspectcosine","aspectsine","Bare","class1","class10","class11","class12","class2"
                                              ,"class3","class4","class5","class6","class7","class8","class9","eviqa1","eviqa2","eviqa3","fpar8qa1","fpar8qa2","fpar8qa3","gppqa1"
                                              ,"gppqa2","gppqa3","lai8qa1","lai8qa2","lai8qa3","ndviqa1","ndviqa2","ndviqa3","dx","dxx","dy","dyy","eastness","geomflat"
                                              ,"geomfootslope","geomhollow","geompeak","geompit","geomridge","geomshoulder","geomslope","geomspur","geomvalley","northness","pcurv","roughness"
                                              ,"Short_vegetation","slope","tcurv","tpi","Tree","tri","vrm","bio1","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18",
                                              "bio2","bio3","bio4","bio5","bio6","bio7","elev","mean.1","mean.2","mean.3","mean.4","mean.5","mean.6","mean.7" )

colnames(resultado$matriz_correlacion_alet)<-c("mean","aetf","30average10","30average30","30average40","alpha","aspectcosine","aspectsine","Bare","class1","class10","class11","class12","class2"
                                               ,"class3","class4","class5","class6","class7","class8","class9","eviqa1","eviqa2","eviqa3","fpar8qa1","fpar8qa2","fpar8qa3","gppqa1"
                                               ,"gppqa2","gppqa3","lai8qa1","lai8qa2","lai8qa3","ndviqa1","ndviqa2","ndviqa3","dx","dxx","dy","dyy","eastness","geomflat"
                                               ,"geomfootslope","geomhollow","geompeak","geompit","geomridge","geomshoulder","geomslope","geomspur","geomvalley","northness","pcurv","roughness"
                                               ,"Short_vegetation","slope","tcurv","tpi","Tree","tri","vrm","bio1","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18",
                                               "bio2","bio3","bio4","bio5","bio6","bio7","elev","mean.1","mean.2","mean.3","mean.4","mean.5","mean.6","mean.7" )

rownames(resultado$matriz_correlacion_alet)<-c("mean","aetf","30average10","30average30","30average40","alpha","aspectcosine","aspectsine","Bare","class1","class10","class11","class12","class2"
                                               ,"class3","class4","class5","class6","class7","class8","class9","eviqa1","eviqa2","eviqa3","fpar8qa1","fpar8qa2","fpar8qa3","gppqa1"
                                               ,"gppqa2","gppqa3","lai8qa1","lai8qa2","lai8qa3","ndviqa1","ndviqa2","ndviqa3","dx","dxx","dy","dyy","eastness","geomflat"
                                               ,"geomfootslope","geomhollow","geompeak","geompit","geomridge","geomshoulder","geomslope","geomspur","geomvalley","northness","pcurv","roughness"
                                               ,"Short_vegetation","slope","tcurv","tpi","Tree","tri","vrm","bio1","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18",
                                               "bio2","bio3","bio4","bio5","bio6","bio7","elev","mean.1","mean.2","mean.3","mean.4","mean.5","mean.6","mean.7" )


## Assign colors for ranges
colores <- c("lightblue", "darkred","lightblue")

## Define cuts for binary graph manually
breaks <- c(-1, -0.7, 0.7, 1)

##Plot correlations
##layerCor correlation matrix
corrplot(resultado$matriz_correlacion_tot, method = "color", order = "original",tl.cex = 0.7)
##cor correlation matrix
corrplot(resultado$matriz_correlacion_alet, method = "color", order = "original", tl.cex = 0.7)

##Plot binary correlations
##layerCor correlation matrix
corrplot(resultado$matriz_correlacion_tot, method = "color", order = "original",col = colores, breaks = breaks,tl.cex = 0.7)
##cor correlation matrix
corrplot(resultado$matriz_correlacion_alet, method = "color", order = "original",col = colores, breaks = breaks, tl.cex = 0.7)
