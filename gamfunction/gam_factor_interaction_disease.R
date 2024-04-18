library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(ecostats)

#### Interaction analysis for linear effects (cognition) by  a categorical variable####
##Function to predict fitted values of a region for a each level of a categorical variable, using a varying coefficients linear-by-category interaction
gam.factor.predict.covariateinteraction <- function(region, dataname, smooth_var, int_var,cog_var, covariates, knots, set_fx = FALSE, increments, stats_only=FALSE){
  #Fit the gam
  gam.data <- get(dataname)
  gam.data[,region] <- as.numeric(gam.data[,region])
  tmp<-gam.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp), tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  tmp<-gam.data[,cog_var]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp), tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  
  parcel <- region
  
  #Fit the gam
  modelformula <- as.formula(sprintf("%1$s ~ %5$s*%6$s +s(%2$s, k=%3$s, fx=%4$s)+ %7$s",cog_var, smooth_var, knots, set_fx, int_var, region, covariates))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  modelformula.null <- as.formula(sprintf("%1$s ~ s(%2$s, k=%3$s, fx=%4$s)+%5$s+ %6$s+ %7$s", cog_var, smooth_var, knots, set_fx, int_var, region, covariates))
  gam.model.null <- gam(modelformula.null, method = "REML", data = gam.data)
  gam.null.results <- summary(gam.model.null)
  
  ##Full versus reduced model anova p-value
  #anova.int.pvalue <- anova(gam.model.null, gam.model,test='Chisq')$`Pr(>Chi)`[2]
  anova.int.pvalue <- anovaPB(gam.model.null, gam.model, n.sim = 1000,test='Chisq')$`Pr(>Chi)`[2]
  gam.int.pvalue <- gam.results$p.table[grep(x=rownames(gam.results$p.table),pattern = ":"),"Pr(>|t|)"][1]
  # interaction effect size
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.model.null$y - gam.model.null$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  ### get predicted data
  ##############################
  modelobj <- gam.model
  #Extract gam input data
  df <- modelobj$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(modelobj$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(modelobj$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(modelobj$terms[[2]]) #the measure to predict
  for (v in c(1:(length(theseVars)-1))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == region) { 
      thisPred[,region] = seq(min(df[,region],na.rm = T),max(df[,region],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "nmatrix.1" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]}, #make predictions based on the modal number
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[levelact]]} #make predictions based on the modal number
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init)
  pred.HC <- pred.Disease <- pred
  pred.HC[,int_var] <- as.factor(0)
  pred.Disease[,int_var] <- as.factor(1)
  
  #Generate fitted (predicted) values based on the gam model and predication data frame
  predicted.smooth.HC <- fitted_values(object = modelobj, data = pred.HC)
  predicted.smooth.HC$fitted.centered <- scale(predicted.smooth.HC$fitted, center=T, scale = F) #subtract the intercept from fitted values
  predicted.smooth.Disease <- fitted_values(object = modelobj, data = pred.Disease)
  predicted.smooth.Disease$fitted.centered <- scale(predicted.smooth.Disease$fitted, center=T, scale = F) #subtract the intercept from fitted values
  predicted.smooth <- rbind(predicted.smooth.HC, predicted.smooth.Disease)
  
  Cog.changed.range.HC <- predicted.smooth.HC$fitted[which.max(predicted.smooth.HC[[region]])]-predicted.smooth.HC$fitted[which.min(predicted.smooth.HC[[region]])]
  Cog.changed.ratio.HC <- predicted.smooth.HC$fitted[which.max(predicted.smooth.HC[[region]])] / predicted.smooth.HC$fitted[which.min(predicted.smooth.HC[[region]])]
  Cog.SCmax.HC <- predicted.smooth.HC$fitted[which.max(predicted.smooth.HC[[region]])]
  Cog.SCmin.HC <- predicted.smooth.HC$fitted[which.min(predicted.smooth.HC[[region]])]
  Cog.changed.range.Disease <- predicted.smooth.Disease$fitted[which.max(predicted.smooth.Disease[[region]])]-predicted.smooth.Disease$fitted[which.min(predicted.smooth.Disease[[region]])]
  Cog.changed.ratio.Disease <- predicted.smooth.Disease$fitted[which.max(predicted.smooth.Disease[[region]])] / predicted.smooth.Disease$fitted[which.min(predicted.smooth.Disease[[region]])]
  Cog.SCmax.Disease <- predicted.smooth.Disease$fitted[which.max(predicted.smooth.Disease[[region]])]
  Cog.SCmin.Disease <- predicted.smooth.Disease$fitted[which.min(predicted.smooth.Disease[[region]])]
  ############################################
  
  stats.reults<-cbind(parcel, int_var, cog_var, anova.int.pvalue,gam.int.pvalue, IntpartialRsq, Cog.changed.range.HC,
                      Cog.changed.ratio.HC, Cog.SCmax.HC, Cog.SCmin.HC, 
                      Cog.changed.range.Disease, Cog.changed.ratio.Disease, Cog.SCmax.Disease, Cog.SCmin.Disease)
  full.results<-list(stats.reults, predicted.smooth)
  if(stats_only == TRUE)
    return(stats.reults)
  if(stats_only == FALSE)
    return(full.results)
}