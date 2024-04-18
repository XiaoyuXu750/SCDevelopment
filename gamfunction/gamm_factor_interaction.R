library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(lme4)
library(gamm4)
library(pbkrtest)

pboot <- function(modelobj, int_var=NA){
  numsims <- 1000
  set.seed(925)
  df <- modelobj$gam$model
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula
  theseVars <- attr(terms(f1),"term.labels")
  if (!is.na(int_var)){theseVars <-  c(theseVars, int_var)}
  if (sum(str_detect(theseVars, "by ="))==1){
    int_var <- str_split_i(theseVars[1], "by = ", 2)
    int_var <- str_split_i(int_var, ", ", 1)
    f2 <- reformulate(c(theseVars[2:(length(theseVars))], int_var),response = thisResp)
  }else{
    f2 <- reformulate(theseVars[2:(length(theseVars))],response = thisResp)
  }
    
  g1 <- gam(f1,data = df)
  g2 <- gam(f2,data = df)
  
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  subID<-df$subID
  y <- df[,thisResp]
  
  m1 <- lmer(y ~ -1 + mat1 + (1|subID))
  m2 <- lmer(y ~ -1 + mat2 + (1|subID))
  refdist <- PBrefdist(m1, m2, nsim=numsims)
  pb <- PBmodcomp(m1, m2, ref = refdist)
  int_pval <- pb$test["PBtest","p.value"]
  
  return(int_pval)
}

#### PREDICT GAM SMOOTH FITTED VALUES FOR A SPECIFIED VALUE OF AN INTERACTING COVARIATE ####
## discrete interaction covariate
##Function to predict fitted values of a region for a given value of a covariate, using a varying coefficients smooth-by-linear covariate interaction
gamm.smooth.predict.interaction <- function(region, dataname, smooth_var, int_var, covariates, knots, set_fx = FALSE, increments, stats_only=TRUE){
  gam.data <- get(dataname)
  tmp<-gam.data[,region]
  outlierindx<-which(tmp<mean(tmp)-3*sd(tmp) | tmp>mean(tmp)+3*sd(tmp))
  if (length(outlierindx)>0){
    gam.data<-gam.data[-outlierindx, ]
  }
  parcel <- region
  
  #Fit the gam
  modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, by=%5$s, k=%3$s, fx=%4$s)+ s(%2$s, k=%3$s, fx=%4$s) + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  gamm.model <- gamm4(modelformula, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.results <- summary(gamm.model$gam)
  
  modelformula.null <- as.formula(sprintf("%1$s ~ %5$s+s(%2$s, k=%3$s, fx=%4$s)+ %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  gamm.model.null <- gamm4(modelformula.null, random=~(1|subID), REML=TRUE, data = gam.data)
  gamm.null.results <- summary(gamm.model.null$gam)
  
  T.disease <- gamm.null.results$p.table[2,3]
  P.disease <- gamm.null.results$p.table[2,4]
  bootstrap.P.disease <- pboot(gamm.model.null)
  # stats
  bootstrap_pvalue<-pboot(gamm.model, int_var)
  gam.int.pvalue <- gamm.results$s.table[grep(x=rownames(gamm.results$s.table),pattern = int_var),"p-value"][1]
  # interaction effect size
  sse.model <- sum((gamm.model$gam$y - gamm.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gamm.model.null$gam$y - gamm.model.null$gam$fitted.values)^2)
  IntpartialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  ### get predicted data
  ##############################
  modelobj <- gamm.model$gam
  #Extract gam input data
  df <- modelobj$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame
  np <- increments #number of predictions to make; predict at np increments of smooth_var
  thisPred <- data.frame(init = rep(0,np)) #initiate a prediction df 
  
  theseVars <- attr(modelobj$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(modelobj$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(modelobj$terms[[2]]) #the measure to predict
  for (v in c(1:length(theseVars))) { #fill the prediction df with data for predictions. These data will be used to predict the output measure (y) at np increments of the smooth_var, holding other model terms constant
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      tab.tmp<-table(df[,thisVar])
      levelact<-which.max(tab.tmp)
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
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
  predicted.smooth.HC <- predicted.smooth.HC %>% dplyr::select(all_of(smooth_var), fitted, se, lower, upper, fitted.centered)
  predicted.smooth.Disease <- fitted_values(object = modelobj, data = pred.Disease)
  predicted.smooth.Disease$fitted.centered <- scale(predicted.smooth.Disease$fitted, center=T, scale = F) #subtract the intercept from fitted values
  predicted.smooth.Disease <- predicted.smooth.Disease %>% dplyr::select(all_of(smooth_var), fitted, se, lower, upper, fitted.centered)
  predicted.smooth <- rbind(predicted.smooth.HC, predicted.smooth.Disease)
  
  changed.range.HC <- predicted.smooth.HC$fitted[which.max(predicted.smooth.HC$age)]-predicted.smooth.HC$fitted[which.min(predicted.smooth.HC$age)]
  changed.ratio.HC <- predicted.smooth.HC$fitted[which.max(predicted.smooth.HC$age)] / predicted.smooth.HC$fitted[which.min(predicted.smooth.HC$age)]
  SCweight.agemax.HC <- predicted.smooth.HC$fitted[which.max(predicted.smooth.HC$age)]
  SCweight.agemin.HC <- predicted.smooth.HC$fitted[which.min(predicted.smooth.HC$age)]
  changed.range.Disease <- predicted.smooth.Disease$fitted[which.max(predicted.smooth.Disease$age)]-predicted.smooth.Disease$fitted[which.min(predicted.smooth.Disease$age)]
  changed.ratio.Disease <- predicted.smooth.Disease$fitted[which.max(predicted.smooth.Disease$age)] / predicted.smooth.Disease$fitted[which.min(predicted.smooth.Disease$age)]
  SCweight.agemax.Disease <- predicted.smooth.Disease$fitted[which.max(predicted.smooth.Disease$age)]
  SCweight.agemin.Disease <- predicted.smooth.Disease$fitted[which.min(predicted.smooth.Disease$age)]
  ############################################
  
  stats.reults<-cbind(parcel, int_var, bootstrap_pvalue, gam.int.pvalue, IntpartialRsq, T.disease, P.disease, bootstrap.P.disease,
                      changed.range.HC, changed.range.Disease, changed.ratio.HC, changed.ratio.Disease,
                      SCweight.agemax.HC, SCweight.agemax.Disease, SCweight.agemin.HC, SCweight.agemin.Disease)
  
  full.results<-list(stats.reults, predicted.smooth)
  if(stats_only == TRUE)
    return(stats.reults)
  if(stats_only == FALSE)
    return(full.results)

}






