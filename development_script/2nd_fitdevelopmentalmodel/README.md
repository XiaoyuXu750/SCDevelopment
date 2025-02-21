# 2nd_fitdevelopmentalmodel
This folder contains codes for developmental analyses that estimate age effects, compute derivatives, and examine correlations between developmental patterns and the S-A connectional axis rank. It includes scripts for generating `Fig. 2, 3b–d, 5a–d,g,h,j,k`. Separate scripts are provided for the ABCD and HCP-D datasets for each step. For sensitivity analyses, users can specify the resolutions of large-scale SC matrices for the first two steps and specify the consistency threshold for all the steps.

## S1st_fitgammodels_SA_ds_sumSCinvnode_*.R
This script fitted developmental models for each edge. The resolutions of the large-scale structural connectivity matrix and the consistency threshold should be specified. Generalized Additive Models (GAM) were used for samples from the HCP-D dataset, and Generalized Additive Mixed Models (GAMM) were used for samples from the ABCD dataset. For each model, the structural connectivity strength per edge was set as the dependent variable, age was set as the smooth term, and sex, and head motion were set as covariates. Knots for the smooth function were set at three, and the restricted maximum likelihood approach was used to estimate smoothing parameters. Sensitivity analyses were conducted by setting different consistency thresholds (`CVthr`) or downsample resolutions (`ds.resolution`).

In this script, statistical parameters were first calculated and saved as `gamresults78_sumSCinvnode_over8_CV75.rds` (12x12 matrix, P75th threshold, can be found in `/wd/interfileFolder_*`). Then, the GAM models were saved as `gammodel78_sumSCinvnode_over8_CV75.rds` (12x12 matrix, P75th threshold). Next, fitted values at 1,000 equispaced values of age were generated and saved as `plotdatasum.df_SA12_sumSCinvnode_CV75.rds`. To avoid the influence of the averaged weight of each edge's strength on derivative analysis or visualization, the structural strength of each edge was divided by its fitted value at the onset age (8 for the HCP-D dataset; 8.9 for the ABCD dataset). The scaled value represents the ratio of structural connectivity strength relative to its initial strength in the interested age span. The scaled data was saved as `SCdata.diw_SA12CV75.rds`. Finally, GAM models, as well as the statistical parameters, were computed and saved as `gamresults78_sumSCinvnode_over8_CV75_scale_TRUE.rds` and `gammodel78_sumSCinvnode_over8_CV75_scale_TRUE.rds` (can be found in `/wd/interfileFolder_*`).

`V1st_fitgammodels_SA_ds_sumSCinvnode_*_covariates.R` are scripts for sensitivity analyses with additional covariates.

## S2nd_calculatederivative_*.R
This script generates derivative and posterior derivative values for scaled GAM or GAMM models. We sampled 1,000 age points at equal intervals from the age range of the two datasets. For each age point, the 1st derivative was computed and saved as `derivative.df78_CV75.rds` (can be found in `/wd/results_*`). Additionally, to ensure the reliability of alignment between derivatives and the connectional axis, we drew 1,000 samples from the posterior distribution of each edge's age smooth function and saved them as `derivative.posterior.df.SA12_CV75.rds` (found in `/wd/results_*`). Scaled GAM or GAMM objects were used to generate derivative and posterior derivative values.

`V2nd_calculatederivative_*_covariates.R` are scripts used to calculate derivatives for models adding additional covariates.

## S3rd_visualizationfitSCcurves_SA12sumSCinvnode_*.R
This script elucidates the developmental trajectories of 78 connections, subsequently aggregating them into deciles based on the S-A connectional axis. First, it illustrates the developmental trajectories by depicting the ratio of structural connectivity strength relative to its initial value at the youngest age within the observed age range (`Fig. 2b`) alongside z-scored fits (`Fig. 2d`) for all 78 connections. Next, to delineate the diverse developmental patterns along the S-A connectional axis, it averages the model fits for connections within each decile of the axis followed by the calculation of z-scores for these averages (`Fig. 3d`).

## S4th_correlationTo_SArank_SA12sumSCinvnode_*.R
This script conducts Spearman correlation analyses between developmental parameters and the S-A connectional axis rank. This script also conducts sensitivity analyses of different consistency thresholds and controling for Euclidean distance.  It performs regression analyses to remove the effect of Euclidean distance from the average second derivatives and the partial R-squared values for the 78 connections. Subsequently, correlation analyses were carried out between the residuals and the connectional axis to further investigate their relationship.

`V4th_correlationTo_SArank_SA_ds_sumSCinvnode_*.R` are scripts for sensitivity analyses of different SC network resolutions.
`V4th_correlationTo_SArank_control_covariates_sumSCinvnode_*.R` are scripts for sensitivity analyses with additional covariates.

## V1st_check_k.R
This script was utilized to select the optimal (k) values for the smoothing functions in GAM and GAMM. We tested (k) values ranging from 3 to 6, and observed that models with a (k) value of 3 consistently exhibited the lowest Akaike Information Criterion (AIC), indicating optimal fit. Consequently, we confirmed the optimality of k=3 for age smooth term in both GAM and GAMM.

## computeEuclidistance.R
Compute Euclidean distance between pairwise systems in 12x12 matrix.

## sup_sigderivative_HCPD.R
This script visualized the age window during which the SC strength develops significantly as shown in `Fig. S5`. The age windows were identified based on the significance of first derivatives.


