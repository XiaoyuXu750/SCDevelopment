# 6th_pfactor
This folder contains scripts to evaluate the associations of a general psychopathological factor (*p*-factor) and SC strength, conduct correlation analysis between the association effect and the S-A connectional axis, and depict the developmental trajectories varying by different *p*-factor levels. For sensitivity analysis, these analyses were repeated using SC strength of the large-scale SC networks, which were filtered based on the P25th coefficient of variation threshold. Additionally, sensitivity analyses for different resolutions of large-scale SC networks, controling for Euclidean distance of systems, and adding additional covariates in GAMMs were conducted. 

## S1st_pfactor_effect_continuous_ABCD.Rmd
This script aims to assess the relationship between *p*-factor and SC strength in the ABCD dataset. GAMMs were employed to control for potential confounders such as sex, head motion, and age smoothing. Subsequently, the GAMM analysis was conducted to investigate the significance of the association between *p*-factor and SC strength. In cases where the association proved significant for at least one edge, further analysis involved computing the correlation coefficient between the effect size of the association, measured by *T* values, and the rank of the S-A connectional axis. `Fig. 5B~E` was created by this script.

`V1st_pfactor_effect_continuous_ABCD_resolution.R`
Sensitivity analyses for different resolutions of large-scale SC matrices.
`V1st_pfactor_effect_continuous_ABCD_add_covariate.R`
Sensitivity analyses for additional covariates.

