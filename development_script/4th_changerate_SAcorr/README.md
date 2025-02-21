# S4th_changerate_SAcorr
This folder contains scripts designed for age-resolved developmental analyses and separate analyses within subdatasets partitioned based on the flip age. This step aims to comprehensively understand the spatial correlations between development and the S-A connectional axis throughout the developmental process. For sensitivity analyses, the procedure in step 1 was replicated by utilizing 7 or 17 resolutions of large-scale SC matrices, controlling Euclidean distance, and adding additional covariates.

## S1st_SAcorr_alongAge_*.Rmd
This script facilitates age-resolved developmental effect analysis. Spearman’s correlations were computed between the rank of the S-A connectional axis and the developmental rate assessed by the first derivative at 1,000 specific age points sampled across age spans. Next, correlation coefficients between the S-A connectional axis rank and the posterior derivative/derivative at each age point were determined. The resulting distribution of posterior correlation coefficients was employed to ascertain the median and 95% credible interval of alignment at each age point. `Fig. 4a~d`, `Fig. 5e,i,l` were generated using this script. 

`V1st_SAcorr_alongAge_*_resolution.R`
Validation script to conduct age-resolved analysis using a large-scale SC network of different resolutions.
`V1st_SAcorr_alongAge_*_control_distance.R`
Validation script to conduct age-resolved analysis while controlling for Euclidean distance.
`V1st_SAcorr_alongAge_*_add_covariates.R`
Validation script to conduct age-resolved analysis after adding additional covariates.

## S2nd_fitgammodels_SA12sumSCinvnode_ageseperate_HCPD.Rmd
At the age point where the alignment of developmental rates with the S-A connectional axis experiences a reversal, the HCP-D dataset is partitioned into two sub-datasets. Next, developmental effects were estimated for each connection within the two subdatasets. Following this, correlations between the age effect size, as measured by partial R-square, and the connectional axis rank were computed. `Fig. 4e,f` were generated.


