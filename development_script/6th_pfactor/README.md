# 6th_pfactor
This folder contains scripts to evaluate the associations of a general psychopathological factor (*p*-factor) and SC strength, conduct correlation analysis between the association effect and the S-A connectional axis, and depict the developmental trajectories varying by different cognition levels. For sensitivity analysis, these analyses were repeated using SC strength of the large-scale SC networks, which were filtered based on the P25th coefficient of variation threshold. Additionally, both 7 and 17 resolutions of large-scale SC networks were employed in the repeated analyses. 

## S1st_pfactor_effect_continuous.R
This script aims to assess the relationship between *p*-factor and SC strength in the ABCD dataset. GAMM were employed to control for potential confounders such as gender, race, handedness, head motion, and age smoothing. Subsequently, the GAMM analysis was conducted to investigate the significance of the association between *p*-factor and SC strength. In cases where the association proved significant for at least one edge, further analysis involved computing the correlation coefficient between the effect size of the association, measured by *T* values, and the rank of the S-A connectional axis. `Fig. 7b~e` was created by this script.

`V1st_pfactor_effect_continuous_resolution.R`
Validation script for different resolutions of large-scale SC matrices.

`V2nd_pFactor_SC_cognition_mediation.R`
This script conducted a linear mediation analysis for the model "*p*-factor --> SC strength --> fluid cognition". `Supplementary Fig. S10` was generated accordingly.



