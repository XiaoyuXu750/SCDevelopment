# 5th_cognition
This folder contains scripts to evaluate the associations between fluid cognition and SC strength, conduct correlation analysis between the association effect and the S-A connectional axis, and depict the developmental trajectories varying by different cognition levels. For sensitivity analyses, these analyses were repeated using SC strength of large-scale SC networks filtered based on the P25th coefficient of variation threshold. Additionally, both 7 and 17 resolutions of large-scale SC networks were employed in the repeated analyses. Furthermore, sensitivity analyses, including adding additional covariates and controlling for Euclidean distance, were also conducted.


## S1st_testCog_correlation_wholesample_*.Rmd
This script aims to assess the relationship between fluid cognition and SC strength in the HCP-D or ABCD dataset. GAMs were employed to control for potential confounders such as sex, head motion, and age smoothing. Subsequently, the GAM analysis was conducted to investigate the significance of the association between fluid cognition and SC strength. In cases where the strength of at least one edge has a significant association with cognition, further analysis involved computing the correlation coefficient between the effect size of the association, measured by *T* values, and the rank of the S-A connectional axis. `Fig. 4B,C` were created by this script.

`V1st_testCog_correlation_wholesample_resolution_ABCD.R`
Script for sensitivity analyses for different resolutions of large-scale SC matrices.
`V1st_testCog_correlation_wholesample_add_covariate_ABCD.R`
Script for sensitivity analyses for adding additional covariates.

## S2nd_compositescorePlot_scatterplot.Rmd
This script generated the schema of fluid cognition components (`Fig. 4A`) and the scatter plots between SC strength and fluid cognition for 3 exemplified connections (`Fig. 4D`). Given that the associations were only significant in the ABCD dataset, the scatter plots were only depicted using the ABCD dataset.

## S3rd_SCdev_vary_by_cognition.Rmd
This script depicts developmental trajectories varying by baseline cognition levels in the ABCD dataset. We fitted an age-by-cognition interaction GAM model for each connection. Using the acquired models, we estimated SC strength by assigning cognitive performance as low and high levels, respectively. To define these levels, we used the 10th percentile of baseline cognitive performance for the low level and the 90th percentile for the high level. We then averaged trajectories for low and high cognition levels independently within deciles of the S-A connectional axis for visualization purposes. This script generates `Fig. 4E`.



