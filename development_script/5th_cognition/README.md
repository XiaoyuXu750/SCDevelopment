# 5th_cognition
This folder contains scripts to evaluate the associations of fluid cognition and SC strength, conduct correlation analysis between the association effect and the S-A connectional axis, and depict the developmental trajectories varying by different cognition levels. For sensitivity analysis, these analyses were repeated using SC strength of the large-scale SC networks, which were filtered based on the P25th coefficient of variation threshold. Additionally, both 7 and 17 resolutions of large-scale SC networks were employed in the repeated analyses.

## S1st_testCog_correlation_wholesample_*.Rmd
This script aims to assess the relationship between fluid cognition and SC strength in the HCP-D or ABCD dataset. GAM were employed to control for potential confounders such as gender, race, handedness, head motion, and age smoothing. Subsequently, the GAM analysis was conducted to investigate the significance of the association between fluid cognition and SC strength. In cases where the association proved significant for at least one edge, further analysis involved computing the correlation coefficient between the effect size of the association, measured by *T* values, and the rank of the S-A connectional axis. `Fig. 6b~c` were created by this script.

`V1st_testCog_correlation_wholesample_resolution_ABCD.R`
Validation script for different resolutions of large-scale SC matrices.

## S2nd_compositescorePlot_scatterplot.Rmd
This script generated the schema of fluid cognition components (`Fig. 6a`) and the scatter plots between SC strength and fluid cognition for 3 exemplified connections (`Fig. 6d`). Given the associations were only significant in the ABCD dataset, the scatter plots were only depicted using the ABCD dataset.

## S3rd_SCdev_vary_by_cognition.Rmd
This script depicts developmental trajectories varying by baseline cognition levels in the ABCD dataset. We fitted an age-by-cognition interaction GAM model for each connection. Using the acquired models, we estimated SC strength by assigning cognitive performance as low and high levels respectively. To define these levels, we used the 10th percentile of baseline cognitive performance for the low level, and the 90th percentile for the high level. We then averaged trajectories for low and high cognition levels independently within deciles of the S-A connectional axis for visualization purposes. This script generates `Fig. 6e`.




