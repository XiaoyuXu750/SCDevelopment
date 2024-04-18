# 1st_dataclean
This folder contains codes to `extract the SC strength of large-scale SC matrices`, `conduct ComBat to harmonize the SC strength from 13 acquisition sites`, `visualize the average SC strength matrix`, and `generate the age distribution plots for the final samples from the ABCD and HCP-D datasets`. Separate scripts are provided for the ABCD and HCP-D datasets for the first three steps. For sensitivity analyses, users can specify the resolutions of large-scale SC matrices and the consistency threshold at the beginning of the scripts for the second and third steps.

## S1st_mergedata_*.R
`This script generated the filter index for the fine-grained SC network based on the thresholds of the 75th percentile (P75th) or 25th percentile (P25th) coefficient of variation (CV).` Firstly, a dataframe was created where each column represents the streamline counts scaled by the volume of nodes in an edge. The normalized streamline counts were extracted from the fine-grained SC network built on the Schaefer-400 atlas. Edges connecting regions in the limbic system were excluded, resulting in 70,786 edges for the 376*376 network. Secondly, CVs were computed for each edge. We applied consistency-based thresholds to the structural connectivity matrices to minimize the influence of false-positive connections. We selected the thresholds of the P75th and P25th based on a previous study, namely removing streamlines from the top quartile or three quartiles of inconsistent connections. The primary results are based on the structural connectivity matrices thresholded at the P75th CV. We also replicated our primary results after applying the P25th threshold and presented them in Supplementary Fig. 2 and 3.

## S2nd_mergedata_SA_ds_sumSC_*.R
This script extract the SC strength of large-scale SC matrices. A dataframe will be generated, in which each column represents the streamline counts scaled by node volumes of a connection in large-scale structural connectivity matrix (e.g., 78 columns for 12\*12 matrix). Before running, ds.resolution should be determined. ds.resolution represents the resolution of large-scale structural connectivity matrix. In this study, we computed the primary results based on the ds.resolution of 12. Additionally, we also replicated our primary results based on the ds.resolution of 7 and 17 (Supplementary Fig. 4 and 5), which are the common resolutions used in previous studies on large-scale networks. The streamlines from the top quartile (P75th CV threshold) or three quartiles (P25th CV threshold) of inconsistent connections in 376*376 matrix were removed. The remaining streamline counts were summed up and divided by the average node volumes they connected.

## S3rd_combat_controlsite_*.R
We conducted ComBat to harmonize the structural connectivity strength from 13 acquisition sites referring to a code improved by Richard Beare [code](https://github.com/PennLINC/Larsen_IronDevelopment/blob/master/combat.R). ComBat was performed separately for observations included in developmental modes (including covariates for age, gender, handedness, head motion, and race), cognitive models (including covariates for fluid cognition, age, gender, handedness, head motion, and race), and p-factor models (including covariates for p-factor, age, gender, handedness, head motion, and race).

## S4th_plotSCdata_SA12_separateage_HCPD.R
This script plots the average structural connectivity matrices of 12*12 at specific ages. The average matrix is used in `Fig. 1f`.

## S5th_demodescrip_plot.R
This script generates the age distribution plots for the final samples from ABCD and HCP-D datasets. The age distribution plots are presented in `Fig. 1a,b`.

# merge_demography_info_and_screen
This folder contains codes to merge demographic, cognitive, and psychopathologic variables and screen data from the HCP-D and ABCD datasets.

