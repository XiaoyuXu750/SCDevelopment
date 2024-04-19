# SCDevelopment
Data and codes for our paper "Structural Connectivity Maturation Spatiotemporally Progresses Along the Sensorimotor-association Connectional Axis in Youth".

Most codes are written in R, and a few in Matlab. The origin data for these analyses are available via the NIMH Data Archive (NDA),  [the Adolescent Brain Cognitive Development, (ABCD)]([https://nda.nih.gov/abcd](https://nda.nih.gov/abcd)), [the Lifespan Human Connectome Project in Development (HCP-D)](https://nda.nih.gov/ccf)

## demopath
This folder contains the demographic, cognitive, and psychopathological characteristics of the participants in the ABCD and HCP-D datasets. `DemodfScreenFinal.csv` is for the ABCD dataset. `HCPD_demo_behav.csv` is for the HCP-D dataset.

## derivatives
This folder contains statistical magnitudes derived from the analyses. Four sub-folders are in it. The abbreviation after the underline indicates an analysis of which dataset the derivatives came from.
* [interdataFolder_ABCD](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/derivatives/interdataFolder_ABCD) : Statistical magnitudes derived from general additive mixed models (GAMM) or general additive models (GAM). 
*  [interdataFolder_HCPD](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/derivatives/interdataFolder_ABCD) : Statistical magnitudes derived from GAM fitted in the HCP-D dataset. 
* [results_ABCD](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/derivatives/results_ABCD) : Correlation coefficients with posterior derivatives at 1,000 specific age points within the age span.
* [results_HCPD](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/derivatives/results_ABCD) : Correlation coefficients with posterior derivatives at 1,000 specific age points within the age span.

## development_script
This folder contains codes for all the analyses in the manuscript, as well as the sensitivity analyses.

* The [1st_dataclean](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/development_script/1st_dataclean) folder contains codes to extract the SC strength of large-scale SC matrics and conduct ComBat analyses.
* The [2nd_fitdevelopmentalmodel](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/development_script/2nd_fitdevelopmentalmodel) folder contains codes to fit GAM or GAMM to evaluate developmental effects and compare the spatial patterns with the sensorimotor-association (S-A) connectional axis. 
* The [3rd_plotConnectionalAxis](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/development_script/3rd_plotConnectionalAxis) folder contains codes to visualize the S-A connectional axis and S-A cortical axis.
* The [4th_changerate_SAcorr](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/development_script/4th_changerate_SAcorr) folder contains codes designed for age-resolved developmental analyses and separate analyses within datasets partitioned based on the age threshold where the alignment of developmental effects undergoes a reversal.
* The [5th_cognition](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/development_script/5th_cognition) folder contains codes to conduct cognitive analyses.
* The [6th_pfactor](https://github.com/XiaoyuXu750/SCDevelopment/tree/main/development_script/6th_pfactor) folder contains codes to conduct psychopathological analyses.

See the `README` in each folder for details of individual scripts.

## gamfunction
This folder contains functions called in the analyses. The functions called in the scripts in the `development_script` folder can be found here.

## GeneralRfunctions
This folder contains the functions for spin tests from https://github.com/frantisekvasa/rotate_parcellation.