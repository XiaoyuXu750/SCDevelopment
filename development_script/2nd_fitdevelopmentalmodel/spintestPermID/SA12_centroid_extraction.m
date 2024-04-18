clc;
clear;
%% schaefer400_SA12 fslr to fsaverage
fslrpath='/Users/xuxiaoyu_work/Cuilab/Atlas/S_A_axis/S-A_ArchetypalAxis/FSLRVertex/';
workbench='/Applications/workbench/bin_macosx64/';
freesurferpath='/Users/xuxiaoyu_work/Cuilab/Atlas/S_A_axis/S-A_ArchetypalAxis/FSaverage5/';

ciftiinput=[fslrpath, 'SensorimotorAssociation_schaefer400_Axis12.dlabel.nii'];
ciftioutput_L=[fslrpath, 'SensorimotorAssociation_schaefer400_Axis12_LH.fslr32k.label.gii'];
ciftioutput_R=[fslrpath, 'SensorimotorAssociation_schaefer400_Axis12_RH.fslr32k.label.gii'];

system([workbench, 'wb_command -cifti-separate ', ciftiinput,...
    ' COLUMN -label CORTEX_LEFT ', ciftioutput_L]);
system([workbench, 'wb_command -cifti-separate ', ciftiinput,...
    ' COLUMN -label CORTEX_RIGHT ', ciftioutput_R]);

currentsphere=[fslrpath, 'fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii'];
newsphere=[freesurferpath, 'fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii'];
metricout_L=[freesurferpath, 'SensorimotorAssociation_schaefer400_Axis12_LH.fsaverage5.label.gii'];
metricout_R=[freesurferpath, 'SensorimotorAssociation_schaefer400_Axis12_RH.fsaverage5.label.gii'];

currentarea=[fslrpath, 'fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii'];
newarea_L=[freesurferpath, 'fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii'];
newarea_R=[freesurferpath, 'fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii'];

system([workbench, 'wb_command -label-resample ', ciftioutput_L, ' ',...
    currentsphere, ' ', newsphere, ' ADAP_BARY_AREA ', metricout_L,...
    ' -area-metrics ', currentarea, ' ', newarea_L]);
system([workbench, 'wb_command -label-resample ', ciftioutput_R, ' ',...
    currentsphere, ' ', newsphere, ' ADAP_BARY_AREA ', metricout_R,...
    ' -area-metrics ', currentarea, ' ', newarea_R]);

%% extract centroid
addpath(genpath('/Users/xuxiaoyu_work/Cuilab/matlab_toolbox/rotate_parcellation/'));
path_annot_lh = [freesurferpath,'lh.SensorimotorAssociation_schaefer400_Axis12_LH.fsaverage5.annot'];
path_annot_rh = [freesurferpath,'rh.SensorimotorAssociation_schaefer400_Axis12_RH.fsaverage5.annot'];
path_sphere_lh = [freesurferpath,'lh.sphere'];
path_sphere_rh = [freesurferpath,'rh.sphere'];

sphere_SA12_lh = centroid_extraction_sphere(path_sphere_lh, path_annot_lh);
sphere_SA12_rh = centroid_extraction_sphere(path_sphere_rh, path_annot_rh);
sphere_SA12_lh(1, :) =[];
sphere_SA12_rh(1, :) =[];

resultFolder = '/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/';
save('sphere_SA12_lh.mat', 'sphere_SA12_lh');
save('sphere_SA12_rh.mat', 'sphere_SA12_rh');
writematrix(sphere_SA12_lh, [resultFolder, 'sphere_SA12_lh.csv']);
writematrix(sphere_SA12_rh, [resultFolder, 'sphere_SA12_rh.csv']);







