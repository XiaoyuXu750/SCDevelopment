clc;
clear;
% this script is to downsample S-A cortex to 17 parcels
%% schaefer 400*7
schaefer400dlabel=cifti_read('/Users/xuxiaoyu_work/Cuilab/Atlas/CBIG-stable_projects-brain_parcellation-Schaefer2018_LocalGlobal/HCP/fslr32k/cifti/Schaefer2018_400Parcels_7Networks_order.dlabel.nii');
schaefer400dlabel_17axis = schaefer400dlabel;
%SA rank
schaefer400_SA=readtable('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD/schaefer400_index_SA.SAorder.delLM.csv');
schaefer400_SA.SAaxisrank = [1:376]';
edges = linspace(1, 376, 18);
schaefer400_SA.SA17node = discretize(schaefer400_SA.SAaxisrank, edges);

schaefer400dlabel_cdata = schaefer400dlabel.cdata;
schaefer400dlabel17_cdata = zeros(length(schaefer400dlabel_cdata),1);
for i = 1:17
    regionidx = schaefer400_SA.index(find(schaefer400_SA.SA17node == i));
    
    % 计算区间的起始索引和结束索引
    Idx = find(ismember(schaefer400dlabel_cdata, regionidx));
    
    % 对区间内的元素赋予编码
    schaefer400dlabel17_cdata(Idx) = i;
end

tabulate(schaefer400dlabel17_cdata)
schaefer400dlabel_17axis.cdata=schaefer400dlabel17_cdata;

tablelabel=schaefer400dlabel.diminfo{1,2}.maps.table;
tablenew=struct;
tablenew(1).name='???'; tablenew(1).key=0; tablenew(1).rgba=tablelabel(1).rgba;
colors = brewermap(17, 'RdBu');
colors = colors([17:-1:1],:);
for i=1:17
    tablenew(i+1).name=['SA', num2str(i)]; 
    tablenew(i+1).key=i; 
    tablenew(i+1).rgba=[colors(i,:), 1];
end
schaefer400dlabel_17axis.diminfo{1,2}.maps.table=tablenew;

SAatlaspath='/Users/xuxiaoyu_work/Cuilab/Atlas/S_A_axis/S-A_ArchetypalAxis/FSLRVertex';
cifti_write(schaefer400dlabel_17axis, [SAatlaspath, '/SensorimotorAssociation_schaefer400_Axis17.dlabel.nii']);

%% schaefer400_SA17 fslr to fsaverage
fslrpath='/Users/xuxiaoyu_work/Cuilab/Atlas/S_A_axis/S-A_ArchetypalAxis/FSLRVertex/';
workbench='/Applications/workbench/bin_macosx64/';
freesurferpath='/Users/xuxiaoyu_work/Cuilab/Atlas/S_A_axis/S-A_ArchetypalAxis/FSaverage5/';

ciftiinput=[fslrpath, 'SensorimotorAssociation_schaefer400_Axis17.dlabel.nii'];
ciftioutput_L=[fslrpath, 'SensorimotorAssociation_schaefer400_Axis17_LH.fslr32k.label.gii'];
ciftioutput_R=[fslrpath, 'SensorimotorAssociation_schaefer400_Axis17_RH.fslr32k.label.gii'];

system([workbench, 'wb_command -cifti-separate ', ciftiinput,...
    ' COLUMN -label CORTEX_LEFT ', ciftioutput_L]);
system([workbench, 'wb_command -cifti-separate ', ciftiinput,...
    ' COLUMN -label CORTEX_RIGHT ', ciftioutput_R]);

currentsphere=[fslrpath, 'fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii'];
newsphere=[freesurferpath, 'fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii'];
metricout_L=[freesurferpath, 'SensorimotorAssociation_schaefer400_Axis17_LH.fsaverage5.label.gii'];
metricout_R=[freesurferpath, 'SensorimotorAssociation_schaefer400_Axis17_RH.fsaverage5.label.gii'];

currentarea=[fslrpath, 'fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii'];
newarea_L=[freesurferpath, 'fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii'];
newarea_R=[freesurferpath, 'fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii'];

system([workbench, 'wb_command -label-resample ', ciftioutput_L, ' ',...
    currentsphere, ' ', newsphere, ' ADAP_BARY_AREA ', metricout_L,...
    ' -area-metrics ', currentarea, ' ', newarea_L]);
system([workbench, 'wb_command -label-resample ', ciftioutput_R, ' ',...
    currentsphere, ' ', newsphere, ' ADAP_BARY_AREA ', metricout_R,...
    ' -area-metrics ', currentarea, ' ', newarea_R]);
%% run in shell
% cp /mnt/nas/xuxiaoyu/Atlas/SensorimotorAssociation_schaefer400_Axis17* /mnt/f/software/freesurfer/subjects/fsaverage5/label/
% cd /mnt/f/software/freesurfer/subjects/fsaverage5/label/
% mris_convert --annot SensorimotorAssociation_schaefer400_Axis17_LH.fsaverage5.label.gii ../surf/lh.sphere SensorimotorAssociation_schaefer400_Axis17_LH.fsaverage5.annot
% mris_convert --annot SensorimotorAssociation_schaefer400_Axis17_RH.fsaverage5.label.gii ../surf/rh.sphere SensorimotorAssociation_schaefer400_Axis17_RH.fsaverage5.annot
% 

%% extract centroid
addpath(genpath('/Users/xuxiaoyu_work/Cuilab/matlab_toolbox/rotate_parcellation/'));
path_annot_lh = [freesurferpath,'lh.SensorimotorAssociation_schaefer400_Axis17_LH.fsaverage5.annot'];
path_annot_rh = [freesurferpath,'rh.SensorimotorAssociation_schaefer400_Axis17_RH.fsaverage5.annot'];
path_sphere_lh = [freesurferpath,'lh.sphere'];
path_sphere_rh = [freesurferpath,'rh.sphere'];

sphere_SA17_lh = centroid_extraction_sphere(path_sphere_lh, path_annot_lh);
sphere_SA17_rh = centroid_extraction_sphere(path_sphere_rh, path_annot_rh);
sphere_SA17_lh(1, :) =[];
sphere_SA17_rh(1, :) =[];

resultFolder = '/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/';
save('sphere_SA17_lh.mat', 'sphere_SA17_lh');
save('sphere_SA17_rh.mat', 'sphere_SA17_rh');
writematrix(sphere_SA17_lh, [resultFolder, 'sphere_SA17_lh.csv']);
writematrix(sphere_SA17_rh, [resultFolder, 'sphere_SA17_rh.csv']);






