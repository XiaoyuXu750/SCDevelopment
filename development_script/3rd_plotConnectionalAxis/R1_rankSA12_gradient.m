clc;
clear;

% Rank the SA 12 systems based on FC gradients
SAAtlaspath = 'D:\xuxiaoyu\Atlas\S_A_axis\S-A_ArchetypalAxis\FSLRVertex\';

% load atlas
SAall = cifti_read([SAAtlaspath, 'SensorimotorAssociation_Axis.dscalar.nii']);
SA12 = cifti_read([SAAtlaspath, 'SensorimotorAssociation_schaefer400_Axis12.dlabel.nii']);
G1fMRI = cifti_read([SAAtlaspath, 'G1.fMRI.dscalar.nii']);
T1T2ratio = cifti_read([SAAtlaspath, 'T1T2ratio.dscalar.nii']);

SA12_L = cifti_struct_dense_extract_surface_data(SA12,'CORTEX_LEFT');
SA12_R = cifti_struct_dense_extract_surface_data(SA12,'CORTEX_RIGHT');
SA12_Label = [SA12_L;SA12_R];

SAall_L = cifti_struct_dense_extract_surface_data(SAall,'CORTEX_LEFT');
SAall_R = cifti_struct_dense_extract_surface_data(SAall,'CORTEX_RIGHT');
SAall_Label = [SAall_L;SAall_R];
SAall_Label = SAall_Label(:, 1);

% G1
G1fMRI_L = cifti_struct_dense_extract_surface_data(G1fMRI,'CORTEX_LEFT');
G1fMRI_R = cifti_struct_dense_extract_surface_data(G1fMRI,'CORTEX_RIGHT');
G1fMRI_Label = [G1fMRI_L;G1fMRI_R];
G1fMRI_Label = G1fMRI_Label(:, 1);

nonzeroidx = find(SA12_Label > 0);

MeanG1 = groupsummary(G1fMRI_Label(nonzeroidx), SA12_Label(nonzeroidx), 'mean');
MedianG1 = groupsummary(G1fMRI_Label(nonzeroidx), SA12_Label(nonzeroidx), 'median');

NetworkLabel = ["SA1", "SA2", "SA3", "SA4", "SA5", "SA6", "SA7", "SA8", "SA9", "SA10", "SA11", "SA12"]';


% T1T2 ratio
T1T2ratio_L = cifti_struct_dense_extract_surface_data(T1T2ratio,'CORTEX_LEFT');
T1T2ratio_R = cifti_struct_dense_extract_surface_data(T1T2ratio,'CORTEX_RIGHT');
T1T2ratio_Label = [T1T2ratio_L;T1T2ratio_R];
T1T2ratio_Label = T1T2ratio_Label(:, 1);

nonzeroidx = find(SA12_Label > 0);

MeanT1T2ratio = groupsummary(T1T2ratio_Label(nonzeroidx), SA12_Label(nonzeroidx), 'mean');
MedianT1T2ratio = groupsummary(T1T2ratio_Label(nonzeroidx), SA12_Label(nonzeroidx), 'median');

% Summarize
SA12_All = table(NetworkLabel, MeanG1, MedianG1, MeanEvolution, MedianEvolution, MeanT1T2ratio, MedianT1T2ratio);
SA12_All.rankG1 = tiedrank(SA12_All.MeanG1);
SA12_All.rankEvolution = tiedrank(SA12_All.MeanEvolution);
SA12_All.rankT1T2ratio = tiedrank(SA12_All.MeanT1T2ratio);

writetable(SA12_All, 'D:\xuxiaoyu\DMRI_network_development\SC_development\interdataFolder_HCPD\SA12_OrtherOrganization.csv');

% Correlation
X = [SAall_Label(nonzeroidx), G1fMRI_Label(nonzeroidx), T1T2ratio_Label(nonzeroidx)];
[rho, pval] = corr(X, ...
    'Type', 'Spearman', ...
    'Rows', 'pairwise');

%% Generate SA12 atlas
% Gradient 1
G1_12 = SA12;

for i = 1:12
    index = find(G1_12.cdata == i);
    G1value = SA12_All.MeanG1(i);
    G1_12.cdata(index) = G1value;
end

G1_12.diminfo{2} = cifti_diminfo_make_scalars(size(G1_12.cdata, 2));
cifti_write(G1_12, [SAAtlaspath '\MeanGradient1_SA12Parcellate.dscalar.nii']);

% T1T2ratio
T1T2ratio_12 = SA12;

for i = 1:12
    index = find(T1T2ratio_12.cdata == i);
    T1T2ratiovalue = SA12_All.MeanT1T2ratio(i);
    T1T2ratio_12.cdata(index) = T1T2ratiovalue;
end

T1T2ratio_12.diminfo{2} = cifti_diminfo_make_scalars(size(T1T2ratio_12.cdata, 2));
cifti_write(T1T2ratio_12, [SAAtlaspath '\MeanT1T2ratio_SA12Parcellate.dscalar.nii']);



