clc;
clear;
% this script is to downsample S-A cortex to 12 parcels
% schaefer 400*7
schaefer400dlabel=cifti_read('/Users/xuxiaoyu_work/Cuilab/Atlas/CBIG-stable_projects-brain_parcellation-Schaefer2018_LocalGlobal/HCP/fslr32k/cifti/Schaefer2018_400Parcels_7Networks_order.dlabel.nii');
schaefer400dlabel_SAaxis = schaefer400dlabel;
rgbID = readtable('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCP/RDBU_rgb.csv');
%SA rank
schaefer400_SA=readtable('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCP/schaefer400_index_SA.SAorder.delLM.csv');
schaefer400dlabel_cdata = schaefer400dlabel.cdata;
schaefer400dlabelSA_cdata = zeros(length(schaefer400dlabel_cdata),1);
for i = 1:376
    SArank = schaefer400_SA.finalrank_wholebrainrank(i);
    regionidx = schaefer400_SA.index(i);
    % 计算区间的起始索引和结束索引
    Idx = find(schaefer400dlabel_cdata==regionidx);
    
    % 对区间内的元素赋予编码
    schaefer400dlabelSA_cdata(Idx) = i;
end
tabulate(schaefer400dlabelSA_cdata)
schaefer400dlabel_SAaxis.cdata=schaefer400dlabelSA_cdata;

tablelabel=schaefer400dlabel.diminfo{1,2}.maps.table;
tablenew=struct;
tablenew(1).name='???'; tablenew(1).key=0; tablenew(1).rgba=tablelabel(1).rgba;
for i=1:376
    tablenew(i+1).name=['SA', num2str(i)]; 
    tablenew(i+1).key=i; 
    tablenew(i+1).rgba=table2array(rgbID(i,:));
end
schaefer400dlabel_SAaxis.diminfo{1,2}.maps.table=tablenew;
SAatlaspath='/Users/xuxiaoyu_work/Cuilab/Atlas/S_A_axis/S-A_ArchetypalAxis/FSLRVertex';
cifti_write(schaefer400dlabel_SAaxis, [SAatlaspath, '/SensorimotorAssociation_schaefer400_RdBu.dlabel.nii']);
