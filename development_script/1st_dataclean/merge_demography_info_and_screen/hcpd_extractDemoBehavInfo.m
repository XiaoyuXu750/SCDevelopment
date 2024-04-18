clc;
clear;

behaviorFolder = '/Users/xuxiaoyu_work/Cuilab/open_dataset_information/HCP/HCPD_behavior';
sublistFolder = '/Users/xuxiaoyu_work/Cuilab/open_dataset_information/HCP/HCPD_SC_info';

% sublistFileID = fopen([sublistFolder, '/sublist.txt']);
% sublist = textscan(sublistFileID, '%s');
% sublist = sublist{1};
% fclose(sublistFileID);
% completeness information for HCP-D
% include subjects with complete dMRI & normal anatomical images
HCPD_subcomp = readtable([sublistFolder, '/HCD_LS_2.0_subject_completeness.csv']);
n = 0;
for i = 1:height(HCPD_subcomp)
    if HCPD_subcomp.dMRI_Compl(i) == 1 && contains(HCPD_subcomp.QC_Issue_Codes{i}, "A") == 0
        n = n + 1;
        subin{n, 1} = HCPD_subcomp.src_subject_id{i};
        subin{n, 2} = HCPD_subcomp.unrelated_subset{i};
        subin{n, 3} = HCPD_subcomp.interview_age(i);
        subin{n, 4} = HCPD_subcomp.sex{i};
    end
end

% add age, gender, handness
hand = readcell([behaviorFolder, '/edinburgh_hand01.txt']);
handscore = cell(653, 2);
for i = 3:length(hand)
    handscore{i-1, 1} = hand{i, 6};
    if sum(ismissing(hand{i, 9})) == 0
        if contains(hand{i, 9}, "r")
            r1 = 1;
            l1 = 0;
        elseif contains(hand{i, 9}, "l")
            r1 = 0;
            l1 = 1;
        else
            r1 = 0;
            l1 = 0;
        end
    else
        r1 = 0;
        l1 = 0;
    end
    if sum(ismissing(hand{i, 11})) == 0
        if contains(hand{i, 11}, "r")
            r2 = 1;
            l2 = 0;
        elseif contains(hand{i, 11}, "l")
            r2 = 0;
            l2 = 1;
        else
            r2 = 0;
            l2 = 0;
        end
    else
        r2 = 0;
        l2 = 0;
    end
    if sum(ismissing(hand{i, 12})) == 0
        if contains(hand{i, 12}, "r")
            r3 = 1;
            l3 = 0;
        elseif contains(hand{i, 12}, "l")
            r3 = 0;
            l3 = 1;
        else
            r3 = 0;
            l3 = 0;
        end
    else
        r3 = 0;
        l3 = 0;
    end 
    if sum(ismissing(hand{i, 14})) == 0
        if contains(hand{i, 14}, "r")
            r4 = 1;
            l4 = 0;
        elseif contains(hand{i, 14}, "l")
            r4 = 0;
            l4 = 1;
        else
            r4 = 0;
            l4 = 0;
        end
    else
        r4 = 0;
        l4 = 0;
    end
    if sum(ismissing(hand{i, 15})) == 0
        if contains(hand{i, 15}, "r")
            r5 = 1;
            l5 = 0;
        elseif contains(hand{i, 15}, "l")
            r5 = 0;
            l5 = 1;
        else
            r5 = 0;
            l5 = 0;
        end
    else
        r5 = 0;
        l5 = 0;
    end
    if sum(ismissing(hand{i, 16})) == 0
        if contains(hand{i, 16}, "r")
            r6 = 1;
            l6 = 0;
        elseif contains(hand{i, 16}, "l")
            r6 = 0;
            l6 = 1;
        else
            r6 = 0;
            l6 = 0;
        end
    else
        r6 = 0;
        l6 = 0;
    end
    if sum(ismissing(hand{i, 23})) == 0
        if contains(hand{i, 23}, "r")
            r7 = 1;
            l7 = 0;
        elseif contains(hand{i, 23}, "l")
            r7 = 0;
            l7 = 1;
        else
            r7 = 0;
            l7 = 0;
        end
    else
        r7 = 0;
        l7 = 0;
    end
    if sum(ismissing(hand{i, 25})) == 0
        if contains(hand{i, 25}, "r")
            r8 = 1;
            l8 = 0;
        elseif contains(hand{i, 25}, "l")
            r8 = 0;
            l8 = 1;
        else
            r8 = 0;
            l8 = 0;
        end
    else
        r8 = 0;
        l8 = 0;
    end
    if sum(ismissing(hand{i, 26})) == 0
        if contains(hand{i, 26}, "r")
            r9 = 1;
            l9 = 0;
        elseif contains(hand{i, 26}, "l")
            r9 = 0;
            l9 = 1;
        else
            r9 = 0;
            l9 = 0;
        end
    else
        r9 = 0;
        l9 = 0;
    end
    if sum(ismissing(hand{i, 29})) == 0
        if contains(hand{i, 29}, "r")
            r10 = 1;
            l10 = 0;
        elseif contains(hand{i, 29}, "l")
            r10 = 0;
            l10 = 1;
        else
            r10 = 0;
            l10 = 0;
        end
    else
        r10 = 0;
        l10 = 0;
    end
    if sum(ismissing(hand{i, 30})) == 0
        if contains(hand{i, 30}, "r")
            r11 = 1;
            l11 = 0;
        elseif contains(hand{i, 30}, "l")
            r11 = 0;
            l11 = 1;
        else
            r11 = 0;
            l11 = 0;
        end
    else
        r11 = 0;
        l11 = 0;
    end
    
    r = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10 + r11;
    l = l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8 + l9 + l10 + l11;
    handscore{i-1, 2} = ((r-l) / (r+l)) * 100;
end
handscore{1,1} = "subID";
handscore{1,2} = "handscore";


for i = 1:length(subin)
    for j = 2:length(handscore)
        if strcmp(handscore{j, 1}, subin{i,1})
            subin{i,5} = handscore{j,2}; % 第5列是利手得分
        end
    end
end

% add EF measures: dccs, flanker, pattern, listsorting
dccs = readcell([behaviorFolder, '/dccs01.txt']);
flanker = readcell([behaviorFolder, '/flanker01.txt']);
pattern = readcell([behaviorFolder, '/pcps01.txt']);
lswm = readcell([behaviorFolder, '/lswmt01.txt']);

for i = 1:length(subin)
    for j = 3:length(dccs)
        if strcmp(dccs{j, 5}, subin{i,1})
            subin{i, 6} = dccs{j, 76}; % 第6列是dccs得分
        end
    end
    for n = 3:length(flanker)
        if strcmp(flanker{n, 5}, subin{i,1})
            subin{i, 7} = flanker{n,55}; % 第7列是flanker得分
        end
    end
    for m = 3:length(pattern)
        if strcmp(pattern{m, 5}, subin{i,1})
            subin{i, 8} = pattern{m,145}; % 第8列是pattern得分
        end
    end
    for x = 3:length(lswm)
        if strcmp(lswm{x, 5}, subin{i,1})
            subin{i, 9} = lswm{x,142}; % 第9列是lswm得分
        end
    end
end

% add puberty developmental scale score
PDS = readcell([behaviorFolder, '/pds01.txt']);

for i = 1:length(subin)
    for j = 3:length(PDS)
        if strcmp(PDS{j, 5}, subin{i,1})
            subin{i, 10} = PDS{j, 22}; % 第10列是男孩PDS得分
            subin{i, 11} = PDS{j, 23}; % 第11列是女孩PDS得分
        end
    end
end
var_name = cell(1, 11);
var_name{1,1} = "subID";
var_name{1,2} = "unrelated_sub";
var_name{1,3} = "age";
var_name{1,4} = "gender";
var_name{1,5} = "handness";
var_name{1,6} = "dccs";
var_name{1,7} = "flanker";
var_name{1,8} = "pattern";
var_name{1,9} = "lswm";
var_name{1,10} = "PDS_boy";
var_name{1,11} = "PDS_girl";

subin = [var_name; subin];
%save('/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.mat', 'subin');
for i = 1:length(subin)
    for j = 1:11
        if ismissing(subin{i,j})
            subin{i,j} = [];
        end
    end
end
subin_tab=cell2table(subin(2:end,:));
subin_tab.Properties.VariableNames=[subin{1,:}];
%load('/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.mat');
emotion = readcell([behaviorFolder, '/er4001.txt']);
for i = 1:size(subin_tab,1)
    for j = 3:length(emotion)
        if strcmp(emotion{j, 5}, subin_tab.subID{i})
            subin_tab.emotion_correct{i} = emotion{j, 9}; % 第12列是emotion recognition正确个数
            
        end
    end
end

motor = readcell([behaviorFolder, '/tlbx_motor01.txt']);
for i = 1:length(motor)
    indx_loco(i) = sum(ismissing(motor{i,32}));
end
indx_loco = find(indx_loco == 0);
locomotion = motor(indx_loco, :);

for i = 1:length(motor)
    indx_Endurance(i) = sum(ismissing(motor{i,48}));
end
indx_Endurance = find(indx_Endurance == 0);
Endurance = motor(indx_Endurance, :);

for i = 1:size(subin_tab,1)
    for j = 2:length(locomotion)
        if strcmp(locomotion{j, 5}, subin_tab.subID{i})
            subin_tab.locomotion{i} = locomotion{j, 32}; % 第13列是locomotion得分
            
        end
    end
end
for i = 1:size(subin_tab,1)
    for j = 2:length(Endurance)
        if strcmp(Endurance{j, 5}, subin_tab.subID{i})
            subin_tab.endurance{i} = Endurance{j, 48}; % 第14列是endurance得分
            
        end
    end
end


%load('/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.mat');
upps01 = readcell([behaviorFolder, '/upps01.txt']);
cbcl01 = readcell([behaviorFolder, '/cbcl01.txt']);
bisbas01 = readcell([behaviorFolder, '/bisbas01.txt']);

for i = 1:size(subin_tab,1)
    for j = 3:length(upps01)
        if strcmp(upps01{j, 5}, subin_tab.subID{i})
            subin_tab.neg_urgent_Impulse{i} = upps01{j,119}; %negative urgent
            subin_tab.pos_urgent_Impulse{i} = upps01{j,120}; %positive urgent
        end
    end
    for j = 3:length(cbcl01)
        if strcmp(cbcl01{j, 5}, subin_tab.subID{i})
            subin_tab.Attention_Tscore{i} = cbcl01{j,205}; %Attention problem T score
            subin_tab.ADHD_Tscore{i} = cbcl01{j,221}; %ADHD T score
        end
    end
    for j = 3:length(bisbas01)
        if strcmp(bisbas01{j, 5}, subin_tab.subID{i})
            subin_tab.BIS_score{i} = bisbas01{j,33}; %Inhibition score
        end
    end
end


%save('/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.mat', 'subin');
%writecell(subin, '/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.csv');

%load('/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.mat');
wisc_v01=readcell([behaviorFolder, '/wisc_v01.txt']);
wais_iv_part101 = readcell([behaviorFolder, '/wais_iv_part101.txt']);
cogcomp01 = readcell([behaviorFolder, '/cogcomp01.txt']);

for i = 1:size(subin_tab,1)
    for j = 3:size(wisc_v01,1)
        if strcmp(wisc_v01{j, 6}, subin_tab.subID{i})
            subin_tab.IQ_matrix_raw{i} = wisc_v01{j,102}; %matrix raw
            subin_tab.IQ_matrix_scale{i} = wisc_v01{j,692}; %scaled matrix
        end
    end
    for j = 3:size(wais_iv_part101,1)
        if strcmp(wais_iv_part101{j, 6}, subin_tab.subID{i})
            subin_tab.IQ_matrix_raw{i} = wais_iv_part101{j,22}; %matrix raw
            subin_tab.IQ_matrix_scale{i} = wais_iv_part101{j,23}; %scaled matrix
        end
    end
    for j = 3:length(cogcomp01)
        if strcmp(cogcomp01{j, 5}, subin_tab.subID{i})
            subin_tab.Ntb_fluid_uncor{i} = cogcomp01{j,14}; %fluid score unadjusted
            subin_tab.Ntb_crystal_uncor{i} = cogcomp01{j,18}; %crystal score unadjusted
            subin_tab.Ntb_ecs_uncor{i} = cogcomp01{j,26}; %early child score unadjusted
            subin_tab.Ntb_total_uncor{i} = cogcomp01{j,30}; %total score unadjusted
            subin_tab.Ntb_fluid_cor{i} = cogcomp01{j,15}; %fluid score adjusted
            subin_tab.Ntb_crystal_cor{i} = cogcomp01{j,19}; %crystal score adjusted
            subin_tab.Ntb_ecs_cor{i} = cogcomp01{j,27}; %early child score adjusted
            subin_tab.Ntb_total_cor{i} = cogcomp01{j,31}; %total score adjusted
        end
    end
    for j = 3:length(dccs)
        if strcmp(dccs{j, 5}, subin_tab.subID{i})
            subin_tab.dccs_cor{i} = dccs{j, 77}; % dccs age corrected
        end
    end
    for n = 3:length(flanker)
        if strcmp(flanker{n, 5}, subin_tab.subID{i})
            subin_tab.flanker_cor{i} = flanker{n,56}; % flanker age corrected
        end
    end
    for m = 3:length(pattern)
        if strcmp(pattern{m, 5}, subin_tab.subID{i})
            subin_tab.pattern_cor{i} = pattern{m,146}; % pattern age corrected
        end
    end
    for x = 3:length(lswm)
        if strcmp(lswm{x, 5}, subin_tab.subID{i})
            subin_tab.lswm_cor{i} = lswm{x,137}; % lswm age corrected
        end
    end
end

save('/Users/xuxiaoyu_work/Cuilab/open_dataset_information/HCP/HCPD_demo_behav.mat', 'subin');
writetable(subin_tab, '/Users/xuxiaoyu_work/Cuilab/open_dataset_information/HCP/HCPD_behavior/HCPD_demo_behav.csv');

%% add preprocessing info
mean_fd=zeros(size(subin_tab,1),1);
max_fd=zeros(size(subin_tab,1),1);
max_rotation=zeros(size(subin_tab,1),1);
max_translation=zeros(size(subin_tab,1),1);
subin_tab=addvars(subin_tab, mean_fd);
subin_tab=addvars(subin_tab, max_fd, max_rotation, max_translation);

for i = 1:size(subin_tab,1)
    qcjson_dir = ['/Users/xuxiaoyu_work/Cuilab/control_energy/HCP_SC/qc_json/sub-', subin_tab.subID{i}, '.json'];
    if exist(qcjson_dir, 'file')
        qcjson = loadjson(qcjson_dir);
        subin_tab.mean_fd(i) = qcjson.subjects.mean_fd;
        subin_tab.max_fd(i) = qcjson.subjects.max_fd;
        subin_tab.max_rotation(i) = qcjson.subjects.max_rotation;
        subin_tab.max_translation(i) = qcjson.subjects.max_translation;
    else
        subin_tab.mean_fd(i) = NaN;
        subin_tab.max_fd(i) = NaN;
        subin_tab.max_rotation(i) = NaN;
        subin_tab.max_translation(i) = NaN;
    end
end

% demography
behavior_path = '/Users/xuxiaoyu_work/Cuilab/open_dataset_information/HCP/HCPD_behavior/';
demograph=readtable([behavior_path, 'socdem01.txt']);
totalincome=zeros(size(subin_tab,1),1);
familymember=zeros(size(subin_tab,1),1);
fatheredu=zeros(size(subin_tab,1),1);
motheredu=zeros(size(subin_tab,1),1);
race=cell(size(subin_tab,1),1);
subin_tab=addvars(subin_tab, totalincome, familymember, fatheredu, motheredu, race);

for i=1:size(subin_tab,1)
    subID=subin_tab.subID{i};
    j=find(strcmp(demograph.src_subject_id, subID));
    subin_tab.totalincome(i)=demograph.annual_fam_inc(j);%43 family annual income
    subin_tab.familymember(i)=demograph.household_number_in_house(j);%44 family number
    subin_tab.fatheredu(i)=demograph.father_edu_cat(j);%45 father highest education
    subin_tab.motheredu(i)=demograph.mother_edu_cat(j);%46 mother highest education
    subin_tab.race{i}=demograph.race(j);%47 race
end

save('/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.mat', 'subin');
writetable(subin_tab, '/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.csv');

subin_tab=readtable('/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.csv');

for i=1:size(subin_tab,1)
    subID=subin_tab.subID{i};
    j=find(strcmp(demograph.src_subject_id, subID));
    subin_tab.interviewdate(i)=demograph.interview_date(j);%43 family annual income
end

save('/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.mat', 'subin_tab');
writetable(subin_tab, '/Users/xuxiaoyu_work/Cuilab/control_energy/demo_behaveinfo/HCPD_demo_behav.csv');




