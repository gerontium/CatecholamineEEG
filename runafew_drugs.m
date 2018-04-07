% Ger Loughnane 10/10/17
% Run a few subjects through the EEG preprocessing stage.
clear
close all
clc
pause(0.5);

eeglab

path_temp = 'E:\Drugs\';
path_temp2 = 'C:\Users\loughnge\Dropbox\Drugs\Drug Analysis\Analysis_12_12_16\';

% define drug conditions
drug_conds = {'MPH','ATM','CIT','PLA'};
% define subject names
for i = 1:40
    if i<10
        subject_folder{i} = ['P0',num2str(i)];
        allsubj{i} = ['0',num2str(i)];
    else
        subject_folder{i} = ['P',num2str(i)];
        allsubj{i} = [num2str(i)];
    end
end

load('drug_orders_good') % loads p3b_drug_order

duds = []; % eliminate bad subjects
single_participants = []; % only run certain subjects

% define noisy/dead channels from previous scripts
allbadchans = {{[],[],[],[]}, ...
    {[15,16,61],[15,16,61],[15,16,61],[15,16,61]}, ...
    {[],[],[],[]}, ...
    {[24],[24],[24],[24]}, ...
    {[10,60],[10,60],[10,60],[10,60]}, ... % 5
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]}, ... % 10
    {[],[],[15],[]}, ...
    {[],[5],[],[]}, ...
    {[],[],[],[]}, ...
    {[15,61],[15],[15],[15]}, ...
    {[63],[],[],[25,27]}, ... % 15
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[20],[20],[20],[20]}, ... % 20
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[62],[62],[62],[62]}, ...
    {[29],[29],[29,64],[29]}, ...
    {[],[],[],[]}, ... % 25
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[9,15],[9,15],[9,15],[9,15]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]}, ... % 30
    {[],[],[],[7]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[21]}, ...
    {[],[40,63],[],[]}, ...
    {[],[],[],[]}, ... % 35
    {[],[],[],[21,22]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]}, ...
    {[],[],[],[]} % 40
    };

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    p3b_drug_order([duds],:) = [];
    allbadchans([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    p3b_drug_order = p3b_drug_order(single_participants,:);
    allbadchans = allbadchans([single_participants]);
end

%% Perform CSD transform, load CSD parameters
E = textread('C:\Program Files\MATLAB\R2016a\toolbox\CSDtoolbox\chans64.asc','%s');
M = ExtractMontage('C:\Program Files\MATLAB\R2016a\toolbox\CSDtoolbox\10-5-System_Mastoids_EGI129.csd',E);  % reading in the montage for the CSD toolbox
[G2,H2] = GetGH(M,5);

for s=1:length(allsubj)

    disp(['Subject: ',num2str(s)])
    disp(['Subject: ',allsubj{s}])
    
    badchans = allbadchans{s};
    
    clear files; k=0;
    for sesh = 1:4
        k=k+1;
        if ~isnan(p3b_drug_order(s,sesh))
            files{k} = [path_temp subject_folder{s} '\' allsubj{s} '_p3b_session' num2str(sesh) '_'  drug_conds{p3b_drug_order(s,sesh)} '.bdf'];
        else
            files{k} = '';
        end
    end
        
    % run the pre-processing script
    p3b_epochs
end
