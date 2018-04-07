% Ger Loughnane 12/12/16
% Plot EEG channel variances in order to find dead or noisy channels
clear
close all
clc

path_temp = 'E:\Drugs\';
path_temp2 = 'C:\Users\loughnge\Dropbox\Drugs\Drug Analysis\Analysis_12_12_16\';

% define the relevant drug conditions
drug_conds = {'ATM','CIT','MPH','PLA'};
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

duds = [];
single_participants = [1];

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    p3b_drug_order([duds]) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    p3b_drug_order = p3b_drug_order(single_participants);
end

nchan = 64;

for s=1:length(subject_folder)
    matfiles_p3b{s} = [path_temp subject_folder{s} '\' allsubj{s} '_p3b_chanvars.mat'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chanlocs = readlocs('cap64.loc');
for s=1:length(subject_folder)
    disp(subject_folder{s})
    load(matfiles_p3b{s})
    chanVar = double(chanVar);
    
    % make a really bad channel equal to its neighbour for clearer view
%     badchans = [3];
%     changechans = [4]; % must be in same order as badchans.
%     chanVar(badchans(1:end),:) = chanVar(changechans(1:end),:);
    
    avVar = mean(chanVar,2); 

    figure;
    topoplot(avVar,chanlocs,'plotchans',[1:64],'electrodes','numbers','maplimits','maxmin');
    title(subject_folder{s})
    
    figure; hold on
    h = plot(chanVar(1:64,:));
    title(subject_folder{s})
    legend(h,'Location','NorthEast');    
    
    figure; hold on
    h = plot(chanVar(1:64,:));
    set(gca,'YLim',[0 100000])
    title(subject_folder{s})
    legend(h,'Location','NorthEast'); 
end