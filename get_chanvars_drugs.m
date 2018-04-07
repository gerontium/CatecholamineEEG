% Ger Loughnane 12/12/16
% Computes and saves the channel variances for each block, so these
% can be plotted to look for bad channels using check4badchans.m.
% Variance calc from FFT spectrum so can avoid certain frequencies.
% list ALL BDF files in a study in a cell array of cell arrays such that
% files{s}{n} is the name of the nth file of session (subject) s.
% also list session IDs (subject initials or whatever) in cell array sessionID.

clear
close all
clc

path_temp = 'E:\Drugs\';
path_temp2 = 'C:\Users\loughnge\Dropbox\Drugs\Drug Analysis\Analysis_12_12_16\';

% define drug conditions
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
single_participants = [];

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
for s = 1:length(subject_folder)
    for sesh = 1:4
        if ~isnan(p3b_drug_order(s,sesh))
            p3b_files{s}{sesh} = [path_temp subject_folder{s} '\' allsubj{s} '_p3b_session' num2str(sesh) '_'  drug_conds{p3b_drug_order(s,sesh)} '.bdf'];
        else
            p3b_files{s}{sesh} = '';
        end
    end
    p3b_matfiles{s} = [path_temp subject_folder{s} '\' allsubj{s} '_p3b_chanvars.mat'];
end

% how much of the spectrum to use?
speclims = [0.5 40];  % Limits in Hz

for s=1:length(allsubj)    
    disp(s), disp(allsubj{s})
    
    % p3b
    files1 = p3b_files{s};
    matfiles1 = p3b_matfiles{s};
        
    clear chanVar
    for b=1:length(files1)
        if isempty(files1{b}), continue, end
        % For the purposes of looking for bad channels, it seems most sensible to leave the BDF referenced as it was recorded.
        % If we average-reference, a bad channel's badness is diluted and may spread to other channels.
        % With a single reference channel, it would be ok, as long as that channel is clean.
        EEG = pop_biosig(files1{b}, 'blockepoch', 'off','channels',[1:nchan]);
        EEG.data = detrend(EEG.data')';
        
        % Fish out the event triggers and times
        clear trigs stimes
        for i=1:length(EEG.event)
            trigs(i)=EEG.event(i).type;
            stimes(i)=EEG.event(i).latency;
        end
        temp = abs(fft(EEG.data(:,stimes(1):stimes(end))'))'; % FFT amplitude spectrum
        tempF = [0:size(temp,2)-1]*EEG.srate/size(temp,2); % Frequency scale
        chanVar(:,b) = mean(temp(:,find(tempF>speclims(1) & tempF<speclims(2))),2);       % ROW of variances        
    end
    save(matfiles1,'chanVar')
end
