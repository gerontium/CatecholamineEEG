% Ger Loughnane 20/10/2017
% Plot ERP data across drug conditions, perform statistical analyses, etc.
clear
close all
clc
pause(0.5)

eeglab

chanlocs = readlocs('cap64.loc');

path_temp = 'E:\Drugs\';
xl_path = 'C:\Users\loughnge\Dropbox\Research\Drugs\Drug Analysis\Analysis_12_12_16\stats_excels\';
save_path = 'C:\Users\loughnge\Dropbox\Research\Drugs\Drug Analysis\Plotting_21_06_17\';

drug_conds = {'MPH','ATM','CIT','PLA'};
cond_switch = [2,3,1,4]; % ATM, CIT, MPH, PLA

load('drug_orders_good') % loads p3b_drug_order
load('VAS_measures')

for i = 1:40
    if i<10
        subject_folder{i} = ['P0',num2str(i)];
        allsubj{i} = ['0',num2str(i)];
    else
        subject_folder{i} = ['P',num2str(i)];
        allsubj{i} = [num2str(i)];
    end
end

duds = []; 
for s = 1:40
    if length(find(ismember(p3b_drug_order(s,:),[1:4])))<length([1:4])
        duds = [duds,s];
    end
end
sick_duds = find(sick_subj); % from VAS measures
jess_duds = [2,11,21,22]; % subjects marked out by Jess to have not reacted well to a particular drug
single_participants = [];

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    p3b_drug_order([duds],:) = [];
    alertness([duds],:,:) = [];
    contentedness([duds],:,:) = [];
    calmness([duds],:,:) = [];
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    p3b_drug_order = p3b_drug_order(single_participants,:);
    alertness = alertness(single_participants,:,:);
    contentedness = contentedness(single_participants,:,:);
    calmness = calmness(single_participants,:,:);
end

switch_order=[];
for s = 1:length(allsubj)
    [~,switch_order(s,:)] = sort(p3b_drug_order(s,:),'ascend');
    alertness(s,:,:) = alertness(s,switch_order(s,:),:);
    contentedness(s,:,:) = contentedness(s,switch_order(s,:),:);
    calmness(s,:,:) = calmness(s,switch_order(s,:),:);
end

%% Channel definitions, triggers, etc.
exclude_chans = [];
plot_chans = [1:64];
left_hemi = [1:23,25:27];
right_hemi = [34:36,39:46,49:60,62:64];
centre_chans = [29:32,38,48];
posterior_chans = [16:32,53:64];
alpha_chans = [25,27,62,64];

targcodes = [50,75];

Fs = 512;
numch = 64;
rtlim = [0.15 0.75]; % [0.1 0.6]

ch_CPP = [31,32];
ch_lr = [60,62;23,25];
ch_rl = [23,25;60,62];

CPP_choice_chans = [31,32,48,20,57,19,56,12,49];
t_CPPr = [-25,0];
howmanychans = 1;

plot_mean = zeros(1,64);
plot_mean(CPP_choice_chans) = 1;
figure
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean) max(plot_mean)], ...
    'electrodes','numbers','plotchans',[1:64]);
pause(1)

components = {'P1','N1','P2'};
ch_P1 = [29];
ch_N1 = [22,23,59,60];
ch_P2 = [29,30];
channel_comps = {ch_P1,ch_N1,ch_P2};

%% Timing
% [-1000,-875,-750,-625,-500,-375,-250,-125,0,125,250,375,500,625,750,875,1000,1125,1250,1375,1500,1625,1750,1875,2000]
% stim-locked erps
ts = -0.125*Fs:0.750*Fs;
t = ts*1000/Fs;

trs = [-.250*Fs:Fs*0];
tr = trs*1000/Fs;

BL_erp = [-125,0];
%% RT binning
no_of_bins = 2;
bin_counter = fliplr(1:no_of_bins);
% plotting parameters
for bin = 1:no_of_bins
    rt_bins_tags{bin} = num2str(bin);
end
rt_bins_tags{1} = [rt_bins_tags{1},' (Fast RT)'];
rt_bins_tags{end} = [rt_bins_tags{end},' (Slow RT)'];

%% EEG preprocessing choices
LPF_choices = [100,35,20,8]; LPF_ind = 3;
CSD_choices = {'_CSD',''}; CSD_ind = 1; % 1 for on, 2 for off
art_thresh = 3; % artifth = [25,50,100];
% ARchan_choice{1} = ARchan_all;
% ARchan_choice{2} = [ARchans_for_blinks,ARchans_for_CPP,ARchans_for_FC,ARchans_for_VEP,ARchans_for_M1];
% ARchan_choice{3} = [ARchans_for_blinks,ARchans_for_CPP];
% ARchan_choice{4} = [ARchans_for_blinks,ARchans_for_FC];
% ARchan_choice{5} = [ARchans_for_blinks,ARchans_for_VEP];
% ARchan_choice{6} = [ARchans_for_blinks,ARchans_for_M1];
archan = 2;
RT_thresh = 5; % RT threshold in z 

%% Subject loop

allslopes = []; allmat = []; allmat_behave = [];
for s=1:length(allsubj)
    load([path_temp subject_folder{s} '\' allsubj{s} '_p3b_epochs'],['erp_',num2str(LPF_choices(LPF_ind)),'Hz',CSD_choices{CSD_ind}], ...
        'allRT','allResp','allTrig','allBlock2','t', ...
        'artifchans','artif_trials');
    eval(['erp = erp_',num2str(LPF_choices(LPF_ind)),'Hz',CSD_choices{CSD_ind},';']);

    allBlock = allBlock2;
    artif_trials2 = squeeze(artif_trials(LPF_ind,archan,art_thresh,:))';
        
    % Baseline erp
    baseline_erp = mean(erp(:,find(t>=BL_erp(1) & t<BL_erp(2)),:),2);
    erp = erp-repmat(baseline_erp,[1,size(erp,2),1]); % baseline full erp

    disp(['Subject loop number ' num2str(s) ': ' allsubj{s} ' number of trials = ' num2str(size(erp,3))])

    % get response-locked ERP
    erpr = NaN(size(erp,1),length(tr),size(erp,3));        
    validrlock = zeros(1,length(allRT)); % length of RTs.
    for n=1:length(allRT)
        if allTrig(n)==targcodes(2)
            [blah,RTsamp] = min(abs(t*Fs/1000-allRT(n))); % get the sample point of the RT.
            if RTsamp+trs(1)>0 & RTsamp+trs(end)<=length(t) & allRT(n)>0 % is the RT larger than 1st stim RT point, smaller than last RT point.
                erpr(:,:,n) = erp(:,RTsamp+trs,n);
                validrlock(n)=1;
            end
        end
    end
    
    % Find extreme RTs
    for drug = 1:4
        RT_temp = allRT(find(allTrig==targcodes(2) & allBlock==drug & allResp==1 & ...
            validrlock==1 & artif_trials2==1 & allRT>rtlim(1)*Fs & allRT<rtlim(2)*Fs)); % in samples!!!
        extreme_RT(s,drug) = RT_thresh*std(RT_temp)+mean(RT_temp);
    end

    for drug = 1:4
        % define which trials belonged to which drug condition
        stan_conds_behave{s,drug} = find(allTrig==targcodes(1) & allBlock==drug & allResp==4);
        dev_conds_behave{s,drug} = find(allTrig==targcodes(2) & allBlock==drug & allResp==1 & ...
            allRT>rtlim(1)*Fs & allRT<rtlim(2)*Fs & allRT<extreme_RT(s,drug));        
        stan_conds{s,drug} = find(allTrig==targcodes(1) & allBlock==drug & allResp==4 & ...
            validrlock==0 & artif_trials2==1);
        dev_conds{s,drug} = find(allTrig==targcodes(2) & allBlock==drug & allResp==1 & ...
            validrlock==1 & artif_trials2==1 & allRT>rtlim(1)*Fs & allRT<rtlim(2)*Fs & allRT<extreme_RT(s,drug));
        RTs{s,drug} = allRT([dev_conds_behave{s,drug}])*1000/Fs;
        RTs_EEG{s,drug} = allRT([dev_conds{s,drug}])*1000/Fs;
        hit_rate_drug(s,drug) = 100*length(find(allTrig==targcodes(2) & allBlock==drug & allResp==1 & ...
            allRT>rtlim(1)*Fs & allRT<rtlim(2)*Fs))/ ...
            (length(find(allTrig==targcodes(2) & allBlock==drug & allResp==1 & ...
            allRT>rtlim(1)*Fs & allRT<rtlim(2)*Fs)) + ...
            length(find(allTrig==targcodes(2) & allBlock==drug & allResp==2))); % hit/hit+miss
        fa_rate_drug(s,drug) = 100*length(find(allTrig==targcodes(1) & allBlock==drug & allResp==3 & ...
            allRT>rtlim(1)*Fs & allRT<rtlim(2)*Fs))/ ...
            (length(find(allTrig==targcodes(1) & allBlock==drug & allResp==3 & ...
            allRT>rtlim(1)*Fs & allRT<rtlim(2)*Fs)) + ...
            length(find(allTrig==targcodes(1) & allBlock==drug & allResp==4))); % fa/fa+cr
            
        % define ERPs per drug condition, per subject
        avERP(s,:,:,drug,1) = squeeze(mean(erp(1:numch,:,[stan_conds{s,drug}]),3));
        avERP(s,:,:,drug,2) = squeeze(mean(erp(1:numch,:,[dev_conds{s,drug}]),3));
        avERP_sub(s,:,:,drug) = avERP(s,:,:,drug,2)-avERP(s,:,:,drug,1);
        avERPr(s,:,:,drug) = squeeze(mean(erpr(1:numch,:,[dev_conds{s,drug}]),3));
        
        % define CPP per drug condition, per subject
        [~,indx] = sort(squeeze(mean(avERPr(s,CPP_choice_chans,find(tr>=t_CPPr(1) & tr<=t_CPPr(2)),drug),3)),'descend');
        CPPrs(s,:,drug) = squeeze(mean(avERPr(s,CPP_choice_chans(indx(1:howmanychans)),:,drug),2));
        CPPs_stim(s,:,drug) = squeeze(mean(avERP(s,CPP_choice_chans(indx(1:howmanychans)),:,drug,2),2));
        CPPs_subs(s,:,drug) = squeeze(mean(avERP_sub(s,CPP_choice_chans(indx(1:howmanychans)),:,drug),2));        
        CPPrs_subj{s,drug} = squeeze(mean(erpr(CPP_choice_chans(indx(1:howmanychans)),:,[dev_conds{s,drug}]),1));
        
        % behavioural variables for modelling
        dev_conds_behave_all{s,drug} = find(allTrig==targcodes(2) & allBlock==drug & ismember(allResp,[1,2]));
        RTs_all{s,drug} = allRT([dev_conds_behave_all{s,drug}])*1000/Fs;
        
        % calculate mean RT, etc. per drug condition, per subject
        valid_trials(s,drug) = length([dev_conds{s,drug}]);
        valid_behave_trials(s,drug) = length([dev_conds_behave{s,drug}]);
        RT_EEG_mean_drug(s,drug) = squeeze(mean([RTs_EEG{s,drug}]));
        RT_mean_drug(s,drug) = squeeze(mean([RTs{s,drug}])); 
        RT_SD_drug(s,drug) = squeeze(std([RTs{s,drug}]));
        RT_CVar_drug(s,drug) = RT_SD_drug(s,drug)/RT_mean_drug(s,drug);      
        
        % grab VAS measures for covariance analysis
        VAS_measures(s,drug,:,1) = alertness(s,drug,:);
        VAS_measures(s,drug,:,2) = contentedness(s,drug,:);
        VAS_measures(s,drug,:,3) = calmness(s,drug,:);
    end
            
    disp(['Subject ',allsubj{s},' Total Valid Trials: ',num2str(length([dev_conds{s,:}])+length([stan_conds{s,:}])),' = ', ...
        num2str(round(100*(length([dev_conds{s,:}])+length([stan_conds{s,:}]))/(length(find(ismember(allBlock,[1:4])))))),'%'])
    disp(['Subject ',allsubj{s},' Total Valid Target Trials: ',num2str(length([dev_conds{s,:}])),' = ', ...
        num2str(round(100*length([dev_conds{s,:}])/(length(find(ismember(allBlock,[1:4]) & allTrig==targcodes(2)))))),'%'])    
    disp(['Subject ',allsubj{s},' Mean RT: ',num2str(mean([RTs{s,:}]))])
end

%% ******************** Behavioural first ***************************
%% Check RT distributions etc.
allrts = [RTs{:,:}];
figure
histogram(allrts);
% prctile(allrts,[5,95])
RT_zscores = zscore(RT_mean_drug,[],1);
figure
plot(RT_zscores)

figure
for s = 1:length(allsubj)
    subplot(5,7,s)
    histogram([RTs{s,1}],7), hold on
    histogram([RTs{s,2}],7)
    set(gca,'xlim',[0,1000],'ylim',[0,20])
    title(allsubj{s})
end

data = valid_trials; % subj x drug
tab = table(data(:,1),data(:,2),data(:,3),data(:,4),'VariableNames',{'MPH','ATM','CIT','PLA'});

% define ANOVA for drug condition
mspec = ['MPH,ATM,CIT,PLA~1'];
within = table([1:4]','VariableNames',{'DrugID'});
rm = fitrm(tab,mspec,'WithinDesign',within);
ranovatbl = ranova(rm)

%% RT x drug

drug_error = std(RT_mean_drug,1)/sqrt(size(RT_mean_drug,1));
figure
h = errorbar(mean(RT_mean_drug,1),drug_error)
h.LineWidth = 2;
set(gca,'FontSize',16,'XTick', [1:length(drug_conds)], 'XTickLabel',drug_conds,'XLim',[0.5 length(drug_conds)+0.5])
xlabel('Drug'); ylabel('RT (ms)')
title('Mean RT x Drug')

data = RT_mean_drug; % subj x drug
tab = table(data(:,1),data(:,2),data(:,3),data(:,4),'VariableNames',{'MPH','ATM','CIT','PLA'});

mspec = ['MPH,ATM,CIT,PLA~1'];
within = table([1:4]','VariableNames',{'DrugID'});
rm = fitrm(tab,mspec,'WithinDesign',within);
ranovatbl = ranova(rm)

disp('MPH vs PLA')
[h,p,ci,stats] = ttest(RT_mean_drug(:,1),RT_mean_drug(:,4))
disp('ATM vs PLA')
[h,p,ci,stats] = ttest(RT_mean_drug(:,2),RT_mean_drug(:,4))
disp('CIT vs PLA')
[h,p,ci,stats] = ttest(RT_mean_drug(:,3),RT_mean_drug(:,4))

%% RT SD x drug

drug_var_error = std(RT_SD_drug,1)/sqrt(size(RT_SD_drug,1));
figure
h = errorbar(mean(RT_SD_drug,1),drug_var_error)
h.LineWidth = 2;
set(gca,'FontSize',16,'XTick', [1:length(drug_conds)], 'XTickLabel',drug_conds,'XLim',[0.5 length(drug_conds)+0.5])
xlabel('Drug'); ylabel('RT (ms)')
title('Mean RT SD x Drug')

data = RT_SD_drug; % subj x drug
tab = table(data(:,1),data(:,2),data(:,3),data(:,4),'VariableNames',{'MPH','ATM','CIT','PLA'});

mspec = ['MPH,ATM,CIT,PLA~1'];
within = table([1:4]','VariableNames',{'DrugID'});
rm = fitrm(tab,mspec,'WithinDesign',within);
ranovatbl = ranova(rm)

disp('MPH vs PLA')
[h,p,ci,stats] = ttest(RT_SD_drug(:,1),RT_SD_drug(:,4))
disp('ATM vs PLA')
[h,p,ci,stats] = ttest(RT_SD_drug(:,2),RT_SD_drug(:,4))
disp('CIT vs PLA')
[h,p,ci,stats] = ttest(RT_SD_drug(:,3),RT_SD_drug(:,4))

%% RT coeff var x drug
data = RT_CVar_drug;

[NormData,SDev,SErr,CInt] = within_subj_summary(data);

figure
h = errorbar(mean(NormData,1),CInt)
h.LineWidth = 2;
set(gca,'FontSize',16,'XTick', [1:length(drug_conds)], 'XTickLabel',drug_conds,'XLim',[0.5 length(drug_conds)+0.5])
xlabel('Drug'); ylabel('RT (ms)')
title('Mean RT Coeff Var x Drug')

tab = table(data(:,1),data(:,2),data(:,3),data(:,4),'VariableNames',{'MPH','ATM','CIT','PLA'});

mspec = ['MPH,ATM,CIT,PLA~1'];
within = table([1:4]','VariableNames',{'DrugID'});
rm = fitrm(tab,mspec,'WithinDesign',within);
ranovatbl = ranova(rm)

disp('MPH vs PLA')
[h,p,ci,stats] = ttest(data(:,1),data(:,4))
disp('ATM vs PLA')
[h,p,ci,stats] = ttest(data(:,2),data(:,4))
disp('CIT vs PLA')
[h,p,ci,stats] = ttest(data(:,3),data(:,4))

figure
plot(zscore(data-repmat(mean(data,2),[1,4])))

%% Hit rate x drug

drug_hr_error = std(hit_rate_drug,1)/sqrt(size(hit_rate_drug,1));
figure
h = errorbar(mean(hit_rate_drug,1),drug_hr_error)
h.LineWidth = 2;
set(gca,'FontSize',16,'XTick', [1:length(drug_conds)], 'XTickLabel',drug_conds,'XLim',[0.5 length(drug_conds)+0.5])
xlabel('Drug'); ylabel('Hit Rate (%)')
title('Mean Hit Rate x Drug')

data = hit_rate_drug; % subj x drug
tab = table(data(:,1),data(:,2),data(:,3),data(:,4),'VariableNames',{'MPH','ATM','CIT','PLA'});

mspec = ['MPH,ATM,CIT,PLA~1'];
within = table([1:4]','VariableNames',{'DrugID'});
rm = fitrm(tab,mspec,'WithinDesign',within);
ranovatbl = ranova(rm)

disp('MPH vs PLA')
[h,p,ci,stats] = ttest(hit_rate_drug(:,1),hit_rate_drug(:,4))
% disp('MPH vs PLA (non-parametric)')
% [p,h,stats] = signrank(hit_rate_drug(:,1),hit_rate_drug(:,4))
disp('ATM vs PLA')
[h,p,ci,stats] = ttest(hit_rate_drug(:,2),hit_rate_drug(:,4))
disp('CIT vs PLA')
[h,p,ci,stats] = ttest(hit_rate_drug(:,3),hit_rate_drug(:,4))

mean_hit_rate = mean(mean(hit_rate_drug,2))
sd_hit_rate = std(mean(hit_rate_drug,2))

%% False alarm rate x drug

drug_fa_error = std(fa_rate_drug,1)/sqrt(size(fa_rate_drug,1));
figure
h = errorbar(mean(fa_rate_drug,1),drug_fa_error)
h.LineWidth = 2;
set(gca,'FontSize',16,'XTick', [1:length(drug_conds)], 'XTickLabel',drug_conds,'XLim',[0.5 length(drug_conds)+0.5])
xlabel('Drug'); ylabel('False Alarm Rate (%)')
title('Mean False Alarm Rate x Drug')

%% ******************** ERP Plots ***************************
%% ERP Plots - plottopos

ERP_group = squeeze(mean(avERP(:,:,:,:,2),1)); % chan x time x drug
figure
plottopo(ERP_group(:,:,:),'chanlocs',chanlocs,'limits',[t(1) t(end) ...
    min(min(min(ERP_group(plot_chans,:,:))))  max(max(max(ERP_group(plot_chans,:,:))))], ...
    'title',['Group ERP x Drug'],'legend',drug_conds,'showleg','on','ydir',1)

clear time_windows
time_windows(1,:) = [0:20:260];
time_windows(2,:) = time_windows(1,:)+20;

figure
for t_temp = 1:size(time_windows,2)
    plot_mean = squeeze(mean(mean(ERP_group(:,find(t>time_windows(1,t_temp) & t<time_windows(2,t_temp)),:),2),3));
    subplot(5,3,t_temp)
%     topoplot(plot_mean,chanlocs,'maplimits', ...
%         [min(min(mean(ERP_group(:,find(t>time_windows(1,1) & t<time_windows(2,end)),:),3))) ...
%         max(max(mean(ERP_group(:,find(t>time_windows(1,1) & t<time_windows(2,end)),:),3)))],'electrodes','off','plotchans',plot_chans);
        topoplot(plot_mean,chanlocs,'maplimits', ...
            [min(plot_mean)  max(plot_mean)],'electrodes','off','plotchans',plot_chans);
    title([num2str(time_windows(1,t_temp)),' ms to ',num2str(time_windows(2,t_temp)),' ms']);
    colorbar
end

t1 = 200; t2 = 220;
plot_mean = squeeze(mean(mean(ERP_group(:,find(t>=t1 & t<=t2),:),2),3));
figure
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean) max(plot_mean)], ...
    'electrodes','numbers','plotchans',plot_chans);

ERPr_group = squeeze(mean(avERPr,1)); % chan x time x drug
figure
plottopo(ERPr_group(:,:,:),'chanlocs',chanlocs,'limits',[tr(1) tr(end) ...
    min(min(min(ERPr_group(plot_chans,:,:))))  max(max(max(ERPr_group(plot_chans,:,:))))], ...
    'title',['Group ERPr x Drug'],'legend',drug_conds,'showleg','on','ydir',1)

ERP_sub_group = squeeze(mean(avERP_sub,1)); % chan x time x drug
figure
plottopo(ERP_sub_group(:,:,:),'chanlocs',chanlocs,'limits',[t(1) t(end) ...
    min(min(min(ERP_sub_group(plot_chans,:,:))))  max(max(max(ERP_sub_group(plot_chans,:,:))))], ...
    'title',['Group ERP Sub x Drug'],'legend',drug_conds,'showleg','on','ydir',1)

t1 = -50; t2 = 0;
plot_mean = squeeze(mean(mean(ERPr_group(:,find(tr>=t1 & tr<=t2),:),2),3));
figure
topoplot(plot_mean,chanlocs,'maplimits', ...
    [min(plot_mean) max(plot_mean)], ...
    'electrodes','numbers','plotchans',plot_chans);
% return

%% Stim-locked CPP x subj
% CPPs_stim(s,:,drug) = squeeze(mean(avERP(s,CPP_choice_chans(indx(1:howmanychans)),:,drug,2),2));
% CPPs_subs(s,:,drug) = squeeze(mean(avERP_sub(s,CPP_choice_chans(indx(1:howmanychans)),:,drug),2));

% CPPs = CPPs_stim;
CPPs = CPPs_subs;
CPPs_mean = squeeze(mean(CPPs,1)); % time x drug
figure
clear h
for drug = 1:size(CPPs_mean,2)
    h(drug) = plot(t,CPPs_mean(:,drug),'LineWidth',2); hold on
end
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','-');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
set(gca,'FontSize',12,'xlim',[t(1),t(end)]), xlabel('Time')
legend(h,drug_conds, ...
    'FontSize',12,'Location','NorthWest');
title('CPP x Drug (Subj chan)')

% get peak latencies
window1 = 37; % half the window in ms.
p3_peak=[];
figure
for drug = 1:4
    for s = 1:length(allsubj)
%         [~,indx] = max(squeeze(CPPs(s,find(t>=100),drug)));
%         p3_peak(s,drug) = t(indx+length(find(t<100)));
        
        CPP_temp = squeeze(CPPs(s,:,drug));
        CPP_smooth=[]; t_smooth=[];
        for tt = 1:length(CPP_temp)-window1
            CPP_smooth(tt) = squeeze(mean(CPP_temp(tt:tt+window1)));
            t_smooth(tt) = mean(t(tt:tt+window1));
        end
        [~,indx] = max(CPP_smooth);
        p3_peak(s,drug) = t_smooth(indx);
    end
    subplot(2,2,drug)
    histogram(p3_peak(:,drug),10)
end

cmap = colormap('lines');
% figure, clear h
% for s = 1:length(allsubj)
%     subplot(5,7,s)
%     for drug = 1:size(CPPs_subs,3)
%         h(drug) = plot(t,squeeze(CPPs_subs(s,:,drug)),'Color',cmap(drug,:),'LineWidth',1); hold on
%         line([p3_peak(s,drug),p3_peak(s,drug)],ylim,'Color',cmap(drug,:),'LineWidth',1,'LineStyle','-');
%     end
%     set(gca,'xlim',[t(1),t(end)],'ylim',[min(min(squeeze(CPPs_subs(s,:,:)))),max(max(squeeze(CPPs_subs(s,:,:))))])
%     line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','-');
%     line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
%     title(allsubj{s})
% end
    
% CPPs = squeeze(mean(CPPs_subs,1)); % time x drug
figure
clear h
for drug = 1:size(CPPs_mean,2)
    h(drug) = plot(t,CPPs_mean(:,drug),'LineWidth',2); hold on
    line([mean(p3_peak(:,drug),1),mean(p3_peak(:,drug),1)],ylim,'Color',cmap(drug,:),'LineWidth',1,'LineStyle','-');
end
line([0,0],ylim,'Color','k','LineWidth',1.5,'LineStyle','-');
line(xlim,[0,0],'Color','k','LineWidth',1.5,'LineStyle','-');
set(gca,'FontSize',12,'xlim',[t(1),t(end)]), ylabel('P3b amplitude(\muV/m^2)'), xlabel('Time (ms)')
legend(h,drug_conds, ...
    'FontSize',12,'Location','NorthWest');
title('Stimulus-locked P3b x Drug')

% t1 = 360-25; t2 = 360+25;
% t1 = 150-50; t2 = 150+50;
% data = squeeze(mean(CPPs_subs(:,find(t>=t1 & t<=t2),:),2)); % subj x drug
data = p3_peak;
tab = table(data(:,1),data(:,2),data(:,3),data(:,4),'VariableNames',{'MPH','ATM','CIT','PLA'});

mspec = ['MPH,ATM,CIT,PLA~1'];
within = table([1:4]','VariableNames',{'DrugID'});
rm = fitrm(tab,mspec,'WithinDesign',within);
ranovatbl = ranova(rm)

disp('MPH vs PLA')
[h,p,ci,stats] = ttest(data(:,1),data(:,4))
disp('ATM vs PLA')
[h,p,ci,stats] = ttest(data(:,2),data(:,4))

% figure, plot(zscore(p3_peak-repmat(mean(p3_peak,2),[1,4])))

% save([save_path,'stim_CPP_data_20Hz.mat'],'CPPs_subs','CPPs_stim','drug_conds','t','window1','p3_peak','avERP','avERP_sub')
