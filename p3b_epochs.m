% Ger Loughnane 10/10/17
% Run one subject through EEG preprocessing stage, called from
% runafew_drugs

%% Trigger definition
std_targcodes = [50];
dev_targcodes = [75];
resp_start_targcode = [100];
resp_targcode = [1];
targcodes = [std_targcodes,dev_targcodes]; % for purposes of time-locking

fs=512; % sample rate
nchan = 64; % number of chans

%% Epoch definition
% for i = 1:length(t)
%     integerTest(i)=~mod(t(i),1);
% end
% theseareintegers = t(find(integerTest==1));
% [-1000,-875,-750,-625,-500,-375,-250,-125,0,125,250,375,500,625,750,875,1000,1125,1250,1375,1500,1625,1750,1875,2000]
% in sample points, the ERP epoch
ts = -0.125*fs:0.750*fs; % -1000ms to 1880ms with 200ms either side.
t = ts*1000/fs;
BLint = [-125 0];   % baseline interval in ms
default_response_time = 450;
ERP_samps = length(t);

%% CSD
% E = textread('C:\Program Files\MATLAB\R2016a\toolbox\CSDtoolbox\chans64.asc','%s');
% M = ExtractMontage('C:\Program Files\MATLAB\R2016a\toolbox\CSDtoolbox\10-5-System_Mastoids_EGI129.csd',E);  % reading in the montage for the CSD toolbox
% [G2,H2] = GetGH(M);

%% Artifact rejection
% frontal channels, occipital channels, and POz and Pz, CPz, CP1, CP2, P1, P2
% ARchans = [1:3,6,7,34:36,41,42, ...
%     21:23,25:27, ...
%     58:60,62:64, ...
%     30,31,32];
% ARchans = [1:3,6,7,34:36,41,42, ...
%     16,17,21:23,25:27, ...
%     53,54,58:60,62:64, ...
%     31];
ARchan_all = [1:64];
ARchan_choice{1} = ARchan_all;

ARchans_for_blinks = [1:3,6,7,34:36,41,42];
ARchans_for_CPP = [30,31,32,20,57,19,56];
ARchans_for_FC = [4,38,39,11,47,46,12,49,48];
ARchans_for_VEP = [16,17,20,21,22,23,25,26,27,29,28,53,54,57,58,59,60,62,63,64];
ARchans_for_M1 = [12,13,14,49,50,51];

ARchan_choice{2} = [ARchans_for_blinks,ARchans_for_CPP,ARchans_for_FC,ARchans_for_VEP,ARchans_for_M1];
ARchan_choice{3} = [ARchans_for_blinks,ARchans_for_CPP];
ARchan_choice{4} = [ARchans_for_blinks,ARchans_for_FC];
ARchan_choice{5} = [ARchans_for_blinks,ARchans_for_VEP];
ARchan_choice{6} = [ARchans_for_blinks,ARchans_for_M1];

% ARchans = [1:23,25:27,29:32,34:36,38:46,48:60,62:64];
artifth = [25,50,100];
artifchans = cell(4,6,3);  % keep track of channels on which the threshold is exceeded, causing trial rejection
artif_trials = [];
hz_tags = [100,35,20,8];

%% Define chanlocs and set up ERP files
chanlocs = readlocs('cap64.loc');
chanlocs = chanlocs(1:nchan)';
  
% figure
% topoplot(zeros(1,64),chanlocs,'maplimits', ...
%     [min(zeros(1,64)) max(zeros(1,64))], ...
%     'electrodes','numbers','plotchans',ARchan_choice{2});

erp_100Hz = []; erp_35Hz = []; erp_20Hz = []; erp_8Hz = []; erp_100Hz_CSD = []; erp_35Hz_CSD = []; erp_20Hz_CSD = []; erp_8Hz_CSD = [];
numtr = 0; % total_numtr = 0;
allRT = []; allResp = []; allTrig = []; allBlock = [];  % note allRT will be in sample points

%% Go through files
for f=1:length(files)
    if isempty(files{f}), continue, end
    disp(f)
    filename=[files{f}];
    EEG = pop_biosig(filename, 'blockepoch', 'off','channels',[1:nchan]);
    EEG = pop_resample(EEG,fs);
    
    numev = length(EEG.event);

    % Fish out the event triggers and times
    clear trigs stimes
    for i=1:numev
        trigs(i)=EEG.event(i).type;
        stimes(i)=round(EEG.event(i).latency);
    end
    
    % Check mains frequencies
%     for i = 1:64
%         ep = squeeze(EEG.data(i,:)); % time
%         nfft = size(ep,2);
%         fftx = abs(fft(ep,[],2))./(nfft/2);
%         fftx = fftx(:,1:ceil((nfft+1)/2));
%         freq_temp = (0:ceil((nfft+1)/2)-1)*fs/nfft;
%         fftx_all(i,:) = fftx;
%     end
%     figure, plot(freq_temp,mean(fftx_all,1))
%     keyboard
    
    % There are peaks at 50,60
    % Get rid of mains frequency
    EEG = pop_eegfiltnew(EEG,49,51,[],1,0,0,0);
    EEG = pop_eegfiltnew(EEG,59,61,[],1,0,0,0);
    % HP Filter
    EEG = pop_eegfiltnew(EEG,0.25,0,[]); % filter to 0.25Hz
    % LP Filter
    EEG_100 = pop_eegfiltnew(EEG,0,100,[]); % 100Hz low pass filter
    EEG_35 = pop_eegfiltnew(EEG,0,35,[]);
    EEG_20 = pop_eegfiltnew(EEG,0,20,[]);
    EEG_8 = pop_eegfiltnew(EEG,0,8,[]);
    clear EEG % don't need the original now
    
    % interpolate bad channels
    if ~isempty(badchans{f})
        EEG_100.chanlocs = chanlocs;
        EEG_100 = eeg_interp(EEG_100,[badchans{f}],'spherical');
        EEG_35.chanlocs = chanlocs;
        EEG_35 = eeg_interp(EEG_35,[badchans{f}],'spherical');
        EEG_20.chanlocs = chanlocs;
        EEG_20 = eeg_interp(EEG_20,[badchans{f}],'spherical');
        EEG_8.chanlocs = chanlocs;
        EEG_8 = eeg_interp(EEG_8,[badchans{f}],'spherical');
    end

    % average-reference the whole continuous data (safe to do this now after interpolation):
    EEG_100.data = EEG_100.data - repmat(mean(EEG_100.data([1:nchan],:),1),[nchan,1]);
    EEG_35.data = EEG_35.data - repmat(mean(EEG_35.data([1:nchan],:),1),[nchan,1]);
    EEG_20.data = EEG_20.data - repmat(mean(EEG_20.data([1:nchan],:),1),[nchan,1]);
    EEG_8.data = EEG_8.data - repmat(mean(EEG_8.data([1:nchan],:),1),[nchan,1]);

%     figure, plot(trigs)
%     keyboard
    
    stim_on = [];
    for n=1:length(trigs)
        if any(targcodes(:)==trigs(n))
            stim_on = [stim_on n];
        end
    end
    
    for n=1:length(stim_on)
        numtr = numtr+1;
        locktime = stimes(stim_on(n));
        
        if trigs(stim_on(n))==75 & n<=length(stim_on)-2
            if trigs(stim_on(n)+2)==1
                response_time = floor((stimes(stim_on(n)+2)-locktime)*1000/fs); % time in ms from beginning of motion to response.
            else
                response_time = default_response_time;
            end
        else
            response_time = default_response_time;
        end
        
        try
            ep_100Hz = EEG_100.data(1:nchan,locktime+ts);   % chop out an epoch
            ep_35Hz = EEG_35.data(1:nchan,locktime+ts);
            ep_20Hz = EEG_20.data(1:nchan,locktime+ts);
            ep_8Hz = EEG_8.data(1:nchan,locktime+ts);
        catch
            disp('EEG ended too soon')
            allTrig(numtr) = NaN;
            allResp(numtr) = 5;
            allRT(numtr) = NaN;
            allBlock(numtr) = p3b_drug_order(s,f);
            
            erp_100Hz(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_35Hz(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_20Hz(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_8Hz(:,:,numtr) = zeros(nchan,ERP_samps);
            
            erp_100Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_35Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_20Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            erp_8Hz_CSD(:,:,numtr) = zeros(nchan,ERP_samps);
            
            continue;
        end

        ep_test = [find(t<response_time)];
        if isempty(ep_test)
            disp('Empty epoch for art rejection')
            keyboard
        end
        
        for hz = 1:4
            ep_temp = eval(['ep_',num2str(hz_tags(hz)),'Hz;']);
            for archan = 1:6
                ARchans = ARchan_choice{archan};
                for art_thresh = 1:3
                    artifchans_thistrial = ARchans(find(max(abs(ep_temp(ARchans,find(t<response_time))),[],2)>artifth(art_thresh) | ...
                        max(abs(ep_temp(ARchans,find(t<default_response_time))),[],2)>artifth(art_thresh)));
                    artifchans{hz,archan,art_thresh} = [artifchans{hz,archan,art_thresh} artifchans_thistrial];
                    if length(artifchans_thistrial) > 0
                        artif_trials(hz,archan,art_thresh,numtr) = 0; % artifact trial
                    else
                        artif_trials(hz,archan,art_thresh,numtr) = 1; % non-artifact trial
                    end   % artifact rejection (threshold test)
                end
            end
        end
        
        erp_100Hz(:,:,numtr) = ep_100Hz;
        erp_35Hz(:,:,numtr) = ep_35Hz;
        erp_20Hz(:,:,numtr) = ep_20Hz;
        erp_8Hz(:,:,numtr) = ep_8Hz;
        
        % CSD
        ep_100Hz_CSD = CSD(ep_100Hz,G2,H2);
        ep_35Hz_CSD = CSD(ep_35Hz,G2,H2);
        ep_20Hz_CSD = CSD(ep_20Hz,G2,H2);
        ep_8Hz_CSD = CSD(ep_8Hz,G2,H2);
        
        erp_100Hz_CSD(:,:,numtr) = ep_100Hz_CSD;
        erp_35Hz_CSD(:,:,numtr) = ep_35Hz_CSD;
        erp_20Hz_CSD(:,:,numtr) = ep_20Hz_CSD;
        erp_8Hz_CSD(:,:,numtr) = ep_8Hz_CSD;

        allTrig(numtr) = trigs(stim_on(n));
        allBlock(numtr) = p3b_drug_order(s,f);
        
        if n<=length(stim_on)-2
            if trigs(stim_on(n))==75
                if trigs(stim_on(n)+2)==1
                    allResp(numtr) = 1; % hit
                    allRT(numtr) = stimes(stim_on(n)+2)-stimes(stim_on(n));
                else
                    allResp(numtr) = 2; % miss
                    allRT(numtr) = NaN;
                end
            elseif trigs(stim_on(n))==50
                if trigs(stim_on(n)+2)==1
                    allResp(numtr) = 3; % false alarm
                    allRT(numtr) = stimes(stim_on(n)+2)-stimes(stim_on(n));
                else
                    allResp(numtr) = 4; % correct rejection
                    allRT(numtr) = NaN;
                end
            end
        else
            allResp(numtr) = 5; % EEG ended too soon
            allRT(numtr) = NaN;
        end
    end
end
%% Save files
% (hz,archan,art_thresh,numtr)
close all
rejected_trials = length(find(artif_trials(2,1,3,:)==0));
figure
subplot(1,2,1)
hist(artifchans{2,1,3},[1:nchan]); title([allsubj{s} ': ' num2str(rejected_trials) ' artifacts = ',num2str(round(100*(rejected_trials/length(allRT)))),'%']) % s from runafew
disp([allsubj{s},' number of trials: ',num2str(length(find(allRT)))])

[counts,centers] = hist(artifchans{2,1,3},[1:nchan]);
subplot(1,2,2)
topoplot(counts,chanlocs,'plotchans',[1:nchan],'electrodes','numbers');
title(subject_folder{s})
saveas(gcf,[path_temp subject_folder{s} '\' allsubj{s} '_p3b_artifact.jpg'],'jpg')
pause(1)

save([path_temp subject_folder{s} '\' allsubj{s} '_p3b_epochs'], ...
    'erp_100Hz','erp_35Hz','erp_20Hz','erp_8Hz','erp_100Hz_CSD','erp_35Hz_CSD','erp_20Hz_CSD','erp_8Hz_CSD', ...
    'allRT','allResp','allTrig','allBlock','t', ...
    'artifchans','artif_trials');
return;