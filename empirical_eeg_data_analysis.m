clear all

addpath(genpath('/Users/rrid9830/Downloads/werk/functions/misc'))
addpath(genpath('/Users/rrid9830/Downloads/werk/functions/analysis'))
addpath('/Users/rrid9830/Library/CloudStorage/OneDrive-TheUniversityofQueensland/MATLAB/toolboxes/eeglab2021.1')
addpath(genpath('/Users/rrid9830/Library/CloudStorage/OneDrive-TheUniversityofQueensland/MATLAB/toolboxes/eeglab2021.1/functions'))

data_folder = 'data/'; % data folder (source)
results_folder = 'results/'; % results folder (destination)
mkdir(results_folder)

trials = 1:3600;  % trials to use (vector integer)
rear_electrodes = [20:31,57:64]; % electrodes to use (vector integer)
num_ori_chans = 6; % number of basis set channels (scalar integer)
num_basis_set_points = 180; % number of basis set points (scalar integer)
num_folds = 10; % number of cross-validation folds (scalar integer)
temporal_window = [-50,450]; % temporal window of analysis around stimulus onset (msec; vector, integer)


% identify data files
files = dir([pwd,filesep,data_folder,'eeg_recordings*.mat']);

for participant = 1:numel(files)
    tic % record analysis duration
    
    % display analysis progress
    fprintf('\nAnalyzing file: %s\n',files(participant).name)
    
    % load data
    load([pwd,filesep,data_folder,files(participant).name])
    Y = eeg.data(rear_electrodes,eeg.times>temporal_window(1) & eeg.times<temporal_window(2),trials); % neural recordings
    X = stimulus.orientations(trials); % stimulus labels
    times = eeg.times(eeg.times>temporal_window(1) & eeg.times<temporal_window(2));
    num_times = size(Y,2); % number of times
    num_trials = size(Y,3); % number of trials
    
    % bin orientations
    chans = linspace(0,pi,num_ori_chans+1);
    chans(end) = [];
    chan_width = round(pi/num_ori_chans,4);
    chans = chans + chan_width/2;
    labels = X*nan;
    for b_idx = 1:num_ori_chans
        labels(round(abs(circ_dist(X*2,chans(b_idx)*2)),4)<=chan_width) = round(chans(b_idx)*180/pi);
    end
    chans = unique(labels);
    
    % define basis set function [rectified cosine]
    funType = @(xx,mu) (cosd(xx-mu)).^(num_ori_chans-mod(num_ori_chans,2));
    
    channel_responses = nan(num_basis_set_points,num_times,num_trials,'single'); % assign placeholder
    for stp_idx = 1:180/num_ori_chans
        
        % generate basis set
        xx = circshift(linspace(1,180,num_basis_set_points),stp_idx-1); % basis set orientations
        basis_set = nan(num_basis_set_points,num_ori_chans);        
        for chan_idx = 1:num_ori_chans
            basis_set(:,chan_idx) = funType(xx,chans(chan_idx));
        end
        
        % generate cross-validation folds
        folds = cell(num_folds,2);
        for chan_idx = 1:num_ori_chans
            % find indices
            index = find(labels(trials) == chans(chan_idx));
            n_index = length(index);
            
            % shuffle
            index = index(randperm(n_index));
            
            % distribute across folds
            group_number = floor((0:(n_index-1))*(num_folds/n_index))+1;
            for fold_idx = 1:num_folds
                folds{fold_idx} = [folds{fold_idx}, index(group_number==fold_idx)'];
            end
        end
        % store fold order
        [~,order] = sort([folds{:}]);
        
        % generate stimulus mask
        stim_mask = zeros(num_trials,length(xx));
        for trial_idx = 1:num_trials
            stim_mask(trial_idx,labels(trial_idx)) = 1;
        end
        
        % generate design matrix
        design = (stim_mask*basis_set)';
        
        part_channel_responses = nan(num_ori_chans,num_times,num_trials,'single'); % assign placeholder
        parfor time_idx = 1:num_times
            temp_channel_responses = [];
            y = squeeze(Y(:,time_idx,:));
            for fold_idx = 1:num_folds
                % define train/test trials
                test_trials = folds{fold_idx};
                train_trials = find(~ismember(trials,test_trials));
                
                % train decoder [forward model]
                decoder = train_beamformerMT([],design(:,train_trials),y(:,train_trials));
                
                % test decoder [inverted model]
                temp_channel_responses = cat(2,temp_channel_responses,decode_beamformer([], decoder,y(:,test_trials)));
            end % for fold_idx = 1:num_folds
            
            part_channel_responses(:,time_idx,:) = temp_channel_responses(:,order);
        end % parfor time_idx = 1:num_times
        
        channel_responses([0:num_ori_chans-1]*180/num_ori_chans+stp_idx,:,:) = part_channel_responses;
    end % for stp_idx = 1:36
    
    % save results
    save([results_folder,files(participant).name(1:end-4),'-results','.mat'],'channel_responses','X','chans','labels','times')
    
    toc % display analysis duration
end % for participant = 1:numel(files)

