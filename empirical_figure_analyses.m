clear all

results_folder = 'results/'; % results folder (destination)

% identify data files
files = dir([pwd,filesep,results_folder,'*-results.mat']);
num_participants = numel(files);

% collate data from across participants
for participant = 1:num_participants
    
    % display analysis progress
    fprintf('Retrieving results from: %s\n',files(participant).name)
    
    % load data
    load([pwd,filesep,results_folder,files(participant).name])
    X = ceil(X*180/pi);
    num_basis_set_points = size(channel_responses,1);
    num_trials = size(channel_responses,2);
    num_chans = numel(chans);
   
    % define circular kernel
    kernel = exp(1i * (chans'*2*pi/180)); % low resolution kernel

    % decode orientations
    Z = reshape(kernel*reshape(channel_responses([1:num_basis_set_points/num_chans:num_basis_set_points],:,:), num_chans, []), num_trials, []);
    theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;
    
    % compute mean and variance of decoded orientations
    for ori_idx = 1:num_chans
        lr.r(:,ori_idx,participant) = mean(exp(1i * (theta(:,labels==chans(ori_idx)) - chans(ori_idx)) * (pi/180) * 2),2);
        lr.sd(:,ori_idx,participant) = circ_std(theta(:,labels==chans(ori_idx)) * (pi/180) * 2,[],[],2)/2;
    end
    
    % define circular kernel
    kernel = exp(1i * ([1:num_basis_set_points]*2*pi/180)); % high resolution kernel
    
    % decode orientations
    channel_responses = circshift(channel_responses,chans(1),1);
    Z = reshape(kernel*reshape(channel_responses, num_basis_set_points, []), num_trials, []);
    theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;
    
    % compute mean and variance of decoded orientations
    for ori_idx = 1:num_basis_set_points
        hr.r(:,ori_idx,participant) = mean(exp(1i * (theta(:,X==ori_idx) - ori_idx) * (pi/180) * 2),2);
        hr.sd(:,ori_idx,participant) = circ_std(theta(:,X==ori_idx) * (pi/180) * 2,[],[],2)/2;
    end
end

%% ACCURACY, PRECSION, AND BIAS AS A FUNCTION OF TIME (Figures 3a-c)
rr = abs(lr.r) .* cos(angle(lr.r));
yl = [-.015,.12]; bl = 0;

figure
color = [67/255,78/255,94/255];
subplot(4,4,1)
shadedErrorBar(times, squeeze(mean(rr,[2,3])), squeeze(std(mean(rr,2),[],3)/sqrt(num_participants-1)),color,1);
xlabel('time from onset (msec)')
ylabel('orientation response (a.u.)')
ylim(yl)
xlim([min(times),max(times)])
yl = ylim;xl = xlim;
line(xl,[bl,bl],'LineStyle','--','color','k')

rr = 1 - (lr.sd/(sqrt(2)/2));
bl = 0;
yl = [bl,.08];

subplot(4,4,3)
shadedErrorBar(times, squeeze(mean(rr,[2,3])), squeeze(std(mean(rr,2),[],3)/sqrt(num_participants-1)),color,1);
xlabel('time from onset (msec)')
ylabel('orientation precision (1/degs)')
ylim(yl)
xlim([min(times),max(times)])
yl = ylim;xl = xlim;
line(xl,[bl,bl],'LineStyle','--','color','k')

theta = double(mod(angle(lr.r), 2*pi));
yl = [-180,180];
bl = 0;
[ym,ul,ll] = circ_mean(theta,[],[2,3]);
ul =  ym-ul;
ll =  ym-ll;
ul(isnan(ul)) = pi;
ll(isnan(ll)) = -pi;

subplot(4,4,6)
shadedErrorBar(times, ym*180/pi, [ul,-ll]*180/pi,color,1);
xlabel('time from onset (msec)')
ylabel('orientation precision (1/degs)')
ylim(yl)
xlim([min(times),max(times)])
yl = ylim;xl = xlim;
line(xl,[bl,bl],'LineStyle','--','color','k')

%% ACCURACY, PRECSION, AND BIAS AS A FUNCTION OF ORIENTATION (Figures 3d-f)
sr = circshift(hr.r,90,2);
ssd = 1-(circshift(hr.sd,90,2)/(sqrt(2)/2));
rr = double(abs(sr) .* cos(angle(sr)));
rsd = double(ssd);
theta = double(mod(angle(sr), 2*pi));
sf = 16; % smoothing window
x = -179:2:179;

figure

subplot(4,4,1)
y = squeeze(mean(rr,[1,3]));
plot(x/2,y,'.k')
hold on
sy{1} = smooth([y,y,y],sf,'moving');sy{1} = sy{1}(181:360);
plot(x/2,sy{1},'-r','lineWidth',2)
xlim([min(x),max(x)]/2)

subplot(4,4,2)
y = squeeze(mean(rsd,[1,3]));
plot(x/2,y,'.k')
hold on
sy{2} = smooth([y,y,y],sf,'moving');sy{2} = sy{2}(181:360);
plot(x/2,sy{2},'-r','lineWidth',2)
xlim([min(x),max(x)]/2)

subplot(4,4,3)
y = squeeze(circ_mean(theta,[],[1,3])/2)'*180/pi;
plot(x/2,y,'.k')
hold on
sy{3} = smooth([y,y,y],sf,'moving');sy{3} = sy{3}(181:360);
plot(x/2,sy{3},'-r','lineWidth',2)
xlim([min(x),max(x)]/2)

