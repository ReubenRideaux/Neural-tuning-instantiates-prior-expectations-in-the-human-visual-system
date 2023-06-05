clear all

% define functions
vmfun = @(xx, kappa, mu) exp(kappa*cosd((xx-mu))); % von Mises function
csfun = @(xx, mu, expo) (cosd(xx-mu)).^expo; % rectified cosine function
dgfun = @(xx, kappa, mu)    diff(exp( kappa*(cosd(xx-mu)-1) )); % differential von Mises function

% define paramters
params.meta.n_subs = 36; % # of particiapnts to simulate (scalar, integer)
params.tcurves.n_ori_prefs  = 16; % # of orientation preferences (scalar, integer)
params.tcurves.n_oris = 180; % # of orientations to test (deg; scalar, integer)
params.tcurves.oris = linspace(1,180,params.tcurves.n_oris); % (deg)
params.tcurves.ori_prefs = linspace(0,180,params.tcurves.n_ori_prefs+1); % (deg)
params.tcurves.ori_prefs(end) = [];
params.tcurves.kappa = 2; % neural tuning kappa value (deg)

params.cardi.oris = [0,90]; % position of cardinal orientations (deg)
params.cardi.sft.kappa = .5; % bias tuning kappa value (deg)
params.cardi.sft.amp = [14,8]; % amplitude of bias at each cardinal (a.u.)

params.neural.n_sensors = 20; % # of EEG sensors to simulate (scalar, integer)
params.neural.sensornoise_sd = 6; % amplitude of sensor noise (a.u.; scalar, float)

params.stim.n_stimuli = 180; % # of different oriented grating to simulate (scalar, integer)
params.stim.n_repeats = 18*2; % # of repeats per grating (scalar, integer)
params.stim.n_trials = params.stim.n_stimuli * params.stim.n_repeats; % total # of trials to simulate
params.stim.stim_lib = repmat(1:params.stim.n_stimuli,1,params.stim.n_repeats);

params.model.n_ori_chans = 6; % # of forward model channels (scalar, integer)
params.model.kernel = exp(1i * (params.tcurves.oris*2*pi/180)); % circular kernel

params.analysis.smoothing = 16; % smoothing window

% simulate neuron tuning curves
upsample_factor = 10000;
upsampled_oris = linspace(0,180,upsample_factor+1);
params.cardi.sft.mod  = dgfun(upsampled_oris*2,params.cardi.sft.kappa,params.cardi.oris(1)*2)*params.cardi.sft.amp(1)+...
    dgfun(upsampled_oris*2,params.cardi.sft.kappa,params.cardi.oris(2)*2)*params.cardi.sft.amp(2);
params.cardi.sft.mod = params.cardi.sft.mod(round(linspace(1,upsample_factor,params.tcurves.n_ori_prefs+1)));
params.cardi.sft.mod(end) = [];

params.cardi.sft.mod = params.cardi.sft.mod/max(abs(params.cardi.sft.mod))*max(abs(params.cardi.sft.amp));
for cc = 1:params.tcurves.n_ori_prefs
    tcurves(:,cc) = vmfun(params.tcurves.oris*2,params.tcurves.kappa,params.tcurves.ori_prefs(cc)*2+params.cardi.sft.mod(cc));
    tcurves(:,cc) = tcurves(:,cc)/max(tcurves(:,cc));
end

% plot tuning curves
figure
subplot(4,4,1)
plot(params.tcurves.oris,tcurves)
ylim([0,1]);
xlim([min(params.tcurves.oris),max(params.tcurves.oris)])

%% begin simulation
warning OFF
for sub = 1:params.meta.n_subs % repeat over multiple subjects
    fprintf('simulating participant %i of %i\n',sub, params.meta.n_subs);
    
    % make a random weighting of neurons on to each sensor
    neuron_to_sensor_weights = rand(params.tcurves.n_ori_prefs,params.neural.n_sensors);
    
    stimuli = Shuffle(params.stim.stim_lib);
    
    sensornoise = randn(params.stim.n_trials,params.neural.n_sensors)*params.neural.sensornoise_sd;
    
    % compute the sensor responses
    sensor_response = (tcurves(stimuli,:))*neuron_to_sensor_weights + sensornoise;
    
    % bin orientations
    phi = stimuli*pi/180;
    chans = linspace(0,pi,params.model.n_ori_chans+1);
    chans(end) = [];
    chan_width = round(pi/params.model.n_ori_chans,4);
    chans = chans + chan_width/2;
    labels = phi*nan;
    for b = 1:params.model.n_ori_chans
        labels(round(abs(circ_dist(phi*2,chans(b)*2)),4)<=chan_width) = round(chans(b)*180/pi);
    end
    chans = unique(labels);
    
    % apply forward encoding to decode simulated data
    train_trials = 1:params.stim.n_trials/2;
    test_trials = params.stim.n_trials/2+1:params.stim.n_trials;
    
    % create stimulus mask
    stim_mask = zeros(params.stim.n_trials,params.tcurves.n_oris);
    for trial = 1:params.stim.n_trials
        stim_mask(trial,labels(trial)) = 1;
    end
    
    for stp = 1:180/params.model.n_ori_chans
        shifted_oris = circshift(params.tcurves.oris,stp-1);
        basis_set = nan(params.tcurves.n_oris,params.model.n_ori_chans);
        
        for cc = 1:params.model.n_ori_chans
            basis_set(:,cc) = csfun(shifted_oris,chans(cc),params.model.n_ori_chans-mod(params.model.n_ori_chans,2));
        end
        
        % compute design matrix
        design = (stim_mask*basis_set);
                    
        decoder = train_beamformerMT([],design(train_trials,:)',sensor_response(train_trials,:)');
        estimatedChannelResponse([0:params.model.n_ori_chans-1]*180/params.model.n_ori_chans+stp,:) = ...
            decode_beamformer([], decoder,sensor_response(test_trials,:)'); 
    end
        
    degs = round(phi(test_trials)*180/pi);
    estimatedChannelResponse = circshift(estimatedChannelResponse,chans(1),1);
    Z = params.model.kernel*estimatedChannelResponse;
    theta = mod(angle(Z), 2*pi) * (180/pi) * 0.5;       % Decoded orientation
    
    for test_o = 1:params.stim.n_stimuli
        r{sub}(test_o) = mean(exp(1i * (theta(:,degs==test_o) - test_o) * (pi/180) * 2),2);
        sd{sub}(test_o) = circ_std(theta(:,degs==test_o) * (pi/180) * 2,[],[],2)/2;
    end
end
warning ON

r = reshape(cell2mat(r),params.stim.n_stimuli,params.meta.n_subs);
sd = reshape(cell2mat(sd),params.stim.n_stimuli,params.meta.n_subs);

% calculate accuracy, precision, and bias
accu = abs(r) .* cos(angle(r)); % (a.u.)
prec = 1./(sd*180/pi); % (deg)
bias = angle(r); % (2pi radians)

% average across subjects
accu = mean(accu,2); % (a.u.)
prec = mean(prec,2); % (deg)
bias = circ_mean(bias,[],2)*180/pi/2; % (deg)

% plot decoding results
x = 1:180;

subplot(4,4,5)
plot(x,accu,'.k')
hold on
sy{1} = smooth([accu,accu,accu],params.analysis.smoothing,'moving');sy{1} = sy{1}(181:360);
plot(x,sy{1},'-r','lineWidth',2)
ylim([min(accu(:)),max(accu(:))]);
xlim([0,180]);

subplot(4,4,6)
plot(x,prec,'.k')
hold on
sy{2} = smooth([prec,prec,prec],params.analysis.smoothing,'moving');sy{2} = sy{2}(181:360);
plot(x,sy{2},'-r','lineWidth',2)
ylim([min(prec(:)),max(prec(:))]);
xlim([0,180]);

subplot(4,4,7)
plot(x,bias,'.k')
hold on
sy{3} = smooth([bias,bias,bias],params.analysis.smoothing,'moving');sy{3} = sy{3}(181:360);
plot(x,sy{3},'-r','lineWidth',2)
ylim([min(bias(:)),max(bias(:))]);
xlim([0,180]);
