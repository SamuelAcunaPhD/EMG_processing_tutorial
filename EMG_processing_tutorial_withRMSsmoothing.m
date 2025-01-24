% Filename: EMG_processing_tutorial_withRMSsmoothing.m
% Author:   Samuel Acuña
% Date:     11 May 2018
% Description:
% Guides you through a typical EMG processing routine.
%
% EMG is more about timing of activation, than the amplitude.
% 
% Make sure you have the accompanying data file:
% EMG_processing_tutorial_DATA.mat  or
% EMG_processing_tutorial_DATA2.mat
%
% Directions: run each section at a time to see steps (usually alt+enter, or command+enter)

%% STEP 1: LOAD EMG DATA
clear; close all; clc;
load('EMG_processing_tutorial_DATA.mat');
%load('EMG_processing_tutorial_DATA2.mat');
% DATA: this is data from the medial gastrocnemius when walking on a treadmill
% DATA2: tibialis anterior, same trial

% Loads matlab structure: trial_data
%   trial_data.emg         : raw EMG data
%   trial_data.label       : muscle where EMG data was collected
%   trial_data.freq        : sampling frequency
%   trial_data.time        : time of each sample of data
%   trial_data.heelStrikes : instances of heel strikes between steps



%% STEP 1: BANDPASS FILTER EMG, 
% the bandpass gets rid of any noise that isnt part of EMG data. Most EMG
% power is between 5-500 Hz
% low cutoff: want to be high enough to get rid of motion artifacts and drift. No higher than 10 hz, per ISEK guidelines
% high cutoff: no lower than 350 Hz, per ISEK guidelines
BP = [10 500]; % bandpass filter parameters in Hz [low cutoff, high cutoff]
[b_BP,a_BP]=butter(4,BP/(trial_data.freq/2)); % bandpass filter coefficients. (4th order butterworth)
emg_BP = filtfilt(b_BP,a_BP,trial_data.emg); % zero phase lag filter

figure(1);
title(['EMG: ' trial_data.label]);
xlabel('time (sec)'); xlim([0 60]);
ylabel('EMG (volts)');
h1 = plot(trial_data.time,emg_BP,'color',[0.6875    0.7656    0.8672]);
hold off;

return
%% STEP 2: Use RMS smoothing on data

windowlength = 250;
overlap = 240;  %windowlength - 1;
emg_RMS = rms_gbiomech(emg_BP,windowlength,overlap);


delta = windowlength - overlap;
indices = 1:delta:length(trial_data.time);
time_RMS = trial_data.time(indices);

figure(1); hold on;
h2 = plot(time_RMS,emg_RMS,'b','LineWidth',2);
hold off;


%% STEP 3: compare to non-RMS smoothed data

emg_ABS = abs(emg_BP); % full wave rectified
LP = 10; % low pass filter for linear envelope, in Hz
[b_LP,a_LP]=butter(4,LP/(trial_data.freq/2),'low'); % linear envelope filter
emg_ENV = filtfilt(b_LP,a_LP,emg_ABS);

figure(1); hold on;
h3 = plot(trial_data.time,emg_ENV,'color',[0  0 0],'LineWidth',2);
hold off;
legend([h1 h2 h3],{'bandpass EMG', 'RMS smoothed EMG', 'non-smooth envelope EMG'});


%% STEP 4: NORMALIZE AMPLITUDE
% scaled to the max value (could also do RMS, or max voluntary contraction)
emg_NORM_RMS = emg_RMS/max(emg_RMS);
emg_NORM = emg_ENV/max(emg_ENV);

figure(2); hold on;
h1 = plot(trial_data.time,emg_NORM,'color',[ 0    0    0],'LineWidth',2);
h2 = plot(time_RMS,emg_NORM_RMS,'b','LineWidth',2);
title(['EMG: ' trial_data.label]);
xlabel('time (sec)'); xlim([0 60]);
ylabel('normalized EMG'); ylim([0 1.5]);
legend([h1 h2],'normalized linear envelope of EMG', 'normalized RMS smoothed EMG');

%% STEP 5: TIME NORMALIZE BY A STRIDE

% divides the EMG signals insto strides
npts = 101; % points per gait cycle
nStrides = length(trial_data.heelStrikes)-1; % number of complete gait cycles
emg_Strides = zeros(101,nStrides); %preallocate
emg_Strides_RMS = zeros(101,nStrides); %preallocate
for j = 1:nStrides
    j1 = find(trial_data.time>trial_data.heelStrikes(j),1); % get index of time at first heel strike
    j2 = find(trial_data.time>trial_data.heelStrikes(j+1),1); % get index of time at second heel strike
    emg_Strides(:,j) = normcycle(emg_NORM(j1:j2),npts); % time normalize
    
    k1 = find(time_RMS>trial_data.heelStrikes(j),1);
    k2 = find(time_RMS>trial_data.heelStrikes(j+1),1);
    emg_Strides_RMS(:,j) = normcycle(emg_NORM_RMS(k1:k2),npts); % time normalize
end


% FIND AVERAGE EMG over a stride
emg_AVG = mean(emg_Strides,2); %average EMG of each stride
emg_STD = std(emg_Strides')'; % standard deviation
emg_AVG_RMS = mean(emg_Strides_RMS,2); %average EMG of each stride
emg_STD_RMS = std(emg_Strides_RMS')'; % standard deviation

% PLOT
figure(3); subplot(2,1,1); 
shadedErrorBar([0:100]',emg_AVG,emg_STD,'k',1);
hold on
shadedErrorBar([0:100]',emg_AVG_RMS,emg_STD_RMS,'b',1);
title(['EMG: ' trial_data.label]);
xlabel('Gait Cycle (0-100%)'); xlim([0 100]);
ylabel('normalized EMG'); ylim([0 1.5]);
legend('EMG (AVG ± STD)','RMS Smoothed EMG (AVG ± STD)')


figure(3); subplot(2,1,2); hold on;
for i = 1:nStrides
        plot([0:100]',emg_Strides(:,i),'k');
        plot([0:100]',emg_Strides_RMS(:,i),'b');
end
hold off;
title(['EMG for every stride']);
xlabel('Gait Cycle (0-100%)'); xlim([0 100]);
ylabel('normalized EMG'); ylim([0 1.5]);





%%%%%%%%%%%%%%

function yf = normcycle(y,n,x)
% yf = normcycle(y,n,x)
% Convert a signal y to n even-spaced data points over a cycle
% Often used for presentation of gait data, default for n is 101 points
% can specify an indpendent variable x (optional)
if ~exist('n','var')
    n=101;
end
[nr,nc]=size(y);
if nc==1 && nr>1
    ny=1;
    nx=nr;
elseif nr==1 && nc>1
    y=y';
    ny=1;
    nx=nc;
elseif nr>1 && nc>1
    ny=nc;
    nx=nr;
else
    disp('normcycle does not work on a scalar value');
    yf=[];
    return
end
if ~exist('x','var')
    x=[0:(nx-1)]/(nx-1);
else
    nx=length(x);
    x=(x-x(1))/(x(end)-x(1));
end
kk=[0:(n-1)]/(n-1);
yf=interp1(x,y,kk,'*pchip');

end

function RMS=rms_gbiomech(signal,windowlength,overlap)
% adapted from: http://g-biomech.blogspot.com/2014/07/emg-signal-processing-smoothing-root.html

% i think the current code assumes an overlapp of windowlength-1

if mod(windowlength,2) % if odd
    error('windowlength must be even integer')
end
    
% zeropad the signal
%if length(signal) - indices(end) + 1 < windowlength

%     % add zeros to the beginning and end of the signal
%     newsignal = [zeros(windowlength/2,1); signal; zeros(windowlength/2,1)];
%     signal = newsignal;

%end

delta = windowlength - overlap;
indices = 1:delta:length(signal);
RMS = zeros(length(indices),1);

signal = signal.^2;

index = 0;

for i = indices
    index = index+1;

    if i <= windowlength/2 % for beginning of signal, look forward in time
        RMS(index) = sqrt(mean(signal(1:windowlength)));
    elseif i > indices(end)-windowlength/2   % for end of signal, look backard in time.
        RMS(index) = sqrt(mean(signal(end-windowlength/2:end)));
    else % for all, look forward and backwards in time
        RMS(index) = sqrt(mean(signal(i-windowlength/2:i+windowlength/2)));
    end
end
end