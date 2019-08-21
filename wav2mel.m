% Convert WAV files into mel spectrograms

clear; clc; close all

%general parameters
plot_flag = 0;
play_flag = 1;
epsilon = 1e-5; % for spectrgram logarithmic rescaling
tau_smooth = 20e-3; % smoothing time constant = 20 msec
th_vec = [0.0625, 0.125, 0.1875, 0.25, 0.3125, ...
    0.375, 0.4375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, ...
    0.875, 0.9375, 1]; % vector of thresholds for different channels
N_th = length(th_vec);

% parameters for Mel spectrogram
WL = 256; % Window Length
OLL = round(0.5*WL); % Overlap Length
FFTL = 2*WL; % FFT Length
NumOfBands = 16; % Like paper of Gutig
f1 = 360; % lower frequency
f2 = 8000; % upper frequency

%DATA_dir = './WAV_Data'; % directory with WAV files
DATA_dir = '/media/shay/DATA/Sound2spk'; % directory with WAV files

% Loop over sound file
N_MAX = 1; %limits the files we process
d = dir(fullfile(DATA_dir, '*.wav'));
for i=1:max(N_MAX,length(d))
    %fname1 = '2193_4qLAn6_xfCY.wav'; % specific file
    %fname = fullfile(DATA_dir, fname1);
    
    f = d(i);
    fname = strcat(f.folder, '/', f.name);
    
    % Load audio file
    [audioIn,fs] = audioread(fname); % read audio file and sampling rate
    i_chan = 1;
    if (size(audioIn, 2) == 2) % in case of stereo - take only one side
        i_chan = randi(2);
    end
    audioIn = audioIn(:,i_chan);
    
    dur = length(audioIn)/fs;
    fprintf('File `%s` duration: %.2f secs ...\n', f.name, dur);
    
    % Cut audio file into 1-sec segments (TODO)
    fs_re = 10000;
    
    %rel = fs/fs_dest;
    %audioIn = ifft(fft(audioIn), length(audioIn)/rel);
    
    %https://www.mathworks.com/help/signal/ref/resample.html
    audioIn_re = resample(audioIn,fs_re, fs);
    
    i_cur = 0;
    for i_start=1:fs_re:length(audioIn_re)-mod(length(audioIn_re), fs_re)
        i_cur = i_cur+1;
        audioIn_cur = audioIn_re(i_start:i_start+fs_re-1);
        assert(length(audioIn_cur)==fs_re)
        
        if (play_flag)
            play(sprintf('%s.%d', f.name, i_cur), audioIn_cur, fs_re);
        end
        
        % Calcualte Mel Spectrogam
        [S,F,T] = melSpectrogram_ex(audioIn_cur,fs_re,plot_flag, ...
            'WindowLength',WL,...
            'OverlapLength',OLL, ...
            'FFTLength',FFTL, ...
            'NumBands',NumOfBands, ...
            'FrequencyRange',[f1 f2]);
        
        % Extract Number of bandpass filters in filterbank and Number of frames in spectrogram
        [numBands,numFrames] = size(S);
        fprintf("Number of bandpass filters in filterbank: %d\n",numBands)
        fprintf("Number of frames in spectrogram: %d\n",numFrames)
        
        % Extract duration of each time frame
        t = T';
        dt = mean(diff(t));
        n_smooth = round(tau_smooth/dt);
        w_smooth = gausswin(n_smooth); % Gaussian window for smoothing
        w_smooth = w_smooth/sum(w_smooth); % normalize window
        fprintf("Duration of each frame: %d\n",dt)
        
        % Display frequency vector
        disp(' ');
        disp('Frequencies = ');
        disp(F);
        
        % Display spectrogram
        % melSpectrogram(audioIn,fs)
        % figure; imagesc(S); colorbar
        
        % Normalize spectrogram (Gutig)
        S = S / max(S(:));
        S = log(S + epsilon) - log(epsilon);
        if plot_flag
            figure; imagesc(flipud(S)); colorbar
        end
        
        % Smooth a single channel using a Gaussian filter
        % y = filter(w_smooth,1,S(1,:));
        % subplot(2,1,1);plot(t,S(1,:));xlabel('t (s)');
        % subplot(2,1,2);plot(t,y);xlabel('t (s)');
        
        nrn_count = 0;
        for k = 1:NumOfBands
            S1 =   S(k,:);
            S1 = filter(w_smooth,1,S1); % smooth signal using a Gaussian filter
            min1 = min(S1);
            S1 = S1 - min1;
            max1 = max(S1);
            S1 = S1 / max(S1);
            S_norm(k,:) = S1;
            
            for m = 1:(N_th-1)
                S1_TH = (S1 >= th_vec(m));
                S1_TH_diff = diff (S1_TH);
                % Onset spikes
                nrn_count = nrn_count + 1;
                spk_times(nrn_count).t = t(S1_TH_diff == 1);
                % Offset spikes
                nrn_count = nrn_count + 1;
                spk_times(nrn_count).t = t(S1_TH_diff == -1);
            end
            
            m = N_th;
            S1_TH = (S1 == th_vec(m));
            S1_TH_diff = diff (S1_TH);
            nrn_count = nrn_count + 1;
            spk_times(nrn_count).t = t(S1_TH_diff == 1);
        end
    end
end


function play(name, audioIn, Fs)
    dur = length(audioIn)/Fs;
    fprintf('Playing `%s` for %.2f secs ...\n', name, dur);
    sound(audioIn, Fs)
    pause(dur);
end
