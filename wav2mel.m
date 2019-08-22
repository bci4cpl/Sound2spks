% Convert WAV files into mel spectrograms

clear; clc; close all
dbstop if error

h_waitbar = waitbar(0,' Starting ...');

%general parameters
plot_flag = 0;
play_flag = 0;
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

%params for matlab conversion
dt_sim = 1/10;

%DATA_dir = './WAV_Data'; % directory with WAV files
DATA_dir = '/media/shay/DATA/Sound2spks'; % directory with WAV files

% Loop over sound file
N_MAX = inf; %100; %TODO: limits the files we process
d = dir(fullfile(DATA_dir, '*.wav'));
N_FILES = min(N_MAX,length(d));

total_failures = 0;
total_spk_times = 0;
all_spk_times = cell(N_FILES*10, 1); %~N_FILES*10
for i_file=1:N_FILES
    %fname1 = '2193_4qLAn6_xfCY.wav'; % specific file
    %fname = fullfile(DATA_dir, fname1);
    
    f = d(i_file);
    fname = strcat(f.folder, '/', f.name);
    
    try
        % Load audio file
        [audioIn,fs] = audioread(fname); % read audio file and sampling rate
        i_chan = 1;
        if (size(audioIn, 2) == 2) % in case of stereo - take only one side
            i_chan = randi(2);
        end
        audioIn = audioIn(:,i_chan);

        dur = length(audioIn)/fs;
        fprintf('File `%s` duration: %.2f secs ...\n', f.name, dur);

        %% Resample at a reasonable delta T
        fs_sim = 1000 / dt_sim;

        %rel = fs/fs_dest;
        %audioIn = ifft(fft(audioIn), length(audioIn)/rel);

        %https://www.mathworks.com/help/signal/ref/resample.html
        audioIn_re = resample(audioIn,fs_sim, fs);

        %% Mel spec trans
        if (play_flag)
            play(sprintf('%s.%d', f.name, i_file), audioIn_re, fs_sim);
        end

        % Calcualte Mel Spectrogam
        [S,F,T] = melSpectrogram_ex(audioIn_re,fs_sim,plot_flag, ...
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
        dt_actual = t(2)-t(1);
        n_smooth = round(tau_smooth/dt_actual);
        w_smooth = gausswin(n_smooth); % Gaussian window for smoothing
        w_smooth = w_smooth/sum(w_smooth); % normalize window
        fprintf("Duration of each frame: %d\n", dt_actual)

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
            figure('Name', 'Normalized spectrogram');
            imagesc(T, F, flipud(S)); colorbar
            xlabel('Time (S)')
            %imagesc(flipud(S)); colorbar
        end

        % Smooth a single channel using a Gaussian filter
        % y = filter(w_smooth,1,S(1,:));
        % subplot(2,1,1);plot(t,S(1,:));xlabel('t (s)');
        % subplot(2,1,2);plot(t,y);xlabel('t (s)');

        i_cur = 0;
        total = length(T);
        cycles = floor(1/dt_actual); %get approx. cycles in one second
        for i_start=1:cycles:total-mod(total, cycles)
            S_cur = S(:,i_start:i_start+cycles-1);
            t_cur = t(:,i_start:i_start+cycles-1)-i_cur; %[0.1]
            t_cur(t_cur<0) = 0;
            assert(all(t_cur>=0) && all(t_cur<=1))

            i_cur = i_cur+1;

            nrn_count = 0;
            spk_times = struct();
            for k = 1:NumOfBands
                S1 = S_cur(k,:);
                S1 = filter(w_smooth,1,S1); % smooth signal using a Gaussian filter
                min1 = min(S1);
                S1 = S1 - min1;
                max1 = max(S1);
                S1 = S1 / max(S1);
                %S_norm(k,:) = S1; TODO: OREN: is it necessary

                for m = 1:(N_th-1)
                    S1_TH = (S1 >= th_vec(m));
                    S1_TH_diff = diff (S1_TH);
                    % Onset spikes
                    nrn_count = nrn_count + 1;
                    spk_times(nrn_count).t = t_cur(S1_TH_diff == 1);
                    % Offset spikes
                    nrn_count = nrn_count + 1;
                    spk_times(nrn_count).t = t_cur(S1_TH_diff == -1);
                end

                m = N_th;
                S1_TH = (S1 == th_vec(m));
                S1_TH_diff = diff (S1_TH);
                nrn_count = nrn_count + 1;
                spk_times(nrn_count).t = t_cur(S1_TH_diff == 1);
            end

            to_msec_idx = @(in_sec) fix(in_sec*fs_sim)+1; %assuming dt=1, fix=floor and to_int(), +1 to start for idx 1
            spk_train = false(nrn_count, fs_sim+1); %NOTE: we add one, so it'd be easier later
            for i_nrn=1:nrn_count
                i_times = to_msec_idx(spk_times(i_nrn).t);
                assert(all(i_times<fs_sim))
                spk_train(i_nrn, i_times) = 1;
            end

            %disp(size(spk_train))

            if plot_flag
                figure('Name', sprintf('spike train #%d', i_cur))

                x = 0:fs_sim;
                for k = 1:nrn_count
                    y = nrn_count-k+1;
                    %y = k;
                    plot(x,y*(spk_train(k,:)'),'.');
                    hold on
                end
                axis([0 fs_sim 0.5 nrn_count+1]);
                xticklabels(xticks()/fs_sim);
                xlabel('Time (S)')
            end

            total_spk_times = total_spk_times+1;
            all_spk_times{total_spk_times} = spk_times;


            % times_tmp2_msec = to_msec_idx(times_tmp2);
            % batchItems(idx_batchitem, k, times_tmp2_msec) = 1;
        end % 1 sec intervals loop
        fprintf('1 sec segments: %d\n', i_cur);
    catch e
        total_failures = total_failures+1;
        fprintf('The processing of `%s` failed.\n', fname)
        fprintf('\tIdentifier: %s\n',e.identifier);
        fprintf('\tMessage: %s\n',e.message);
    end
    
    if ~mod(i_file,1000)
        waitbar(i_file/N_FILES, h_waitbar, ' Processing ...');
    end
end %N_FILES
close(h_waitbar);
fprintf('\n');
fprintf('all_spk_trains size: %d\n', length(all_spk_times))
fprintf('total_failures: %d\n', total_failures);

all_spk_times = all_spk_times(1:total_spk_times); %only valid
DEST = fullfile(DATA_dir, 'wav_dataset.mat');
save(DEST,'all_spk_times');
fprintf('Saved results to `%s`.', DEST)

function play(name, audioIn, Fs)
    dur = length(audioIn)/Fs;
    fprintf('Playing `%s` for %.2f secs ...\n', name, dur);
    sound(audioIn, Fs)
    pause(dur);
end

