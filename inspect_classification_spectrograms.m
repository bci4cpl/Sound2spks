% This script will open the sound2spks eval csv for tempotron
% load associated WAV files, and inspect their spectrograms
clc; close all; clear;

plot_flag = 1;
play_flag = 1;
FILES_TO_DISP = 4;

%params for binary spike trains conversion
dt_sim = 1/10;
fs_sim = 1000 / dt_sim;

%% loading meta-data
T_all = readtable('/datasets/spiking/Sound2spks/tempotron/eval_dataset_for_tempotron.v2.csv');
WAV_DIR = '/datasets/spiking/Sound2spks/wav_data';
MAT_DIR = '/datasets/spiking/Sound2spks/result_mats';

assert(isfolder(WAV_DIR));


%% main loop
for label=[0,1]
    ME = [];
    %h_waitbar = waitbar(0,' Starting ...');
    try
        T = T_all(T_all.label == label, :);
        
        N_FILES = size(T, 1);
        MAX_FILES = N_FILES;
        N_FILES = min(MAX_FILES, N_FILES);
        for i_file=1:N_FILES
            if (~mod(i_file,ceil(N_FILES/20)) || i_file == 1)
                %waitbar(i_file/N_FILES, h_waitbar, ' Processing ...');
                fprintf('%.1f%% Processing ...\n', i_file/N_FILES*100);

            end

            suffix = sprintf('%d/%s', T.label(i_file), T.file{i_file});
            [~,name,ext] = fileparts(T.file{i_file});
            wav_full_name = fullfile(WAV_DIR, sprintf('%s.%s', name, 'wav'));
            assert(isfile(wav_full_name))
            
            mat_full_name = fullfile(MAT_DIR, sprintf('%s.%s', name, 'mat'));
            assert(isfile(mat_full_name))

            %% Part 1: load wav

            % Load audio file
            [audioIn,fs] = audioread(wav_full_name); % read audio file and sampling rate
            i_chan = 1;
            if (size(audioIn, 2) == 2) % in case of stereo - take only one side
                i_chan = randi(2);
            end
            audioIn = audioIn(:,i_chan);

            dur_in_sec = length(audioIn)/fs;
            fprintf('[%d/%d] File `%s`, [%d]`%s`, duration: %.2f secs ...\n', ...
                i_file, N_FILES, ...
                name, T.label(i_file), T.annotation{i_file}, dur_in_sec);

            %% Part 2: plot spects
            
            % parameters for Mel spectrogram
            WL = 256; % Window Length
            OLL = round(0.5*WL); % Overlap Length
            FFTL = 2*WL; % FFT Length
            NumOfBands = 16; % Like paper of Gutig
            f1 = 80; %was 360; lower frequency
            f2 = 8000; % upper frequency

            if plot_flag
                % Calcualte Mel Spectrogam
                figure('Name', suffix);

                %spectrogram(audioIn, WL);
                ax = subplot(2,1,1); 
                
            end
            
            [S,F,T_] = melSpectrogram_ex(audioIn,fs,plot_flag, ...
                'WindowLength',WL,...
                'OverlapLength',OLL, ...
                'FFTLength',FFTL, ...
                'NumBands',NumOfBands, ...
                'FrequencyRange',[f1 f2]);
            
            M = hz2mel(F)/1000;
            for b=1:NumOfBands
                %yline(b*((f2-f1)/16)/1000, '-', 'Color', [.5 .5 .5], 'Alpha', .2);
                yline(M(b), '-', 'Color', [.5 .5 .5], 'Alpha', .5);
                hold on;
            end
            
            if plot_flag
                xticks(ax, 0:ceil(length(audioIn)/fs))
            end
            
            %% Part 3: plot spikes
            if plot_flag
                ax = subplot(2,1,2);
                
                load(mat_full_name);
                [n_batches, nrn_count, n_iterations] = size(batch_items);
                i_batch_item = 5; %FUTURE: can we choose it wiser?
                %figure('Name', sprintf('spike train #%d', i_batch_item))
                
                spk_train = squeeze(batch_items(i_batch_item, :, :));
                
                for b=[1:NumOfBands:nrn_count]
                    yline(b, '-', 'Color', [.5 .5 .5], 'Alpha', .2);
                    hold on;
                end

                x = 0:fs_sim;
                for k = 1:nrn_count
                    y = k;
                    plot(x,y*(spk_train(k,:)'),'.');
                    hold on
                end
                axis([0 fs_sim 0.5 nrn_count+1]);
                xticklabels(xticks()/fs_sim);
                
                ax.YDir = 'reverse'; %order neurons from top to bottom on y axis - we align them with spect
                
                xlabel('Time (S)')
                
                
                N_LOWEST = 2;
                fprintf('\tTaking %d lowest-freq neurons amounts to the range of %.1f-%.1f Hz.\n', ...
                    N_LOWEST, f1, f1+sum(F(1:N_LOWEST)));
            end
            
            %% Part 4: Play wav
            if (play_flag)
                % play(sprintf('%s.%d', f.name, i_file), audioIn, fs_sim);
                p = audioplayer(audioIn, fs);
                playblocking(p);
            end

            %%
            if i_file==FILES_TO_DISP
                break
            end
            
        end
    catch ME
    end

    %waitbar(1, h_waitbar, ' DONE!');
    %pause(1);
    %close(h_waitbar);

    if ~isempty(ME), rethrow(ME); end
end