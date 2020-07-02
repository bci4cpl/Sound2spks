% Go over all .mat files in result_mats dir
% remove the fs_sim and add t_final with corrected values

clear; clc; close all
dbstop if error

h_waitbar = waitbar(0,' Starting ...');

%% Parameters
%general parameters
t_final = 1000;

if isunix()
    %DELL
    DATA_DIR = '/datasets/spiking/Sound2spks/result_mats'; % directory with MAT files
    %DATA_DIR = '~/Temp/Sound2spks/PoC';
    %DATA_DIR = '~/Infomax/Tempotrons/PoC';
    DEST_DIR = DATA_DIR;
end


assert(isdir(DATA_DIR));
assert(isdir(DEST_DIR));

% Loop over sound file
N_MAX = inf; %100; %NOTE: limits the files we process
d = dir(fullfile(DATA_DIR, '*.mat'));
N_FILES = min(N_MAX,length(d));

%% Wav to mel -> spikes according to THs -> binary spike trains
total_failures = 0;
for i_file=1:N_FILES
    f = d(i_file);
    mat_full_name = strcat(f.folder, '/', f.name);
    n_batch_items = 0;

    [~, baseFileName, extension] = fileparts(mat_full_name);
    DEST_FILE = fullfile(DEST_DIR, strcat(baseFileName, '.mat'));
    
    try
        %% Part 1: load mat 
        load(mat_full_name);
        fprintf('[%d/%d] File `%s` ...\n', ...
            i_file, N_FILES, ...
            f.name);
        
        if ~exist('fs_sim','var') %already modified
            disp('   skipping');
            continue
        end

        %% Part 2: Adjust fields 
        t_final = fs_sim*dt_sim; %in msec

        %% Part 3: persist to mat
        save(DEST_FILE, 'dur_in_sec', ...
            'F', 'numBands' ,'numFrames', 'dt_frame', ...
            'dt_sim', 't_final', ...
            'batch_items' ...
            ,'-v7.3','-nocompression');
    catch e
        total_failures = total_failures+1;
        fprintf('The processing of `%s` failed.\n', mat_full_name)
        fprintf('\tIdentifier: %s\n',e.identifier);
        fprintf('\tMessage: %s\n',e.message);
    end

    if ~mod(i_file,1000)
        waitbar(i_file/N_FILES, h_waitbar, ' Processing ...');
    end
end %N_FILES
close(h_waitbar);
fprintf('\n');
fprintf('total_failures: %d\n', total_failures);

