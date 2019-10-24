% Converts MEL spectrogram spike times Mat files into CSVs

clear; clc; close all
dbstop if error

DATA_DIR = '/media/shay/DATA/Sound2spks/result_mats';
DEST_DIR = '/media/shay/DATA/Sound2spks/result_csvs';

assert(isdir(DATA_DIR));
assert(isdir(DEST_DIR));

%Empirical values
N_MAX_TIMES = 27;
N_NRNS = 496;
    
% Loop over results
N_MAX = 100; %inf %100; %TODO: limits the files we process
d = dir(fullfile(DATA_DIR, '*.mat'));
N_FILES = min(N_MAX,length(d));

total_failures = 0;
total_spk_times = 0;
max_nrn_spk_times = 0;

for i_file=1:N_FILES
            
    f = d(i_file);
    fname = strcat(f.folder, '/', f.name);
    data = load(fname);
    
    fprintf('[%d] File `%s` duration: %.2f secs ...\n', i_file, data.wav_filename, data.dur_in_sec);
    
    [~, baseFileName, extension] = fileparts(data.wav_filename);
       
    
    n_rows = N_NRNS*data.n_file_spk_times;
    st = nan(n_rows, N_MAX_TIMES, 'single');
    nrns = zeros(n_rows, 1, 'uint16');
    cycles = zeros(n_rows, 1, 'uint16');
    
    i_row = 0;
    for i_cycle=1:data.n_file_spk_times
        batch_item = data.file_spk_times{i_cycle};
        assert(isequal(size(batch_item), [1 N_NRNS]))

        for i_nrn=1:N_NRNS
            i_row = i_row +1;
            L = size(batch_item(i_nrn).t, 2);
            st(i_row, 1:L) = batch_item(i_nrn).t(1:L);


            nrns(i_row) = i_nrn;
            cycles(i_row) = i_cycle;
        end
        
        
    end
    
    %persist
    DEST_FILE = fullfile(DEST_DIR, strcat(baseFileName, '.csv'));
    
    T = array2table(st);

    T.index = (1:n_rows)';
    T.nrn = nrns;
    T.cycle = cycles;

    %reorder
    T = [T(:, ~contains(T.Properties.VariableNames, 'st')) ...
        T(:, contains(T.Properties.VariableNames, 'st'))];

    %make sure we write blanks instead of `NaN`
    %https://stackoverflow.com/a/45805330/1640414
    for name = T.Properties.VariableNames  % Loop over variable names
      temp = num2cell(T.(name{1}));        % Convert numeric array to cell array
      temp(cellfun(@isnan, temp)) = {[]};  % Set cells with NaN to empty
      T.(name{1}) = temp;                  % Place back into table
    end
    
    %~5 secs
    writetable(T, DEST_FILE)
    
    fprintf('Saved results to `%s`.\n', DEST_FILE)
end

