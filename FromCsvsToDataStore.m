% Converts MEL spectrogram spike times CSVs into A Data Store single CSV

clear; clc; close all
dbstop if error

DATA_DIR = '/media/shay/DATA/Sound2spks/result_csvs';
DEST_DIR = '/media/shay/DATA/Sound2spks/result_datastore';

assert(isdir(DATA_DIR));
assert(isdir(DEST_DIR));



ds = tabularTextDatastore(DATA_DIR,'FileExtensions','.csv');
ds.ReadSize = 'file';
%T = readall(ds);
%size(T)

%T.Properties

%print(T.Properties.VariableNames)

%preview(ds)

data = [];
ix = 0;
n_files = 0;
last_cycle = 0;
while hasdata(ds)
    [T, info] = read(ds);
    
    [~, baseFileName, extension] = fileparts(info.Filename);
    fprintf('Processing `%s` ...\n', baseFileName)
    n_rows = size(T,1);
    iFrom = ix+1; iTo = ix+n_rows;
    T.index = (iFrom:iTo)';
    T.filename = repmat({baseFileName}, [n_rows, 1]);
    T.file_ix = repmat(n_files+1, [n_rows, 1]);
    T.cycle = last_cycle + T.cycle;
    
    ix = iTo;    
    last_cycle = T.cycle(end);
    
    %reorganize
    T = [T(:, ~contains(T.Properties.VariableNames, 'st')) ...
        T(:, contains(T.Properties.VariableNames, 'st'))];
    
    data = [data; T];
    %size(data)
    
    n_files = n_files+1;
    
    %if (n_files>1)
    %    break
    %end

end

size(data)
%head(data)
%tail(data)

fprintf('Writing to final datastore ...\n')
    
DEST_FILE = fullfile(DEST_DIR, 'datastore.csv');
writetable(data, DEST_FILE, 'QuoteStrings', true);

DEST_FILE
