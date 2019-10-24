% Loads the datastore CSV full of spktrains of neurons and transforms them
% into EntropyMax compatible batch items of binary spike trains

clear; clc; close all
dbstop if error

DATASTORE = '/media/shay/DATA/Sound2spks/result_datastore/datastore.csv';

N_NRNS = 496;
N_MAX_TIMES = 27;
BATCH_SIZE = 10;

%params for matlab conversion
dt_sim = 1; %p.dt = castAsUnderlineType(p, 1); % time step (1/(how many slots per msec))
fs_sim = 1000 / dt_sim;

assert(isfile(DATASTORE));


ds = tabularTextDatastore(DATASTORE);
%ds.ReadSize = 'file';
T = readall(ds);
size(T)

%T.Properties

%head(T)
tail(T)

N_HEADER_COLS = sum(cellfun('isempty', regexp(T.Properties.VariableNames, 'st\d+')));
ST_COLS = N_HEADER_COLS+1:size(T,2);

[n_rows, n_cols] = size(T);
ixs_of_batch_items = 1:496:n_rows;

to_msec_idx = @(in_sec) fix(in_sec*fs_sim)+1; %assuming dt=1, fix=floor and to_int(), +1 to start for idx 1

ixs_sampled = randsample(ixs_of_batch_items, BATCH_SIZE);

sz = ...
    [BATCH_SIZE, ... %TOTAL BATCH ITEMS
    N_NRNS, ... %INPUT NEURONS
    fs_sim+1]; %TIME OF SIMULATION IN DT UNITS
%NOTE: we add one, so it'd be easier later

batchItems = false(sz);
            
for i_sample=1:length(ixs_sampled)
    ix = ixs_sampled(i_sample);
    iFrom = ix; iTo = iFrom+N_NRNS-1;
        
    %assert(T{iFrom, 'nrn'} == 1);
    %assert(T{iTo, 'nrn'} == N_NRNS);
    mat = T{iFrom:iTo, ST_COLS};
    assert(isequal(size(mat), [N_NRNS, N_MAX_TIMES]));
    
    spk_train = false(N_NRNS, fs_sim+1); %NOTE: we add one, so it'd be easier later
    for i_nrn=1:N_NRNS
        spk_times = mat(i_nrn, :);
        spk_times = spk_times(~isnan(spk_times));
        i_times = to_msec_idx(spk_times);
        assert(all(i_times<fs_sim))
        batchItems(i_sample, i_nrn, i_times) = 1;
    end
    
    figure;
    ax = gca;
    for iNeuron = 1:N_NRNS
        sp = find(batchItems(i_sample,iNeuron,:));
        spy = iNeuron+0*sp; %put neuron index in spike times of neuron
        plot(ax, sp,spy,'.');
        hold(ax, 'on');
    end
    axis(ax, [0 fs_sim+1 0 N_NRNS+1]);
    title(ax, sprintf('Input layer raster (BatchItem: %d)', i_sample));
    ylabel(ax, '# Neuron')
    xticklabels(ax, []);
end

    
