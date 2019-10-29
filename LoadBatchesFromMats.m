% Loads the datastore CSV full of spktrains of neurons and transforms them
% into EntropyMax compatible batch items of binary spike trains

clear; clc; close all
dbstop if error

DATA_DIR = '/media/shay/DATA/Sound2spks/result_mats';
% DATA_DIR = '~/Temp/Sound2spks';
assert(isfolder(DATA_DIR));

%Determintic psedu-random generator
RNG = RandStream.create('mrg32k3a','NumStreams',1);

%Empirical values
N_NRNS = 496;

%Script Params
BATCH_SIZE = 3;

PLOT = 0;

%params for binary spike trains conversion
dt_sim = 1/10;
fs_sim = 1000 / dt_sim;

d = dir(fullfile(DATA_DIR, '*.mat'));
N_FILES = length(d);

files = randsample(RNG, 1:N_FILES, BATCH_SIZE);

for i_file=1:length(files)
    ix = files(i_file);
    f = d(ix);
    
    fname = strcat(f.folder, '/', f.name);
    data = load(fname);
    
    fprintf('[%d] File `%s` duration: %.2f secs ...\n', i_file, f.name, data.dur_in_sec);
    
    assert(data.dt_sim == dt_sim, 'Incosistent dt_sim')
    assert(data.fs_sim == fs_sim, 'Incosistent fs_sim')
    
    n_batch_items = size(data.batch_items, 1);
    %sample a batch within the file
    i_sample = randi(RNG, n_batch_items);
    spk_trains = squeeze(data.batch_items(i_sample,:,:)); %3d to 2d
    
    %% plot
    if (PLOT)
        
        figure;
        ax = gca;
        for i_nrn = 1:N_NRNS
            sp = find(spk_trains(i_nrn, :));
            spy = i_nrn+0*sp; %put neuron index in spike times of neuron
            plot(ax, sp,spy,'.');
            hold(ax, 'on');
        end
        axis(ax, [0 fs_sim+1 0 N_NRNS+1]);
        title(ax, sprintf('Input layer raster (BatchItem: %d)', i_sample));
        ylabel(ax, '# Neuron')
        xticklabels(ax, []);    

    end
end