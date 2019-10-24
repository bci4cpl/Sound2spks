DATA_FILE = '/media/shay/DATA/Sound2spks/result/wav_dataset.mat';
ix = 1;

% DATA_FILE = '~/Sound2spks/wav_dataset.mat';
% 
% tic
% data_load = load(DATA_FILE);
% 
% %data_load.all_spk_times
% sample_load = data_load.all_spk_times{ix};
% assert(isequal(size(sample_load), [1 496]))
% toc


tic
data = matfile(DATA_FILE,'Writable', false);
cell_with_sample = data.all_spk_times(ix,1);
sample = cell_with_sample{1}; %access the cell array
assert(isequal(size(sample), [1 496]))
toc

%assert(isequal(sample_load, sample))
assert(isequal(sample1_load, sample))
%assert(isequal(sample2_load, sample))