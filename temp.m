clear; clc; close all
dbstop if error

% h_waitbar = waitbar(0,' Starting ...');

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

%DATA_DIR = './WAV_Data'; % directory with WAV files
DATA_DIR = '/media/shay/DATA/Sound2spks'; % directory with WAV files
%TODO: DEST_DIR = '/media/shay/DATA/Sound2spks/results'
DEST_DIR = '~/Sound2spks';

assert(isdir(DATA_DIR));
assert(isdir(DEST_DIR));

% Loop over sound file
N_MAX = inf %100; %TODO: limits the files we process
d = dir(fullfile(DATA_DIR, '*.wav'));
N_FILES = min(N_MAX,length(d));
N_FILES
