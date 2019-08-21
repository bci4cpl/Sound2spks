
function varargout = melSpectrogram_ex(x,fs,plot_flag,varargin)
%MELSPECTROGRAM Compute mel spectrogram.
%   S = MELSPECTROGRAM(audioIn,fs) returns the mel spectrogram of the audio
%   input, audioIn. fs is the input sample rate, in Hz. Columns of the
%   input are treated as individual channels. S is returned as an
%   L-by-M-by-N array:
%       L - Number of frequency points in the spectrogram. This is
%           determined by the NUMBANDS property.
%       M - Number of frames the audio signal is partitioned into.
%           This is determined by the WINDOWLENGTH and OVERLAPLENGTH
%           properties.
%       N - Number of channels.
%
%   S = MELSPECTROGRAM(...,'WindowLength',WINDOWLENGTH) specifies the
%   analysis window length used to calculate the spectrogram. Specify the
%   window length in samples as a positive scalar. If unspecified,
%   WINDOWLENGTH defaults to round(0.030 * fs).
%
%   S = MELSPECTROGRAM(...,'OverlapLength',OVERLAPLENGTH) specifies the
%   number of samples overlap between adjacent windows. Specify the overlap
%   length as a positive scalar integer smaller than the window length. If
%   unspecified, OVERLAPLENGTH defaults to round(0.020 * fs).
%
%   S = MELSPECTROGRAM(...,'FFTLength',FFTLENGTH) specifies the FFT length.
%   Specify the FFT length as a positive scalar integer larger than or
%   equal to the window length. If unspecified, FFTLENGTH defaults to
%   WINDOWLENGTH.
%
%   S = MELSPECTROGRAM(...,'NumBands',NUMBANDS) specifies the number of
%   bands in the mel filter bank. If unspecified, NUMBANDS defaults to 32.
%
%   S = MELSPECTROGRAM(...,'FrequencyRange',FREQUENCYRANGE) specifies the
%   frequency range over which to compute the spectrogram. If unspecified,
%   FREQUENCYRANGE defaults to [0, fs/2].
%
%   S = MELSPECTROGRAM(...,'SpectrumType',SPECTRUMTYPE) specifies the type
%   of the spectrogram as either 'power' or 'magnitude'. If unspecified,
%   SPECTRUMTYPE defaults to 'power'.
%
%   [S,F,T] = MELSPECTROGRAM(...) returns the center frequencies of the
%   bands (in Hz) and the location (in seconds) of each window of data. The
%   location corresponds to the center of each window.
%
%    MELSPECTROGRAM(...) with no output arguments plots the mel spectrogram
%    on a surface in the current figure.
%
%   EXAMPLE 1: Visualize the mel spectrogram for entire speech file
%     % Use the default settings to compute and visualize the mel
%     % spectrogram from a speech file.
%
%       [audioIn,fs] = audioread('Counting-16-44p1-mono-15secs.wav');
%       melSpectrogram(audioIn,fs)
%
%   EXAMPLE 2: Compute the mel spectrogram for entire audio file
%     % Use non-default settings to compute the mel spectrogram
%     % from an audio file.
%
%       [audioIn,fs] = audioread('FunkyDrums-44p1-stereo-25secs.mp3');
%       [S,F,T] = melSpectrogram(audioIn,fs,'WindowLength',1024,...
%                                'OverlapLength',512,'NumBands',64,...
%                                'FFTLength',2048);
%
% See also MFCC, CEPSTRALFEATUREEXTRACTOR, SPECTROGRAM

% Copyright 2018 The Mathworks, Inc.

%#codegen

% Parse and validate inputs
validateRequiredInputs(x, fs)

% Parse and validate optional inputs
params =  audio.internal.MelSpectrogramValidator(fs,size(x,1),varargin{:});

range          = hz2mel(cast(params.FrequencyRange,class(x)));
bandEdges      = mel2hz(linspace(range(1),range(end),params.NumBands+2));
[filterBank,F,FFTLengthTooSmall] = audio.internal.designMelFilterBank( ...
    fs,bandEdges,params.FFTLength,2,'Bandwidth','Hz',class(x));

if isempty(coder.target) && FFTLengthTooSmall
    coder.internal.warning('audio:melspectrogam:FFTLengthTooSmall');
end

[~,~,Z] = audio.internal.spectralDescriptors.stft(x,fs,params);

Z = abs(Z);
if strcmp(params.SpectrumType , 'power')
    Z = Z.^2;
end

P = filterBank.' * Z;

S = reshape(P , size(P,1) , size(P,2) / size(x,2) , size(x,2));

N            = params.WindowLength;
hopLength    = N - params.OverlapLength;
nRow         = size(x,1);
numHops      = floor((nRow-N)/hopLength) + 1;
T            = cast(((0:(numHops-1))*hopLength + N/2)' / fs,'like',x);

if nargout > 0
    varargout{1} = S;
    
    if nargout > 1
        varargout{2} = F;
    end
    
    if nargout > 2
        varargout{3} = T;
    end
end

if plot_flag
    plotopts.isFsnormalized = false;
    plotopts.freqlocation   = 'yaxis';
    if strcmp(params.SpectrumType,'power')
        S =  10*log10(S(:,:,1)+eps(class(x)));
    else
        S =  20*log10(S(:,:,1)+eps(class(x)));
    end
    
    ax = signalwavelet.internal.convenienceplot.plotTFR(T,hz2mel(F), S,plotopts);
    yticks(ax , hz2mel(F)/1000);
    
    numF   = length(F);
    % Include 5 ticks
    if numF > 5
        indices = fix(linspace(1, numF,5));
    else
        indices = 1:numF;
    end
    labels = cell(1,length(F));
    for index = indices
        labels{index} = num2str(F(index)/1000 , 3);
    end
    yticklabels(ax , labels);
end

end

% -------------------------------------------------------------------------
% Validate required inputs
% -------------------------------------------------------------------------
function validateRequiredInputs(x,fs)
validateattributes(x,{'single','double'},...
    {'nonempty','2d','real'}, ...
    'melSpectrogram','audioIn')
validateattributes(fs,{'single','double'}, ...
    {'nonempty','positive','real','scalar','nonnan','finite'}, ...
    'melSpectrogram','fs');
end
