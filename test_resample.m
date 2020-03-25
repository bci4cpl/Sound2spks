% Oren Tested if resample() produces artifacts at certain freq bands

clear
close all

fs = 44100;
dt = 1/fs;
fs_sim = 10000;
dt_dim = 1/fs_sim;

T = 15; % sec
N = T*fs;
t = (1:N)*dt;

% white noise
s1 = randn(1,N);
[pxx1,f1] = pwelch(s1,fs/2,fs/4,fs/2,fs);
s1r = resample(s1, fs_sim, fs);
[pxx1r,f1r] = pwelch(s1r,fs_sim/2,fs_sim/4,fs_sim/2,fs_sim);

figure('Color','w')
subplot(2,1,1)
loglog(f1,pxx1)
set(gca,'Xlim',[1 20000]);
xlabel('f (Hz)');
ylabel('Power')
subplot(2,1,2)
loglog(f1r,pxx1r)
set(gca,'Xlim',[1 20000]);
xlabel('f (Hz)');
ylabel('Power')

% sine wave
f_stim = 15000; % Hz
s2 = sin(2*pi*f_stim*t);
[pxx2,f2] = pwelch(s2,fs/2,fs/4,fs/2,fs);
s2r = resample(s2, fs_sim, fs);
[pxx2r,f2r] = pwelch(s2r,fs_sim/2,fs_sim/4,fs_sim/2,fs_sim);

figure('Color','w')
subplot(2,1,1)
loglog(f2,pxx2)
set(gca,'Xlim',[1 20000]);
xlabel('f (Hz)');
ylabel('Power')
subplot(2,1,2)
loglog(f2r,pxx2r)
set(gca,'Xlim',[1 20000]);
xlabel('f (Hz)');
ylabel('Power')


