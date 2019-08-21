d = dir('WAV_Data/*.wav');
for i=1:min(1,length(d))
    f = d(i);
    wav_path = strcat(f.folder, '/', f.name);
    [audioIn,fs] = audioread(wav_path);
    
    %1 sec
    i_chan = 1;
    audioIn = audioIn(:,i_chan); 
    audioIn = audioIn(1:fs); %first sec
        
    dur = length(audioIn)/fs;
    fprintf('Playing `%s` for %.2f secs ...\n', f.name, dur);
    sound(audioIn, fs)
    pause(dur);
        
    fprintf('Playing resample ...\n', f.name, dur);
    fs_dest = 1000;
    rel = fs/fs_dest;
    %audioIn = ifft(fft(audioIn), length(audioIn)/rel);
    audioIn_re = resample(audioIn,fs_dest, fs);
    sound(audioIn_re, fs_dest)
    pause(dur);
end