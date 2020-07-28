clear; clc; close all
dbstop if error

%we are going to multiply this by the "fs_sim"
%this will also be each batch_item size
sim_dur_in_sec = 1;

dt_sim = 1/10;
fs_sim = sim_dur_in_sec*1000 / dt_sim;

fs_sim

to_msec_idx = @(in_sec) floor(in_sec*1000/dt_sim)+1;
to_msec_idx(1)-1
to_msec_idx(2.5)-1