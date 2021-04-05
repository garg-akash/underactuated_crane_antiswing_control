clc;
clear all;

T = 20.0;
Ts = 0.001;
t = 0:Ts:T;
numSteps = size(t,2);


% Get Travelling and Hoisting Trajectories
rx_wpts = [0 1.5 1.5];
[rx, drx, ddrx] = trapveltraj(rx_wpts, T/Ts + 1, 'EndTime', [12.5 20.0], 'Acceleration', [0.04 0.05]);

rl_wpts = [1.5 0.6 0.6 1.5 1.5];
[rl, drl, ddrl] = trapveltraj(rl_wpts, T/Ts + 1, 'EndTime', [5.s5 3.5 4.5 20.0], 'Acceleration', [0.12 0.18 0.18 0.05]);

td = (mean(rl(1:4000)) + min(rl(1:4000))) / 2;
ta = (mean(rl(5500:8000)) + min(rl(5500:8000))) / 2;

aa = max(drx)/ta;
ad = max(drx)/td;