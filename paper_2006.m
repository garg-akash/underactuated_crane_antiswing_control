clc;
clear all;

T = 20.0;
Ts = 0.001;
t = 0:Ts:T;
numSteps = size(t,2);

rx = zeros(numSteps,1); drx = zeros(numSteps,1); ddrx = zeros(numSteps,1);
rl = zeros(numSteps,1); drl = zeros(numSteps,1); ddrl = zeros(numSteps,1);
rth = zeros(numSteps,1); drth = zeros(numSteps,1); ddrth = zeros(numSteps,1);

x = zeros(numSteps,1); dx = zeros(numSteps,1); ddx = zeros(numSteps,1);
l = zeros(numSteps,1); dl = zeros(numSteps,1); ddl = zeros(numSteps,1);
th = zeros(numSteps,1); dth = zeros(numSteps,1); ddth = zeros(numSteps,1);

% Get Travelling and Hoisting Trajectories
rx_wpts = [0 1.5 1.5];
[rx, drx, ddrx] = trapveltraj(rx_wpts, T/Ts + 1, 'EndTime', [12.5 20.0], 'Acceleration', [0.04 0.05]);

rl_wpts = [1.5 0.6 0.6 1.5 1.5];
[rl, drl, ddrl] = trapveltraj(rl_wpts, T/Ts + 1, 'EndTime', [5.5 3.5 4.5 20.0], 'Acceleration', [0.12 0.18 0.18 0.05]);

%% Plots

grid on
hold on
plot(t,rl)
ylabel('X desired');
xlabel('Time [sec]');
hold on
