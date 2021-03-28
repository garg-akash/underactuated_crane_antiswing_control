clc
clear all

T = 20;
Ts = 0.001;
t = 0:Ts:T;
numSteps = size(t,2);

xd = zeros(numSteps,1); dxd = zeros(numSteps,1); ddxd = zeros(numSteps,1);
ld = zeros(numSteps,1); dld = zeros(numSteps,1); ddld = zeros(numSteps,1);
x = zeros(numSteps,1); dx = zeros(numSteps,1); ddx = zeros(numSteps,1);
l = zeros(numSteps,1); dl = zeros(numSteps,1); ddl = zeros(numSteps,1);
th = zeros(numSteps,1); dth = zeros(numSteps,1); ddth = zeros(numSteps,1);

ex = zeros(numSteps,1); dex = zeros(numSteps,1); ddex = zeros(numSteps,1);
el = zeros(numSteps,1); del = zeros(numSteps,1); ddel = zeros(numSteps,1);
eth = zeros(numSteps,1); deth = zeros(numSteps,1); ddeth = zeros(numSteps,1);

q = zeros(3,numSteps); dq = zeros(3,numSteps); ddq = zeros(3,numSteps);
ux = zeros(numSteps,1); ul = zeros(numSteps,1);
U = zeros(3,numSteps); dU = zeros(3,numSteps); ddU = zeros(3,numSteps);

global M m g
global Frx_0 krx epsilon dL dTH
M = 6.5; m=1; g=9.8;
pdx = 0.6; pdl_dash = 0.25;
kvx = 0.4; kvl = 0.4; kax = 0.4; kal = 0.4;
epsilonX = 3.5; epsilonL = 3.5;
psiX = 0.02; psiL = 0.015;
kpx = 240; kdx = 100; kpl = 45; kdl = 10;
lambdaX = 0.1; lambdaL = 0.1; 
phiX_bar = 0; phiL_bar = 0;%phiX_bar = 2; phiL_bar = 5; %TODO
Frx_0 = 4.4; 
krx = -0.5; 
epsilon = 0.01;
dL = 6.5;
l(1) = 0.5;
pdl = l(1) + pdl_dash;
dTH = 0; %TODO

xd_prev = 0;
ld_prev = 0;
dxd_prev = 0;
dld_prev = 0;
for k=1:numSteps
    xd(k) = pdx/2 + kvx^2*log(cosh(2*kax*t(k)/kvx-epsilonX)/cosh(2*kax*t(k)/kvx-epsilonX-2*kax*pdx/kvx^2))/4*kax;
    ld(k) = (pdl_dash+2*l(1))/2 + kvl^2*log(cosh(2*kal*t(k)/kvl-epsilonL)/cosh(2*kal*t(k)/kvl-epsilonL-2*kal*pdl_dash/kvl^2))/4*kal;
%     dxd(k) = (xd(k)-xd_prev)/Ts;
%     dld(k) = (ld(k)-ld_prev)/Ts;
%     ddxd(k) = (dxd(k)-dxd_prev)/Ts;
%     ddld(k) = (dld(k)-dld_prev)/Ts;
%     xd_prev = xd(k);
%     ld_prev = ld(k);
%     dxd_prev = dxd(k);
%     dld_prev = dld(k);    
end
