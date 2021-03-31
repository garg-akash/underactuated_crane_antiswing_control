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
M = 6.5; m=1; g=9.8;
pdx = 0.6; pdl_dash = 0.25;
kvx = 0.4; kvl = 0.4; kax = 0.4; kal = 0.4;
epsilonX = 3.5; epsilonL = 3.5;
psiX = 0.02; psiL = 0.015;
kpx = 240; kdx = 100; kpl = 45; kdl = 10;
lambdaX = 0.1; lambdaL = 0.1; 
phiX_bar = 2; phiL_bar = 5;
%phiX_bar = 2; phiL_bar = 5; %TODO

Frx_0 = 4.4; 
krx = -0.5; 
epsilon = 0.01;
dL = 6.5;
l(1) = 0.5;
pdl = l(1) + pdl_dash;
dTH = 0; %TODO

xd_prev = 0;
ld_prev = l(1);
dxd_prev = 0;
dld_prev = 0;

for k=1:numSteps
    xd(k) = pdx/2 + kvx^2*log(cosh(2*kax*t(k)/kvx-epsilonX)/cosh(2*kax*t(k)/kvx-epsilonX-2*kax*pdx/kvx^2))/(4*kax);
    ld(k) = (pdl_dash+2*l(1))/2 + kvl^2*log(cosh(2*kal*t(k)/kvl-epsilonL)/cosh(2*kal*t(k)/kvl-epsilonL-2*kal*pdl_dash/kvl^2))/(4*kal);
    dxd(k) = (xd(k)-xd_prev)/Ts;
    dld(k) = (ld(k)-ld_prev)/Ts;
    ddxd(k) = (dxd(k)-dxd_prev)/Ts;
    ddld(k) = (dld(k)-dld_prev)/Ts;
    
    ex(k) = x(k)-xd(k);
    el(k) = l(k)-ld(k);
    eth(k) = th(k);
    
    dex(k) = dx(k)-dxd(k);
    del(k) = dl(k)-dld(k);
    deth(k) = dth(k);
    
%     ddex(k) = ddx(k)-ddxd(k);
%     ddel(k) = ddl(k)-ddld(k);
%     ddeth(k) = ddth(k);

    phiX = 5.5 * dx(k); %TODO
    phiL = 21.0 * dl(k);
    phiX_store(k) = phiX;
    phiL_store(k) = phiL;
    
    q(:,k) = [x(k);l(k);th(k)];
    dq(:,k) = [dx(k);dl(k);dth(k)];
    Frx = (Frx_0 * tanh(dx(k)/epsilon)) - (krx * abs(dx(k)) * dx(k));
    [Mc,C,G] = computeMCG(q(:,k),dq(:,k));
    ux(k) = -kpx*ex(k) - 2*lambdaX*psiX^2*ex(k)/(psiX^2-ex(k)^2)^2 - kdx*dex(k) + (M+m)*ddxd(k) + m*ddld(k)*sin(th(k)) + m*dld(k)*dth(k)*cos(th(k)) + Frx - phiX_bar*sign(dex(k));
    ul(k) = -kpl*el(k) - 2*lambdaL*psiL^2*el(k)/(psiL^2-el(k)^2)^2 - kdl*del(k) + m*ddxd(k)*sin(th(k)) + m*ddld(k) - m*g + dL*dl(k) - phiL_bar*sign(del(k));
    U(:,k) = [ux(k);ul(k);0];
    Fd = [-Frx+phiX; -dL*dl(k)+phiL; -dTH*dth(k)];
    ddq(:,k) = inv(Mc)*(U(:,k) + Fd - G - C*dq(:,k));
    dq(:,k+1) = dq(:,k) + ddq(:,k)*Ts; 
    q(:,k+1) = q(:,k) + dq(:,k)*Ts + 0.5*ddq(:,k)*Ts.^2; 
    q(3,k+1) = wrapToPi(q(3,k+1));
    
    x(k+1) = q(1,k+1);
    l(k+1) = q(2,k+1);
    th(k+1) = q(3,k+1);
    dx(k+1) = dq(1,k+1);
    dl(k+1) = dq(2,k+1);
    dth(k+1) = dq(3,k+1);
    ddx(k) = ddq(1,k);
    ddl(k) = ddq(2,k);
    ddth(k) = ddq(3,k);
    
    xd_prev = xd(k);
    ld_prev = ld(k);
    dxd_prev = dxd(k);
    dld_prev = dld(k);
    
end



%% Functions

function [Mc,C,G] = computeMCG(q,dq)
    global M m g 
    Mc = [M+m           m*sin(q(3))  m*q(2)*cos(q(3));
          m*sin(q(3))   m            0;
          m*q(2)*cos(q(3)) 0         m*q(2)^2];
    
    C = [0      2*m*cos(q(3))*dq(3)     -m*q(2)*sin(q(3))*dq(3);
         0      0                       -m*q(2)*dq(3);
         0      2*m*q(2)*dq(3)           0];
     
    G = [0 -m*g*cos(q(3)) m*g*q(2)*sin(q(3))]';
end

