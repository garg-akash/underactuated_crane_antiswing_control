clc;
clear all;
T = 20;
Ts = 0.001;
t = 0:Ts:T;
numSteps = size(t,2);

xd = zeros(numSteps,1); dxd = zeros(numSteps,1); ddxd = zeros(numSteps,1);
ld = zeros(numSteps,1); dld = zeros(numSteps,1); ddld = zeros(numSteps,1);
x = zeros(numSteps,1); dx = zeros(numSteps,1); ddx = zeros(numSteps,1);
l = zeros(numSteps,1); dl = zeros(numSteps,1); ddl = zeros(numSteps,1);
th = zeros(numSteps,1); dth = zeros(numSteps,1); ddth = zeros(numSteps,1);
noise_amp_ = 5; %(in degrees)
noise_start = 6; noise_end = 8;
freq = 1/(noise_end-noise_start);
noise_amp = [zeros(size(0:Ts:noise_start,2),1);deg2rad(noise_amp_).*ones(size(noise_start:Ts:noise_end,2)-1,1);zeros(size(noise_end:Ts:T,2)-1,1)];
nn = noise_amp.*sin(2*pi*freq*t');

ex = zeros(numSteps,1); dex = zeros(numSteps,1); ddex = zeros(numSteps,1);
el = zeros(numSteps,1); del = zeros(numSteps,1); ddel = zeros(numSteps,1);
eth = zeros(numSteps,1); deth = zeros(numSteps,1); ddeth = zeros(numSteps,1);

q = zeros(3,numSteps); dq = zeros(3,numSteps); ddq = zeros(3,numSteps);
ux = zeros(numSteps,1); ul = zeros(numSteps,1);
U = zeros(3,numSteps); dU = zeros(3,numSteps); ddU = zeros(3,numSteps);

global M m ml g
M = 6.5; m=1; g=9.8; ml=0;
pdx = 0.6; pdl_dash = 0.25;
kvx = 0.4; kvl = 0.4; kax = 0.4; kal = 0.4;
epsilonX = 3.5; epsilonL = 3.5;
psiX = 0.02; psiL = 0.015;
kpx = 240; kdx = 100; kpl = 45; kdl = 10;
lambdaX = 0.1; lambdaL = 0.1; 
% phiX_bar = 2; phiL_bar = 5;
% %phiX_bar = 0.01; phiL_bar = 0.1; %TODO

% Frx_0 = 3.5*4.4; 
% krx = -0.5; 
epsilon = 0.01;
% dL = 4*6.5; dX= 0;
l(1) = 0.5;
pdl = l(1) + pdl_dash;
dTH = 4.5; %TODO

% TODO: These params are from previous paper. Currently set to zero. Has to
% be tuned
kpx = 3.5;          kpl = 2;
kas = 0.6;          kli = 4;
eta_x = 2;          eta_l = 2;
dX= 8*2.5;          dL = 0.2;          

el_int = 0;
for k=1:numSteps
    xd(k) = pdx/2 + kvx^2*log(cosh(2*kax*t(k)/kvx-epsilonX)/cosh(2*kax*t(k)/kvx-epsilonX-2*kax*pdx/kvx^2))/(4*kax);
    ld(k) = (pdl_dash+2*l(1))/2 + kvl^2*log(cosh(2*kal*t(k)/kvl-epsilonL)/cosh(2*kal*t(k)/kvl-epsilonL-2*kal*pdl_dash/kvl^2))/(4*kal);
    dxd(k) = (kvx*sinh((2*kax*pdx)/kvx^2))/(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2));
    dld(k) = (kvl*sinh((2*kal*pdl_dash)/kvl^2))/(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2));
    ddxd(k) = -(kax*(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)^2 - 2*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2)^2))/(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)^2*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2)^2);
    ddld(k) = -(kal*(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)^2 - 2*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2)^2))/(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)^2*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2)^2);
%     ex(k) = x(k)-xd(k);
%     el(k) = l(k)-ld(k);
%     eth(k) = th(k);
%     
%     dex(k) = dx(k)-dxd(k);
%     del(k) = dl(k)-dld(k);
%     deth(k) = dth(k);
    ex(k) = -x(k)+xd(k);
    el(k) = -l(k)+ld(k);
    eth(k) = -th(k);
    
    dex(k) = -dx(k)+dxd(k);
    del(k) = -dl(k)+dld(k);
    deth(k) = -dth(k);
    
%     ddex(k) = ddx(k)-ddxd(k);
%     ddel(k) = ddl(k)-ddld(k);
%     ddeth(k) = ddth(k);

%     %phiX = 0.6 * dx(k); %TODO
%     %phiL = 21.0 * dl(k);
%     %phiX =  min(max(dx(k)/phiX_bar, -1), 1);
%     phiX =  min(max(dx(k), -phiX_bar), phiX_bar);
%     %phiL = min(max(dl(k)/phiL_bar, -1), 1);
%     phiL =  min(max(dl(k), -phiL_bar), phiL_bar);
%     phiX_store(k) = phiX;
%     phiL_store(k) = phiL;
    
    q(:,k) = [x(k);l(k);th(k)];
    dq(:,k) = [dx(k);dl(k);dth(k)];
%     Frx = (Frx_0 * tanh(dx(k)/epsilon)) - (krx * abs(dx(k)) * dx(k));
    [Mc,D,Cg] = computeMCG(q(:,k),dq(:,k),dX,dL);
    %ux(k) = -kpx*ex(k) - 2*lambdaX*psiX^2*ex(k)/(psiX^2-ex(k)^2)^2 - kdx*dex(k) + (M+m)*ddxd(k) + m*ddld(k)*sin(th(k)) + m*dld(k)*dth(k)*cos(th(k)) + Frx - phiX_bar*tanh(10*dex(k));%*sign(dex(k));
    %ul(k) = -kpl*el(k) - 2*lambdaL*psiL^2*el(k)/(psiL^2-el(k)^2)^2 - kdl*del(k) + m*ddxd(k)*sin(th(k)) + m*ddld(k) - m*g + dL*dl(k) - phiL_bar*tanh(10*del(k));%*sign(del(k));
    
    ux(k) = ddxd(k) + kpx*dex(k) - kas*deth(k);
    ul(k) = ddld(k) + kpl*del(k) + kli*el(k);
    
    ddz(:,k) = [ux(k);ul(k)];
    el_int = el_int + el(k)*Ts;
    fun = @(x) el(k);
    sx(k) = dex(k) + kpx*ex(k) - kas*eth(k);
    sl(k) = del(k) + kpl*el(k) + kli*el_int;
    
    eta_sgn = [eta_x*tanh(10*sx(k)) eta_l*tanh(10*sl(k))]';
    
    U(1:2,k) = Mc*(ddz(:, k) + eta_sgn) + D*dq(1:2,k) + Cg;

    %Fd = [-Frx+phiX; -dL*dl(k)+phiL; -dTH*dth(k)];
    
    ddq(1:2,k) = inv(Mc)*(U(1:2,k) - Cg - D*dq(1:2,k));
    ddq(3,k) = -(cos(q(3,k))*ddq(1,k) + 2*dq(2,k)*dq(3,k) + g*sin(dq(3,k)))/(q(2,k));
    
    dq(:,k+1) = dq(:,k) + ddq(:,k)*Ts; 
    q(:,k+1) = q(:,k) + dq(:,k)*Ts + 0.5*ddq(:,k)*Ts.^2; 
% %     q(3,k+1) = q(3,k+1) + noise_amp(k)*(-sin(2*pi*freq*t(k)));%noise_amp(k).*rand(1);
%     if (t(k)>=6 && t(k)<=8)
%         q(3,k+1) = noise_amp(k)*(-sin(2*pi*freq*t(k)));%noise_amp(k).*rand(1);
%     else
%         q(3,k+1) = q(3,k+1);
%     end
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
end
eth = rad2deg(eth);

%% Plot
figure
grid on
hold on
plot(t,xd)
hold on
plot(t,x(1:end-1))
ylabel('X desired and actual');
xlabel('Time [sec]');

figure
grid on
hold on
plot(t,ld)
hold on
plot(t,l(1:end-1))
ylabel('l desired and actual');
xlabel('Time [sec]');

figure
grid minor;
line(t, ex', 'Color', 'blue', 'LineWidth', 1);
hold on
plot(t,psiX*ones(size(t)),'--','Color','g')
hold on
plot(t,-psiX*ones(size(t)),'--','Color','g')
ylabel('Error in x [meters]');
xlabel('Time [sec]');
ylim([-0.1 0.1])

figure
grid minor;
line(t, el', 'Color', 'red', 'LineWidth', 1);
hold on
plot(t,psiL*ones(size(t)),'--','Color','g')
hold on
plot(t,-psiL*ones(size(t)),'--','Color','g')
ylabel('Error in l [meters]');
xlabel('Time [sec]');
ylim([-0.04 0.04])

figure
grid minor;
plot(t, eth', '-.r');
ylabel('Error in theta [deg]');
xlabel('Time [sec]');
ylim([-5 5])

figure
grid minor;
line(t, U(1,:)', 'Color', 'blue', 'LineWidth', 1);
ylabel('Ux(t) [N]');
xlabel('Time [sec]');
ylim([-5 20])

figure
grid minor;
line(t, U(2,:)', 'Color', 'blue', 'LineWidth', 1);
ylabel('Ul(t) [N]');
xlabel('Time [sec]');
ylim([-12 -4])

%% Functions
function [Mc,D,Cg] = computeMCG(q,dq,dX,dL)
    global M m ml g
    Mc = [M+m*(sin(q(3))^2) m*sin(q(3));
            m*sin(q(3))     ml+m];
%     Mc = [M+m           m*sin(q(3))  m*q(2)*cos(q(3));
%           m*sin(q(3))   m            0;
%           m*q(2)*cos(q(3)) 0         m*q(2)^2];
    
    D = [dX 0;
         0   dL];
     
%     C = [0      2*m*cos(q(3))*dq(3)     -m*q(2)*sin(q(3))*dq(3);
%          0      0                       -m*q(2)*dq(3);
%          0      2*m*q(2)*dq(3)           0];
     
    Cg = [-sin(q(3))*(m*q(2)*dq(3)^2-m*g*cos(q(3)));
           -m*q(2)*dq(3)^2-m*g*cos(q(3))]; 
%     G = [0 -m*g*cos(q(3)) m*g*q(2)*sin(q(3))]';
end


