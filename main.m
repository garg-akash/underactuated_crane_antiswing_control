clc;
clear all;
T = 20;
Ts = 0.001;
t = 0:Ts:T;
numSteps = size(t,2);

disturbance_function = @(L,A) A*(0.5+0.5*(sin(linspace(-pi/2,pi*3/2,L/0.001))));
disturbance_amplitude = 5.0;    % In degrees
disturbance_length = 3.0;       % In seconds
disturbance_start = 6.0;        % In seconds
disturbance = zeros(numSteps,1);
disturbance(disturbance_start/Ts : (disturbance_start+disturbance_length)/Ts-1) = disturbance_function(disturbance_length, disturbance_amplitude);
disturbance = deg2rad(disturbance);

xd = zeros(numSteps,1); dxd = zeros(numSteps,1); ddxd = zeros(numSteps,1);
ld = zeros(numSteps,1); dld = zeros(numSteps,1); ddld = zeros(numSteps,1);
x = zeros(numSteps,1); dx = zeros(numSteps,1); ddx = zeros(numSteps,1);
l = zeros(numSteps,1); dl = zeros(numSteps,1); ddl = zeros(numSteps,1);
th = zeros(numSteps,1); dth = zeros(numSteps,1); ddth = zeros(numSteps,1);

x_dist = zeros(numSteps,1); dx_dist = zeros(numSteps,1); ddx_dist = zeros(numSteps,1);
l_dist = zeros(numSteps,1); dl_dist = zeros(numSteps,1); ddl_dist = zeros(numSteps,1);
th_dist = zeros(numSteps,1); dth_dist = zeros(numSteps,1); ddth_dist = zeros(numSteps,1);

noise_amp_ = 5; %(in degrees)
noise_start = 6; noise_end = 8;
freq = 1/(noise_end-noise_start);
noise_amp = [zeros(size(0:Ts:noise_start,2),1);deg2rad(noise_amp_).*ones(size(noise_start:Ts:noise_end,2)-1,1);zeros(size(noise_end:Ts:T,2)-1,1)];
nn = noise_amp.*sin(2*pi*freq*t');

ex = zeros(numSteps,1); dex = zeros(numSteps,1); ddex = zeros(numSteps,1);
el = zeros(numSteps,1); del = zeros(numSteps,1); ddel = zeros(numSteps,1);
eth = zeros(numSteps,1); deth = zeros(numSteps,1); ddeth = zeros(numSteps,1);

ex_dist = zeros(numSteps,1); dex_dist = zeros(numSteps,1); ddex_dist = zeros(numSteps,1);
el_dist = zeros(numSteps,1); del_dist = zeros(numSteps,1); ddel_dist = zeros(numSteps,1);
eth_dist = zeros(numSteps,1); deth_dist = zeros(numSteps,1); ddeth_dist = zeros(numSteps,1);

q = zeros(3,numSteps); dq = zeros(3,numSteps); ddq = zeros(3,numSteps);
ux = zeros(numSteps,1); ul = zeros(numSteps,1);
U = zeros(3,numSteps); dU = zeros(3,numSteps); ddU = zeros(3,numSteps);

q_dist = zeros(3,numSteps); dq_dist = zeros(3,numSteps); ddq_dist = zeros(3,numSteps);
ux_dist = zeros(numSteps,1); ul_dist = zeros(numSteps,1);
U_dist = zeros(3,numSteps); dU_dist = zeros(3,numSteps); ddU_dist = zeros(3,numSteps);

global M m g
M = 6.5; m=1; g=9.8;
pdx = 0.6; pdl_dash = 0.25;
kvx = 0.4; kvl = 0.4; kax = 0.4; kal = 0.4;
epsilonX = 3.5; epsilonL = 3.5;
psiX = 0.02; psiL = 0.015;
kpx = 240; kdx = 100; kpl = 45; kdl = 10;
lambdaX = 0.1; lambdaL = 0.1; 
phiX_bar = 2; phiL_bar = 5;
%phiX_bar = 0.01; phiL_bar = 0.1; %TODO

Frx_0 = 3.5*4.4; 
krx = -0.5; 
epsilon = 0.01;
dL = 4*6.5;
l(1) = 0.5;
l_dist(1) = 0.5;
pdl = l(1) + pdl_dash;
dTH = 4.5; %TODO

for k=1:numSteps
    xd(k) = pdx/2 + kvx^2*log(cosh(2*kax*t(k)/kvx-epsilonX)/cosh(2*kax*t(k)/kvx-epsilonX-2*kax*pdx/kvx^2))/(4*kax);
    ld(k) = (pdl_dash+2*l(1))/2 + kvl^2*log(cosh(2*kal*t(k)/kvl-epsilonL)/cosh(2*kal*t(k)/kvl-epsilonL-2*kal*pdl_dash/kvl^2))/(4*kal);
    dxd(k) = (kvx*sinh((2*kax*pdx)/kvx^2))/(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2));
    dld(k) = (kvl*sinh((2*kal*pdl_dash)/kvl^2))/(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2));
    ddxd(k) = -(kax*(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)^2 - 2*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2)^2))/(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)^2*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2)^2);
    ddld(k) = -(kal*(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)^2 - 2*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2)^2))/(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)^2*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2)^2);
    
    % Error
    ex(k) = x(k)-xd(k);
    el(k) = l(k)-ld(k);
    eth(k) = th(k);
    
    dex(k) = dx(k)-dxd(k);
    del(k) = dl(k)-dld(k);
    deth(k) = dth(k);

    % Error with disturbance
    ex_dist(k) = x_dist(k)-xd(k);
    el_dist(k) = l_dist(k)-ld(k);
    eth_dist(k) = th_dist(k);
    
    dex_dist(k) = dx_dist(k)-dxd(k);
    del_dist(k) = dl_dist(k)-dld(k);
    deth_dist(k) = dth_dist(k);

    phiX =  min(max(dx(k), -phiX_bar), phiX_bar);
    phiL =  min(max(dl(k), -phiL_bar), phiL_bar);
    phiX_store(k) = phiX;
    phiL_store(k) = phiL;

    phiX_dist =  min(max(dx_dist(k), -phiX_bar), phiX_bar);
    phiL_dist =  min(max(dl_dist(k), -phiL_bar), phiL_bar);
    phiX_store_dist(k) = phiX_dist;
    phiL_store_dist(k) = phiL_dist;

    % Set current q
    q(:,k) = [x(k); l(k); th(k)];
    dq(:,k) = [dx(k); dl(k); dth(k)];
    Frx = (Frx_0 * tanh(dx(k)/epsilon)) - (krx * abs(dx(k)) * dx(k));
    
    % Set current q with disturbance
    q_dist(:,k) = [x_dist(k); l_dist(k); th_dist(k)];
    dq_dist(:,k) = [dx_dist(k); dl_dist(k); dth_dist(k)];
    Frx_dist = (Frx_0 * tanh(dx_dist(k)/epsilon)) - (krx * abs(dx_dist(k)) * dx_dist(k));
    
    [Mc,C,G] = computeMCG(q(:,k),dq(:,k));
    [Mc_dist, C_dist, G_dist] = computeMCG(q_dist(:,k),dq_dist(:,k));

    % Control 
    ux(k) = -kpx*ex(k) - 2*lambdaX*psiX^2*ex(k)/(psiX^2-ex(k)^2)^2 - kdx*dex(k) + (M+m)*ddxd(k) + m*ddld(k)*sin(th(k)) + m*dld(k)*dth(k)*cos(th(k)) + Frx - phiX_bar*tanh(10*dex(k));
    ul(k) = -kpl*el(k) - 2*lambdaL*psiL^2*el(k)/(psiL^2-el(k)^2)^2 - kdl*del(k) + m*ddxd(k)*sin(th(k)) + m*ddld(k) - m*g + dL*dl(k) - phiL_bar*tanh(10*del(k));
    U(:,k) = [ux(k);ul(k);0];

    % Control with disturbance
    ux_dist(k) = -kpx*ex_dist(k) - 2*lambdaX*psiX^2*ex_dist(k)/(psiX^2-ex_dist(k)^2)^2 - kdx*dex_dist(k) + (M+m)*ddxd(k) + m*ddld(k)*sin(th_dist(k)) + m*dld(k)*dth(k)*cos(th_dist(k)) + Frx_dist - phiX_bar*tanh(10*dex_dist(k));
    ul_dist(k) = -kpl*el_dist(k) - 2*lambdaL*psiL^2*el_dist(k)/(psiL^2-el_dist(k)^2)^2 - kdl*del_dist(k) + m*ddxd(k)*sin(th_dist(k)) + m*ddld(k) - m*g + dL*dl(k) - phiL_bar*tanh(10*del_dist(k));
    U_dist(:,k) = [ux_dist(k);ul_dist(k);0];
    
    Fd = [-Frx+phiX; -dL*dl(k)+phiL; -dTH*dth(k)];
    Fd_dist = [-Frx_dist+phiX_dist; -dL*dl_dist(k)+phiL_dist; -dTH*dth_dist(k)];

    ddq(:,k) = inv(Mc)*(U(:,k) + Fd - G - C*dq(:,k));
    dq(:,k+1) = dq(:,k) + ddq(:,k)*Ts; 
    q(:,k+1) = q(:,k) + dq(:,k)*Ts + 0.5*ddq(:,k)*Ts.^2; 

    ddq_dist(:,k) = inv(Mc_dist)*(U_dist(:,k) + Fd_dist - G_dist - C_dist*dq_dist(:,k));
    dq_dist(:,k+1) = dq_dist(:,k) + ddq_dist(:,k)*Ts; 
    q_dist(:,k+1) = q_dist(:,k) + dq_dist(:,k)*Ts + 0.5*ddq_dist(:,k)*Ts.^2;
%     q(3,k+1) = q(3,k+1) + noise_amp(k)*(-sin(2*pi*freq*t(k)));%noise_amp(k).*rand(1);
    
    
    
    if (t(k)>=6 && t(k)<=8)
        q_dist(3,k+1) = noise_amp(k)*(-sin(2*pi*freq*t(k)));%noise_amp(k).*rand(1);
    else
        q_dist(3,k+1) = q_dist(3,k+1);
    end
    
    q(3,k+1) = wrapToPi(q(3,k+1));
    q_dist(3,k+1) = wrapToPi(q_dist(3,k+1));
    
    
    % Update states
    x(k+1) = q(1,k+1);
    l(k+1) = q(2,k+1);
    th(k+1) = q(3,k+1);
    dx(k+1) = dq(1,k+1);
    dl(k+1) = dq(2,k+1);
    dth(k+1) = dq(3,k+1);
    ddx(k) = ddq(1,k);
    ddl(k) = ddq(2,k);
    ddth(k) = ddq(3,k);    

    % Update states with disturbance
    x_dist(k+1) = q_dist(1,k+1);
    l_dist(k+1) = q_dist(2,k+1);
    th_dist(k+1) = q_dist(3,k+1);
    dx_dist(k+1) = dq_dist(1,k+1);
    dl_dist(k+1) = dq_dist(2,k+1);
    dth_dist(k+1) = dq_dist(3,k+1);
    ddx_dist(k) = ddq_dist(1,k);
    ddl_dist(k) = ddq_dist(2,k);
    ddth_dist(k) = ddq_dist(3,k);

end
eth = rad2deg(eth);
eth_dist = rad2deg(eth_dist);

%% Plots
close all

desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
myGroup = desktop.addGroup('Plots');
desktop.setGroupDocked('Plots', 0);
myDim   = java.awt.Dimension(3, 3);   % 3 columns, 2 rows
desktop.setDocumentArrangement('Plots', 2, myDim)
figH    = gobjects(1, 9);
bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
% Error Plots

clc
figH(1) = figure('WindowStyle', 'docked', 'Name', sprintf('X desired and X actual Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(1)), 'javaframe'), 'GroupName', 'Plots');
grid on
hold on
plot(t,xd)
hold on
plot(t,x(1:end-1))
ylabel('X desired and actual');
xlabel('Time [sec]');
saveas(gcf, 'plots/X_desired_actual.png');
hold on

figH(2) = figure('WindowStyle', 'docked', 'Name', sprintf('l desired and l actual Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(2)), 'javaframe'), 'GroupName', 'Plots');
grid on
hold on
plot(t,ld)
hold on
plot(t,l(1:end-1))
ylabel('l desired and actual');
xlabel('Time [sec]');
saveas(gcf, 'plots/l_desired_actual.png');
hold on

figH(3) = figure('WindowStyle', 'docked', 'Name', sprintf('X Error Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(3)), 'javaframe'), 'GroupName', 'Plots');
grid minor;
line(t, ex', 'Color', 'blue', 'LineWidth', 1.5);
hold on
line(t, ex_dist', 'Color', 'red', 'LineWidth', 1.5, 'LineStyle','--');
hold on
plot(t,psiX*ones(size(t)),'--','Color','g', 'LineWidth', 1.5)
hold on
plot(t,-psiX*ones(size(t)),'--','Color','g', 'LineWidth', 1.5)
ylabel('Error in x [meters]');
xlabel('Time [sec]');
ylim([-0.1 0.1])
saveas(gcf, 'plots/X_error.png');
hold on

figH(4) = figure('WindowStyle', 'docked', 'Name', sprintf('L Error Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(4)), 'javaframe'), 'GroupName', 'Plots');
grid minor;
line(t, el', 'Color', 'blue', 'LineWidth', 1.5);
hold on
line(t, el_dist', 'Color', 'red', 'LineWidth', 1.5, 'LineStyle','--');
hold on
plot(t,psiL*ones(size(t)),'--','Color','g', 'LineWidth', 1.5)
hold on
plot(t,-psiL*ones(size(t)),'--','Color','g', 'LineWidth', 1.5)
ylabel('Error in l [meters]');
xlabel('Time [sec]');
ylim([-0.04 0.04])
saveas(gcf, 'plots/L_error.png');
hold on

figH(5) = figure('WindowStyle', 'docked', 'Name', sprintf('Theta Error Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(5)), 'javaframe'), 'GroupName', 'Plots');
grid minor;
plot(t, eth', 'b', 'LineWidth', 1.5);
hold on
plot(t, eth_dist', '-.r', 'LineWidth', 1.5);
ylabel('Error in theta [deg]');
xlabel('Time [sec]');
ylim([-5 5])
saveas(gcf, 'plots/Theta_error.png');
hold on

figH(6) = figure('WindowStyle', 'docked', 'Name', sprintf('Ux Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(6)), 'javaframe'), 'GroupName', 'Plots');
grid minor;
line(t, ux', 'Color', 'blue', 'LineWidth', 1.5);
hold on
line(t, ux_dist', 'Color', 'red', 'LineWidth', 1.5, 'LineStyle','--');
ylabel('Ux(t) [N]');
xlabel('Time [sec]');
ylim([-5 20])
saveas(gcf, 'plots/Ux.png');
hold on

figH(7) = figure('WindowStyle', 'docked', 'Name', sprintf('Ul Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(7)), 'javaframe'), 'GroupName', 'Plots');
grid minor;
line(t, ul', 'Color', 'blue', 'LineWidth', 1.5);
hold on
line(t, ul_dist', 'Color', 'red', 'LineWidth', 1.5, 'LineStyle','--');
ylabel('Ul(t) [N]');
xlabel('Time [sec]');
ylim([-12 -4])
saveas(gcf, 'plots/Ul.png');
hold on

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


