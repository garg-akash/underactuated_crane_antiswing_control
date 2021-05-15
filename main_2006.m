clc;
clear all;
T = 20;
Ts = 0.001;
t = 0:Ts:T;
numSteps = size(t,2);

disturbance_function = @(L,A) A*(0.5+0.5*(sin(linspace(-pi/2,pi*3/2,L/0.001))));
disturbance_amplitude = 4.0;    % In degrees
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

noise_amp_ = 5*2; %(in degrees)
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
l_dist(1) = 0.5;
pdl = l(1) + pdl_dash;
dTH = 0.0; %TODO

% TODO: These params are from previous paper. Currently set to zero. Has to
% be tuned
kpx = 4.5;          kpl = 2;
kas = 4.6;          kli = 4;
eta_x = 2;          eta_l = 2;
dX= 20.0;          dL = 0.2;          
delta_sx = 0.1;    delta_sl = 0.005;

el_int = 0.02;
el_int_dist = 0.02;

for k=1:numSteps
    xd(k) = pdx/2 + kvx^2*log(cosh(2*kax*t(k)/kvx-epsilonX)/cosh(2*kax*t(k)/kvx-epsilonX-2*kax*pdx/kvx^2))/(4*kax);
    ld(k) = (pdl_dash+2*l(1))/2 + kvl^2*log(cosh(2*kal*t(k)/kvl-epsilonL)/cosh(2*kal*t(k)/kvl-epsilonL-2*kal*pdl_dash/kvl^2))/(4*kal);
    dxd(k) = (kvx*sinh((2*kax*pdx)/kvx^2))/(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2));
    dld(k) = (kvl*sinh((2*kal*pdl_dash)/kvl^2))/(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2));
    ddxd(k) = -(kax*(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)^2 - 2*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2)^2))/(2*cosh((epsilonX*kvx - 2*kax*t(k))/kvx)^2*cosh((epsilonX*kvx^2 - 2*kax*t(k)*kvx + 2*kax*pdx)/kvx^2)^2);
    ddld(k) = -(kal*(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)^2 - 2*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2)^2))/(2*cosh((epsilonL*kvl - 2*kal*t(k))/kvl)^2*cosh((epsilonL*kvl^2 - 2*kal*t(k)*kvl + 2*kal*pdl_dash)/kvl^2)^2);

    ex(k) = -x(k)+xd(k);
    el(k) = -l(k)+ld(k);
    eth(k) = -th(k);
    
    dex(k) = -dx(k)+dxd(k);
    del(k) = -dl(k)+dld(k);
    deth(k) = -dth(k);

    ex_dist(k) = -x_dist(k)+xd(k);
    el_dist(k) = -l_dist(k)+ld(k);
    eth_dist(k) = -th_dist(k);
    
    dex_dist(k) = -dx_dist(k)+dxd(k);
    del_dist(k) = -dl_dist(k)+dld(k);
    deth_dist(k) = -dth_dist(k);
    
    q(:,k) = [x(k);l(k);th(k)];
    dq(:,k) = [dx(k);dl(k);dth(k)];

    q_dist(:,k) = [x_dist(k); l_dist(k); th_dist(k)];
    dq_dist(:,k) = [dx_dist(k); dl_dist(k); dth_dist(k)];
    
    [Mc,D,Cg] = computeMCG(q(:,k),dq(:,k),dX,dL);
    [Mc_dist,D_dist,Cg_dist] = computeMCG(q_dist(:,k),dq_dist(:,k),dX,dL);

    ux(k) = ddxd(k) + kpx*dex(k) - kas*deth(k);
    ul(k) = ddld(k) + kpl*del(k) + kli*el(k);

    ux_dist(k) = ddxd(k) + kpx*dex_dist(k) - kas*deth_dist(k);
    ul_dist(k) = ddld(k) + kpl*del_dist(k) + kli*el_dist(k);

    ddz(:,k) = [ux(k);ul(k)];
    el_int = el_int + el(k)*Ts;
    fun = @(x) el(k);
    sx(k) = dex(k) + kpx*ex(k) - kas*eth(k);
    sl(k) = del(k) + kpl*el(k) + kli*el_int;

    ddz_dist(:,k) = [ux_dist(k);ul_dist(k)];
    el_int_dist = el_int_dist + el_dist(k)*Ts;
    fun_dist = @(x) el_dist(k);
    sx_dist(k) = dex_dist(k) + kpx*ex_dist(k) - kas*eth_dist(k);
    sl_dist(k) = del_dist(k) + kpl*el_dist(k) + kli*el_int_dist;
    
    if sx(k)/delta_sx <= 1
        sat_sx = sx(k)/delta_sx;
    else
        sat_sx = sign(sx(k)/delta_sx);
    end
    if sl(k)/delta_sl <= 1
        sat_sl = sl(k)/delta_sl;
    else
        sat_sl = sign(sl(k)/delta_sl);
    end

    if sx_dist(k)/delta_sx <= 1
        sat_sx_dist = sx_dist(k)/delta_sx;
    else
        sat_sx_dist = sign(sx_dist(k)/delta_sx);
    end
    if sl_dist(k)/delta_sl <= 1
        sat_sl_dist = sl_dist(k)/delta_sl;
    else
        sat_sl_dist = sign(sl_dist(k)/delta_sl);
    end

    eta_sgn = [eta_x*sat_sx eta_l*sat_sl]';
    eta_sgn_dist = [eta_x*sat_sx_dist eta_l*sat_sl_dist]';
    
    U(1:2,k) = Mc*(ddz(:, k) + eta_sgn) + D*dq(1:2,k) + Cg;
    U_dist(1:2,k) = Mc_dist*(ddz_dist(:, k) + eta_sgn) + D_dist*dq_dist(1:2,k) + Cg_dist;

    %Fd = [-Frx+phiX; -dL*dl(k)+phiL; -dTH*dth(k)];
    
    ddq(1:2,k) = inv(Mc)*(U(1:2,k) - Cg - D*dq(1:2,k));
    ddq(3,k) = -(cos(q(3,k))*ddq(1,k) + 2*dq(2,k)*dq(3,k) + g*sin(dq(3,k)))/(q(2,k));
    dq(:,k+1) = dq(:,k) + ddq(:,k)*Ts; 
    q(:,k+1) = q(:,k) + dq(:,k)*Ts + 0.5*ddq(:,k)*Ts.^2; 

    ddq_dist(1:2,k) = inv(Mc_dist)*(U_dist(1:2,k) - Cg_dist - D_dist*dq_dist(1:2,k));
    ddq_dist(3,k) = -(cos(q_dist(3,k))*ddq_dist(1,k) + 2*dq_dist(2,k)*dq_dist(3,k) + g*sin(dq_dist(3,k)))/(q_dist(2,k));
    dq_dist(:,k+1) = dq_dist(:,k) + ddq_dist(:,k)*Ts; 
    q_dist(:,k+1) = q_dist(:,k) + dq_dist(:,k)*Ts + 0.5*ddq_dist(:,k)*Ts.^2; 
    

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

%% Plot
close all

desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
myGroup = desktop.addGroup('Plots');
desktop.setGroupDocked('Plots', 0);
myDim   = java.awt.Dimension(4, 4);   % 4 columns, 4 rows
desktop.setDocumentArrangement('Plots', 2, myDim)
figH    = gobjects(1, 8);
bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

clc
figH(1) = figure('WindowStyle', 'docked', 'Name', sprintf('X desired and actual'), 'NumberTitle', 'off');
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
saveas(gcf, 'plots_2006/x_desired_actual_2006.png');
hold on

figH(2) = figure('WindowStyle', 'docked', 'Name', sprintf('l desired and actual'), 'NumberTitle', 'off');
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
saveas(gcf, 'plots_2006/l_desired_actual_2006.png');
hold on

figH(3) = figure('WindowStyle', 'docked', 'Name', sprintf('Error in x [meters]'), 'NumberTitle', 'off');
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
saveas(gcf, 'plots_2006/X_error_2006.png');
hold on

figH(4) = figure('WindowStyle', 'docked', 'Name', sprintf('Error in l [meters]'), 'NumberTitle', 'off');
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
saveas(gcf, 'plots_2006/L_error_2006.png');
hold on

figH(5) = figure('WindowStyle', 'docked', 'Name', sprintf('Error in theta [deg]'), 'NumberTitle', 'off');
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
saveas(gcf, 'plots_2006/Theta_error_2006.png');
hold on

figH(6) = figure('WindowStyle', 'docked', 'Name', sprintf('Ux [N]'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(6)), 'javaframe'), 'GroupName', 'Plots');
line(t, U(1,:)', 'Color', 'blue', 'LineWidth', 1.5);
hold on
line(t, U_dist(1,:)', 'Color', 'red', 'LineWidth', 1, 'LineStyle','--');
ylabel('Ux(t) [N]');
xlabel('Time [sec]');
ylim([-5 20])
saveas(gcf, 'plots_2006/Ux_2006.png');
hold on

figH(7) = figure('WindowStyle', 'docked', 'Name', sprintf('Ul [N]'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(7)), 'javaframe'), 'GroupName', 'Plots');
grid minor;
line(t, U(2,:)', 'Color', 'blue', 'LineWidth', 1.5);
hold on
line(t, U_dist(2,:)', 'Color', 'red', 'LineWidth', 1, 'LineStyle','--');
ylabel('Ul(t) [N]');
xlabel('Time [sec]');
ylim([-12 -4])
saveas(gcf, 'plots_2006/Ul_2006.png');
hold on

figH(8) = figure('WindowStyle', 'docked', 'Name', sprintf('theta'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(8)), 'javaframe'), 'GroupName', 'Plots');
grid on
hold on
plot(t,disturbance)
hold on
ylabel('theta');
xlabel('Time [sec]');
hold on


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


