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

ex = zeros(numSteps,1); dex = zeros(numSteps,1); ddex = zeros(numSteps,1);
el = zeros(numSteps,1); del = zeros(numSteps,1); ddel = zeros(numSteps,1);
eth = zeros(numSteps,1); deth = zeros(numSteps,1); ddeth = zeros(numSteps,1);

q = zeros(3,numSteps); dq = zeros(3,numSteps); ddq = zeros(3,numSteps);
ux = zeros(numSteps,1); ul = zeros(numSteps,1);

sx = zeros(numSteps,1); sl = zeros(numSteps,1);
ddz = zeros(3,numSteps); f = zeros(3,numSteps); ddU = zeros(3,numSteps);


% Control Parameters
kx = 0.6;           kl = 200;
kas = 1.3;          kli = 300;
eta_x = 0.4;        eta_l = 8;
delta_sx = 0.01;    delta_sl = 0.005;

% Prototype Parameters
g=9.8;
m = 15.6;
mx = 36.2;
ml = 1.7;
M_tot = mx + ml;
dvx = 2.5;
dvl = 0.2;


% Get Travelling and Hoisting Trajectories
rx_wpts = [0 1.5 1.5];
[rx, drx, ddrx] = trapveltraj(rx_wpts, T/Ts + 1, 'EndTime', [12.5 20.0], 'Acceleration', [0.04 0.05]);

rl_wpts = [1.5 0.6 0.6 1.5 1.5];
[rl, drl, ddrl] = trapveltraj(rl_wpts, T/Ts + 1, 'EndTime', [5.5 3.5 4.5 20.0], 'Acceleration', [0.12 0.18 0.18 0.05]);


l(1) = 1.5;

for k=1:numSteps
    ex(k) = rx(k) - x(k);
    el(k) = rl(k) - l(k);
    eth(k) = th(k);
    
    dex(k) = drx(k)-dx(k);
    del(k) = drl(k)-dl(k);
    deth(k) = dth(k);
    
    q(:,k) =  [x(k); l(k); th(k)];
    dq(:,k) = [dx(k); dl(k); dth(k)];
    ddq(:,k) = [ddx(k); ddl(k); ddth(k)];
    
%     M = [mx + m*sin(q(3))^2          m*sin(q(3));
%          m*sin(q(3))                 ml + m];
%     
%     D = [dvx,     0;
%          0 ,      dvl];
%      
%     cg = [-sin(q(3)) * ( (m*q(2)*dq(3)^2) - (m*g*cos(q(3))) ); 
%           -m*q(2)*dq(3)^2 - m*g*cos(q(3))]; 
      
    M = [M_tot+m                 m*sin(q(3,k))  m*q(2,k)*cos(q(3,k));
         m*sin(q(3,k))          m            0;
         m*q(2,k)*cos(q(3,k))   0         m*q(2,k)^2];
    
    D = [0      2*m*cos(q(3,k))*dq(3,k)     -m*q(2,k)*sin(q(3,k))*dq(3,k);
         0      0                       -m*q(2,k)*dq(3,k);
         0      2*m*q(2,k)*dq(3,k)           0];
     
    cg = [0 -m*g*cos(q(3,k)) m*g*q(2,k)*sin(q(3,k))]';
    
    ux(k) = ddrx(k) + kx*dex(k) - kas*deth(k);
    ul(k) = ddrl(k) + kl*del(k) + kli*el(k);
    
    ddz(:, k) = [ux(k); ul(k); 0];
    
    fun = @(x) el(k);
    sx(k) = dex(k) + kx*ex(k) - kas*eth(k);
    sl(k) = del(k) + kl*el(k) + kli*integral(fun, 0, t(k), 'ArrayValued', true);
    
    eta_sgn = [eta_x*sign(sx(k)) eta_l*sign(sl(k)) 0]';
    
    f(:,k) = M*(ddz(:, k) + eta_sgn) + D*dq(:,k) + cg;
    
    ddq(:,k) = inv(M)*(f(:,k) - D*dq(:,k) - cg);
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
    
end

%% Plots

desktop = com.mathworks.mde.desk.MLDesktop.getInstance;
myGroup = desktop.addGroup('Plots');
desktop.setGroupDocked('Plots', 0);
myDim   = java.awt.Dimension(2, 2);   % 3 columns, 2 rows
desktop.setDocumentArrangement('Plots', 2, myDim)
figH    = gobjects(1, 4);
bakWarn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
% Error Plots

clc
figH(1) = figure('WindowStyle', 'docked', 'Name', sprintf('X desired and X actual Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(1)), 'javaframe'), 'GroupName', 'Plots');
grid on
hold on
plot(t,rx)
hold on
plot(t,x(1:end-1))
ylabel('X desired and actual');
xlabel('Time [sec]');
hold on

figH(2) = figure('WindowStyle', 'docked', 'Name', sprintf('l desired and l actual Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(2)), 'javaframe'), 'GroupName', 'Plots');
grid on
hold on
plot(t,rl)
hold on
plot(t,l(1:end-1))
ylabel('l desired and actual');
xlabel('Time [sec]');
hold on

figH(3) = figure('WindowStyle', 'docked', 'Name', sprintf('th desired and th actual Plot'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(3)), 'javaframe'), 'GroupName', 'Plots');
grid on
hold on
plot(t,rth)
hold on
plot(t,th(1:end-1))
ylabel('l desired and actual');
xlabel('Time [sec]');
hold on
    
figH(4) = figure('WindowStyle', 'docked', 'Name', sprintf('Solid and Hoisting forces'), 'NumberTitle', 'off');
drawnow;
pause(0.02);
set(get(handle(figH(4)), 'javaframe'), 'GroupName', 'Plots');
grid on
hold on
plot(t,f(1,:)*10)
hold on
plot(t,f(2,:))
ylabel('l desired and actual');
xlabel('Time [sec]');
hold on



