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
ddz = zeros(2,numSteps); f = zeros(2,numSteps); ddU = zeros(3,numSteps);


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
dvx = 2.5;
dvl = 0.2;


% Get Travelling and Hoisting Trajectories
rx_wpts = [0 1.5 1.5];
[rx, drx, ddrx] = trapveltraj(rx_wpts, T/Ts + 1, 'EndTime', [12.5 20.0], 'Acceleration', [0.04 0.05]);

rl_wpts = [1.5 0.6 0.6 1.5 1.5];
[rl, drl, ddrl] = trapveltraj(rl_wpts, T/Ts + 1, 'EndTime', [5.5 3.5 4.5 20.0], 'Acceleration', [0.12 0.18 0.18 0.05]);

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
    
    M = [mx + m*sin(q(3))^2          m*sin(q(3));
         m*sin(q(3))                 ml + m];
    
    D = [dvx,     0;
         0 ,      dvl];
     
    cg = [-sin(q(3)) * ( (m*q(2)*dq(3)^2) - (m*g*cos(q(3))) ); 
          -m*q(2)*dq(3)^2 - m*g*cos(q(3))]; 
    
    ux(k) = ddrx(k) + kx*dex(k) - kas*deth(k);
    ul(k) = ddrl(k) + kl*del(k) + kli*el(k);
    
    ddz(:, k) = [ux(k); ul(k)];
    
    fun = @(x) el(k);
    sx(k) = dex(k) + kx*ex(k) - kas*eth(k);
    sl(k) = del(k) + kl*el(k) + kli*integral(fun, 0, t(k), 'ArrayValued', true);
    
    eta_sgn = [eta_x*sign(sx(k)) eta_l*sign(sl(k))]';
    
    f(:,k) = M*(ddz(:, k) + eta_sgn) + D*dq(1:2,k) + cg;
    
    ddq(1:2,k) = inv(M)*(f(:,k) - D*dq(1:2,k) - cg);
    dq(1:2,k+1) = dq(1:2,k) + ddq(1:2,k)*Ts; 
    q(1:2,k+1) = q(1:2,k) + dq(1:2,k)*Ts + 0.5*ddq(1:2,k)*Ts.^2; 
    
    % TODO: Theta not in the dynamic model
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

grid on
hold on
plot(t,rl)
ylabel('X desired');
xlabel('Time [sec]');
hold on
