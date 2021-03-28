clc
clear all

T = 20;
Ts = 0.001;
t = 0:Ts:T;

numSteps = size(t,2);

xd = zeros(numSteps,1); 
ld = zeros(numSteps,1); 

pdx = 0.6;

kvx = 0.4; kvl = 0.4;
kax = 0.4; kal = 0.4;

epsilonX = 3.5; epsilonL = 3.5;

pdl_dash = 0.25;

l(1) = 0.5;

for k=1:numSteps
    term1 = pdx / 2;
    term2 = (kvx^2) / (4*kax);
    term3 = (2*kax*t(k)) / kvx;
    term4 = (2*kax*pdx) / (kvx^2);
    xd(k) = term1 + (term2 * log(cosh(term3 - epsilonX) / cosh(term3 - epsilonX - term4)));    

    term1 = (pdl_dash + (2*l(1)) ) / 2;
    term2 = (kvl^2) / (4*kal);
    term3 = (2*kal*t(k)) / kvl;
    term4 = (2*kal*pdl_dash) / (kvl^2);
    ld(k) = term1 + (term2 * log(cosh(term3 - epsilonL) / cosh(term3 - epsilonL - term4)));    
    
end
