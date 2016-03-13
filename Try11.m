%ODE Solver for simple Pseudo-Promoter Aptamer 
clc;
close all;
Ktx1 = 83;
Ktx2 = 83;
Kdm = 0.03;
Kdp = 0.03;
Ka1 = 0.04;
Kd1 = 0.02;
Ka2 = 1.04; 
Kd2 = 0.04;
g = 5;
p = 0.01;
a = 0.5;
Ktr = 2;

Np1 = 5;
Np2 = 30;





f = @(t,x) [Ktx1 - Kdm*x(1)-p*x(2);
    Ka2*x(1)-Kd2*x(2);
    Ktx2*a*(x(1)/(0.04+x(1)))-Kdm*x(3);
    Ktr*x(3) - Kdp*x(4)];

[t,xa] = ode45(f,[0 250],[0 0 0 0]);

plot(t,xa)
%Fix the titles later
title('Model of RNA Aptamer Regulatory Element')
xlabel('Time (Min)'), ylabel('Concentration (nM)')