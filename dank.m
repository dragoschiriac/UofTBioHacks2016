%ODE Solver for interfering mRNA & Pseudo-Promoter Aptamer 

clc;
close all;
Ktx1 = 83;
Ktx2 = 83;
Kdm = 0.03;
Kdp = 0.03;
Ka1 = 1.02;
Kd1 = 0.04;
Ka2 = 1.04;  
Kd2 = 0.04;
g = 5;
p =0.01;
a = 0.5;
Ktr = 2;

Np1 = 5;
Np2 = 30;





f = @(t,x) [Ktx1-Kdm*x(1)-g*x(2);
    Ka1*x(4)*x(1)-Kd1*x(2);
    Ka2*x(4)-Kd2*x(3);
    Ktx1 - Kdm*x(4)-g*x(2)-p*x(3);
    Ktx2*a*(x(4)/(0.04+x(4)))-Kdm*x(5);
    Ktr*x(5) - Kdp*x(6)];

[t,xa] = ode45(f,[0 250],[0 0 0 0 0 0]);

plot(t,xa)
%Fix the titles later
title('Model of RNA Aptamer Regulatory Element')
xlabel('Time (Min)'), ylabel('Concentration (nM)')