%Dank ODE Solver 

Ktx = 0;
Kdm = 0;
Kdp = 0;
w = 0;
Ka1 = 0;
Kd1 = 0;
Ka2 = 0;
Kd2 = 0;
g = 0;
p = 0;
a = 0;
Ktr = 0;

Np1 = 0;
Np2 = 0;

Koc = 0;
Rpol = 0;




f = @(t,x) [Ktx-Kdm*x(1)-w*x(2);
    Ka1*x(4)*x(1)-Kd1*x(2);
    Ka2*x(4)*Np2-Kd2*x(3);
    Ktx - Kdm*x(4)-g*x(2)-p*x(3);
    Koc*Np2*a*x(3)-Kdm*x(5);
    Ktr*x(5) - Kdp*x(6)];

[t,xa] = ode45(f,[0 100],[0 1/2 3]);

plot(t,xa)
%Fix the titles later
title('y(t)')
xlabel('t'), ylabel('y')