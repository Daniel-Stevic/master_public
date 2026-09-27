clc; clear;
syms L T R U1 d1 d2 Uc2 U2 real
% mit fallender Flanke delta I gerechnet
eqn1 = (Uc2-U1)*d2^2+(3*U1*d1-3*Uc2*d1)*d2+2*Uc2*(d1^2)-2*U1*(d1^2)-(L*(2*Uc2-U1)/(R*T))== 0;
d2_sol1 = solve(eqn1,d2)

d2_pos1 = (d2_sol1(1))

%mit steigender Flanke delta I gerechnet

eqn2 = (2*Uc2-U1)/R==Uc2*d1*T/L*(d2-2*d1);
d2_sol2 = solve(eqn2,d2)

% Nummerische Werte aus Simulation
L=47e-6;
R=1000;
T=10e-6;
d1=0.053;
U1=12;
Uc2=17.673;

double(subs(d2_pos1))

double(subs(d2_sol2))