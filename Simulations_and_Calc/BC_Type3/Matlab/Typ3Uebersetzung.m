clc; clear;
syms L T R U1 d1 Uc U2 real

% DLB Übersetzung

eqn = 2*L*Uc^2-(3*U1*L+d1^2*T*U1*R)*Uc+U1^2*L == 0;
Uc_sol = solve(eqn,Uc);

%Uc=(3*U1*L+(d1^2)*T*U1*R+sqrt( (3*U1*L+(d1^2)*T*U1*R)^2 -8*U1^2*L^2))/(4*L)
%pretty(simplify(Uc))

% Select one solution branch (usually the '+' root)
Uc = Uc_sol(1);

% Numeric values
L=47e-6;
R=500;
T=10e-6;

U2 = 2*Uc-U1;

% latex Ausgabe
DLB_Uebersetzung = simplify(U2/U1)
latex(DLB_Uebersetzung)

% Transfer ratio
ratio = simplify(U2/U1);

disp('U2/U1 =')
pretty(ratio)

% Convert symbolic expression to MATLAB function
f = matlabFunction(subs(ratio), 'Vars', d1);

% d1 range
d1_vals = linspace(0, 0.45, 1000);

% Evaluate
y = f(d1_vals);

% Plot
figure;
plot(d1_vals, y, 'LineWidth', 2);
grid on;

xlabel('d1');
ylabel('U2/U1');
title('Übersetzung U2/U1 über d1');

