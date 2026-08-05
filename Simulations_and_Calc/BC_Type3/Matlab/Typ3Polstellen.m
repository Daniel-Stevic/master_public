clc; clear; close all;
%% ── Symbolische Polberechnung des Nennerpolynoms 4. Ordnung ─────────────────
syms s L1 L2 C1 C2 R D0 real

%% ── Koeffizienten des Nennerpolynoms ────────────────────────────────────────
% a4*s^4 + a3*s^3 + a2*s^2 + a1*s + a0

a4 = 1;
a3 = 2 / (C2*R);
a2 = (1-D0)^2/(C1*L2) + (1-D0)^2/(C2*L1) + D0^2/(C1*L1) + D0^2/(C2*L2);
a1 = 2*(1-D0)^2/(C1*C2*L2*R) + 2*D0^2/(C1*C2*L1*R);
a0 = (1-2*D0)^2 / (C1*C2*L1*L2);

%% ── Polynom aufstellen und Nullstellen (Pole) berechnen ─────────────────────
p = a4*s^4 + a3*s^3 + a2*s^2 + a1*s + a0;

disp('════════════════════════════════════════════════════')
disp('  Nennerpolynom')
disp('════════════════════════════════════════════════════')
pretty(collect(expand(p), s))

%% ── Symbolische Polberechnung via solve ─────────────────────────────────────
disp(' ')
disp('════════════════════════════════════════════════════')
disp('  Symbolische Pole (solve)')
disp('════════════════════════════════════════════════════')

poles_sym = solve(p == 0, s);

for k = 1:length(poles_sym)
    fprintf('\ns_%d =\n', k)
    pretty(simplify(poles_sym(k)))
    fprintf('LaTeX: %s\n', latex(simplify(poles_sym(k))))
end

%% ── Numerische Probe mit Arbeitspunktwerten ─────────────────────────────────
L1_val = 47e-6;
L2_val = 47e-6;
C1_val = 330e-6;
C2_val = 330e-6;
R_val  = 20;
D0_val = 0.25;

subs_vars = [L1,    L2,    C1,    C2,    R,     D0   ];
subs_vals = [L1_val,L2_val,C1_val,C2_val,R_val, D0_val];

disp(' ')
disp('════════════════════════════════════════════════════')
disp('  Numerische Pole am Arbeitspunkt')
disp('════════════════════════════════════════════════════')

% Koeffizientenvektor numerisch
a4_n = double(subs(a4, subs_vars, subs_vals));
a3_n = double(subs(a3, subs_vars, subs_vals));
a2_n = double(subs(a2, subs_vars, subs_vals));
a1_n = double(subs(a1, subs_vars, subs_vals));
a0_n = double(subs(a0, subs_vars, subs_vals));

poles_num = roots([a4_n, a3_n, a2_n, a1_n, a0_n]);

for k = 1:length(poles_num)
    fprintf('s_%d = %.4f + %.4fj\n', k, real(poles_num(k)), imag(poles_num(k)))
end

%% ── Stabilitätsprüfung (Hurwitz) ────────────────────────────────────────────
disp(' ')
disp('════════════════════════════════════════════════════')
disp('  Hurwitz-Stabilitätsbedingungen (symbolisch)')
disp('════════════════════════════════════════════════════')

% Für n=4: alle Koeffizienten > 0  UND  a3*a2 > a4*a1  UND  a1*(a3*a2-a4*a1) > a0*a3^2
H1 = simplify(a3);
H2 = simplify(a3*a2 - a4*a1);
H3 = simplify(a1*(a3*a2 - a4*a1) - a0*a3^2);

disp('H1 = a3 > 0:')
pretty(H1)

disp('H2 = a3*a2 - a4*a1 > 0:')
pretty(simplify(H2))

disp('H3 = a1*H2 - a0*a3^2 > 0:')
pretty(simplify(H3))

% Numerische Auswertung der Hurwitz-Bedingungen
H1_n = double(subs(H1, subs_vars, subs_vals));
H2_n = double(subs(H2, subs_vars, subs_vals));
H3_n = double(subs(H3, subs_vars, subs_vals));

fprintf('\nNumerisch am Arbeitspunkt:\n')
fprintf('H1 = %.6f  →  %s\n', H1_n, tern(H1_n>0,'stabil','INSTABIL'))
fprintf('H2 = %.6f  →  %s\n', H2_n, tern(H2_n>0,'stabil','INSTABIL'))
fprintf('H3 = %.6f  →  %s\n', H3_n, tern(H3_n>0,'stabil','INSTABIL'))

%% ── Hilfsfunktion ───────────────────────────────────────────────────────────
function out = tern(cond, a, b)
    if cond; out = a; else; out = b; end
end




%% ── Toleranz- und Parameteranalyse der Pollagen ────────────────────────────
% Untersucht den Einfluss von:
%   - Bauteiltoleranz: L1, L2, C1, C2  je +/-20%
%   - Lastwiderstand:  R   von 10 bis 50 Ohm
%   - Tastverhältnis:  D0  von 0.05 bis 0.45

%% ── Nominalwerte ────────────────────────────────────────────────────────────
L1_nom = 47e-6;
L2_nom = 47e-6;
C1_nom = 330e-6;
C2_nom = 330e-6;
R_nom  = 20;
D0_nom = 0.25;

tol    = 0.20;   % +/- 20%
N_mc   = 2000;   % Monte-Carlo Stichproben für Toleranzanalyse
N_sw   = 200;    % Punkte für R- und D0-Sweep

%% ── Hilfsfunktion: Pole aus Koeffizienten ───────────────────────────────────
function p = get_poles(L1, L2, C1, C2, R, D0)
    a4 =  1;
    a3 =  2 / (C2*R);
    a2 =  (1-D0)^2/(C1*L2) + (1-D0)^2/(C2*L1) ...
        + D0^2/(C1*L1)      + D0^2/(C2*L2);
    a1 =  2*(1-D0)^2/(C1*C2*L2*R) + 2*D0^2/(C1*C2*L1*R);
    a0 =  (1-2*D0)^2 / (C1*C2*L1*L2);
    p  = roots([a4, a3, a2, a1, a0]);
end

%% ════════════════════════════════════════════════════════════════════════════
%% TEIL 1 – Monte-Carlo Toleranzanalyse (L1,L2,C1,C2 je +/-20%)
%% ════════════════════════════════════════════════════════════════════════════
poles_mc = zeros(N_mc, 4);

for k = 1:N_mc
    L1_k = L1_nom * (1 + tol*(2*rand-1));
    L2_k = L2_nom * (1 + tol*(2*rand-1));
    C1_k = C1_nom * (1 + tol*(2*rand-1));
    C2_k = C2_nom * (1 + tol*(2*rand-1));
    p_k  = get_poles(L1_k, L2_k, C1_k, C2_k, R_nom, D0_nom);
    poles_mc(k,:) = sort(p_k.').';   % nach Imaginärteil sortieren
end

% Nominalpole
p_nom = get_poles(L1_nom, L2_nom, C1_nom, C2_nom, R_nom, D0_nom);

figure('Name','Monte-Carlo Toleranzanalyse','Position',[100 100 1200 500]);

subplot(1,2,1)
scatter(real(poles_mc(:)), imag(poles_mc(:)), 4, ...
    repmat([0.6 0.6 0.6],N_mc*4,1), 'filled'); hold on
scatter(real(p_nom), imag(p_nom), 80, 'r', 'x', 'LineWidth', 2)
xline(0,'k--','LineWidth',1)
grid on; xlabel('Re(s)'); ylabel('Im(s)')
title(sprintf('Polwolke: L,C je \\pm%d%% (N=%d)', tol*100, N_mc))
legend('Monte-Carlo','Nominalpole','Location','best')

subplot(1,2,2)
histogram(real(poles_mc(:)), 60, 'FaceColor',[0.2 0.4 0.8],'EdgeColor','none')
xline(0,'r--','LineWidth',1.5)
grid on; xlabel('Re(s_k)'); ylabel('Häufigkeit')
title('Verteilung Realteil aller Pole')
sgtitle('Toleranzanalyse ±20% (L1,L2,C1,C2)')

%% ════════════════════════════════════════════════════════════════════════════
%% TEIL 2 – R-Sweep (10 … 50 Ohm), nominale L,C,D0
%% ════════════════════════════════════════════════════════════════════════════
R_vec      = linspace(10, 50, N_sw);
poles_R    = zeros(N_sw, 4);

for k = 1:N_sw
    p_k = get_poles(L1_nom, L2_nom, C1_nom, C2_nom, R_vec(k), D0_nom);
    poles_R(k,:) = sort(p_k.').';
end

figure('Name','R-Sweep','Position',[100 100 1200 800]);

subplot(2,2,1)
plot(R_vec, real(poles_R), 'LineWidth', 1.5); grid on
xlabel('R [\Omega]'); ylabel('Re(s_k)')
title('Realteil der Pole vs R')
legend('s_1','s_2','s_3','s_4','Location','best')
yline(0,'k--')

subplot(2,2,2)
plot(R_vec, imag(poles_R), 'LineWidth', 1.5); grid on
xlabel('R [\Omega]'); ylabel('Im(s_k)')
title('Imaginärteil der Pole vs R')
legend('s_1','s_2','s_3','s_4','Location','best')

subplot(2,2,[3 4])
for k = 1:4
    plot(real(poles_R(:,k)), imag(poles_R(:,k)), 'LineWidth', 1.5); hold on
end
scatter(real(poles_R(1,:)),   imag(poles_R(1,:)),   60,'gs','filled')   % R=10
scatter(real(poles_R(end,:)), imag(poles_R(end,:)), 60,'r^','filled')   % R=50
xline(0,'k--'); grid on
xlabel('Re(s)'); ylabel('Im(s)')
title('Wurzelortskurve: R = 10 \rightarrow 50 \Omega')
legend('s_1','s_2','s_3','s_4','R=10\Omega','R=50\Omega','Location','best')
sgtitle('Parametervariation Lastwiderstand R')

%% ════════════════════════════════════════════════════════════════════════════
%% TEIL 3 – D0-Sweep (0.05 … 0.45), nominale L,C,R
%% ════════════════════════════════════════════════════════════════════════════
D0_vec     = linspace(0.05, 0.45, N_sw);
poles_D0   = zeros(N_sw, 4);

for k = 1:N_sw
    p_k = get_poles(L1_nom, L2_nom, C1_nom, C2_nom, R_nom, D0_vec(k));
    poles_D0(k,:) = sort(p_k.').';
end

figure('Name','D0-Sweep','Position',[100 100 1200 800]);

subplot(2,2,1)
plot(D0_vec, real(poles_D0), 'LineWidth', 1.5); grid on
xlabel('D_0'); ylabel('Re(s_k)')
title('Realteil der Pole vs D_0')
legend('s_1','s_2','s_3','s_4','Location','best')
yline(0,'k--')

subplot(2,2,2)
plot(D0_vec, imag(poles_D0), 'LineWidth', 1.5); grid on
xlabel('D_0'); ylabel('Im(s_k)')
title('Imaginärteil der Pole vs D_0')
legend('s_1','s_2','s_3','s_4','Location','best')

subplot(2,2,[3 4])
for k = 1:4
    plot(real(poles_D0(:,k)), imag(poles_D0(:,k)), 'LineWidth', 1.5); hold on
end
scatter(real(poles_D0(1,:)),   imag(poles_D0(1,:)),   60,'gs','filled')   % D0=0.05
scatter(real(poles_D0(end,:)), imag(poles_D0(end,:)), 60,'r^','filled')   % D0=0.45
xline(0,'k--'); grid on
xlabel('Re(s)'); ylabel('Im(s)')
title('Wurzelortskurve: D_0 = 0.05 \rightarrow 0.45')
legend('s_1','s_2','s_3','s_4','D_0=0.05','D_0=0.45','Location','best')
sgtitle('Parametervariation Tastverhältnis D_0')

%% ════════════════════════════════════════════════════════════════════════════
%% TEIL 4 – Stabilitätskarte: R × D0 Gitter
%% ════════════════════════════════════════════════════════════════════════════
R_grid  = linspace(10, 50,   80);
D0_grid = linspace(0.05, 0.45, 80);
stabil  = zeros(length(D0_grid), length(R_grid));

for ir = 1:length(R_grid)
    for id = 1:length(D0_grid)
        p_k = get_poles(L1_nom, L2_nom, C1_nom, C2_nom, ...
                        R_grid(ir), D0_grid(id));
        stabil(id,ir) = all(real(p_k) < 0);
    end
end

figure('Name','Stabilitätskarte','Position',[100 100 700 500]);
imagesc(R_grid, D0_grid, stabil)
set(gca,'YDir','normal')
colormap([0.85 0.2 0.2; 0.2 0.75 0.2])   % rot=instabil, grün=stabil
colorbar('Ticks',[0.25 0.75],'TickLabels',{'instabil','stabil'})
xlabel('R [\Omega]'); ylabel('D_0')
title('Stabilitätskarte: nominale L,C — alle Pole Re(s)<0 ?')
grid on
sgtitle('Stabilitätsgebiet im R-D_0-Parameterraum')

%% ── Konsolenausgabe Nominalpole ─────────────────────────────────────────────
fprintf('\n════════════════════════════════════════════════════\n')
fprintf('  Nominalpole (L,C,R,D0 nominal)\n')
fprintf('════════════════════════════════════════════════════\n')
for k = 1:4
    fprintf('  s_%d = %+.2f %+.2f j\n', k, real(p_nom(k)), imag(p_nom(k)))
end
fprintf('\n  System ist: %s\n', ternary(all(real(p_nom)<0),'STABIL','INSTABIL'))

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end