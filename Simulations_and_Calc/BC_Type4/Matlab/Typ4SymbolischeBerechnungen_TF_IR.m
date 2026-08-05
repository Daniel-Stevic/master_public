clc; clear; close all;
syms s L1 L2 C1 C2 D0 R real
syms UC10 UC20 U10 IL10 IL20 real

%% -- Systemmatrizen ----------------------------------------------------------
% Zustandsmatrix A (zeitbasiert)
A_mat = [ 0,          0,          D0/L1,      (D0-1)/L1;
          0,          0,         (D0-1)/L2,    D0/L2;
         -D0/C1,      (1-D0)/C1,  0,           0;
          (1-D0)/C2, -D0/C2,      0,           0      ];

% Eingangsmatrix B  (Spalten: U1, IR, D)
B_mat = [ 1/L1,   0,        (UC10+UC20)/L1;
          0,   0,        (UC10+UC20)/L2;
          0,           0,       -(IL10+IL20)/C1;
          0,          -1/C2,    -(IL10+IL20)/C2      ];

% Ausgangsgleichung: u2 = 2*UC2 - u1
% u2 = C_out*x + D_out*[u1; iR]
% C_out = [0 0 0 2],  D_out = [-1 0]

%% -- Cramersche Regel --------------------------------------------------------
function TF = cramer(M, B, zustands_idx, eingangs_idx)
    M_mod = M;
    M_mod(:, zustands_idx) = B(:, eingangs_idx);
    TF = simplify(det(M_mod) / det(M));
end

% (sI - A) für offenes System
M_open = s*eye(4) - A_mat;

% Zustandsbezeichnungen für Ausgabe
zustands_namen = {'IL1', 'IL2', 'UC1', 'UC2'};
eingangs_namen = {'U1', 'IR', 'D'};

%% -- Alle offenen TFs berechnen (4 Zustände × 3 Eingänge) -------------------
disp('====================================================')
disp('===         OFFENE ÜBERTRAGUNGSFUNKTIONEN         ===')
disp('====================================================')

G = sym(zeros(4,3));
for i = 1:4
    for j = 1:3
        G(i,j) = cramer(M_open, B_mat, i, j);
        fprintf('\n=== %s / %s (offen) ===\n', zustands_namen{i}, eingangs_namen{j})
        [num, den] = numden(G(i,j));
        pretty(collect(expand(num), s))
        disp('---------------------------------------------')
        pretty(collect(expand(den), s))
    end
end

%% -- Lastschluss algebraisch: iR = (2*UC2 - U1) / R ------------------------
%
% Offenes System:
% Xi = G(i,1)*U1 + G(i,2)*IR + G(i,3)*D
%
% Last:
% IR = (2*UC2 - U1)/R
%
% Zuerst UC2 schließen:
%
% UC2/U1 = (G(4,1) - G(4,2)/R) / (1 - 2*G(4,2)/R)
% UC2/D  =  G(4,3)             / (1 - 2*G(4,2)/R)
%
% Danach alle übrigen Zustände:
%
% Xi/U1 = (G(i,1)-G(i,2)/R) + (2*G(i,2)/R)*UC2/U1
% Xi/D  =  G(i,3)             + (2*G(i,2)/R)*UC2/D

Nenner_closed = simplify(1 - 2*G(4,2)/R);

disp(' ')
disp('====================================================')
disp('===      GESCHLOSSENE ÜBERTRAGUNGSFUNKTIONEN      ===')
disp('===         (Lastschluss: iR = u2/R)              ===')
disp('====================================================')

G_closed = sym(zeros(4,2));   % Spalten: U1, D
eingangs_namen_closed = {'U1', 'D'};

%% -- Schritt 1: UC2 schließen ----------------------------------------------

G_closed(4,1) = simplify( ...
    (G(4,1) - G(4,2)/R) / Nenner_closed );

G_closed(4,2) = simplify( ...
    G(4,3) / Nenner_closed );

%% -- Schritt 2: übrige Zustände über UC2 rückführen -----------------------

for i = 1:3

    % Xi / U1
    G_closed(i,1) = simplify( ...
        (G(i,1) - G(i,2)/R) ...
        + (2*G(i,2)/R)*G_closed(4,1) );

    % Xi / D
    G_closed(i,2) = simplify( ...
        G(i,3) ...
        + (2*G(i,2)/R)*G_closed(4,2) );

end

%% -- Ausgabe ---------------------------------------------------------------

for i = 1:4
    for j = 1:2

        fprintf('\n=== %s / %s (geschlossen) ===\n', ...
            zustands_namen{i}, ...
            eingangs_namen_closed{j})

        [num, den] = numden(G_closed(i,j));

        pretty(collect(expand(num), s))
        disp('---------------------------------------------')
        pretty(collect(expand(den), s))

    end
end

%% -- Ausgangs-TFs: u2 = 2*UC2 - u1 -----------------------------------------
disp(' ')
disp('====================================================')
disp('===          AUSGANGS-TFs  (u2 = 2*UC2 - u1)     ===')
disp('====================================================')

TF_u2_U1 = simplify( 2*G_closed(4,1) - 1 );
TF_u2_D  = simplify( 2*G_closed(4,2)     );

fprintf('\n=== u2 / U1 ===\n')
[num, den] = numden(TF_u2_U1);
pretty(collect(expand(num), s)), disp('-----')
pretty(collect(expand(den), s))

fprintf('\n=== u2 / D ===\n')
[num, den] = numden(TF_u2_D);
pretty(collect(expand(num), s)), disp('-----')
pretty(collect(expand(den), s))

%% -- Arbeitspunkt einsetzen --------------------------------------------------
UC10_val = 18;        % [V]
UC20_val = 18;        % [V]
U10_val  = 12;        % [V]
IL10_val = 1.77;      % [A]
IL20_val = 0.59;      % [A]
D0_val   = 0.25;      % [-]
L1_val   = 47e-6;     % [H]
L2_val   = 47e-6;     % [H]
C1_val   = 330e-6;    % [F]
C2_val   = 330e-6;    % [F]
R_val    = 20;        % [Ohm]

subs_vars = [UC10, UC20, U10, IL10, IL20, D0, L1, L2, C1, C2, R];
subs_vals = [UC10_val, UC20_val, U10_val, IL10_val, IL20_val, ...
             D0_val, L1_val, L2_val, C1_val, C2_val, R_val];

%% -- Hilfsfunktion: symbolische TF → MATLAB tf Objekt -----------------------
function sys = TF_to_sys(TF_sym, s, subs_vars, subs_vals)
    % 1. Numerische Werte einsetzen
    TF_num = subs(TF_sym, subs_vars, subs_vals);
    
    % 2. Zähler und Nenner trennen
    [num, den] = numden(TF_num);
    
    % 3. Symbolic polynomial → numerischer Koeffizientenvektor
    %    sym2poly gibt direkt absteigende Potenzen zurück (kompatibel mit tf)
    %    vpa erhöht Präzision bei großen Koeffizientenunterschieden (z.B. 1 bis 1e15)
    n = double(sym2poly(vpa(expand(num))));
    d = double(sym2poly(vpa(expand(den))));
    
    sys = tf(n, d);
end

%% -- Plots: alle geschlossenen TFs ------------------------------------------
figure('Name','Sprungantworten (geschlossen)');
titel_closed = {'IL1/U1','IL2/U1','UC1/U1','UC2/U1'; ...
                'IL1/D', 'IL2/D', 'UC1/D', 'UC2/D'};
for j = 1:2
    for i = 1:4
        subplot(2,4, (j-1)*4 + i)
        sys = TF_to_sys(G_closed(i,j), s, subs_vars, subs_vals);
        step(sys); grid on;
        title(titel_closed{j,i});
    end
end
sgtitle('Sprungantworten – geschlossenes System');

figure('Name','Bodediagramme (geschlossen)');
for j = 1:2
    for i = 1:4
        subplot(2,4, (j-1)*4 + i)
        sys = TF_to_sys(G_closed(i,j), s, subs_vars, subs_vals);
        bode(sys); grid on;
        title(titel_closed{j,i});
    end
end
sgtitle('Bodediagramme – geschlossenes System');

%% -- Plots: Ausgangs-TFs u2 --------------------------------------------------
sys_u2_U1 = TF_to_sys(TF_u2_U1, s, subs_vars, subs_vals);
sys_u2_D  = TF_to_sys(TF_u2_D,  s, subs_vars, subs_vals);

figure('Name','u2 Sprungantworten');
subplot(1,2,1); step(sys_u2_U1); grid on; title('Sprungantwort: u2/U1');
subplot(1,2,2); step(sys_u2_D);  grid on; title('Sprungantwort: u2/D');

figure('Name','u2 Bodediagramme');
subplot(1,2,1); bode(sys_u2_U1); grid on; title('Bode: u2/U1');
subplot(1,2,2); bode(sys_u2_D);  grid on; title('Bode: u2/D');


%% -- Pole / Nullstellen ----------------------------------------------------
disp(' ')
disp('=== Pole u2/U1 ===')
disp(pole(sys_u2_U1))

disp('=== Nullstellen u2/U1 ===')
disp(zero(sys_u2_U1))


%% -- Normierung & LaTeX-Ausgabe aller geschlossenen TFs ---------------------
function print_normalized_sym(TF_sym, s, label, norm_factor)

    % symbolische Zerlegung
    [num, den] = numden(simplify(TF_sym));

    % Normierung
    num_norm = simplifyFraction(num / norm_factor);
    den_norm = simplifyFraction(den / norm_factor);

    % Nach s sammeln und ausmultiplizieren
    num_c = collect(expand(num_norm), s);
    den_c = collect(expand(den_norm), s);

    fprintf('\n=== %s ===\n', label)
    disp('-- Zähler (pretty) --')
    pretty(num_c)
    fprintf('LaTeX: %s\n', latex(num_c))

    disp('-- Nenner (pretty) --')
    pretty(den_c)
    fprintf('LaTeX: %s\n', latex(den_c))
    disp('---------------------------------------------')
end

norm_factor = C1*C2*L1*L2*R;

disp(' ')
disp('====================================================')
disp('===   NORMIERTE ZÄHLER / NENNER  (÷ C1·C2·L1·L2·R) ===')
disp('====================================================')

labels_closed = {'IL1/U1','IL2/U1','UC1/U1','UC2/U1'; ...
                 'IL1/D', 'IL2/D', 'UC1/D', 'UC2/D'};

for j = 1:2
    for i = 1:4
        print_normalized_sym(G_closed(i,j), s, labels_closed{j,i}, norm_factor);
    end
end

disp(' ')
disp('====================================================')
disp('===        AUSGANGS-TFs  u2 (normiert)            ===')
disp('====================================================')

print_normalized_sym(TF_u2_U1, s, 'u2/U1', norm_factor);
print_normalized_sym(TF_u2_D,  s, 'u2/D',  norm_factor);

% [num, den] = numden(G_closed(4,1));
% 
% 
% den_norm = simplify(den / (C1*C2*L1*L2*R));
% num_norm = simplify(num / (C1*C2*L1*L2*R));
% 
% latex(collect(expand(num_norm), s))
% 
% pretty(collect(expand(num_norm), s))
% disp('---------------------------------------------')
% pretty(collect(expand(den_norm), s))
% 
% latex(collect(expand(den_norm), s))

%% -- Textreport öffnen ------------------------------------------------------
fid = fopen('Transferfunktionen_Report.txt','w');

fprintf(fid,'=====================================================\n');
fprintf(fid,' SYMBOLISCHE ÜBERTRAGUNGSFUNKTIONEN\n');
fprintf(fid,'=====================================================\n\n');

%% -- Hilfsfunktion ----------------------------------------------------------
function write_tf_report(fid, TF_sym, s, label, norm_factor)

    [num, den] = numden(simplify(TF_sym));

    num_norm = simplify(num / norm_factor);
    den_norm = simplify(den / norm_factor);

    num_c = collect(expand(num_norm), s);
    den_c = collect(expand(den_norm), s);

    fprintf(fid,'=====================================================\n');
    fprintf(fid,'%s\n', label);
    fprintf(fid,'=====================================================\n\n');

    fprintf(fid,'ZAEHLER\n');
    fprintf(fid,'-----------------------------------------------------\n');
    fprintf(fid,'%s\n\n', char(num_c));

    fprintf(fid,'LATEX ZAEHLER\n');
    fprintf(fid,'-----------------------------------------------------\n');
    fprintf(fid,'%s\n\n', latex(num_c));

    fprintf(fid,'NENNER\n');
    fprintf(fid,'-----------------------------------------------------\n');
    fprintf(fid,'%s\n\n', char(den_c));

    fprintf(fid,'LATEX NENNER\n');
    fprintf(fid,'-----------------------------------------------------\n');
    fprintf(fid,'%s\n\n', latex(den_c));

    fprintf(fid,'LATEX GESAMTE TF\n');
    fprintf(fid,'-----------------------------------------------------\n');
    fprintf(fid,'\\begin{equation}\n');
    fprintf(fid,'G(s)=\\frac{%s}{%s}\n', ...
        latex(num_c), ...
        latex(den_c));
    fprintf(fid,'\\end{equation}\n\n\n');

end

norm_factor = C1*C2*L1*L2*R;

labels_closed = {'IL1/U1','IL2/U1','UC1/U1','UC2/U1'; ...
                 'IL1/D', 'IL2/D', 'UC1/D', 'UC2/D'};

for j = 1:2
    for i = 1:4
        write_tf_report( ...
            fid, ...
            G_closed(i,j), ...
            s, ...
            labels_closed{j,i}, ...
            norm_factor);
    end
end

write_tf_report(fid,TF_u2_U1,s,'u2/U1',norm_factor);
write_tf_report(fid,TF_u2_D ,s,'u2/D' ,norm_factor);

fclose(fid);

