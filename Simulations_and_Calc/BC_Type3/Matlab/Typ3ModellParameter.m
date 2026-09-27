clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter des Typ3 Hochsetzstellers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dutyCycleStart  = 0;    % [1]
dutyCycleRun    = 0.25; % [1]
dutystepTime    = 15e-3;% [s]
u1stepTime      = 35e-3;% [s]
simTime         = 50e-3;% [s]

% Umgebungsparameter
R   = 20;           % [ohm]
u1ini = 12;
u1  = 15;           % [V]
fs  = 100e3;        % [Hz]

% Stufe 1
L11 = 47e-6;        % [H]
L21 = 47e-6;        % [H]
C11 = 330e-6;       % [F]
C21 = 330e-6;       % [F]
UD1 = 0.01;         % [V]

RC11 = 0.090;       % [ohm]
RC21 = 0.090;       % [ohm]
RL11 = 0.00295;     % [ohm]
RL21 = 0.00295;     % [ohm]
RD1  = 0.0048;      % [ohm]
RS1  = 0.0048;      % [ohm]
RDS1 = RS1;         % [ohm]
UDS1 = 1.2          % [V]

% Stufe 2
L12 = 47e-6;        % [H]
L22 = 47e-6;        % [H]
C12 = 330e-6;       % [F]
C22 = 330e-6;       % [F]
UD2 = 0.01;         % [V]

RC12 = 0.090;       % [ohm]
RC22 = 0.090;       % [ohm]
RL12 = 0.00295;     % [ohm]
RL22 = 0.00295;     % [ohm]
RD2  = 0.0048;      % [ohm]
RS2  = 0.0048;      % [ohm]
RDS2 = RS2;         % [ohm]
UDS2 = 1.2          % [V]




