%% ======================================================================%%
%  ------------------------ MICROWAVE AMPLIFIER ------------------------- %
% Author: Nedal M. Benelmekki                                             %
% Date (DD/MM/YYYY): 18/05/2025                                           %

%% ======================================================================%%
%  ----------------------- TRANSISTOR PARAMETERS ------------------------ %
% Transistor MPN: AT-41411                                                %
% Operating Frequency: 3 GHz                                              %
S11_mag = 0.529;                                                          %
S11_phase = 139.153; % Degrees                                            %
S21_mag = 2.896;                                                          %
S21_phase = 55.280; % Degrees                                             %
S12_mag = 0.098;                                                          %
S12_phase = 60.012; % Degrees                                             %
S22_mag = 0.335;                                                          %
S22_phase = -27.192; % Degrees                                            %

%% ======================================================================%%
%  ------------------------- HELPER FUNCTIONS --------------------------- %

function [z] = polarToCartesian(magnitude, phaseDeg)
    phase = deg2rad(phaseDeg);                 % Convert to radians
    z = magnitude * exp(1i * phase);           % Polar to complex
end

function [x] = linearTodB(linear)
    x = 10 * log10(linear);
end

function [delta] = calculateDelta(S11, S12, S21, S22)
    delta = S11 * S22 - S12 * S21;
end

function [k] = calculateK(S11, S22, delta, S12, S21)
    k = (1 - abs(S11)^2 - abs(S22)^2 + abs(delta)^2) / (2 * abs(S12 * S21));
end

function [B1] = calculateB1(S11, S22, delta)
    B1 = 1 + abs(S11)^2 - abs(S22)^2 - abs(delta)^2;
end

function [B2] = calculateB2(S11, S22, delta)
    B2 = 1 + abs(S22)^2 - abs(S11)^2 - abs(delta)^2;
end

function [C1] = calculateC1(S11, S22, delta)
    C1 = S11 - delta * conj(S22);
end

function [C2] = calculateC2(S11, S22, delta)
    C2 = S22 - delta * conj(S11);
end

function [rho_s] = calculateRhoS(B1, C1)
    discriminant = sqrt(B1^2 - 4 * abs(C1)^2);
    rho1 = (B1 + discriminant) / (2 * C1);
    rho2 = (B1 - discriminant) / (2 * C1);
    rho_s = rho1; if abs(rho1) > 1, rho_s = rho2; end
end

function [rho_l] = calculateRhoL(B2, C2)
    discriminant = sqrt(B2^2 - 4 * abs(C2)^2);
    rho1 = (B2 + discriminant) / (2 * C2);
    rho2 = (B2 - discriminant) / (2 * C2);
    rho_l = rho1; if abs(rho1) > 1, rho_l = rho2; end
end

function [Gs_max] = calculateGSmax(rho_s)
    Gs_max = 1 / (1 - abs(rho_s)^2);
end

function [Go] = calculateGo(S21)
    Go = abs(S21)^2;
end

function [Gl_max] = calculateGLmax(rho_l, S22)
    Gl_max = (1 - abs(rho_l)^2)/(abs(1-S22*rho_l)^2);
end

function [Gt_max] = calculateGTmax(Gs_max, Go, Gl_max)
    Gt_max = Gs_max * Go * Gl_max;
end

%% ======================================================================%%
%  --------------------------- CALCULATIONS ----------------------------- %

% Convert S-parameters to complex
S11 = polarToCartesian(S11_mag, S11_phase);
S12 = polarToCartesian(S12_mag, S12_phase);
S21 = polarToCartesian(S21_mag, S21_phase);
S22 = polarToCartesian(S22_mag, S22_phase);

% Delta
delta = calculateDelta(S11, S12, S21, S22);

% B and C terms
B1 = calculateB1(S11, S22, delta);
B2 = calculateB2(S11, S22, delta);
C1 = calculateC1(S11, S22, delta);
C2 = calculateC2(S11, S22, delta);

% Reflection coefficients
Rho_s = calculateRhoS(B1, C1);
Rho_l = calculateRhoL(B2, C2);

% Gain terms
Gs_max = calculateGSmax(Rho_s);
Go = calculateGo(S21);
Gl_max = calculateGLmax(Rho_l, S22);
Gt = linearTodB(calculateGTmax(Gs_max, Go, Gl_max));

% Stability
K = calculateK(S11, S22, delta, S12, S21);

%% ======================================================================%%
%  ---------------------------- RESULTS --------------------------------- %

fprintf('\n================= RESULTS =================\n');
fprintf('Δ (Delta)              : %.4f ∠ %.2f°\n', abs(delta), rad2deg(angle(delta)));
fprintf('Gt (Max Gain)          : %.2f dB\n', Gt);
fprintf('K (Stability Factor)   : %.4f\n', K);
fprintf('Input Match |Γs|       : %.4f\n', abs(Rho_s));
fprintf('Output Match |Γl|      : %.4f\n', abs(Rho_l));
fprintf('============================================\n');


