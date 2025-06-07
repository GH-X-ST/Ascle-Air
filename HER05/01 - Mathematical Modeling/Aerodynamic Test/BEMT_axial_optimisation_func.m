function [FM4, eta, Ct_u_final, Ct_l_final] = BEMT_axial_optimisation_func(theta_tw_u, theta_tw_l, TR_u, TR_l, Vc, polar)

% Get all constants
constants = getConstants();  % or use setupConstants() directly if you prefer

theta_tw_u = rad2deg(theta_tw_u);
theta_tw_l = rad2deg(theta_tw_l);

% Access individual constants using dot notation
GTOW    = constants.GTOW;      % kg
Nb      = constants.Nb;
AR      = constants.AR;
Vtip    = constants.Vtip;      % m/s
R       = constants.R;         % m
r_0     = constants.r0;        % ND root cutout
A       = constants.A;         % m^2
Ae      = constants.Ae;        % m^2 (effective area)
cbar    = constants.cbar;      % mean chord
Cd0     = constants.Cd0;       % profile drag coefficient
sigma_e = constants.sigma_e;   % thrust-weighted equivalent solidity
rho     = constants.rho;       % kg/m^3
rho_cr  = constants.rho_cr;    % cruise density
Cw      = constants.Cw;        % weight coefficient

% start by thrust trim - equal thrust sharing
Ct_u_req = Cw/2;
Ct_l_req = Cw/2;

d = 0;

count  = 1;
error = 1;
tol = 1e-6;
while abs(error) > tol
% for ii = 1:30
    % thrust trim
    [r, Ct_u, Cp_u, ~, dCp_u, lambda_u_induced, lambda_u, ~, ~, ~,Cp_induced_u, Cp_profile_u] = BEMT_axial_UpperRotor(theta_tw_u, TR_u, Vc, Ct_u_req,polar);
    [~, Ct_l, Cp_l, ~, dCp_l, lambda_l_induced, lambda_l, ~, ~, ~,Cp_induced_l, Cp_profile_l] = BEMT_axial_LowerRotor(theta_tw_l, TR_l, lambda_u, d, Vc, Ct_l_req,polar);

    error = Cp_u-Cp_l;
    eps(count) = error;

    if Cp_u < Cp_l
        Ct_u_req = Ct_u_req - error;
        Ct_l_req = Ct_l_req + error;

    elseif Cp_u > Cp_l
        Ct_u_req = Ct_u_req + error;
        Ct_l_req = Ct_l_req - error;
    end

    ct_u(count) = Ct_u_req;
    ct_l(count) = Ct_l_req;

    cp_u(count) = Cp_u;
    cp_l(count) = Cp_l;

    dCp_u_dr = dCp_u/(r(2)-r(1));
    dCp_l_dr = dCp_l/(r(2)-r(1));

    count = count +1;
    
end

Ct_u_final = ct_u(end-1);
Ct_l_final = ct_l(end-1);

Cp_u_final = Cp_u;
Cp_l_final = Cp_l;

% Figure of Merit calculation
k_int = 1.219; % torque trimmed coaxial interference factor
k = 1.1; % induced loss factor

FM1 = Ct_l_final^(3/2)/sqrt(2) * ((Ct_u_final/Ct_l_final)^(3/2)+1) / ...
    (k_int*k*(Ct_l_final^(3/2)/sqrt(2) * ((Ct_u_final/Ct_l_final)^(3/2)+1)) + 2*sigma_e*Cd0/9);

FM2 = Cw^(3/2)/sqrt(2)/(k_int*k*Cw^(3/2)/sqrt(2) + 2*sigma_e*Cd0/9);

% FM from BEMT
k_int = 1.219; % torque trimmed coaxial interference factor
k = 1.1; % induced loss factor
Cp_ideal = sum(Cp_induced_u + Cp_induced_l);
Cp_profile = sum(Cp_profile_u + Cp_profile_l);
FM3 = Cp_ideal/(k_int*k*Cp_ideal + Cp_profile);
FM4 = Cw^(3/2)/sqrt(2) / (k_int*k*Cp_ideal + Cp_profile);

% propulsive efficiency evaluated at 6 m/s
lambda_c = 6/Vtip;
[~, Ct_u_climb, ~, ~, ~, lambda_u_induced_climb, lambda_u_climb, ~, ~, ~, ~, ~] = BEMT_axial_UpperRotor(theta_tw_u, TR_u, 6, Ct_u_final, polar);
[~, Ct_l_climb, ~, ~, ~, lambda_l_induced_climb, ~, ~, ~, ~, ~, ~] = BEMT_axial_LowerRotor(theta_tw_l, TR_l, lambda_u_climb, d, 6, Ct_l_final, polar);
eta_u = Ct_u_climb*lambda_c/(Ct_u_climb*lambda_c + k_int*k*Ct_u_final*lambda_u_induced_climb + sigma_e*Cd0/9);
eta_l = Ct_l_climb*lambda_c/(Ct_l_climb*lambda_c + k_int*k*Ct_l_final*lambda_l_induced_climb + sigma_e*Cd0/9);
eta = 0.5*(eta_u+eta_l);

end