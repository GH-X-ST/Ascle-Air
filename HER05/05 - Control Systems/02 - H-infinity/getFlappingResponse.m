function [Beta, Beta_dot, beta_1c_rad, beta_1s_rad] = getFlappingResponse(Vx, rho_input, m_input, CT, theta0_u_rad, theta0_l_rad, theta_1c_rad, theta_1s_rad, rotor)

% All input and output angles are radians
    
% ===== Physical parameters =====

Cl_alpha = 6.2305;         % Lift curve slope (1/rad)

% Get all constants
constants = getConstants();  % or use setupConstants() directly if you prefer

% Access individual constants using dot notation
% GTOW    = constants.GTOW;      % kg
GTOW    = m_input;
e       = constants.e;
Nb      = constants.Nb;
AR      = constants.AR;
V_tip    = constants.Vtip;      % m/s
R       = constants.R;         % m
r_0     = constants.r0;        % ND root cutout
A       = constants.A;         % m^2
Ae      = constants.Ae;        % m^2 (effective area)
c    = constants.cbar;      % mean chord
Cd0     = constants.Cd0;       % profile drag coefficient
sigma_e = constants.sigma_e;   % thrust-weighted equivalent solidity
rho     = constants.rho;       % kg/m^3
rho_cruise  = constants.rho_cr;    % cruise density
Cd_f    = constants.Cd_f;      % fuselage drag coeff
S_f     = constants.S_f;       % fuselage frontal area
Omega   = constants.Omega;     % angular speed

theta_tw_u_rad = constants.theta_tw_u;
theta_tw_l_rad = constants.theta_tw_l;

% Flexbeam parameters (from structure team)

b = 0.2;
L = 0.14* 5.3;
H = 0.06* 5.3;
h = 0.025;

yield_strength = 680e6;

E1 = 139.2e9;

I1 = b.*h.^3./12;

b1 = 0.015;
b2 = 0.13;
h1 = 0.03;
h2 = 0.016;

Ix1 = (1/12) * b1 * h1^3;
Iy1 = (1/12) * h1 * b1^3;

Ix2_local = (1/12) * b2 * h2^3;
Iy2 = (1/12) * h2 * b2^3;

d = (h1/2) - (h2/2);
A2 = b2 * h2;
Ix2 = Ix2_local + A2 * d^2;

Ix = Ix1 + Ix2;
Iy = Iy1 + Iy2;

k_beta =(3*E1*Ix)/(0.14*R);

m = 8.88;            % Uniform mass distribution (kg/m) --- based on blade material

% === Mass moment of inertia===
Ib = (1/3)*m*R^3*(1-e)^3;

% === Lock Number ===
gamma = rho_input * c * Cl_alpha * R^4/Ib;

% === Rotating flap frequency ===
nu_flap = sqrt(1 + (3*e)/(2*(1-e)) + k_beta/(Ib*Omega^2));

% === Advance Ratio ===
V_forward = Vx;   % forward speed in m/s
mu = V_forward / V_tip;      % advance ratio
beta_p_deg = 2.5;  % Pre-cone angle in deg
beta_p_rad = deg2rad(beta_p_deg);


% === Lambda ===
lambda = sqrt(0.5*CT + 0.25*mu.^2 ) - 0.5*mu;

% 

if rotor == "upper"

    beta0_rad = (gamma/nu_flap^2) * ...
            ( theta0_u_rad/8 * (1+mu^2) ...
            + theta_tw_u_rad/10 * (1+5*mu^2/6) ...
            + mu/6 * theta_1s_rad ...
            - lambda/6 ) ...
            + k_beta * beta_p_rad /(Ib*Omega^2*nu_flap^2);

    eqs = @(X) [
    X(2) - (gamma/(nu_flap^2-1))*...
           ( 1/8*(theta_1s_rad + X(1))*(1-0.5*mu^2) ...
           + mu*theta0_u_rad/3 ...
           - mu*lambda/4 ...
           + mu^2*theta_1s_rad/4 ...
           + mu*theta_tw_u_rad/4 );
    X(1) - (gamma/(nu_flap^2-1))*...
           ( 1/8*(theta_1c_rad - X(2))*(1+0.5*mu^2) ...
           - mu*beta0_rad/6 );
           ];

    opts = optimoptions('fsolve', 'Display', 'off'); 

    sol = fsolve(eqs, [0 0], opts);  % Initial guess
    beta_1c_rad = sol(1);
    beta_1s_rad = sol(2);

elseif rotor == "lower"
    
    beta0_rad = (gamma/nu_flap^2) * ...
            ( theta0_l_rad/8 * (1+mu^2) ...
            + theta_tw_l_rad/10 * (1+5*mu^2/6) ...
            + mu/6 * theta_1s_rad ...
            - lambda/6 ) ...
            + k_beta * beta_p_rad /(Ib*Omega^2*nu_flap^2);

    eqs = @(X) [
    X(2) - (gamma/(nu_flap^2-1))*...
           ( 1/8*(theta_1s_rad + X(1))*(1-0.5*mu^2) ...
           + mu*theta0_l_rad/3 ...
           - mu*lambda/4 ...
           + mu^2*theta_1s_rad/4 ...
           + mu*theta_tw_l_rad/4 );
    X(1) - (gamma/(nu_flap^2-1))*...
           ( 1/8*(theta_1c_rad - X(2))*(1+0.5*mu^2) ...
           - mu*beta0_rad/6 );
           ];

    
    opts = optimoptions('fsolve', 'Display', 'off'); 

    sol = fsolve(eqs, [0 0], opts);  % Initial guess
    beta_1c_rad = sol(1);
    beta_1s_rad = sol(2);
   
end

Beta = @(psi) beta0_rad + beta_1c_rad*cos(psi) + beta_1s_rad*sin(psi);

Beta_dot = @(psi) (-beta_1c_rad*sin(psi) + beta_1s_rad*cos(psi))*Omega;

end