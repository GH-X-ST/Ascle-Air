function [beta_1c_rad, beta_1s_rad] = getFlappingForwardResponse(e,CT,theta0_deg,theta_tw_deg,theta_1c_deg,theta_1s_deg)

    % All output beta angles are radians !!
    % All input theta angles are degrees !!
   

    
    % ===== Physical parameters =====

Cl_alpha = 5;         % Lift curve slope (1/rad)

% Get all constants
constants = getConstants();  % or use setupConstants() directly if you prefer

% Access individual constants using dot notation
GTOW    = constants.GTOW;      % kg
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


% theta0 & theta_twist values

theta0_rad = theta0_deg * (pi/180);  % in radian
theta_tw_rad = theta_tw_deg * (pi/180); % in radian

% theta_1s & theta_1c values 
theta_1s_rad = theta_1s_deg * (pi/180); % in radian

theta_1c_rad = theta_1c_deg * (pi/180); % in radian

% Flexbeam parameters
E = 140e9;            % Material flexbeam
I = (0.2*0.025^3)/12 ;  % Dimension of flexbeam ---- need to update!!
k_beta = (E*I)/(0.14*R);  % Equivalent spring stiffness(Nm/rad)  ---- need to update!!

m = 8.88;            % Uniform mass distribution (kg/m) --- based on blade material

% === Mass moment of inertia===
Ib = (1/3)*m*R^3*(1-e)^3;

% === Lock Number ===
gamma = rho_cruise * c * Cl_alpha * R^4/Ib;

% === Rotating flap frequency ===
nu_flap = sqrt(1 + (3*e)/(2*(1-e)) + k_beta/(Ib*Omega^2));

% === Advance Ratio ===
V_forward = 140 / 1.94384;   % forward speed in m/s
mu = V_forward / V_tip;      % advance ratio


% Get outputs
beta0_rad = (gamma/(nu_flap^2)) * ((theta0_rad/8)*(1+mu^2) + (theta_tw_rad/10)*(1+(5/6)*mu^2) + (mu/6)*theta_1s_rad - (CT/(12*mu)));

% x = beta_1c, y = beta_1s
syms x y

eq1 = y == (gamma/(nu_flap^2 - 1)) * ((1/8)*(theta_1s_rad + x)*(1 - 0.5*mu^2) + mu*theta0_rad/3 - CT/8 + mu^2*theta_1s_rad/4 + mu*theta_tw_rad/4);

eq2 = x == (gamma/(nu_flap^2 - 1)) * ((1/8)*(theta_1c_rad - y)*(1 + 0.5*mu^2) - mu*beta0_rad/6);


sol = solve([eq1, eq2], [x, y]);

beta_1c_rad = double(sol.x);

beta_1s_rad = double(sol.y);


Beta = @(psi) beta0_rad + beta_1c_rad*cos(psi) + beta_1s_rad*sin(psi);

Beta_dot = @(psi) (-beta_1c_rad*sin(psi) + beta_1s_rad*cos(psi))*Omega;