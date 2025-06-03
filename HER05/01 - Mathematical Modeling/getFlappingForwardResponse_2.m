function [Beta, Beta_dot] = getFlappingForwardResponse(Vx,e,CT,theta0_deg,theta_tw_deg,theta_1c_deg,theta_1s_deg,rotor)

    % All output beta angles are radians !!
    % All input theta angles are degrees !!
   

    
    % ===== Physical parameters =====

Cl_alpha = 6.2305;         % Lift curve slope (1/rad)

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
gamma = rho_cruise * c * Cl_alpha * R^4/Ib;

% === Rotating flap frequency ===
nu_flap = sqrt(1 + (3*e)/(2*(1-e)) + k_beta/(Ib*Omega^2));

% === Advance Ratio ===
V_forward = Vx;   % forward speed in m/s
mu = V_forward / V_tip;      % advance ratio
beta_p_deg = 2.5;  % Pre-cone angle in deg
beta_p_rad = deg2rad(beta_p_deg);

% 
% % Get outputs
% beta0_rad = (gamma/(nu_flap^2)) * ((theta0_rad/8)*(1+mu^2) + (theta_tw_rad/10)*(1+(5/6)*mu^2) + (mu/6)*theta_1s_rad - (CT/(12*mu)));
% 
% % x = beta_1c, y = beta_1s
% syms x y
% 
% eq1 = y == (gamma/(nu_flap^2 - 1)) * ((1/8)*(theta_1s_rad + x)*(1 - 0.5*mu^2) + mu*theta0_rad/3 - CT/8 + mu^2*theta_1s_rad/4 + mu*theta_tw_rad/4);
% 
% eq2 = x == (gamma/(nu_flap^2 - 1)) * ((1/8)*(theta_1c_rad - y)*(1 + 0.5*mu^2) - mu*beta0_rad/6);
% 
% 
% sol = solve([eq1, eq2], [x, y]);
% 
% beta_1c_rad = double(sol.x);
% 
% beta_1s_rad = double(sol.y);
% 
% 
% Beta = @(psi) beta0_rad + beta_1c_rad*cos(psi) + beta_1s_rad*sin(psi);
% 
% Beta_dot = @(psi) (-beta_1c_rad*sin(psi) + beta_1s_rad*cos(psi))*Omega;

if rotor == "upper"

    theta0_u = 11.61; % from BEMT
    theta0_u_rad = deg2rad(theta0_u);
    theta_tw_u = -11;
    theta_tw_u_rad = deg2rad(theta_tw_u);

    beta0_rad = (gamma/(nu_flap^2)) * ((theta0_u_rad/8)*(1+mu^2) + (theta_tw_u_rad/10)*(1+(5/6)*mu^2) + (mu/6)*theta_1s_rad - (CT/(12*mu))) + k_beta*beta_p_rad/(Ib*Omega^2*nu_flap^2);
    
    
    % x = beta_1c, y = beta_1s
    syms x1 y1
    
    eq1 = y1 == (gamma/(nu_flap^2 - 1)) * ((1/8)*(theta_1s_rad + x1)*(1 - 0.5*mu^2) + mu*theta0_u_rad/3 - CT/8 + mu^2*theta_1s_rad/4 + mu*theta_tw_u_rad/4);
    
    eq2 = x1 == (gamma/(nu_flap^2 - 1)) * ((1/8)*(theta_1c_rad - y1)*(1 + 0.5*mu^2) - mu*beta0_rad/6);
    
    
    sol = solve([eq1, eq2], [x1, y1]);
    
    beta_1c_rad = double(sol.x1);
   
    beta_1s_rad = double(sol.y1);



elseif rotor == "lower"
    
    theta0_l = 7.82; 
    theta0_l_rad = deg2rad(theta0_l);
    theta_tw_l = -6;
    theta_tw_l_rad = deg2rad(theta_tw_l);
    
    
    
    
    beta0_rad = (gamma/(nu_flap^2)) * ((theta0_l_rad/8)*(1+mu^2) + (theta_tw_l_rad/10)*(1+(5/6)*mu^2) + (mu/6)*theta_1s_rad - (CT/(12*mu))) + k_beta*beta_p_rad/(Ib*Omega^2*nu_flap^2);
    
    
    % x = beta_1c, y = beta_1s
    syms x2 y2
    
    eq1 = y2 == (gamma/(nu_flap^2 - 1)) * ((1/8)*(theta_1s_rad + x2)*(1 - 0.5*mu^2) + mu*theta0_l_rad/3 - CT/8 + mu^2*theta_1s_rad/4 + mu*theta_tw_l_rad/4);
    
    eq2 = x2 == (gamma/(nu_flap^2 - 1)) * ((1/8)*(theta_1c_rad - y2)*(1 + 0.5*mu^2) - mu*beta0_rad/6);
    
    
    sol = solve([eq1, eq2], [x2, y2]);
    
    beta_1c_rad = double(sol.x2);
    
    beta_1s_rad = double(sol.y2);
   
end

Beta = @(psi) beta0_rad + beta_1c_rad*cos(psi) + beta_1s_rad*sin(psi);

Beta_dot = @(psi) (-beta_1c_rad*sin(psi) + beta_1s_rad*cos(psi))*Omega;