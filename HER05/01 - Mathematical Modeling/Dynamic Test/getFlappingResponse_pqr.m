function [Beta, Beta_dot, beta_1c_rad, beta_1s_rad] = getFlappingResponse_pqr(Vx, rho_input, m_input, CT, theta0_u_rad, theta0_l_rad, theta_1c_rad, theta_1s_rad, rotor, Vy, p, q)

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
beta_p_deg = 2.5;  % Pre-cone angle in deg
beta_p_rad = deg2rad(beta_p_deg);

K_beta  = k_beta;     %  rename to match symbolic expression
beta_P  = beta_p_rad; %  rename to match symbolic expression
v_beta  = nu_flap;    %  rename to match symbolic expression
I_b     = Ib;         %  rename to match symbolic expression

%% Kinematics
% Normalised velocities (body axis)
u_bar = Vx / V_tip;               % forward speed ratio
v_bar = Vy / V_tip;               % lateral speed ratio

% Normalised body rates
p_hat_body =  p / Omega;
q_hat_body =  q / Omega;

p_hat = -p_hat_body;
q_hat = -q_hat_body;

% Advance ratio for inflow model (|μ|) and inflow λ (same as before)
mu     = sqrt(u_bar.^2 + v_bar.^2);
lambda = sqrt(0.5*CT + 0.25*mu.^2) - 0.5*mu;


%% Rotor‑specific collective & twist
if strcmpi(rotor,"upper")

    theta_0  = theta0_u_rad;
    theta_tw = theta_tw_u_rad;

else

    theta_0  = theta0_l_rad;
    theta_tw = theta_tw_l_rad;

end

%% Analytic flap solution

theta_1c = theta_1c_rad;
theta_1s = theta_1s_rad;

beta_0 = (120*K_beta*beta_P - 20*I_b*Omega^2*gamma*lambda + 15*I_b*Omega^2*gamma*theta_0 + 12*I_b*Omega^2*gamma*theta_tw + 15*I_b*Omega^2*gamma*theta_0*u_bar^2 + 10*I_b*Omega^2*gamma*theta_tw*u_bar^2 + 15*I_b*Omega^2*gamma*theta_0*v_bar^2 + 10*I_b*Omega^2*gamma*theta_tw*v_bar^2 + 10*I_b*Omega^2*gamma*p_hat*u_bar + 10*I_b*Omega^2*gamma*q_hat*v_bar + 20*I_b*Omega^2*gamma*theta_1s*u_bar + 20*I_b*Omega^2*gamma*theta_1c*v_bar)/(120*(I_b*Omega^2*v_beta^2 + K_beta));
beta_1c = -(I_b*Omega^2*(1920*K_beta^2*beta_P*gamma*u_bar - 1440*K_beta^2*gamma*q_hat - 1440*K_beta^2*gamma*theta_1c - 23040*K_beta^2*p_hat*v_beta^2 + 2880*K_beta^2*gamma*lambda*v_bar - 3840*K_beta^2*gamma*theta_0*v_bar - 2880*K_beta^2*gamma*theta_tw*v_bar - 720*K_beta^2*gamma*theta_1c*u_bar^2 - 2160*K_beta^2*gamma*theta_1c*v_bar^2 + 23040*I_b^2*Omega^4*p_hat*v_beta^4 - 23040*I_b^2*Omega^4*p_hat*v_beta^6 + 320*I_b^2*Omega^4*gamma^2*lambda*u_bar - 40*I_b^2*Omega^4*gamma^3*lambda*v_bar + 1440*I_b^2*Omega^4*gamma*q_hat*v_beta^2 - 4320*I_b^2*Omega^4*gamma*q_hat*v_beta^4 - 240*I_b^2*Omega^4*gamma^2*theta_0*u_bar - 192*I_b^2*Omega^4*gamma^2*theta_tw*u_bar + 30*I_b^2*Omega^4*gamma^3*theta_0*v_bar + 1440*I_b^2*Omega^4*gamma*theta_1c*v_beta^2 - 1440*I_b^2*Omega^4*gamma*theta_1c*v_beta^4 + 24*I_b^2*Omega^4*gamma^3*theta_tw*v_bar + 1440*I_b*K_beta*Omega^2*gamma*q_hat + 1440*I_b*K_beta*Omega^2*gamma*theta_1c + 20*I_b^2*Omega^4*gamma^3*lambda*v_bar^3 - 160*I_b^2*Omega^4*gamma^2*p_hat*u_bar^2 + 180*I_b^2*Omega^4*gamma^2*p_hat*v_beta^2 + 20*I_b^2*Omega^4*gamma^3*q_hat*v_bar^2 - 10*I_b^2*Omega^4*gamma^3*q_hat*v_bar^4 - 240*I_b^2*Omega^4*gamma^2*theta_0*u_bar^3 - 320*I_b^2*Omega^4*gamma^2*theta_1s*u_bar^2 - 160*I_b^2*Omega^4*gamma^2*theta_tw*u_bar^3 + 15*I_b^2*Omega^4*gamma^3*theta_0*v_bar^3 - 15*I_b^2*Omega^4*gamma^3*theta_0*v_bar^5 + 40*I_b^2*Omega^4*gamma^3*theta_1c*v_bar^2 - 20*I_b^2*Omega^4*gamma^3*theta_1c*v_bar^4 + 180*I_b^2*Omega^4*gamma^2*theta_1s*v_beta^2 + 8*I_b^2*Omega^4*gamma^3*theta_tw*v_bar^3 - 10*I_b^2*Omega^4*gamma^3*theta_tw*v_bar^5 - 1440*K_beta^2*gamma*theta_1s*u_bar*v_bar + 180*I_b*K_beta*Omega^2*gamma^2*p_hat + 180*I_b*K_beta*Omega^2*gamma^2*theta_1s + 23040*I_b*K_beta*Omega^2*p_hat*v_beta^2 - 46080*I_b*K_beta*Omega^2*p_hat*v_beta^4 + 20*I_b^2*Omega^4*gamma^3*lambda*u_bar^2*v_bar - 680*I_b^2*Omega^4*gamma^2*lambda*u_bar*v_beta^2 - 10*I_b^2*Omega^4*gamma^3*p_hat*u_bar*v_bar^3 - 10*I_b^2*Omega^4*gamma^3*p_hat*u_bar^3*v_bar - 1440*I_b^2*Omega^4*gamma*q_hat*u_bar^2*v_beta^4 + 1440*I_b^2*Omega^4*gamma*q_hat*v_bar^2*v_beta^4 - 240*I_b^2*Omega^4*gamma^2*theta_0*u_bar*v_bar^2 + 15*I_b^2*Omega^4*gamma^3*theta_0*u_bar^2*v_bar - 15*I_b^2*Omega^4*gamma^3*theta_0*u_bar^4*v_bar + 720*I_b^2*Omega^4*gamma^2*theta_0*u_bar*v_beta^2 - 20*I_b^2*Omega^4*gamma^3*theta_1s*u_bar*v_bar^3 - 20*I_b^2*Omega^4*gamma^3*theta_1s*u_bar^3*v_bar + 720*I_b^2*Omega^4*gamma*theta_1c*u_bar^2*v_beta^2 - 720*I_b^2*Omega^4*gamma*theta_1c*u_bar^2*v_beta^4 - 160*I_b^2*Omega^4*gamma^2*theta_tw*u_bar*v_bar^2 + 8*I_b^2*Omega^4*gamma^3*theta_tw*u_bar^2*v_bar - 10*I_b^2*Omega^4*gamma^3*theta_tw*u_bar^4*v_bar + 552*I_b^2*Omega^4*gamma^2*theta_tw*u_bar*v_beta^2 + 2160*I_b^2*Omega^4*gamma*theta_1c*v_bar^2*v_beta^2 - 2160*I_b^2*Omega^4*gamma*theta_1c*v_bar^2*v_beta^4 + 240*I_b*K_beta*Omega^2*beta_P*gamma^2*v_bar - 680*I_b*K_beta*Omega^2*gamma^2*lambda*u_bar - 5760*I_b*K_beta*Omega^2*gamma*q_hat*v_beta^2 + 720*I_b*K_beta*Omega^2*gamma^2*theta_0*u_bar + 720*I_b*K_beta*Omega^2*gamma*theta_1c*u_bar^2 + 552*I_b*K_beta*Omega^2*gamma^2*theta_tw*u_bar + 2160*I_b*K_beta*Omega^2*gamma*theta_1c*v_bar^2 - 2880*I_b*K_beta*Omega^2*gamma*theta_1c*v_beta^2 - 180*I_b^2*Omega^4*gamma^2*lambda*u_bar^3*v_beta^2 + 250*I_b^2*Omega^4*gamma^2*p_hat*u_bar^2*v_beta^2 - 10*I_b^2*Omega^4*gamma^3*q_hat*u_bar^2*v_bar^2 - 90*I_b^2*Omega^4*gamma^2*p_hat*v_bar^2*v_beta^2 - 30*I_b^2*Omega^4*gamma^3*theta_0*u_bar^2*v_bar^3 - 20*I_b^2*Omega^4*gamma^3*theta_1c*u_bar^2*v_bar^2 + 480*I_b^2*Omega^4*gamma^2*theta_0*u_bar^3*v_beta^2 + 680*I_b^2*Omega^4*gamma^2*theta_1s*u_bar^2*v_beta^2 + 135*I_b^2*Omega^4*gamma^2*theta_1s*u_bar^4*v_beta^2 - 20*I_b^2*Omega^4*gamma^3*theta_tw*u_bar^2*v_bar^3 + 340*I_b^2*Omega^4*gamma^2*theta_tw*u_bar^3*v_beta^2 - 45*I_b^2*Omega^4*gamma^2*theta_1s*v_bar^4*v_beta^2 - 120*I_b*K_beta*Omega^2*beta_P*gamma^2*v_bar^3 - 180*I_b*K_beta*Omega^2*gamma^2*lambda*u_bar^3 + 250*I_b*K_beta*Omega^2*gamma^2*p_hat*u_bar^2 - 90*I_b*K_beta*Omega^2*gamma^2*p_hat*v_bar^2 + 480*I_b*K_beta*Omega^2*gamma^2*theta_0*u_bar^3 + 680*I_b*K_beta*Omega^2*gamma^2*theta_1s*u_bar^2 + 135*I_b*K_beta*Omega^2*gamma^2*theta_1s*u_bar^4 + 340*I_b*K_beta*Omega^2*gamma^2*theta_tw*u_bar^3 - 45*I_b*K_beta*Omega^2*gamma^2*theta_1s*v_bar^4 - 2880*I_b^2*Omega^4*gamma*lambda*v_bar*v_beta^2 + 2880*I_b^2*Omega^4*gamma*lambda*v_bar*v_beta^4 + 20*I_b^2*Omega^4*gamma^3*p_hat*u_bar*v_bar - 160*I_b^2*Omega^4*gamma^2*q_hat*u_bar*v_bar - 320*I_b^2*Omega^4*gamma^2*theta_1c*u_bar*v_bar + 40*I_b^2*Omega^4*gamma^3*theta_1s*u_bar*v_bar + 3840*I_b^2*Omega^4*gamma*theta_0*v_bar*v_beta^2 - 3840*I_b^2*Omega^4*gamma*theta_0*v_bar*v_beta^4 + 2880*I_b^2*Omega^4*gamma*theta_tw*v_bar*v_beta^2 - 2880*I_b^2*Omega^4*gamma*theta_tw*v_bar*v_beta^4 - 1920*I_b*K_beta*Omega^2*beta_P*gamma*u_bar - 2880*I_b*K_beta*Omega^2*gamma*lambda*v_bar + 3840*I_b*K_beta*Omega^2*gamma*theta_0*v_bar + 2880*I_b*K_beta*Omega^2*gamma*theta_tw*v_bar - 180*I_b^2*Omega^4*gamma^2*lambda*u_bar*v_bar^2*v_beta^2 + 480*I_b^2*Omega^4*gamma^2*theta_0*u_bar*v_bar^2*v_beta^2 + 180*I_b^2*Omega^4*gamma^2*theta_1c*u_bar*v_bar^3*v_beta^2 + 180*I_b^2*Omega^4*gamma^2*theta_1c*u_bar^3*v_bar*v_beta^2 + 340*I_b^2*Omega^4*gamma^2*theta_tw*u_bar*v_bar^2*v_beta^2 - 120*I_b*K_beta*Omega^2*beta_P*gamma^2*u_bar^2*v_bar - 180*I_b*K_beta*Omega^2*gamma^2*lambda*u_bar*v_bar^2 - 1440*I_b*K_beta*Omega^2*gamma*q_hat*u_bar^2*v_beta^2 + 1440*I_b*K_beta*Omega^2*gamma*q_hat*v_bar^2*v_beta^2 + 480*I_b*K_beta*Omega^2*gamma^2*theta_0*u_bar*v_bar^2 + 180*I_b*K_beta*Omega^2*gamma^2*theta_1c*u_bar*v_bar^3 + 180*I_b*K_beta*Omega^2*gamma^2*theta_1c*u_bar^3*v_bar - 1440*I_b*K_beta*Omega^2*gamma*theta_1c*u_bar^2*v_beta^2 + 340*I_b*K_beta*Omega^2*gamma^2*theta_tw*u_bar*v_bar^2 - 4320*I_b*K_beta*Omega^2*gamma*theta_1c*v_bar^2*v_beta^2 + 2880*I_b^2*Omega^4*gamma*p_hat*u_bar*v_bar*v_beta^4 + 1440*I_b^2*Omega^4*gamma*theta_1s*u_bar*v_bar*v_beta^2 - 1440*I_b^2*Omega^4*gamma*theta_1s*u_bar*v_bar*v_beta^4 + 1440*I_b*K_beta*Omega^2*gamma*theta_1s*u_bar*v_bar + 90*I_b^2*Omega^4*gamma^2*theta_1s*u_bar^2*v_bar^2*v_beta^2 + 90*I_b*K_beta*Omega^2*gamma^2*theta_1s*u_bar^2*v_bar^2 + 340*I_b^2*Omega^4*gamma^2*q_hat*u_bar*v_bar*v_beta^2 + 680*I_b^2*Omega^4*gamma^2*theta_1c*u_bar*v_bar*v_beta^2 + 1920*I_b*K_beta*Omega^2*beta_P*gamma*u_bar*v_beta^2 + 5760*I_b*K_beta*Omega^2*gamma*lambda*v_bar*v_beta^2 + 340*I_b*K_beta*Omega^2*gamma^2*q_hat*u_bar*v_bar + 680*I_b*K_beta*Omega^2*gamma^2*theta_1c*u_bar*v_bar - 7680*I_b*K_beta*Omega^2*gamma*theta_0*v_bar*v_beta^2 - 5760*I_b*K_beta*Omega^2*gamma*theta_tw*v_bar*v_beta^2 + 2880*I_b*K_beta*Omega^2*gamma*p_hat*u_bar*v_bar*v_beta^2 - 2880*I_b*K_beta*Omega^2*gamma*theta_1s*u_bar*v_bar*v_beta^2))/(45*(I_b*Omega^2*v_beta^2 + K_beta)*(- I_b^2*Omega^4*gamma^2*u_bar^4 - 2*I_b^2*Omega^4*gamma^2*u_bar^2*v_bar^2 - I_b^2*Omega^4*gamma^2*v_bar^4 + 4*I_b^2*Omega^4*gamma^2 + 256*I_b^2*Omega^4*v_beta^4 - 512*I_b^2*Omega^4*v_beta^2 + 256*I_b^2*Omega^4 + 512*I_b*K_beta*Omega^2*v_beta^2 - 512*I_b*K_beta*Omega^2 + 256*K_beta^2));
beta_1s = (I_b*Omega^2*(1440*K_beta^2*gamma*p_hat - 23040*K_beta^2*q_hat*v_beta^2 + 1440*K_beta^2*gamma*theta_1s + 1920*K_beta^2*beta_P*gamma*v_bar - 2880*K_beta^2*gamma*lambda*u_bar + 3840*K_beta^2*gamma*theta_0*u_bar + 2880*K_beta^2*gamma*theta_tw*u_bar + 2160*K_beta^2*gamma*theta_1s*u_bar^2 + 720*K_beta^2*gamma*theta_1s*v_bar^2 + 23040*I_b^2*Omega^4*q_hat*v_beta^4 - 23040*I_b^2*Omega^4*q_hat*v_beta^6 + 40*I_b^2*Omega^4*gamma^3*lambda*u_bar + 320*I_b^2*Omega^4*gamma^2*lambda*v_bar - 1440*I_b^2*Omega^4*gamma*p_hat*v_beta^2 + 4320*I_b^2*Omega^4*gamma*p_hat*v_beta^4 - 30*I_b^2*Omega^4*gamma^3*theta_0*u_bar - 24*I_b^2*Omega^4*gamma^3*theta_tw*u_bar - 240*I_b^2*Omega^4*gamma^2*theta_0*v_bar - 1440*I_b^2*Omega^4*gamma*theta_1s*v_beta^2 + 1440*I_b^2*Omega^4*gamma*theta_1s*v_beta^4 - 192*I_b^2*Omega^4*gamma^2*theta_tw*v_bar - 1440*I_b*K_beta*Omega^2*gamma*p_hat - 1440*I_b*K_beta*Omega^2*gamma*theta_1s - 20*I_b^2*Omega^4*gamma^3*lambda*u_bar^3 - 20*I_b^2*Omega^4*gamma^3*p_hat*u_bar^2 + 10*I_b^2*Omega^4*gamma^3*p_hat*u_bar^4 - 160*I_b^2*Omega^4*gamma^2*q_hat*v_bar^2 + 180*I_b^2*Omega^4*gamma^2*q_hat*v_beta^2 - 15*I_b^2*Omega^4*gamma^3*theta_0*u_bar^3 + 15*I_b^2*Omega^4*gamma^3*theta_0*u_bar^5 - 40*I_b^2*Omega^4*gamma^3*theta_1s*u_bar^2 + 20*I_b^2*Omega^4*gamma^3*theta_1s*u_bar^4 - 8*I_b^2*Omega^4*gamma^3*theta_tw*u_bar^3 + 10*I_b^2*Omega^4*gamma^3*theta_tw*u_bar^5 - 240*I_b^2*Omega^4*gamma^2*theta_0*v_bar^3 - 320*I_b^2*Omega^4*gamma^2*theta_1c*v_bar^2 + 180*I_b^2*Omega^4*gamma^2*theta_1c*v_beta^2 - 160*I_b^2*Omega^4*gamma^2*theta_tw*v_bar^3 + 1440*K_beta^2*gamma*theta_1c*u_bar*v_bar + 180*I_b*K_beta*Omega^2*gamma^2*q_hat + 180*I_b*K_beta*Omega^2*gamma^2*theta_1c + 23040*I_b*K_beta*Omega^2*q_hat*v_beta^2 - 46080*I_b*K_beta*Omega^2*q_hat*v_beta^4 - 20*I_b^2*Omega^4*gamma^3*lambda*u_bar*v_bar^2 - 680*I_b^2*Omega^4*gamma^2*lambda*v_bar*v_beta^2 - 1440*I_b^2*Omega^4*gamma*p_hat*u_bar^2*v_beta^4 + 10*I_b^2*Omega^4*gamma^3*q_hat*u_bar*v_bar^3 + 10*I_b^2*Omega^4*gamma^3*q_hat*u_bar^3*v_bar + 1440*I_b^2*Omega^4*gamma*p_hat*v_bar^2*v_beta^4 - 240*I_b^2*Omega^4*gamma^2*theta_0*u_bar^2*v_bar - 15*I_b^2*Omega^4*gamma^3*theta_0*u_bar*v_bar^2 + 15*I_b^2*Omega^4*gamma^3*theta_0*u_bar*v_bar^4 + 20*I_b^2*Omega^4*gamma^3*theta_1c*u_bar*v_bar^3 + 20*I_b^2*Omega^4*gamma^3*theta_1c*u_bar^3*v_bar - 2160*I_b^2*Omega^4*gamma*theta_1s*u_bar^2*v_beta^2 + 2160*I_b^2*Omega^4*gamma*theta_1s*u_bar^2*v_beta^4 - 160*I_b^2*Omega^4*gamma^2*theta_tw*u_bar^2*v_bar - 8*I_b^2*Omega^4*gamma^3*theta_tw*u_bar*v_bar^2 + 10*I_b^2*Omega^4*gamma^3*theta_tw*u_bar*v_bar^4 + 720*I_b^2*Omega^4*gamma^2*theta_0*v_bar*v_beta^2 - 720*I_b^2*Omega^4*gamma*theta_1s*v_bar^2*v_beta^2 + 720*I_b^2*Omega^4*gamma*theta_1s*v_bar^2*v_beta^4 + 552*I_b^2*Omega^4*gamma^2*theta_tw*v_bar*v_beta^2 - 240*I_b*K_beta*Omega^2*beta_P*gamma^2*u_bar - 680*I_b*K_beta*Omega^2*gamma^2*lambda*v_bar + 5760*I_b*K_beta*Omega^2*gamma*p_hat*v_beta^2 - 2160*I_b*K_beta*Omega^2*gamma*theta_1s*u_bar^2 + 720*I_b*K_beta*Omega^2*gamma^2*theta_0*v_bar - 720*I_b*K_beta*Omega^2*gamma*theta_1s*v_bar^2 + 2880*I_b*K_beta*Omega^2*gamma*theta_1s*v_beta^2 + 552*I_b*K_beta*Omega^2*gamma^2*theta_tw*v_bar - 180*I_b^2*Omega^4*gamma^2*lambda*v_bar^3*v_beta^2 + 10*I_b^2*Omega^4*gamma^3*p_hat*u_bar^2*v_bar^2 - 90*I_b^2*Omega^4*gamma^2*q_hat*u_bar^2*v_beta^2 + 250*I_b^2*Omega^4*gamma^2*q_hat*v_bar^2*v_beta^2 + 30*I_b^2*Omega^4*gamma^3*theta_0*u_bar^3*v_bar^2 + 20*I_b^2*Omega^4*gamma^3*theta_1s*u_bar^2*v_bar^2 - 45*I_b^2*Omega^4*gamma^2*theta_1c*u_bar^4*v_beta^2 + 20*I_b^2*Omega^4*gamma^3*theta_tw*u_bar^3*v_bar^2 + 480*I_b^2*Omega^4*gamma^2*theta_0*v_bar^3*v_beta^2 + 680*I_b^2*Omega^4*gamma^2*theta_1c*v_bar^2*v_beta^2 + 135*I_b^2*Omega^4*gamma^2*theta_1c*v_bar^4*v_beta^2 + 340*I_b^2*Omega^4*gamma^2*theta_tw*v_bar^3*v_beta^2 + 120*I_b*K_beta*Omega^2*beta_P*gamma^2*u_bar^3 - 180*I_b*K_beta*Omega^2*gamma^2*lambda*v_bar^3 - 90*I_b*K_beta*Omega^2*gamma^2*q_hat*u_bar^2 + 250*I_b*K_beta*Omega^2*gamma^2*q_hat*v_bar^2 - 45*I_b*K_beta*Omega^2*gamma^2*theta_1c*u_bar^4 + 480*I_b*K_beta*Omega^2*gamma^2*theta_0*v_bar^3 + 680*I_b*K_beta*Omega^2*gamma^2*theta_1c*v_bar^2 + 135*I_b*K_beta*Omega^2*gamma^2*theta_1c*v_bar^4 + 340*I_b*K_beta*Omega^2*gamma^2*theta_tw*v_bar^3 + 2880*I_b^2*Omega^4*gamma*lambda*u_bar*v_beta^2 - 2880*I_b^2*Omega^4*gamma*lambda*u_bar*v_beta^4 - 160*I_b^2*Omega^4*gamma^2*p_hat*u_bar*v_bar - 20*I_b^2*Omega^4*gamma^3*q_hat*u_bar*v_bar - 40*I_b^2*Omega^4*gamma^3*theta_1c*u_bar*v_bar - 3840*I_b^2*Omega^4*gamma*theta_0*u_bar*v_beta^2 + 3840*I_b^2*Omega^4*gamma*theta_0*u_bar*v_beta^4 - 320*I_b^2*Omega^4*gamma^2*theta_1s*u_bar*v_bar - 2880*I_b^2*Omega^4*gamma*theta_tw*u_bar*v_beta^2 + 2880*I_b^2*Omega^4*gamma*theta_tw*u_bar*v_beta^4 - 1920*I_b*K_beta*Omega^2*beta_P*gamma*v_bar + 2880*I_b*K_beta*Omega^2*gamma*lambda*u_bar - 3840*I_b*K_beta*Omega^2*gamma*theta_0*u_bar - 2880*I_b*K_beta*Omega^2*gamma*theta_tw*u_bar - 180*I_b^2*Omega^4*gamma^2*lambda*u_bar^2*v_bar*v_beta^2 + 480*I_b^2*Omega^4*gamma^2*theta_0*u_bar^2*v_bar*v_beta^2 + 180*I_b^2*Omega^4*gamma^2*theta_1s*u_bar*v_bar^3*v_beta^2 + 180*I_b^2*Omega^4*gamma^2*theta_1s*u_bar^3*v_bar*v_beta^2 + 340*I_b^2*Omega^4*gamma^2*theta_tw*u_bar^2*v_bar*v_beta^2 + 120*I_b*K_beta*Omega^2*beta_P*gamma^2*u_bar*v_bar^2 - 180*I_b*K_beta*Omega^2*gamma^2*lambda*u_bar^2*v_bar - 1440*I_b*K_beta*Omega^2*gamma*p_hat*u_bar^2*v_beta^2 + 1440*I_b*K_beta*Omega^2*gamma*p_hat*v_bar^2*v_beta^2 + 480*I_b*K_beta*Omega^2*gamma^2*theta_0*u_bar^2*v_bar + 180*I_b*K_beta*Omega^2*gamma^2*theta_1s*u_bar*v_bar^3 + 180*I_b*K_beta*Omega^2*gamma^2*theta_1s*u_bar^3*v_bar + 4320*I_b*K_beta*Omega^2*gamma*theta_1s*u_bar^2*v_beta^2 + 340*I_b*K_beta*Omega^2*gamma^2*theta_tw*u_bar^2*v_bar + 1440*I_b*K_beta*Omega^2*gamma*theta_1s*v_bar^2*v_beta^2 - 2880*I_b^2*Omega^4*gamma*q_hat*u_bar*v_bar*v_beta^4 - 1440*I_b^2*Omega^4*gamma*theta_1c*u_bar*v_bar*v_beta^2 + 1440*I_b^2*Omega^4*gamma*theta_1c*u_bar*v_bar*v_beta^4 - 1440*I_b*K_beta*Omega^2*gamma*theta_1c*u_bar*v_bar + 90*I_b^2*Omega^4*gamma^2*theta_1c*u_bar^2*v_bar^2*v_beta^2 + 90*I_b*K_beta*Omega^2*gamma^2*theta_1c*u_bar^2*v_bar^2 + 340*I_b^2*Omega^4*gamma^2*p_hat*u_bar*v_bar*v_beta^2 + 680*I_b^2*Omega^4*gamma^2*theta_1s*u_bar*v_bar*v_beta^2 + 1920*I_b*K_beta*Omega^2*beta_P*gamma*v_bar*v_beta^2 - 5760*I_b*K_beta*Omega^2*gamma*lambda*u_bar*v_beta^2 + 340*I_b*K_beta*Omega^2*gamma^2*p_hat*u_bar*v_bar + 7680*I_b*K_beta*Omega^2*gamma*theta_0*u_bar*v_beta^2 + 680*I_b*K_beta*Omega^2*gamma^2*theta_1s*u_bar*v_bar + 5760*I_b*K_beta*Omega^2*gamma*theta_tw*u_bar*v_beta^2 - 2880*I_b*K_beta*Omega^2*gamma*q_hat*u_bar*v_bar*v_beta^2 + 2880*I_b*K_beta*Omega^2*gamma*theta_1c*u_bar*v_bar*v_beta^2))/(45*(I_b*Omega^2*v_beta^2 + K_beta)*(- I_b^2*Omega^4*gamma^2*u_bar^4 - 2*I_b^2*Omega^4*gamma^2*u_bar^2*v_bar^2 - I_b^2*Omega^4*gamma^2*v_bar^4 + 4*I_b^2*Omega^4*gamma^2 + 256*I_b^2*Omega^4*v_beta^4 - 512*I_b^2*Omega^4*v_beta^2 + 256*I_b^2*Omega^4 + 512*I_b*K_beta*Omega^2*v_beta^2 - 512*I_b*K_beta*Omega^2 + 256*K_beta^2));

beta0_rad = beta_0;
beta_1c_rad = beta_1c;
beta_1s_rad = beta_1s;

Beta = @(psi) beta0_rad + beta_1c_rad*cos(psi) + beta_1s_rad*sin(psi);

Beta_dot = @(psi) (-beta_1c_rad*sin(psi) + beta_1s_rad*cos(psi))*Omega;

end