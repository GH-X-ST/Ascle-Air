function [X, Y, Z, L, M, N, I_xx, I_yy, I_zz, I_xz] = FM(state, ctrl, params)
% Introduction:
%   This function computes the aerodynamic forces and moments acting on the
%   helicopter in body axes, considering contributions from two coaxial
%   rotors, gravity and fuselage.
%
% Author:
%   Hanchen Li (hl3422@ic.ac.uk)
%
% References:
%   \bibitem{Padfield} Rotor blade element and flapping dynamics 
%   \bibitem{Coleman}  Coaxial rotor aerodynamic research
%   \bibitem{Leishman} BET/BEMT for rotor thrust/torque
%   \bibitem{Bramwell} Helicopter dynamics
%
% Inputs:
%   state     - [u, v, w, p, q, r, phi, theta, psi]
%   u         - Forward velocity (m/s)
%   v         - Lateral velocity (m/s)
%   w         - Vertical velocity (m/s)
%   p         - Roll rate (rad/s)
%   q         - Pitch rate (rad/s)
%   r         - Yaw rate (rad/s)
%   phi       - Euler roll angle (rad)
%   theta     - Euler pitch angle (rad)
%   psi       - Euler yaw angle (rad)
%
%   ctrl      - [theta_LR, theta_UR, theta_1s, theta_1c, delta_E, delta_R]
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_1s  - Longitudinal cyclic (rad)
%   theta_1c  - Lateral cyclic (rad)
%   delta_E   - Elevator deflection (rad)
%   delta_R   - Rudder deflection (rad)
%
%   params    - [h, m, g, z_R, z_HT, z_VT, x_HT, x_VT]
%   h         - Altitude (m)
%   m         - Mass (kg)
%   g         - Acceleration of gravity (m s^-2)
%   z_R       - Vertical distance between upper rotor and lower rotor (m)
%   z_HT      - Vertical distance from horizontal tailplane aerodynamic
%               centre to the lower rotor hub (+ve below) (m)
%   z_VT      - Vertical distance from vertical tailplane aerodynamic
%               centre to the lower rotor hub (+ve below) (m)
%   x_HT      - Horizontal distance from horizontal tailplane
%               aerodynamic centre to shaft (m)
%   x_VT      - Horizontal distance from vertical tailplane
%               aerodynamic centre to shaft (m)
%
% Outputs:
%   X         - Force along x body axis (+ve forward) (N)
%   Y         - Force along y body axis (+ve to starboard) (N)
%   Z         - Force along z body axis (+ve downward) (N)
%   L         - Rolling moment (+ve port side moves up) (Nm)
%   M         - Pitching moment (+ve nose up) (Nm)
%   N         - Yawing moment (+ve nose turns to the starboard) (Nm)  
%   I_xx      - Roll second moment of inertia (kgm^2)
%   I_yy      - Pitch second moment of inertia (kgm^2)
%   I_zz      - Yaw second moment of inertia (kgm^2)
%   I_xz      - Cross second moment of inertia (kgm^2)

% Conventions:
%   Body axes - x forward, y starboard, z downward
%   Rotors    - Upper rotor rotates anticlockwise when viewed
%               from above, Lower rotor clockwise

%% 0 Basic Parameters

% 0.1 Extract parameters
% state - [u, v, w, p, q, r, phi, theta, psi]
u = state(1);
v = state(2);
w = state(3);
p = state(4);
q = state(5);
r = state(6);
phi = state(7);
theta = state(8);
psi = state(9);
% ctrl - [theta_LR, theta_UR, theta_1s, theta_1c, delta_E, delta_R]
theta_LR = ctrl(1);
theta_UR = ctrl(2);
theta_1s = ctrl(3);
theta_1c = ctrl(4);
delta_E = ctrl(5);
delta_R = ctrl(6);
% params - [h, m, g, z_R, z_HT, z_VT, x_HT, x_VT]
h = params(1);
m = params(2);
g = params(3);
z_R = params(4);
z_HT = params(5);
z_VT = params(6);
x_HT = params(7);
x_VT = params(8);

% 0.2 Initialize outputs
% Forces
X = 0;
Y = 0;
Z = 0;
% Moments
L = 0;
M = 0;
N = 0;
% Second moment of inertia
I_xx = 0;
I_yy = 0;
I_zz = 0;
I_xz = 0;
    
% 0.3 Initialize rotor outputs placeholders
T_UR = 0;      % Upper rotor thrust (N)
T_LR = 0;      % Lower rotor thrust (N)
H_UR = 0;      % Upper rotor drag (N)
H_LR = 0;      % Lower rotor drag (N)
Y_UR = 0;      % Upper rotor side force (N)
Y_LR = 0;      % Lower rotor side force (N)
% Q_diff = 0;    % Q_LR - Q_UR
Q_UR = 0;      % Upper rotor torque (+ve clockwise) (Nm)
Q_LR = 0;      % Lower rotor torque (+ve anticlockwise) (Nm)
beta1c_UR = 0; % Upper rotor longitudinal flapping
               % (+ve disk tilt fore-aft) (rad)
beta1s_UR = 0; % Lower rotor longitudinal flapping
               % (+ve disk tilt fore-aft) (rad)
beta1c_LR = 0; % Upper rotor lateral flapping
               % (+ve port side down) (rad)
beta1s_LR = 0; % Lower rotor lateral flapping
               % (+ve port side down) (rad)
Mx_UR = 0;     % Upper rotor rolling moment (Nm)      
Mx_LR = 0;     % Lower rotor rolling moment (Nm)
My_UR = 0;     % Upper rotor rolling moment (Nm)
My_LR = 0;     % Lower rotor pitching moment (Nm)

% 0.4 Initialize fuselage outputs placeholders
D = 0;         % Fuselage drag (N)
Y_F = 0;       % Fuselage side force (N)
Mx_F = 0;      % Fuselage rolling moment (Nm)
My_F = 0;      % Fuselage pitching moment (Nm)

% 0.5 Initialize empennage outputs placeholders
L_HT = 0;      % Horzizontal tailplane lift (+ve upward) (N)
L_VT = 0;      % Vertical tailplane lift (+ve to right) (N)
D_HT = 0;      % Horzizontal tailplane drag (N)
D_VT = 0;      % Vertical tailplane drag (N)

% 0.6 Initialize centre of gravity
x_CG = 0;
y_CG = 0;
z_CG = 0;

%% 1 Force arms

theta_FP = tan(w / u);

alpha_s = theta_FP + theta;

h_UR = z_R + z_CG;

h_LR = z_CG;

h_VT = z_CG - z_VT;

h_HT = z_CG - z_HT;

l_VT = x_CG + x_VT;

l_HT = x_CG + x_HT;

%% 2 International Standard Atmosphere model
[~, a, ~, rho, ~, mu] = atmosisa(h);

%% 3 System & Structure Model

%% 4 Rotor Aerodynamics & Dynamics Model

% 4.1 Load airfoil data
addpath('Airfoil');
polar = loadPolarData('xf-rc410-il-1000000.txt');

if u == 0
    
    % 4.2 Hover and climb aerodynamics & dynamics

    [T_UR, T_LR, ~, ~, ~, ~] = Stability_hover_getT(u, v, w, a, rho, mu, m, g, polar, theta_UR, theta_LR);

    % [~, ~, beta1c_UR, beta1s_UR] = getFlappingResponse(u, rho, m, CT, theta_UR, theta_LR, theta_1c, theta_1s, 'upper');
    % [~, ~, beta1c_LR, beta1s_LR] = getFlappingResponse(u, rho, m, CT, theta_UR, theta_LR, theta_1c, theta_1s, 'lower');

    % It will return empty value as u = 0

else
    
    % 4.4 Forward flight aerodynamics & dynamics

    [T_UR, T_LR, H_UR, H_LR, Y_UR, Y_LR, Q_UR, Q_LR, Mx_UR, My_UR, Mx_LR, My_LR, ~, ~, beta1c_UR, beta1s_UR, beta1c_LR, beta1s_LR] = Stability_forward(u, v, w, a, rho, mu, m, g, polar, theta_UR, theta_LR, theta_1c, theta_1s);

    Q_UR = -Q_UR; % to reaction torque

    Q_LR = -Q_LR; % to reaction torque

    beta1s_LR = -beta1s_LR; % clockwise rotation

end

%% 5 Fuselage Aerodynamics Model

%% 6 Empennage Aerodynamics Model
% 6.1 Horizontal tailplane
[L_HT, D_HT] = FHT(u, v, w, rho, delta_E);

% 6.1 Vertical tailplane
[L_VT, D_VT] = FVT(u, v, w, rho, delta_R);

%% 7 Forces and Moments
% 7.1 X, Logitudinal Force
X = -D * cos(alpha_s) - L * sin(alpha_s) - 2 * D_VT - D_HT ...
    - (H_LR * cos(beta1c_LR) + H_UR * cos(beta1c_UR)) ...
    + (T_LR * sin(beta1c_LR) * cos(beta1s_LR) ...
    + T_UR * sin(beta1c_UR) * cos(beta1s_UR));

% 7.2 Y, Lateral Force
Y = Y_F - 2 * L_VT ...
    + (Y_LR * cos(beta1s_LR) + Y_UR * cos(beta1s_UR)) ...
    - (T_LR * cos(beta1c_LR) * sin(beta1s_LR) ...
    + T_UR * cos(beta1c_UR) * sin(beta1s_UR));

% 7.3 Z, Lateral Force
Z = D * sin(alpha_s) - L * cos(alpha_s) - L_HT ...
    - (H_LR * sin(beta1c_LR) + H_UR * sin(beta1c_UR)) ...
    - (Y_LR * sin(beta1s_LR) + Y_UR * sin(beta1s_UR)) ...
    - (T_LR * cos(beta1c_LR) * cos(beta1s_LR) ...
    + T_UR * cos(beta1c_UR) * cos(beta1s_UR));

% 7.4 L, Rolling moment
L = Mx_F + Mx_LR + Mx_UR + 2 * L_VT * h_VT + L_HT * y_CG ...
    + ((T_LR * cos(beta1c_LR) * cos(beta1s_LR) ...
    + T_UR * cos(beta1c_UR) * cos(beta1s_UR)) ...
    + (H_LR * sin(beta1c_LR) + H_UR * sin(beta1c_UR)) ...
    + (Y_LR * sin(beta1s_LR) + Y_UR * sin(beta1s_UR))) * y_CG ...
    + (Y_LR * cos(beta1s_LR) - T_LR * cos(beta1c_LR) * sin(beta1s_LR)) * h_LR ...
    + (Y_UR * cos(beta1s_UR) - T_UR * cos(beta1c_UR) * sin(beta1s_UR)) * h_UR;

% 7.5 M, Pitching moment
M = My_F + My_LR + My_UR - L_HT * l_HT + 2 * D_VT * h_VT + D_HT * h_HT ...
    - ((T_LR * cos(beta1c_LR) * cos(beta1s_LR) ...
    + T_UR * cos(beta1c_UR) * cos(beta1s_UR)) ...
    + (H_LR * sin(beta1c_LR) + H_UR * sin(beta1c_UR)) ...
    + (Y_LR * sin(beta1s_LR) + Y_UR * sin(beta1s_UR))) * x_CG ...
    + (H_LR * cos(beta1c_LR) - T_LR * sin(beta1c_LR) * cos(beta1s_LR)) * h_LR ...
    + (H_UR * cos(beta1c_UR) - T_UR * sin(beta1c_UR) * cos(beta1s_UR)) * h_UR;

% 7.6 N, Yawing moment
N = Q_UR - Q_LR - L_VT * l_VT - (2 * D_VT + D_HT) * y_CG ...
    + ((T_LR * cos(beta1c_LR) * sin(beta1s_LR) ...
    + T_UR * cos(beta1c_UR) * sin(beta1s_UR)) ...
    - (Y_LR * cos(beta1s_LR) + Y_UR * cos(beta1s_UR))) * x_CG ...
    + ((T_LR * sin(beta1c_LR) * cos(beta1s_LR) ...
    + T_UR * sin(beta1c_UR) * cos(beta1s_UR)) ...
    - (H_LR * cos(beta1c_LR) + H_UR * cos(beta1c_UR))) * y_CG;

end