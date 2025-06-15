function [A, B, A_lon, B_lon, A_lat, B_lat] = Linearisation(state, ctrl, params, stepsize)
% Introduction:
%   This is the provisional linearisation code.
%
% Author:
%   Hanchen Li (hl3422@ic.ac.uk)
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
%   stepsize  - finite-difference step size (tunable for accuracy)
%
% Outputs:
%   A         - System matrix (9x9)
%   B         - Control matrix (6x9)

%% 0 Basic Parameters

% 0.0 Setup
delta = stepsize; % finite-difference step size (tunable for accuracy)
nx = 9;
nu = 6;
A = zeros(nx, nx);
B = zeros(nx, nu);
state  = state(:);

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

% for linearisation
theta_0 = 0.5 * (theta_LR + theta_UR);
theta_diff = 0.5 * (theta_UR - theta_LR);
ctrl_trans = [theta_0, theta_diff, theta_1s, theta_1c, delta_E, delta_R];

% params - [h, m, g, z_R, z_HT, z_VT, x_HT, x_VT]
h = params(1);
m = params(2);
g = params(3);
z_R = params(4);
z_HT = params(5);
z_VT = params(6);
x_HT = params(7);
x_VT = params(8);

%% 1 Linearisation

% 1.1  Derivative at the nominal point (optional, kept for reference)
f0 = stateDerivative(state, ctrl_trans, params);

% 1.2  Build A  (perturb each state component)
for i = 1:nx
    dx        = zeros(nx,1);   % 9 × 1 column – same orientation as ‘state’
    dx(i)     = delta;
    f_up      = stateDerivative(state + dx, ctrl_trans, params);
    f_down    = stateDerivative(state - dx, ctrl_trans, params);
    A(:,i)    = (f_up - f_down) / (2 * delta);   % dx(i) == delta
end

% 1.3  Build B  (perturb each control input)
for j = 1:nu
    du        = zeros(1,nu);   % 1 × 6 row – matches ‘ctrl_trans’
    du(j)     = delta;
    f_up      = stateDerivative(state, ctrl_trans + du, params);
    f_down    = stateDerivative(state, ctrl_trans - du, params);
    B(:,j)    = (f_up - f_down) / (2 * delta);   % du(j) == delta
end

% 1.4 Extract matrices

idxLon = [1  3  5  8];
idxLat = [2  4  6  7  9];
idxCLon = [1  3  5];
idxCLat = [2  4  6];

A_lon = A(idxLon, idxLon);
A_lat = A(idxLat, idxLat);

B_lon = B(idxLon, idxCLon);
B_lat = B(idxLat, idxCLat);

end

