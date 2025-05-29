function [c, ceq] = Nonlinear_H_C(x)
% Introduction:
% This function is used to define the non-linear constraint values.
%
% Author:
%   Hanchen Li (hl3422@ic.ac.uk)
%
% Inputs:
%   x         - (theta_LR, theta_UR, theta_1s, theta_1c, phi, theta)
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_1s  - Longitudinal cyclic (rad)
%   theta_1c  - Lateral cyclic (rad)
%   phi       - Euler roll angle (rad)
%   theta     - Euler pitch angle (rad)
%
% Outputs:
%   c         - Inequality constriant
%   ceq       - Equality constraint

%% 0 Basic Parameters
%   x         - (theta_LR, theta_UR, theta_1s, theta_1c, phi, theta)
theta_LR  = x(1);
theta_UR  = x(2);
theta_1s  = x(3);
theta_1c  = x(4);
phi       = x(5);
theta     = x(6);

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
state = [0, 0, 0, 0, 0, 0, phi, theta, 0];

%   ctrl      - [theta_LR, theta_UR, theta_1s, theta_1c, delta_E, delta_R]
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_1s  - Longitudinal cyclic (rad)
%   theta_1c  - Lateral cyclic (rad)
%   delta_E   - Elevator deflection (rad)
%   delta_R   - Rudder deflection (rad)
ctrl = [theta_LR, theta_UR, theta_1s, theta_1c, 0, 0];

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
params = [1524, 3500, 9.80665, 1, 0, 0, 5.5, 5.5];

%% 1 Calling function
[X, Y, Z, L, M, N, ~, ~, ~, ~] = FM_N(state, ctrl, params);

%% 3 Define the inequality constriant

% 3.1 Define the inequality constriant
% Torlerence
tol = 1e-4;

c = [
    abs(X) - tol;
    abs(Y) - tol;
    abs(Z) - tol;
    abs(L) - tol;
    abs(M) - tol;
    abs(N) - tol
];

% 3.2 Define the equality constraint
ceq = [];

end