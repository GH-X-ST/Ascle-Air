function cost = Trim_H_C(theta_LR, theta_UR)
% Author:
%   Hanchen Li (hl3422@ic.ac.uk)
%
% Inputs:
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%
% Outputs:
%   X         - Force along x body axis (+ve forward) (N)
%   Y         - Force along y body axis (+ve to starboard) (N)
%   Z         - Force along z body axis (+ve downward) (N)
%   L         - Rolling moment (+ve port side moves up) (Nm)
%   M         - Pitching moment (+ve nose up) (Nm)
%   N         - Yawing moment (+ve nose turns to the starboard) (Nm)

%% 0 Basic Parameters
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
state = [0, 0, 0, 0, 0, 0, 0, 0, 0];

%   ctrl      - [theta_LR, theta_UR, theta_1s, theta_1c, delta_E, delta_R]
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_1s  - Longitudinal cyclic (rad)
%   theta_1c  - Lateral cyclic (rad)
%   delta_E   - Elevator deflection (rad)
%   delta_R   - Rudder deflection (rad)
ctrl = [theta_LR, theta_UR, 0, 0, 0, 0];

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
params = [1524, 3491, 9.80665, 1, 0, 0, 5.5, 5.5];

%% 1 Calling function
[X, Y, Z, L, M, N, ~, ~, ~, ~] = FM_N(state, ctrl, params);

%% 2 Cost function
Weighting = diag([1, 1, 1, 1, 1, 1]);

cost = sqrt([X; Y; Z; L; M; N]' * Weighting * [X; Y; Z; L; M; N]);

end