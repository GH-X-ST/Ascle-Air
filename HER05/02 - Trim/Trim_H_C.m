function [X, Y, Z, L, M, N] = Trim_H_C(theta_LR, theta_UR, theta_ls, theta_lc, phi, theta)
% Introduction:
%   This function computes the net forces and moments acting on the
%   helicopter in body axes during hover and forward flight, considering
%   contributions from two coaxial rotors, gravity and fuselage.
%
% References:
%   \bibitem{Padfield} Rotor blade element and flapping dynamics 
%   \bibitem{Coleman}  Coaxial rotor aerodynamic research
%   \bibitem{Leishman} BET/BEMT for rotor thrust/torque
%   \bibitem{Bramwell} Helicopter dynamics
%
% Inputs:
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_ls  - Longitudinal cyclic (rad)
%   theta_lc  - Lateral cyclic (rad)
%   delta_E   - Elevator deflection (rad)
%   delta_R   - Rudder deflection (rad)
%   phi       - Euler roll angle (rad)
%   theta     - Euler pitch angle (rad)
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
states = [0, 0, 0, 0, 0, 0, 0, 0, 0];

%   ctrl      - [theta_LR, theta_UR, theta_ls, theta_lc, delta_E, delta_R]
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_ls  - Longitudinal cyclic (rad)
%   theta_lc  - Lateral cyclic (rad)
%   delta_E   - Elevator deflection (rad)
%   delta_R   - Rudder deflection (rad)
ctrl = [0, 0, 0, 0, 0, 0];

%   % params  - [h, m, g, z_R, z_HT, z_VT, x_HT, x_VT]
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
params = [0, 600, 9.80665, 0, 0, 0, 0, 0];


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
[~, ~, ~, rho, ~, ~] = atmosisa(h);

%% 3 System & Structure Model

%% 4 Rotor Aerodynamic Model
if u == 0

    if w == 0 % 2.1 Hover

        T_UR = 0;

    else      % 2.2 Climb

        T_UR = 0;

    end

else          % 1.3 Forward flight

    T_UR = 0;

end

%% 5 Fuselage Aerodynamic Model

%% 6 Empennage Aerodynamic Model

%% 7 Forces and Moments
% 7.1 X, Logitudinal Force
X = -m * g * sin(theta) - D * cos(alpha_s) - 2 * D_VT - D_HT ...
    - (H_LR * cos(beta1c_LR) + H_UR * cos(beta1c_UR)) ...
    + (T_LR * sin(beta1c_LR) * cos(beta1s_LR) ...
    + T_UR * sin(beta1c_UR) * cos(beta1s_UR));

% 7.2 Y, Lateral Force
Y = m * g * sin(phi) * cos(theta) + Y_F - 2 * L_VT ...
    + (Y_LR * cos(beta1s_LR) + Y_UR * cos(beta1s_UR)) ...
    - (T_LR * cos(beta1c_LR) * sin(beta1s_LR) ...
    + T_UR * cos(beta1c_UR) * sin(beta1s_UR));

% 7.3 Z, Lateral Force
Z = m * g * cos(phi) * cos(theta) + D * sin(alpha_s) - L_HT ...
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
