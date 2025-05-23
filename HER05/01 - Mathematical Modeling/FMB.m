function [F, M] = FMB(state, ctrl, params)
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
%   ctrl      - [theta_LR, theta_UR, theta_ls, theta_lc, delta_E, delta_R]
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_ls  - Longitudinal cyclic (rad)
%   theta_lc  - Lateral cyclic (rad)
%   delta_E   - Elevator deflection (rad)
%   delta_R   - Rudder deflection (rad)
%
%   params    - [h, m, g, h_LR, h_UR, h_HT, h_VT, x_HT, x_VT]
%   h         - Altitude (m)
%   m         - Mass (kg)
%   g         - Acceleration of gravity (m s^-2)
%   h_LR      - Vertical distance from lower rotor to z_CG (+ve above) (m)
%   h_UR      - Vertical distance from upper rotor to z_CG (+ve above) (m)
%   h_HT      - Vertical distance from horizontal tailplane
%               to z_CG (+ve above) (m)
%   h_VT      - Vertical distance from vertical tailplane
%               to z_CG (+ve above) (m)
%   x_HT      - Horizontal distance from horizontal tailplane
%               aerodynamic centre to shaft (m)
%   x_VT      - Horizontal distance from vertical tailplane
%               aerodynamic centre to shaft (m)
%
% Outputs:
%   F = [X; Y; Z]   - Total forces in body axes (N)
%   M = [L; M; N]   - Total moments in body axes (Nm)
%
% Conventions:
%   Body axes       - X forward, Y starboard, Z downward
%   Rotors          - Top rotor rotates counter-clockwise when viewed from
%                     above, bottom rotor clockwise
%   Flapping angles - β1c = longitudinal flapping (disk tilt fore-aft)
%                     Positive β1c tilts the rotor disk back
%                     (gives nose-up moment)
%                     β1s = lateral flapping (disk tilt side-to-side)
%                     Positive β1s tilts the disk to the left 
%                     (port side down)
%   Moments         - L (roll) positive for right-wing-down
%                     M (pitch) positive nose-up                          
%                     N (yaw) positive nose                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             -right

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
% ctrl - [theta_avg, theta_diff, theta_ls, theta_lc]
theta_avg = ctrl(1);
theta_diff = ctrl(2);
theta_ls = ctrl(3);
theta_lc = ctrl(4);
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
F = [0; 0; 0]; % Forces (X, Y, Z)
M = [0; 0; 0]; % Moments (L, M, N)
    
% 0.3 Initialize rotor outputs placeholders
T_UR = 0;
T_LR = 0;  % Rotor thrust (N)
H_UR = 0;
H_LR = 0;  % Rotor drag (N)
Y_UR = 0;
Y_LR = 0;  % Rotor side force (N)
Q_UR = 0;
Q_LR = 0;  % Rotor torque (Nm)
Mx_UR = 0;
Mx_LR = 0; % Rotor rolling moment (Nm)
My_UR = 0;
My_LR = 0; % Rotor pitching moment (Nm)
D = 0;     % Fuselage drag (N)
Y_F = 0;   % Fuselage side force (N)
Mx_F = 0;  % Fuselage rolling moment (Nm)
My_F = 0;  % Fuselage pitching moment (Nm)
L_HT = 0;  % Horzizontal tailplane lift (+ve upward) (N)
L_VT = 0;  % Vertical tailplane lift (+ve to right) (N)
x_CG = 0;  % Horizontal centre of gravity (+ve aft of shaft) (kg)
y_CG = 0;  % Lateral centre of gravity (+ve right of shaft) (kg)
z_CG = 0;  % vertical centre of gravity (+ve right of shaft) (kg)


%% 1 Input Functions
% 1.1 Hover
if config(1) == 1

end

% 1.2 Forward flight
if config(1) == 2

end

%% 2. Basic parameters
 
%% 3. Forces
% 3.1 X, Logitudinal Force

(T_upper + )


end
