function xdot = Non_linear(state, ctrl, params)
% Introduction:
%   This is the provisional non-linear model.
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
% Outputs:
%   xdot      - state derivative

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

%% 1 Calculate force and moments

[X, Y, Z, L, M, N, I_xx, I_yy, I_zz, I_xz] = FM(state, ctrl, params);

%% 2 Calculate derivative terms

% 2.1 Calculate xdot

udot = -(w * q - v * r) + X / m - g * sin(theta);
vdot = -(u * r - w * p) + Y / m + g * cos(theta) * sin(phi);
wdot = -(v * p - u * q) + Z / m + g * cos(theta) * cos(phi);

I = [I_xx, 0,   -I_xz;
     0,   I_yy,  0;
    -I_xz, 0,    I_zz];

omega = [p; q; r];
Mvec  = [L; M; N];

omega_dot = I \ ( Mvec - cross(omega, I * omega) );

pdot = omega_dot(1);
q_dot = omega_dot(2);
rdot = omega_dot(3);

eulerdot = EulerRatesSafe(phi, theta, p, q, r);
phi_dot   = eulerdot(1);
theta_dot = eulerdot(2);
psi_dot   = eulerdot(3);

% 2.2 Cast into vector
xdot = [udot; vdot; wdot; pdot; q_dot; rdot; phi_dot; theta_dot; psi_dot];

%% A Avoid singularity

function eulerdot = EulerRatesSafe(phi, theta, p, q, r)
% Compute eulerdot without blowing up at 90 degrees
%
% Inputs:
%   phi, theta – current Euler angles (rad)
%   p, q, r    – body-axis angular rates (rad/s)
%
% Output:
%   eulerdot   - [phidot, thetadot, psidot]

tol      = 1e-8; % switch to quaternion method if |cosθ| < tol
singular = abs(cos(theta)) < tol;

if ~singular

    phi_dot = p + q * sin(phi) * tan(theta) + r *  cos(phi) * tan(theta);
    theta_dot = q * cos(phi) - r * sin(phi);
    psi_dot = q * sin(phi) * sec(theta) + r * cos(phi) * sec(theta);

    eulerdot  = [phi_dot; theta_dot; psi_dot];

else

    dt  = 1e-9;                    % tiny step for finite difference

    % Build a quaternion with zero yaw
    q0  = eul2quat([0 theta phi]); % keep psi = 0 because only derivatives matter

    % body-rate matrix
    Omega = 0.5 * [   0,  -p,  -q,  -r;
                      p,   0,   r,  -q;
                      q,  -r,   0,   p;
                      r,   q,  -p,   0 ];

    q_dot = Omega * q0(:);

    q1   = q0(:) + q_dot * dt;

    q1   = q1 / norm(q1);              % normalise

    eul0 = quat2eul(q0(:).');          % back to [phi theta psi]

    eul1 = quat2eul(q1(:).');

    eulerdot = (eul1 - eul0).' / dt;   % finite-difference derivative

end

end

end

