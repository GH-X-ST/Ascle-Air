function [F, M] = FM_h(state, ctrl, params)
% Introduction:
%   This function computes the net forces and moments acting on the
%   helicopter in body axes, considering contributions from two coaxial
%   rotors, gravity and fuselage.
%   It is tailored for hover and low-speed forward flight (advance ratio
%   μ small, 0–40 knots).
%
% References:
%   \bibitem{Padfield} Rotor blade element and flapping dynamics 
%   \bibitem{Coleman}  Coaxial rotor aerodynamic research
%   \bibitem{Leishman} BET/BEMT for rotor thrust/torque
%   \bibitem{Bramwell} Helicopter dynamics
%
% Inputs
%   state.U         - Body-axis velocities (m/s) along X (forward)
%   state.V         - Body-axis velocities (m/s) along Y (starboard)
%   state.W         - Body-axis velocities (m/s) along Z (down)
%
%   ctrl.coll_avg   - Average collective pitch (rad) for both rotors
%   ctrl.coll_diff  - Differential collective (rad)
%                     (adds to top rotor, subtracts from bottom rotor)
%   ctrl.cyclic_lon - Longitudinal cyclic (rad)
%   ctrl.cyclic_lat - Lateral cyclic (rad)
%   ctrl.r          - Rudder deflection (rad)
%
%   params.R        - Rotor radius (m)
%   params.Nb       - Number of blades per rotor
%   params.theta_tw - Blade linear twist (rad)
%                     (negative for outwash, assume θ_tw < 0)
%   params.Omega    - Rotor angular speed (rad/s)
%   params.rho      - Air density (kg/m^3)
%   params.I_beta   - Blade flapping moment of inertia about hinge (kg·m^2)
%   params.e        - Hinge offset or effective flapping hinge offset (m)
%   params.mass     - Helicopter mass (kg)
%   params.g        - Gravitational acceleration (m/s^2)
%   params.Cd0_f    - Fuselage drag coefficient (frontal area reference)
%   params.A_f      - Fuselage frontal area (m^2)
%   params.S_vt     - Vertical tail (rudder) area (m^2)
%   params.l_vt     - Level arm from CG to vertical tail (m)
%   params.Cn_delta - Rudder yaw moment coefficient (per rad deflection)   
%
% (array):
%   params.c        - Rotor radius (m)
%   params.a0       - Lift-curve slope of blade airfoil (1/rad)
%   params.Cd0      - Blade section profile drag coefficient
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
%                     N (yaw) positive nose-right
%
% Note:
%   For simplicity, this model assumes small angles (β1c, β1s << 1 rad)
%   and uniform inflow across the rotor disk.

%% 0 Basic Parameters
% Extract parameters
% (double):
R        = params.R;
Nb       = params.Nb;
theta_tw = params.theta_tw;
Omega    = params.Omega;
rho      = params.rho;
I_beta   = params.I_beta;
e        = params.e;
m        = params.mass;
g        = params.g;
Cd0_f    = params.Cd0_f;
A_f      = params.A_f;
S_vt     = params.S_vt;
l_vt     = params.l_vt;
Cn_delta = params.Cn_delta;
% (linear):
c        = params.c;
a0       = params.a0;
Cd0      = params.Cd0;
N        = length(c); % find the length of the array

% Derived parameters
% (double):
A        = pi * R^2;                          % Rotor disc area (m^2)
% (array):
sigma    = (Nb .* c) ./ (pi .* R);            % Rotor solidity
gamma    = (rho .* c .* a0 .* R^4) ./ I_beta; % Lock number
    
% Control inputs (collective for each rotor)
theta0_top  = ctrl.coll_avg + 0.5 * ctrl.coll_diff;
% collective pitch top rotor (rad)
theta0_bot  = ctrl.coll_avg - 0.5 * ctrl.coll_diff;
% collective pitch bottom rotor (rad)
    
% Body-axis velocities (for low-speed inflow adjustment)
U = state.U; % forward velocity (m/s)
V = state.V; % sideward velocity (starboard positive) (m/s)
W = state.W; % vertical velocity (down positive) (m/s)
    
% Advance ratios (non-dimensional horizontal velocity components)
mu_x = U / (Omega * R);        % longitudinal advance ratio
mu_y = V / (Omega * R);        % lateral advance ratio
mu   = sqrt(mu_x^2 + mu_y^2);  % total advance ratio

% Initialize outputs
F = [0; 0; 0];   % forces (X, Y, Z)
M = [0; 0; 0];   % moments (L, M, N)
    
% Initialize rotor outputs placeholders
T_top = 0;
T_bottom = 0;       % thrust (N)
Q_top = 0;
Q_bottom = 0;       % torque (Nm)
beta1c_top = 0;
beta1s_top = 0;     % top rotor flapping angles (rad)
beta1c_bottom = 0;
beta1s_bottom = 0;  % bottom rotor flapping angles (rad)

%% 1. Compute uniform inflow for each rotor




end
