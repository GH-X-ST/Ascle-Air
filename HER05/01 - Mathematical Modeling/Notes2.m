function [F, M] = forces_and_moments(state, controls, params)
% state: structure containing flight state (e.g., forward velocity V, air density rho, etc.)
% controls: structure with pilot inputs (collective, cyclic, rudder deflection)
% params: structure with helicopter parameters (rotor geometry, airframe data, etc.)
%
% Outputs:
%   F = [X; Y; Z] total forces in body axes (N)
%   M = [L; M; N] total moments about body axes (N·m)
%
% References: Padfield (2007):contentReference[oaicite:6]{index=6}, NASA TP-3675 (Coleman, 1997):contentReference[oaicite:7]{index=7}:contentReference[oaicite:8]{index=8},
%             blade element theory:contentReference[oaicite:9]{index=9}, and other sources as cited below.

    % Unpack state and control variables for clarity
    V     = state.V;           % forward speed (m/s) – assume small (μ << 1)
    rho   = state.rho;         % air density (kg/m^3)
    % No side velocity or wind assumed (hover/low-speed forward flight)

    delta0   = controls.collective;   % collective pitch (rad)
    delta_lon = controls.lon_cyclic;  % longitudinal cyclic input (rad)
    delta_lat = controls.lat_cyclic;  % lateral cyclic input (rad)
    delta_r   = controls.rudder;      % rudder deflection (rad)

    % Unpack key parameters
    R      = params.rotor.R;        % rotor radius (m)
    A      = pi * R^2;              % rotor disc area (m^2)
    Omega  = params.rotor.Omega;    % rotor angular speed (rad/s)
    a      = params.rotor.a;        % lift-curve slope of blade airfoil (1/rad)
    sigma  = params.rotor.sigma;    % rotor solidity (blade area / disc area)
    Cd0    = params.fuselage.Cd0;   % fuselage drag coefficient (frontal area reference)
    A_fp   = params.fuselage.A_front; % fuselage frontal area (m^2)
    S_vt   = params.empennage.S;    % vertical tail (rudder) area (m^2)
    l_vt   = params.empennage.l;    % lever arm from CG to vertical tail (m)
    Cn_delta = params.empennage.Cn_delta; % rudder yaw moment coefficient (per rad deflection)

    % -- Main Rotor Thrust and Torque (Blade Element Theory, uniform inflow) --
    % Compute thrust coefficient CT using blade element theory (uniform inflow assumed):contentReference[oaicite:10]{index=10}.
    % Small-angle approximation: CT ≈ (sigma * a / 2) * (theta0 - theta_i), where theta_i is induced angle:contentReference[oaicite:11]{index=11}.
    % For simplicity, use effective collective (delta0) and assume uniform inflow factor k_i.
    theta_i = params.rotor.theta_i0;         % uniform inflow angle (rad) or initial guess
    CT = (sigma * a / 2) * (delta0 - theta_i);  % thrust coefficient (linearized)
    T_single = CT * rho * A * (Omega * R)^2;    % Thrust of one rotor (N):contentReference[oaicite:12]{index=12}

    % Compute torque coefficient CQ similarly (from profile + induced drag).
    % Approximate CQ = CT * ct_to_cq + profile contribution (simplified).
    % (ct_to_cq is a conversion factor from thrust coeff to torque coeff, depends on blade drag; assume given or computed)
    ct_to_cq = params.rotor.ct_to_cq;        % torque coefficient per unit CT (empirical or from theory)
    CQ = CT * ct_to_cq;                      % simple torque coefficient estimate
    Q_single = CQ * rho * A * (Omega * R)^2 * R;  % Torque of one rotor (N·m)

    % For coaxial rotors, assume two rotors share the load equally (symmetric):contentReference[oaicite:13]{index=13}.
    T_total = 2 * T_single;    % total thrust (upper + lower rotor)
    Q_upper = Q_single; 
    Q_lower = Q_single; 
    % Counter-rotation causes torques to oppose each other, canceling net yaw moment in hover:contentReference[oaicite:14]{index=14}:
    Q_total = Q_upper - Q_lower;   % net torque about yaw axis (should be ~0 if equal)

    % -- Rotor Flapping Dynamics (First Harmonic) --
    % Compute first-harmonic flapping angles β1c (longitudinal) and β1s (lateral).
    % Model as linear functions of cyclic inputs and advance ratio μ (μ = V/(Omega*R), assumed ≪ 1):contentReference[oaicite:15]{index=15}.
    mu = V / (Omega * R);  % advance ratio (non-dimensional forward speed)
    % Padfield (2007) suggests β1c and β1s vary approximately with pilot cyclic and μ in low-speed flight:contentReference[oaicite:16]{index=16}.
    Kmu   = params.rotor.Kmu;        % flapping sensitivity to forward speed (1)
    K_lon = params.rotor.K_lon;      % flapping gain for longitudinal cyclic (1)
    K_lat = params.rotor.K_lat;      % flapping gain for lateral cyclic (1)
    beta1c = Kmu * mu + K_lon * delta_lon;  % longitudinal flapping angle (rad)
    beta1s = K_lat * delta_lat;             % lateral flapping angle (rad)
    % (β1c tilts rotor tip-path plane fore-aft, β1s tilts it laterally:contentReference[oaicite:17]{index=17})

    % -- Rotor Force Components due to Tilted Thrust (Flapping) --
    % The rotor thrust vector tilts with the flapping angles, producing horizontal force components.
    % Small-angle assumption: sin(β) ~ β, cos(β) ~ 1.
    % Thrust in body Z-axis (vertical): mostly upward. We take Z positive down (so upward thrust is negative Z).
    Z_rotor = -T_total;  % thrust acts upward against gravity (negative Z by convention)

    % X-axis force from rotor thrust tilt (forward/backward).
    % β1c > 0 means rotor disk tilting backward (nose-up), causing thrust to have a rearward component.
    % X_rotor = -T_total * sin(β1c) ≈ -T_total * β1c (N, positive X forward) 
    X_rotor = -T_total * beta1c;  % forward component of thrust (acts opposite to direction of disk tilt)

    % Y-axis force from rotor thrust tilt (left/right).
    % β1s > 0 tilts disk to the right, yielding thrust component to the right.
    % Y_rotor = T_total * sin(β1s) ≈ T_total * β1s 
    Y_rotor = T_total * beta1s;   % lateral force (positive to right)

    % Note: If rotor hub is above CG, X_rotor or Y_rotor could create pitch/roll moments. 
    % Here we assume rotor forces act through the CG (no extra moment from hub offset).

    % -- Fuselage Aerodynamic Drag (Parasite Drag) --
    % Use flat-plate drag model: D = 0.5 * rho * Cd0 * A_front * V^2:contentReference[oaicite:18]{index=18}.
    D = 0.5 * rho * Cd0 * A_fp * V^2;    % parasitic drag (N)
    % Drag opposes forward motion (X-axis), so acts negative in body X when V is positive forward.
    X_drag = -D;
    Y_drag = 0;  % assume no side slip, so no lateral drag
    Z_drag = 0;  % assume negligible vertical drag component at low speeds

    % -- Empennage Rudder Yaw Moment --
    % The vertical tail/rudder produces a yaw moment N_rudder = 0.5 * rho * V^2 * S_vt * l_vt * Cn_delta * delta_r.
    % This is analogous to a flat-plate lift generating a side-force times lever arm for yaw:contentReference[oaicite:19]{index=19}.
    Q_dyn = 0.5 * rho * V^2;   % dynamic pressure
    N_rudder = Q_dyn * S_vt * l_vt * Cn_delta * delta_r;  % yaw moment from rudder (N·m)
    % (Positive delta_r yields yaw moment in the positive N direction by convention)

    % The rudder side force is F_y_tail = Q_dyn * S_vt * Cn_delta * delta_r (if needed for completeness).
    Y_rudder = 0;  % (Assume net side-force from rudder is small or cancels with tail fin; included in moment via lever arm)

    % -- Summation of Forces and Moments --
    % Sum forces from rotor thrust tilt, fuselage drag, and rudder (if any side-force considered).
    X_total = X_rotor + X_drag;      % forward force (trimmed to ~0 in steady flight)
    Y_total = Y_rotor + Y_drag + Y_rudder;
    Z_total = Z_rotor + Z_drag;      % vertical force (will balance weight in hover/level flight)

    % Sum moments: roll (L), pitch (M), yaw (N).
    % Roll moment L and pitch moment M are primarily from rotor flapping in advanced models.
    % Here we assume symmetric rotor system and no hub moment (articulated rotor), so L, M from rotors ≈ 0.
    L_rotor = 0;
    M_rotor = 0;
    % If considering rotor hub offset h (m) above CG: 
    %   L_rotor ≈ Y_rotor * h (roll moment due to lateral force), 
    %   M_rotor ≈ -X_rotor * h (pitch moment due to forward force).
    % For simplicity, we neglect these hub moment contributions.

    L_total = L_rotor;  % (no other significant roll moments in hover/low-speed)
    M_total = M_rotor;  % (assuming fuselage has negligible pitch moment in hover trim)
    N_total = N_rudder + Q_total;  % yaw moment from rudder plus any net rotor torque (Q_total ~ 0 if symmetric)

    % Package outputs
    F = [X_total; Y_total; Z_total];
    M = [L_total; M_total; N_total];
end
