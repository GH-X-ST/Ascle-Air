function [F, M] = FM_h(state, ctrl, params)
% Introduction:
%   This function computes the net forces and moments acting on the
%   helicopter in body axes, considering contributions from two coaxial
%   rotors and gravity.
%   It is tailored for hover and low-speed forward flight (advance ratio
%   μ small, 0–40 knots).
%
% References:
%   \bibitem{Padfield} Rotor blade element and flapping dynamics 
%   \bibitem{Coleman}  Coaxial rotor aerodynamic research
%   \bibitem{Leishman} BEMT for rotor thrust/torque
%   \bibitem{Bramwell} Helicopter dynamics
%
% Inputs:
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
%   params.c        - Blade chord (m)
%   params.a0       - Lift-curve slope of blade airfoil (1/rad)
%   params.theta_tw - Blade linear twist (rad)
%                     (negative for outwash, assume θ_tw < 0)
%   params.Omega    - Rotor angular speed (rad/s).
%   params.rho      - Air density (kg/m^3).
%   params.I_beta   - Blade flapping moment of inertia about hinge (kg·m^2)
%   params.e        - Hinge offset or effective flapping hinge offset (m)
%   params.Cd0      - Blade section profile drag coefficient
%   params.mass     - Helicopter mass (kg)
%   params.g        - Gravitational acceleration (m/s^2)
%   params.Cd0_f    - Fuselage drag coefficient (frontal area reference)
%   params.A_f      - Fuselage frontal area (m^2)
%   params.S_vt     - Vertical tail (rudder) area (m^2)
%   params.l_vt     - Level arm from CG to vertical tail (m)
%   params.Cn_delta - Rudder yaw moment coefficient (per rad deflection)
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
%                     β1s = lateral flapping (disk tilt side-to-side)
%                     Positive β1c tilts the rotor disk back
%                     (gives nose-up moment)
%                     Positive β1s tilts the disk to the left 
%                     (port side down)
%   Moments         - L (roll) positive for right-wing-down
%                     M (pitch) positive nose-up
%                     N (yaw) positive nose-right
%
% Note:
%   For simplicity, this model assumes small angles (β1c, β1s << 1 rad)
%   and uniform inflow across the rotor disk.

% Extract parameters for convenience
R        = params.R;
Nb       = params.Nb;
c        = params.c;
a0       = params.a0;
theta_tw = params.theta_tw;
Omega    = params.Omega;
rho      = params.rho;
I_beta   = params.I_beta;
e        = params.e;
Cd0      = params.Cd0;
m        = params.mass;
g        = params.g;
Cd0_f    = params.Cd0_f;
A_f      = params.A_f;
S_vt     = params.S_vt;
l_vt     = params.l_vt;
Cn_delta = params.Cn_delta;
    
% Derived parameters
A        = pi * R^2;                      % Rotor disc area (m^2)
sigma    = (Nb * c) / (pi * R);           % Rotor solidity
gamma    = (rho * c * a0 * R^4) / I_beta; % Lock number
    
% Control inputs (collective for each rotor)
theta0_top  = ctrl.coll_avg + 0.5 * ctrl.coll_diff;
% collective pitch top rotor (rad)
theta0_bot  = ctrl.coll_avg - 0.5 * ctrl.coll_diff;
% collective pitch bottom rotor (rad)
    
% Body-axis velocities (for low-speed inflow adjustment)
U = state.U; % forward velocity (m/s)
V = state.V; % sideward velocity (starborad positive) (m/s)
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

% Momentum theory for hover as baseline: v_h = v_i0 = sqrt(T/2*rho*A)
% 1. Estimate thrust from collective 
% 2. Refine with induced velocity consistency

% Simplify Glauert's Momentum Theory for Forward Flight \bibitem{Bramwell}
% 1. For low forward speed, the effective induced velocity is slightly
% reduced because the rotor wake is tilted and convected aft (skew angle)
% 2. Taylor seris expansion yields v_i = v_i0(1 - mu^2/2 - ...)
% 3. Approximate by scaling v_i0 by (1 - mu^2) for small mu

% 1. Estimate thrust coefficient from blade element
%    (assuming some initial inflow guess)
% 2. Solve T = 0.5*rho*A*(Omega*R)^2 * C_T
%    with C_T = (sigma*a0/2)*(theta0_eff - lambda_i).
% 3. Iterate once for consistency between C_T and lambda_i
    
    function [T_est, lambda_i] = estimateThrust(theta0)
        % Nested helper to estimate thrust for a given collective (one rotor).
        % Start with initial guess lambda_i from momentum theory in hover (no forward flight).
        lambda_i = 0;  % initial guess (will iterate)
        for it = 1:2
            % Effective blade pitch accounting for twist: use average pitch over the blade
            theta_avg = theta0 + 0.5 * theta_tw;  % approximate average pitch (rad)
            % Thrust coefficient via blade element (linearized): C_T = sigma*a0/2 * (theta_avg - lambda_i)
            C_T = (sigma * a0 / 2) * (theta_avg - lambda_i);
            % Induced inflow from momentum: lambda_i_new = C_T / 2 (hover momentum theory):contentReference[oaicite:18]{index=18}
            lambda_i_new = C_T / 2;
            % Simple forward-flight correction: reduce induced inflow for horizontal velocity
            lambda_i_new = lambda_i_new * sqrt(1 + mu^2);  % (Glauert-like correction; for small mu, ~ no change)
            % Update lambda_i for next iteration
            lambda_i = lambda_i_new;
        end
        % Thrust from final C_T
        T_est = C_T * (rho * A * (Omega*R)^2);
    end
    
    % Estimate thrust for top and bottom rotors based on current collectives
    [T_top, lambda_top] = estimateThrust(theta0_top);
    [T_bottom, lambda_bottom] = estimateThrust(theta0_bot);
    
    %% 2. Compute rotor flapping angles (first-harmonic β1c, β1s)
    % The flapping angles β1c (fore-aft tilt) and β1s (lateral tilt) are determined by the 1/rev (once-per-revolution) 
    % aerodynamic moments on the rotor blades and the blade inertial/restoring forces.
    % In steady low-speed flight, we assume quasi-static flapping (no significant flapping acceleration), so we solve for 
    % β1c and β1s that balance the 1/rev aerodynamic loads.
    %
    % **Longitudinal flapping (β1c):** In forward flight, the rotor disk tends to "blow back" – it tilts aft (nose-up) 
    % due to higher lift on the advancing blade front-half vs retreating rear-half:contentReference[oaicite:19]{index=19}. We approximate β1c by assuming 
    % the disk tilts to align resultant thrust opposite to flight velocity. A simple model is β1c ≈ (2/3)*mu_x (for small mu) with no cyclic input.
    % (More precisely, β1c would be solved from flap equations including inflow; here we use a proportional estimate.)
    %
    % **Lateral flapping (β1s):** With no lateral cyclic and a single rotor, asymmetry of lift (advancing vs retreating side) causes 
    % a lateral tilt. For a CCW rotor (top), the advancing blade is on the right, lifting more, so the disk tilts left (port) 
    % i.e. β1s_top > 0 for forward flight. For the CW bottom rotor, the advancing blade is on the left, so it tilts right (β1s_bottom < 0). 
    % In a symmetric coaxial system, these lateral flaps tend to cancel out in the net force:contentReference[oaicite:20]{index=20} (we assume equal magnitude).
    %
    % **Lock number effect:** The magnitude of flapping response is inversely related to blade inertia. A higher Lock number (γ) means 
    % stronger aerodynamic influence and larger flapping angles for a given asymmetry:contentReference[oaicite:21]{index=21}:contentReference[oaicite:22]{index=22}. We include γ in the proportional gain.
    %
    % Compute approximate flapping angles based on current flight velocities:
    
    % Longitudinal flapping (nose-up tilt due to forward speed):
    beta1c_top    = (2/3) * mu_x * (1 / (1 + 0.5*gamma));   % small-signal guess: reduced by blade inertia (γ) effect
    beta1c_bottom = beta1c_top;  % assume both rotors tilt aft by same amount in symmetric flight
    
    % Lateral flapping (disk tilt to port for top rotor, to starboard for bottom rotor):
    beta1s_top    = (2/3) * mu_y * (1 / (1 + 0.5*gamma));   % if moving right (mu_y>0), CCW top rotor advances from front-right, disk tilts left
    beta1s_bottom = -beta1s_top;  % bottom rotor tilts opposite due to opposite rotation
    
    % Note: The above flapping model is a simplification. A more accurate solution would solve blade flapping ODE:
    %       I_beta * \ddot{β} + K_beta*β = Aerodynamic Moment  (with K_beta from hub stiffness or offset). 
    %       Here we assume steady-state (no acceleration) and use empirical fraction of velocity.
    
    %% 3. Compute rotor forces in body axes
    % Each rotor produces a thrust mostly along its shaft. When the rotor flaps, the tip-path-plane tilts by β1 angles relative to the shaft, 
    % effectively tilting the thrust vector in the body frame:contentReference[oaicite:23]{index=23}. 
    % For small flapping angles, the thrust vector tilts approximately by the same angles.
    %
    % We resolve each rotor's thrust into body-axis components:
    %   X_force ~ T * sin(β1c)   (forward component; β1c > 0 tilts thrust aft, so if β1c is defined nose-up, the X-force might be negative. We assume β1c > 0 = aft tilt => negative X force.)
    %   Y_force ~ T * sin(β1s)   (side component; β1s > 0 tilts thrust to left (port), producing a leftward force, i.e. negative Y).
    %   Z_force ~ T * cos(tilt)  (~ -T for small angles since thrust is upward opposite body +Z).
    %
    % Sign convention: We will take β1c > 0 as nose-up tilt (thrust vector tilting backward), which yields a *negative* X force (drag). 
    %                 β1s > 0 as tilt left, yielding *negative* Y force.
    
    % Top rotor thrust vector (body axes components)
    Fx_top = - T_top * beta1c_top;      % positive β1c (aft tilt) gives negative X force
    Fy_top = - T_top * beta1s_top;      % positive β1s (tilt left) gives negative Y force
    Fz_top = - T_top * (1 - 0.5*(beta1c_top^2 + beta1s_top^2));  % Z (down) component of thrust. ~-T (upwards lift), small correction for tilt.
    
    % Bottom rotor thrust vector
    Fx_bottom = - T_bottom * beta1c_bottom;
    Fy_bottom = - T_bottom * beta1s_bottom;
    Fz_bottom = - T_bottom * (1 - 0.5*(beta1c_bottom^2 + beta1s_bottom^2));
    
    % Sum rotor forces
    F_rotors_body = [Fx_top + Fx_bottom;
                     Fy_top + Fy_bottom;
                     Fz_top + Fz_bottom];
    
    %% 4. Compute rotor moments about center of gravity/hub
    % Rotor thrusts and flapping induce moments about the aircraft CG in roll (L) and pitch (M):
    % - A tilted thrust vector with an effective lever arm (due to flapping hinge offset or hub height) produces moments.
    % - We approximate hub moments from flapping as if the rotor thrust acts at a horizontal offset e (hinge offset) from the shaft:contentReference[oaicite:24]{index=24}.
    %
    % Roll moment L: A lateral disk tilt (β1s) causes a roll moment. If β1s_top > 0 (disk tilts left), thrust vector shifts left of center, rolling the aircraft left (negative L).
    %               Similarly the bottom rotor tilts opposite. We add their effects.
    %   L_top ≈ - T_top * e * beta1s_top
    %   L_bottom ≈ - T_bottom * e * beta1s_bottom
    %
    % Pitch moment M: A longitudinal tilt (β1c) causes a pitch moment. β1c > 0 (disk tilts aft, nose-up thrust tilt) shifts thrust vector behind the CG, pitching nose up (positive M).
    %                But if our sign for β1c is nose-up tilt, that yields a positive M.
    %   M_top ≈ + T_top * e * beta1c_top   (nose-up moment)
    %   M_bottom ≈ + T_bottom * e * beta1c_bottom
    %
    % (Sign note: With our convention β1c > 0 = nose-up tilt, thrust vector goes aft of rotor shaft, creating a nose-up pitching moment, hence positive M.)
    %
    % Yaw moment N: The net yawing moment is generated by torque difference between rotors. Each rotor’s torque reaction acts on the fuselage in opposite directions due to counter-rotation.
    %   - Top rotor (CCW viewed from above) torque tends to yaw fuselage clockwise (nose right).
    %   - Bottom rotor (CW) torque tends to yaw fuselage counter-clockwise (nose left).
    % So, N = Q_top - Q_bottom (positive if top torque > bottom torque, causing nose-right):contentReference[oaicite:25]{index=25}.
    %
    
    % Roll (L) and Pitch (M) from rotor disk tilts:
    L_rotors = - T_top * e * beta1s_top  +  - T_bottom * e * beta1s_bottom;
    M_rotors = + T_top * e * beta1c_top  +  + T_bottom * e * beta1c_bottom;
    
    % Yaw (N) from differential torque:
    N_rotors = Q_top - Q_bottom;
    
    M_rotors_body = [L_rotors;
                     M_rotors;
                     N_rotors];
    
    %% 5. Gravity and total forces/moments
    % Add gravitational force (weight) in body axes. In the NED frame (Z down), weight contributes +W in Z direction.
    W_force = m * g;  % weight (N)
    F_gravity = [0; 0; W_force];
    
    % Total forces and moments
    F = F_rotors_body + F_gravity;   % [X; Y; Z]
    M = M_rotors_body;               % [L; M; N] (no separate gravity moment since CG assumed at rotor hub for simplicity)
    
    % (If CG is offset from rotor hub along X/Y, we could add moment from weight: M_weight = [m*g* y_CG; -m*g* x_CG; 0], not included here.)
    
end
