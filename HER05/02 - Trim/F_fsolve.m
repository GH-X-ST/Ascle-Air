function [theta_LR, theta_UR, theta_1s, theta_1c, phi, theta, delta_E, delta_R] = F_fsolve(u)
% Author:
%   Hanchen Li (hl3422@ic.ac.uk)
%
% Inputs:
%   u         - Forward flight velocity (m/s)
%
% Outputs:
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_1s  - Longitudinal cyclic (rad)
%   theta_1c  - Lateral cyclic (rad)
%   phi       - Euler roll angle (rad)
%   theta     - Euler pitch angle (rad)
%   delta_E   - Elevator deflection (rad)
%   delta_R   - Rudder deflection (rad)
%

%% 0 Basic Parameters

% 0.1 Initial guess
x0 = [0.2788; 0.2323; 0; 0; 0; 0; 0; 0];

% fsolve options
opts = optimoptions('fsolve', ...
    'Display','iter', ...
    'FunctionTolerance',1e-8, ...
    'StepTolerance',1e-8, ...
    'MaxIterations',200);

%% 2. Residual function

function F = trimResiduals(x)
    
    params = [1524, 3491, 9.80665, 0.9, 0.682, 0.682, 5.5, 5.5];

    % Unpack
    theta_LR  = x(1);
    theta_UR  = x(2);
    theta_1s  = x(3);
    theta_1c  = x(4);
    phi       = x(5);
    theta     = x(6);
    delta_E   = x(7);
    delta_R   = x(8);

    % Build ctrl & state
    ctrl  = [theta_LR, theta_UR, theta_1s, theta_1c, delta_E, delta_R];
    state = [u, 0, 0, 0, 0, 0, 0, 0, 0];
    state(7:8) = [phi, theta];  % insert Euler angles into state

    % Call your force-moment model
    [X, Y, Z, L, M, N, ~, ~, ~, ~] = FM_N(state, ctrl, params);

    % Residuals: forces & moments must be zero in steady trim
    F = [X; Y; Z; L; M; N];
end

%% 3. Solve with fsolve
[x_trim, ~, exitflag, output] = fsolve(@trimResiduals, x0, opts);

if exitflag <= 0

    warning('Trim did not converge (flag = %d)', exitflag);

else
    
    fprintf('Converged in %d iterations.\n', output.iterations);
end

%% 4. Inspect solution
theta_LR = x_trim(1);
theta_UR = x_trim(2);
theta_1s  = x_trim(3);
theta_1c  = x_trim(4);
phi       = x_trim(5);
theta     = x_trim(6);
delta_E   = x_trim(7);
delta_R   = x_trim(8);

end