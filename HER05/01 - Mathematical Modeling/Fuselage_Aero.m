function [L_F, D_F, Y_F, Mx_F, My_F] = Fuselage_Aero(u, v, w, rho)
% Introduction:
% This function is developed from the Jim Liu's CFD result for 140 kts 
% forward flight at 5000 ft, based on very-rough “velocity-squared” scaling
% and simple geometry calculation
%
% Author:
%   Hanchen Li (hl3422@ic.ac.uk)
%
% Inputs:
%   u         - Forward velocity (m/s)
%   v         - Lateral velocity (m/s)
%   w         - Vertical velocity (m/s)
%   rho       - Density (kg/m^3)
%
% Outputs:
%   L_F       - Fuselage lift (N)
%   D_F       - Fuselage drag (N)
%   Y_F       - Fuselage side force (N)
%   Mx_F      - Fuselage rolling moment (Nm)
%   My_F      - Fuselage pitching moment (Nm)

[~, ~, ~, rho_test, ~, ~] = atmosisa(1524); % International Standard Atmosphere model at 5000 ft

S = 4.095;    % Fuselage reference area (m^2)

% forward flight
CL_fwd  = -0.00716;
CD_fwd  =  0.15;

% backward flight
CL_back = 0.059;
CD_back = -0.208;

% velocity direction
V = sqrt(u^2 + v^2 + w^2);
v_hat = [u, v, w] / V;

if V < 1e-6
    L_F  = 0;
    D_F  = 0;
    Y_F  = 0;
    Mx_F = 0;
    My_F = 0;
    return
end

% Cosine of angle between velocity and +X-body
cos_theta = u / V;

% Interpolate CL and CD
CL = 0.5 * (1 + cos_theta) * CL_fwd + 0.5 * (1 - cos_theta) * CL_back;
CD = 0.5 * (1 + cos_theta) * CD_fwd + 0.5 * (1 - cos_theta) * CD_back;

% Drag acts opposite to velocity
e_D = -v_hat;

% Lift is perpendicular to velocity in X–Z plane (assume no side-slip)
if abs(u) < 1e-3 && abs(w) < 1e-3

    e_L = [0, 0, 0];

else

    e_L = [-w, 0, u];
    e_L = e_L / norm(e_L);

end

% Dynamic pressure
q = 0.5 * rho * V^2;

% Force vectors
F_D = q * S * CD * e_D;
F_L = q * S * CL * e_L;

% Total aerodynamic force
F_total = F_D + F_L;

% Extract components in body axes
D_F = F_total(1);   % Drag (positive aft)
Y_F = F_total(2);   % Side force (positive right)
L_F = F_total(3);   % Lift (positive up)

% Moments using velocity-squared scaling
Mx_F = 104  * (rho * V^2) / (rho_test * (74^2));
My_F = 1700 * (rho * V^2) / (rho_test * (74^2));

end