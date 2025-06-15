function [lambda_i, lambda] = LinearInflow_ff(mu_x, mu_z, alpha_s, C_T, r, psi)

% Linear skewed-wake inflow (Drees + Pitt/Peters first harmonic)
% All angles in rad, r nondimensional radius, psi azimuth (rad)

% 1. Hub induced inflow (exact momentum theory)
lambda_0 = sqrt( 0.5*C_T + 0.25*mu_x^2 ) - 0.5*mu_x;

% 2. Wake skew angle
chi = atan2( mu_x , lambda_0 + mu_x*tan(alpha_s) );   % always finite

% 3. Drees/Pitt constants
k_x = 4/3 * ( 1 - 1.8*mu_x^2 ) / max( sin(chi) , 1e-6 );  % protect /0
k_y = -2 * mu_x;

% 4. Induced inflow distribution
lambda_i = lambda_0 .* ( 1 + k_x*r.*cos(psi) + k_y*r.*sin(psi) );

% 5. Total axial inflow
lambda   = lambda_i + mu_z + mu_x*tan(alpha_s);
end