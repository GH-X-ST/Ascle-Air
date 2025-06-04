function [lambda_i, lambda] = LinearInflow_ff(mu_x, mu_z, alpha_s, C_T, r, psi)

% inflow calculation using iteration in forward flight
% Define the function whose root we want to find
f = @(lambda) lambda - (mu_x * tan(alpha_s) + C_T / (2 * sqrt(mu_x^2 + lambda.^2)));
lambda0 = sqrt(C_T/2); % Initial guess from hover value
lambda_ff = fzero(f, lambda0); % Use fzero to find total inflow 
lambda_0 = C_T/(2*sqrt(mu_x^2 + lambda_ff^2));

% overall mu
mu = sqrt(mu_x^2 + mu_z^2);

% wake angle
chi = atan(mu_x/(mu_x*tan(alpha_s)+lambda_0));

% dree's constants
k_x = 4/3 * (1 - cos(chi) - 1.8*mu^2) / sin(chi);
k_y = -2*mu;

% linear inflow model
lambda_i = mu_x*tan(alpha_s) + lambda_0*(1+k_x*r*cos(psi) + k_y*r*sin(psi));

% total inflow = forward flight inflow + climb inflow + induced inflow 
lambda = mu_x*tan(alpha_s) + lambda_i + mu_z;

end