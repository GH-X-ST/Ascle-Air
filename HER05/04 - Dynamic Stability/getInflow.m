function [lambda_tot,lambda_0,lambda_1c,lambda_1s] = getInflow(u,v,w,CT,geometry)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Extract Parameters
V_tip = geometry.Vtip;

% Advance ratio
mu = sqrt(u^2 + v^2) / V_tip;

% Total inflow coefficient
lambdafun = @(lambda) w/V_tip + .5*CT/sqrt(lambda^2 + mu^2) - lambda;
lambdainit = sqrt(CT/2);
options = optimoptions('fsolve','Display','off');
lambda_tot = real(fsolve(lambdafun,lambdainit,options));

% Decomposed Inflow coefficient
lambdai = .5*CT/sqrt(lambda_tot^2 + mu^2);

% Wake Angle
chi = atan2(mu,lambda_tot);

% Dree's constant;
Kc = 4/3 * (1 - 1.8*mu^2) * tan(chi/2);
Ks = -2*mu;

% Rotate to flight path direction
rotate = atan2(u,v);
lambda_0 = lambdai;
lambda_1c = lambdai*(Kc*cos(rotate) - Ks*sin(rotate));
lambda_1s = lambdai*(Ks*cos(rotate) + Kc*sin(rotate));

end