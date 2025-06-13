function [T_u, T_l, Ct_j_u, Ct_j_l, Cp_j, Cp_j_l, Q_UR, Q_LR] = Stability_hover_getT(u, v, w, a_input, rho_input, mu_input, m_input, g_input, polar, theta_0_u, theta_0_l)

% Get all constants
constants = getConstants();  % or use setupConstants() directly if you prefer

% Access individual constants using dot notation
% GTOW    = constants.GTOW;    % kg
GTOW    = m_input;             % kg
Nb      = constants.Nb;
AR      = constants.AR;
Vtip    = constants.Vtip;      % m/s
R       = constants.R;         % m
r_0     = constants.r0;        % ND root cutout
A       = constants.A;         % m^2
Ae      = constants.Ae;        % m^2 (effective area)
cbar    = constants.cbar;      % mean chord
Cd0     = constants.Cd0;       % profile drag coefficient
sigma_e = constants.sigma_e;   % thrust-weighted equivalent solidity
rho     = constants.rho;       % kg/m^3
rho_cr  = constants.rho_cr;    % cruise density
Omega   = constants.Omega;

% load airfoil data
% addpath('Airfoil');
% polar = loadPolarData('xf-rc410-il-1000000.txt');

% Parameter Input
theta_tw_u = rad2deg(constants.theta_tw_u); % root to tip, degree
theta_tw_u = theta_tw_u/R*pi/180*R;
theta_tw_l = rad2deg(constants.theta_tw_l);
theta_tw_l = theta_tw_l/R*pi/180*R;

TR_u = 1;
TR_l = 1; 
d = 0;

mu_x =  u / Vtip;
mu_y =  v / Vtip;
mu   =  hypot(mu_x,mu_y);

Vc = -w;

% [~, ~, Ct_u_final, Ct_l_final] = BEMT_axial_optimisation_func(theta_tw_u,theta_tw_l, TR_u, TR_l, Vc, polar);

% Chord distribution from Taper Ratio
c0_u = cbar*(3*((TR_u-1)/4 + 1/3))^-1;
c0_l = cbar*(3*((TR_l-1)/4 + 1/3))^-1;

c_u = @(r) c0_u* ((TR_u - 1)*r + 1); 
c_l = @(r) c0_l* ((TR_l - 1)*r + 1);

% Solidity 
sigma_u = @(r) Nb*c_u(r)/(pi*R);
sigma_l = @(r) Nb*c_l(r)/(pi*R);

% Airfoil data
Cl_alpha0 = 6.2305; % lift curve slope
AoA_zl = -1*pi/180; % zero-lift AoA

% Non dimensional terms
lambda_c = Vc/Vtip;

Cw = GTOW*g_input/(rho_input*Ae*Vtip^2); % weight coeff

% Equal thrust sharing
% Ct_u_req = Ct_u_final; % thrust coeff required by upper rotor (Ct_u/Ct_l = 1.2461)

% Non dimensional induced velocity
lambda_hover = sqrt(Cw/2/2); % hover induced velocity
lambda_induced = -lambda_c/2 + lambda_hover*sqrt((0.5*lambda_c/lambda_hover)^2 + 1); % induced velocity at climb (Momentum Theory)
% lambda_total = lambda_induced+lambda_c; % lambda_total = lambda_induced + lambda_climb


% theta_0_u_initial = ((2*Ct_u_req*pi*R/Cl_alpha0/Nb) ...
                    % - theta_tw_u*c0_u/4 * (4/5*TR_u+1/5) ...
                    % + lambda_total*c0_u/2 * (2/3*TR_u+1/3)) * 3*(c0_u*(3*TR_u/4+1/4))^-1 ;

% Discretisation
N = 100;
tol = 1e-5;
r = linspace(r_0,1,N);
dr = r(2)-r(1);
% theta_0_u = theta_0_u_initial;
% lambda_induced_j = lambda_induced;
eps = 1;

% disc discretisation - for forward flight analysis
NN = 90;
azimuth = linspace(0,2*pi,NN);
dr = r(2)-r(1);
dpsi = azimuth(2)-azimuth(1);

% initialise arrays
lambda_j = zeros(1,length(r));
dCt_j = zeros(1,length(r));
sigma = zeros(1,length(r));
pitch = zeros(1,length(r));
inflow_angle = zeros(1,length(r));
prandlt_loss = zeros(1,length(r));
Mach = zeros(1,length(r));
dCp_j = zeros(1,length(r));
AoA = zeros(1,length(r));
AoA_l = zeros(1,length(r));

count = 1;

% while abs(eps) > tol
% % for i = 1:10

for i = 1:length(azimuth)
    psi = azimuth(i);

    % for each discretised blade element
    for j = 2:length(r)-1
       

        sigma_j = sigma_u(r(j)); % 
        pitch_j = theta_0_u + theta_tw_u*(r(j));
 
        % calculate blade sectional velocity
        U_T = Omega*r(j)*R + Vtip*( mu_x * sin(psi)  -  mu_y * cos(psi) );       
        % U_P = (lambda_c + lambda_induced_j)*Vtip;
        % U_t(i,j) = U_T;
        % U_p(i,j) = U_P;
        % U_R = mu_x*Vtip*cos(psi);
        U_P = (lambda_c + lambda_induced)*Vtip;
        U = sqrt(U_T^2 + U_P^2);
        Mach(j) = U/a_input;

        % Prandlt tip loss
        phi = atan(U_P/U_T);
        f = Nb/2* (1-r(j))/(r(j)*phi);
        F = 2/pi * acos(exp(-f));
        % 
        % % local angle of attack
        AoA(j) = pitch_j - phi - AoA_zl;
        
        % compressibility correction
        Cl_alpha = Cl_alpha0/sqrt(1-Mach(j)^2);

        % CL and CD from airfoil data (compressibility corrected)
        [Cl, Cd] = airfoilCoeffs(AoA(j)*180/pi, polar, Mach(j));

        % record
        sigma(j) = sigma_j;
        pitch(j) = pitch_j;
        % inflow_angle(j) = phi;
        % prandlt_loss(j) = F;

        % calculate inflow at each blade element
        lambda_j(j) = ((sigma_j*Cl_alpha/16/F - lambda_c/2)^2 ...
                    + sigma_j*Cl_alpha*pitch_j*r(j)/8/F)^(1/2) ...
                    - (sigma_j*Cl_alpha/16/F - lambda_c/2);
        % lambda_j3(j) = sigma_j*Cl_alpha/16 * (sqrt(1 + 32/(sigma_j*Cl_alpha) * pitch_j*r(j)) - 1);
    
        
        % incremental thrust
        % dCt_j(j) = sigma_j*Cl_alpha/2 * (pitch_j - phi) * r(j)^2 * dr; % without small angle approx
        dCt_j2(i,j) = sigma_j*Cl_alpha/2 * ((pitch_j) * r(j)^2 - lambda_j(j)*r(j)) * dr; % with small angle approx
        dCt_j3(i,j) = 4*F*lambda_j(j)*(lambda_j(j) - lambda_c)*r(j)*dr;
        dCl_j(i,j) = Cl_alpha*(pitch_j - phi - AoA_zl);

        % dCt_j4(j) = Nb*0.5*rho*U^2*c_u(r(j))*Cl*dr*R/(rho_input*Ae*Vtip^2);
        dT_j(i,j) = dCt_j2(i,j)*(rho_input*Ae*Vtip^2);

        dCp_j(i,j) = dCt_j3(i,j)*(lambda_j(j) - lambda_c);

        % incremental power
  
        % dCp_induced_j(j) = dCt_j4(j)*lambda_j(j);
        % dCp_profile_j(j) = Nb*(0.5*rho*U^2*c_u(r(j))*Cd*dr*R)*U_T/(rho*Ae*Vtip^3); % profile power

    end

    Ct_psi_u(i) = sum(dCt_j2(i ,:));
    Cp_psi_u(i) = sum(dCp_j(i ,:));

end

    % Sum sectional thrust to get overall Ct
    Ct_j_u = trapz(azimuth, Ct_psi_u) / (2*pi);
    Cp_j   = trapz(azimuth, Cp_psi_u) / (2*pi);
    T_u    = Ct_j_u * rho_input * Ae * Vtip^2;

    Q_UR = Cp_j * rho_input * Ae * Vtip^3 / Omega;
    
    % Cp_induced_j = sum(dCp_induced_j);
    % Cp_profile_j = sum(dCp_profile_j);

    % % Ct_j2 = sum(dCt_j2)
    % record_Ct(count) = Ct_j;
    % 
    % % Recalculate lambdas
    % lambda_hover_j = sqrt(Ct_j/2); % hover induced velocity
    % lambda_induced_j = -lambda_c/2 + lambda_hover_j*sqrt((0.5*lambda_c/lambda_hover_j)^2 + 1); % induced velocity at climb (Momentum Theory)
    % lambda_total_j = lambda_induced_j+lambda_c; % lambda_total = lambda_induced + lambda_climb\
    % 
    % theta_0_u_new = theta_0_u + ((2*(Ct_u_req - Ct_j)*pi*R/Cl_alpha0/Nb) ...
    %                 + (lambda_total - lambda_total_j)*c0_u/2 * (2/3*TR_u+1/3)) * 3*(c0_u*(3*TR_u/4+1/4))^-1;

%     % Calculate residual
%     eps(count) = theta_0_u_new - theta_0_u;
%     eps2(count) = (Ct_u_req - Ct_j);
% 
%     % update 
%     theta_0_u = theta_0_u_new;
%     % lambda_total = lambda_total_j;
% 
%     record(count) = theta_0_u;
%     count = count+1;
% end

% --------------- LOW ROTOR ------------------- %
% Equal thrust sharing
% Ct_l_req = Ct_l_final; % thrust coeff required by upper rotor (Ct_u/Ct_l = 1.2461)

% Non dimensional induced velocity
lambda_hover = sqrt(Cw/2/2); % hover induced velocity
lambda_induced = -lambda_c/2 + lambda_hover*sqrt((0.5*lambda_c/lambda_hover)^2 + 1); % induced velocity at climb (Momentum Theory)
% lambda_total = lambda_induced+lambda_c; % lambda_total = lambda_induced + lambda_climb

% initial estimate of collective pitch
% theta_0_l_initial = ((2*Ct_l_req*pi*R/Cl_alpha0/Nb) ...
                    % - theta_tw_l*c0_l/4 * (4/5*TR_l+1/5) ...
                    % + lambda_total*c0_l/2 * (2/3*TR_l+1/3)) * 3*(c0_u*(3*TR_l/4+1/4))^-1; 

% theta_0_l = theta_0_l_initial;
lambda_u = lambda_j;

A_Ac = 2;
r_c = 1/sqrt(2);
count = 1;
eps_l = 1;

% while abs(eps_l) > tol
% % for i = 1:10

for i = 1:length(azimuth)
    psi = azimuth(i);

    % for each discretised blade element
    for j = 2:length(r)-1
        


        sigma_j = sigma_l(r(j)); % 
        pitch_j = theta_0_l + theta_tw_l*(r(j));

        % calculate blade sectional velocity
        U_T = Omega*r(j)*R + Vtip*( mu_x * sin(psi)  -  mu_y * cos(psi) );

        % if r(j) <= r_c
        %     U_P = (lambda_c + lambda_induced+A_Ac*lambda_u(j))*Vtip;
        % elseif r(j) > r_c
        %     U_P = (lambda_c + lambda_induced)*Vtip;
        % end
        U_P = (lambda_c + lambda_induced)*Vtip;
        % U_t(i,j) = U_T;
        % U_p(i,j) = U_P;
        % U_R = mu_x*Vtip*cos(psi);
        U = sqrt(U_T^2 + U_P^2);
        Mach(j) = U/a_input;

        % compressibility correction
        Cl_alpha = Cl_alpha0/sqrt(1-Mach(j)^2);

        % Prandlt tip loss
        % phi = atan((lambda_c + lambda_induced_j)/r(j));
        phi = atan(U_P/U_T);
        f = Nb/2* (1-r(j))/(r(j)*phi);
        F = 2/pi * acos(exp(-f));

        % local angle of attack
        AoA_l(j) = pitch_j - phi - AoA_zl;

        % CL and CD from airfoil data (compressibility corrected)
        [Cl, Cd] = airfoilCoeffs(AoA_l(j)*180/pi, polar, Mach(j));

        % record
        sigma(j) = sigma_j;
        pitch_l(j) = pitch_j;
        inflow_angle_l(j) = phi;
        prandlt_loss(j) = F;

        if r(j) <= r_c
            % calculate inflow at each blade element
            lambda_j_l(j) = ((sigma_j*Cl_alpha/16/F - (lambda_c + A_Ac*lambda_u(j))/2)^2 ...
                        + sigma_j*Cl_alpha*pitch_j*r(j)/8/F)^(1/2) ...
                        - (sigma_j*Cl_alpha/16/F - (lambda_c + A_Ac*lambda_u(j))/2); 
        elseif r(j) > r_c
            lambda_j_l(j) = ((sigma_j*Cl_alpha/16/F - (lambda_c)/2)^2 ...
                            + sigma_j*Cl_alpha*pitch_j*r(j)/8/F)^(1/2) ...
                            - (sigma_j*Cl_alpha/16/F - (lambda_c)/2); 
        end


        % incremental thrust
        % dCt_j(j) = sigma_j*Cl_alpha/2 * (pitch_j - phi) * r(j)^2 * dr; % without small angle approx
        dCt_j_l(i,j) = sigma_j*Cl_alpha/2 * ((pitch_j) * r(j)^2 - lambda_j_l(j)*r(j)) * dr; % with small angle approx
        % dCt_j2(j) = 4*F*lambda_j(j)*(lambda_j(j) - lambda_c)*r(j)*dr;
        dCl_j_l(i,j) = Cl_alpha*(pitch_j - phi - AoA_zl);

        % dCt_j2_l(j) = Nb*0.5*rho*U^2*c_l(r(j))*Cl*dr*R/(rho*Ae*Vtip^2);
        dT_j_l(i,j) = dCt_j_l(i,j)*(rho_input*Ae*Vtip^2);

        % incremental power
        dCp_j_l(i,j) = dCt_j_l(i,j)*(lambda_j_l(j) - lambda_c);

        % dCp_induced_j_l(j) = dCt_j2_l(j)*lambda_j_l(j);
        % dCp_profile_j_l(j) = Nb*(0.5*rho*U^2*c_l(r(j))*Cd*dr*R)*U_T/(rho*Ae*Vtip^3); % profile power

    end

    Ct_psi_l(i) = sum(dCt_j_l(i ,:));
    Cp_psi_l(i) = sum(dCp_j_l(i ,:));

end

    % Sum sectional thrust to get overall Ct
    Ct_j_l = trapz(azimuth, Ct_psi_l) / (2*pi);
    Cp_j_l = trapz(azimuth, Cp_psi_l) / (2*pi);
    T_l    = Ct_j_l * rho_input * Ae * Vtip^2;

    % Sum sectional thrust to get overall Ct
    % Ct_j_l = sum(dCt_j_l);
    % Cp_j_l = sum(dCp_j_l);

    % T_l = sum(dT_j_l);

    Q_LR = Cp_j_l * rho_input * Ae * Vtip^3 / Omega;

    % Cp_induced_j_l = sum(dCp_induced_j_l);
    % Cp_profile_j_l = sum(dCp_profile_j_l);

    % % Ct_j2 = sum(dCt_j2)
    % record_Ct(count) = Ct_j_l;
    % 
    % % Recalculate lambdas
    % lambda_hover_j = sqrt(Ct_j_l/2); % hover induced velocity
    % lambda_induced_j = -lambda_c/2 + lambda_hover_j*sqrt((0.5*lambda_c/lambda_hover_j)^2 + 1); % induced velocity at climb (Momentum Theory)
    % lambda_total_j = lambda_induced_j+lambda_c; % lambda_total = lambda_induced + lambda_climb
    % 
%     theta_0_l_new = theta_0_l + ((2*(Ct_l_req - Ct_j_l)*pi*R/Cl_alpha0/Nb) ...
%                     + (lambda_total - lambda_total_j)*c0_l/2 * (2/3*TR_l+1/3)) * 3*(c0_l*(3*TR_l/4+1/4))^-1;
% 
%     % Calculate residual
%     eps_l(count) = theta_0_l_new - theta_0_l;
%     eps2_l(count) = (Ct_l_req - Ct_j_l);
% 
%     % update 
%     theta_0_l = theta_0_l_new;
%     % lambda_total = lambda_total_j;
% 
%     record(count) = theta_0_l;
%     count = count+1;
% end

% theta_0_l = theta_0_l*180/pi;


% FM calcuation
% k_int = 1.219; % torque trimmed coaxial interference factor
% k = 1.1; % induced loss factor
% Cp_ideal = sum(Cp_induced_j + Cp_induced_j_l);
% Cp_profile = sum(Cp_profile_j + Cp_profile_j_l);
% FM = Cp_ideal/(k_int*k*Cp_ideal + Cp_profile);

end