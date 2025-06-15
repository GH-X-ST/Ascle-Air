function [r, Ct_j, Cp_j_l, dCt_j, dCp_j_l, lambda_induced_j, lambda_j, dCl_j, AoA, theta_0_l, Cp_induced_j, Cp_profile_j] = BEMT_axial_LowerRotor(theta_tw_l, TR_l, lambda_u, d, Vc, Ct_l_req, polar)

% Get all constants
constants = getConstants();  % or use setupConstants() directly if you prefer

% Access individual constants using dot notation
GTOW    = constants.GTOW;      % kg
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

% % load airfoil data
% addpath('Airfoil');
% polar = loadPolarData('xf-rc410-il-1000000.txt');

% Parameter Input
theta_tw_l = theta_tw_l/R*pi/180*R;

% BladeGeometryInput = [theta_tw_u,theta_tw_l, TR_u, TR_l, d];

% Chord distribution from Taper Ratio
c0_l = cbar*(3*((TR_l-1)/4 + 1/3))^-1;
c_l = @(r) c0_l* ((TR_l - 1)*r + 1); 

% Solidity 
sigma_l = @(r) Nb*c_l(r)/(pi*R);

% Airfoil data
Cl_alpha0 = 6.2305; % lift curve slope
AoA_zl = -1*pi/180; % zero-lift AoA

% Non dimensional terms
lambda_c = Vc/Vtip;

% Non dimensional induced velocity
lambda_hover = sqrt(Ct_l_req/2); % hover induced velocity
lambda_induced = -lambda_c/2 + lambda_hover*sqrt((0.5*lambda_c/lambda_hover)^2 + 1); % induced velocity at climb (Momentum Theory)
lambda_total = lambda_induced+lambda_c; % lambda_total = lambda_induced + lambda_climb

% initial estimate of collective pitch
theta_0_l_initial = ((2*Ct_l_req*pi*R/Cl_alpha0/Nb) ...
                    - theta_tw_l*c0_l/4 * (4/5*TR_l+1/5) ...
                    + lambda_total*c0_l/2 * (2/3*TR_l+1/3)) * 3*(c0_l*(3*TR_l/4+1/4))^-1; 

% Calculate the rotor spacing and area contraction ratio
d; % how to use?
A_Ac = 2;
r_c = 1/sqrt(2);


% Discretisation
N = 100;
tol = 1e-5;
r = linspace(r_0,1,N);
dr = r(2)-r(1);
theta_0_l = theta_0_l_initial;
lambda_induced_j = lambda_induced;
eps = 1;

% initialise arrays
lambda_j = zeros(1,length(r));
dCt_j = zeros(1,length(r));
dCl_j = zeros(1,length(r));
sigma = zeros(1,length(r));
pitch = zeros(1,length(r)); 
AoA = zeros(1,length(r));
inflow_angle = zeros(1,length(r));
prandlt_loss = zeros(1,length(r));
Mach= zeros(1,length(r));
dCp_j_l= zeros(1,length(r));
dCt_j2 = zeros(1,length(r));
dCp_induced_j = zeros(1,length(r));
dCp_profile_j = zeros(1,length(r));

count = 1;

while abs(eps) > tol
% for i = 1:1

    % for each discretised blade element
    for j = 2:length(r)-1
        
        sigma_j = sigma_l(r(j)); % 
        pitch_j = theta_0_l + theta_tw_l*(r(j));

        % calculate blade sectional velocity
        U_T = Omega*r(j)*R;

        if r(j) <= r_c
            U_P = (lambda_c + lambda_induced+A_Ac*lambda_u(j))*Vtip;
        elseif r(j) > r_c
            U_P = (lambda_c + lambda_induced)*Vtip;
        end
        % U_t(i,j) = U_T;
        % U_p(i,j) = U_P;
        % U_R = mu_x*Vtip*cos(psi);
        U = sqrt(U_T^2 + U_P^2);
        Mach(j) = U/334.3;
        % Reynolds(j) = rho_cr*sqrt(U_T^2 + U_P^2)*c_u(r(j))/1.628e-5;

        % Prandlt tip loss
        phi = atan(U_P/U_T);
        f = Nb/2* (1-r(j))/(r(j)*phi);
        F = 2/pi * acos(exp(-f));

        % local angle of attack
        AoA(j) = pitch_j - phi - AoA_zl;

        % record
        sigma(j) = sigma_j;
        pitch(j) = pitch_j;
        inflow_angle(j) = phi;
        prandlt_loss(j) = F;

        % compressibility correction
        Cl_alpha = Cl_alpha0/sqrt(1-Mach(j)^2);

        % CL and CD from airfoil data (compressibility corrected)
        [Cl, Cd] = airfoilCoeffs(AoA(j)*180/pi, polar, Mach(j));

        if r(j) <= r_c
            % calculate inflow at each blade element
            lambda_j(j) = ((sigma_j*Cl_alpha/16/F - (lambda_c + A_Ac*lambda_u(j))/2)^2 ...
                        + sigma_j*Cl_alpha*pitch_j*r(j)/8/F)^(1/2) ...
                        - (sigma_j*Cl_alpha/16/F - (lambda_c + A_Ac*lambda_u(j))/2); 
        elseif r(j) > r_c
            lambda_j(j) = ((sigma_j*Cl_alpha/16/F - (lambda_c)/2)^2 ...
                            + sigma_j*Cl_alpha*pitch_j*r(j)/8/F)^(1/2) ...
                            - (sigma_j*Cl_alpha/16/F - (lambda_c)/2); 
        end


        % incremental thrust
        % dCt_j(j) = sigma_j*Cl_alpha/2 * (pitch_j - phi) * r(j)^2 * dr; % without small angle approx
        dCt_j(j) = sigma_j*Cl_alpha/2 * ((pitch_j) * r(j)^2 - lambda_j(j)*r(j)) * dr; % with small angle approx
        % dCt_j2(j) = 4*F*lambda_j(j)*(lambda_j(j) - lambda_c)*r(j)*dr;
        dCl_j(j) = Cl_alpha*(pitch_j - phi - AoA_zl);
        dCt_j2(j) = Nb*0.5*rho*U^2*c_l(r(j))*Cl*dr*R/(rho*Ae*Vtip^2);

        % incremental power
        dCp_j_l(j) = dCt_j(j)*lambda_j(j);

        dCp_induced_j(j) = dCt_j2(j)*lambda_j(j);
        dCp_profile_j(j) = Nb*(0.5*rho*U^2*c_l(r(j))*Cd*dr*R)*U_T/(rho*Ae*Vtip^3); % profile power

    end

    % Sum sectional thrust to get overall Ct
    Ct_j = sum(dCt_j);
    Cp_j_l = sum(dCp_j_l);

    Cp_induced_j = sum(dCp_induced_j);
    Cp_profile_j = sum(dCp_profile_j);

    % Ct_j2 = sum(dCt_j2)
    record_Ct(count) = Ct_j;

    % Recalculate lambdas
    lambda_hover_j = sqrt(Ct_j/2); % hover induced velocity
    lambda_induced_j = -lambda_c/2 + lambda_hover_j*sqrt((0.5*lambda_c/lambda_hover_j)^2 + 1); % induced velocity at climb (Momentum Theory)
    lambda_total_j = lambda_induced_j+lambda_c; % lambda_total = lambda_induced + lambda_climb
    
    theta_0_l_new = theta_0_l + ((2*(Ct_l_req - Ct_j)*pi*R/Cl_alpha0/Nb) ...
                    + (lambda_total - lambda_total_j)*c0_l/2 * (2/3*TR_l+1/3)) * 3*(c0_l*(3*TR_l/4+1/4))^-1;

    % Calculate residual
    eps(count) = theta_0_l_new - theta_0_l;
    eps2(count) = (Ct_l_req - Ct_j);
    
    % update 
    theta_0_l = theta_0_l_new;
    % lambda_total = lambda_total_j;

    record(count) = theta_0_l;
    count = count+1;
end

theta_0_l = theta_0_l*180/pi;

% disp([' ----------- Lower rotor -----------']);
% disp(['Collective pitch angle initial estimate at Vc = ', num2str(Vc), 'm/s is: ', num2str(theta_0_l_initial*180/pi), ' degree.']);
% disp(['Converged value: ', num2str(180*theta_0_l_new/pi), ' degree.']);
% disp(['CT: ', num2str(Ct_j),'; Required CT: ',num2str(Ct_l_req)]);
% disp(['CP/CQ: ', num2str(Cp_j_l)]);
% disp(['Number of iteration: ', num2str(count-1)]);
% disp([' -----------------------------------']);
