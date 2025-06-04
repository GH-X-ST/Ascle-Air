function [T, H, Y, Q, Mx, My, Cq_i, beta1c_UR, beta1s_UR] = Stability_forward_upper(u, v, w, a_input, rho_input, mu_input, m_input, g_input, polar, theta_0_u, theta_0_l, theta_1c, theta_1s)

% input in degrees

% Get all constants
constants = getConstants();  % or use setupConstants() directly if you prefer

% Access individual constants using dot notation
% GTOW    = constants.GTOW;      % kg
GTOW    = m_input;
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
Cd_f    = constants.Cd_f;      % fuselage drag coeff
S_f     = constants.S_f;       % fuselage frontal area
Omega   = constants.Omega;     % angular speed

% e = constants.e; % flapping hinge

% theta input is radian now (Hanchen)

% theta_1c = constants.theta_1c; % deg;
% theta_1s = constants.theta_1s; % deg;

% addpath('Airfoil');
% polar = loadPolarData('xf-rc410-il-1000000.txt');

theta_tw_u = rad2deg(constants.theta_tw_u);
TR_u = 1;
% v  = 0;
Vc = w;
Vx = u;
Cw = (g_input*GTOW/(rho_input*Ae*Vtip^2));

% Parameter Input
theta_tw_u = theta_tw_u/R*pi/180*R;
theta_0_u = theta_0_u*pi/180;

% BladeGeometryInput = [theta_tw_u,theta_tw_l, TR_u, TR_l, d];

% Chord distribution from Taper Ratio
c0_u = cbar*(3*((TR_u-1)/4 + 1/3))^-1;
c_u = @(r) c0_u* ((TR_u - 1)*r + 1); 

% Solidity 
sigma_u = @(r) Nb*c_u(r)/(pi*R);

% Airfoil data
Cl_alpha0 = 6.2305; % lift curve slope at Mach = 0
AoA_zl = -1*pi/180; % zero-lift AoA

% Non dimensional terms
lambda_c = Vc/Vtip;
mu_x = Vx/Vtip;
k_hover = 1.1; % induced loss factor
k_cr = k_hover*cosh(7.5*mu_x^2);

% --------- Initial Inflow Calculation ----------- %
% drag of helicopter in forward flight
D = 0.5*rho_input*Vx^2*S_f*Cd_f; % fuselage drag
alpha_s = atan(Vc/Vx) + atan(D/(g_input*GTOW)); % Rotor angle of attacked


Ct_u_req = Cw/2/cos(alpha_s);

% inflow calculation using iteration in forward flight
% Define the function whose root we want to find
f = @(lambda) lambda - (mu_x * tan(alpha_s) + Ct_u_req / (2 * sqrt(mu_x^2 + lambda.^2)));
lambda0 = sqrt(Ct_u_req/2); % Initial guess from hover value
lambda_ff = fzero(f, lambda0); % Use fzero to find total inflow 

% induced inflow
lambda_induced_ff = Ct_u_req/(2*sqrt(mu_x^2 + lambda_ff^2));

% initial estimate of collective pitch - simplified CT expression in
% forward flight
% theta_0_u_initial = ((2*Ct_u_req*pi*R/Cl_alpha0/Nb) ...
%                     - theta_tw_u*c0_u/4 * (4/5*TR_u+1/5)*(1+mu_x^2) ...
%                     + lambda_ff*c0_u/2 * (2/3*TR_u+1/3)) * 3*(1+3/2*mu_x^2)^(-1)*(c0_u*(3*TR_u/4+1/4))^(-1);

% lambda_induced_j = lambda_induced_ff;


% Discretisation
N = 30;
tol = 1e-5;

% radial discretisation
r = linspace(r_0,1,N);

% disc discretisation - for forward flight analysis
NN = 90;
azimuth = linspace(0,2*pi,NN);

dr = r(2)-r(1);
dpsi = azimuth(2)-azimuth(1);

% theta_0_u = theta_0_u_initial;

eps = 1;

% initialise arrays
lambda_j = zeros(1,length(r));
dCl_j = zeros(1,length(r));
sigma = zeros(1,length(r));
pitch = zeros(1,length(r)); 
AoA = zeros(length(azimuth),length(r));
inflow_angle = zeros(1,length(r));
prandlt_loss = zeros(1,length(r));


dT = zeros(length(azimuth),length(r));
dH = zeros(length(azimuth),length(r));
dQ = zeros(length(azimuth),length(r));
dY = zeros(length(azimuth),length(r));
dCt_j = zeros(length(azimuth),length(r));
dCt_j2 = zeros(length(azimuth),length(r));
dCq_j = zeros(length(azimuth),length(r));
dCp_j = zeros(length(azimuth),length(r));
dCp_induced_j = zeros(length(azimuth),length(r));
dCp_profile_j = zeros(length(azimuth),length(r));
dCq_j2 = zeros(length(azimuth),length(r));
dCH_j = zeros(length(azimuth),length(r));
dCY_j = zeros(length(azimuth),length(r));
Cq_j = zeros(1,length(azimuth));
Cp_j = zeros(1,length(azimuth));
Cp_induced_j = zeros(1,length(azimuth));
Cp_profile_j = zeros(1,length(azimuth));
Ct_j = zeros(1,length(azimuth));
CH_j = zeros(1,length(azimuth));
CY_j = zeros(1,length(azimuth));
lambda = zeros(length(azimuth),length(r));
lambda_induced = zeros(length(azimuth),length(r));
U_t = zeros(length(azimuth),length(r));
U_p = zeros(length(azimuth),length(r));
Mach = zeros(length(azimuth),length(r));
Reynolds = zeros(length(azimuth),length(r));

count = 1;

% while abs(eps) > tol
% for ii = 1:5

    % flapping respons - at 140 knots, 
    % assuming e = 0.1; theta_1s_deg = -3; % from paper -- flight test; theta_1c_deg = 0.2; % from paper -- flight test
    % [beta,beta_dot] = getFlappingForwardResponse(Vx, e,Ct_u_req,theta_0_u*180/pi,theta_tw_u*180/pi,theta_1c,theta_1s, 'upper'); % inputs are degree

    % [Beta, Beta_dot, beta_1c_rad, beta_1s_rad] = getFlappingResponse(Vx, rho_input, m_input, CT, theta0_u_rad, theta0_l_rad, theta_1c_rad, theta_1s_rad, rotor)
    [beta, beta_dot, beta1c_UR, beta1s_UR] = getFlappingResponse(Vx, rho_input, m_input, Ct_u_req, theta_0_u, theta_0_l, theta_1c, theta_1s, 'upper'); % inputs are radian

    % at each discretised azimuth location dpsi
    for i = 1:length(azimuth)
        % for each discretised blade element
        psi = azimuth(i);
        for j = 1:length(r)
                       
            sigma_j = sigma_u(r(j)); % 
            pitch_j = theta_0_u + theta_tw_u*(r(j)) + theta_1c*cos(psi) + theta_1s*sin(psi); % involve pilot cyclic input

            % linear inflow model at each blade element
            [lambda_i,lambda_j] = LinearInflow_ff(mu_x, lambda_c, alpha_s, Ct_u_req, r(j), psi);
Mach
            lambda(i,j) = lambda_j;
            lambda_induced(i,j) = lambda_i;

            % calculate blade sectional velocity
            U_T = Omega*r(j)*R + mu_x*Vtip*sin(psi);
            U_P = (lambda_c + lambda_i)*Vtip + r(j)*R*beta_dot(psi) + Vx*beta(psi)*cos(psi); % ignored flapping terms
            U_t(i,j) = U_T;
            U_p(i,j) = U_P;
            % U_R = mu_x*Vtip*cos(psi);
            U = sqrt(U_T^2 + U_P^2);
            Mach(i,j) = U/a_input;
            Reynolds(i,j) = rho_input*sqrt(U_T^2 + U_P^2)*c_u(r(j))/mu_input;
    

            % Prandlt tip loss
            phi = atan((U_P/U_T));
            InflowAngle(i,j) = phi;
            f = Nb/2* (1-r(j))/(r(j)*phi);
            F = 2/pi * acos(exp(-f));
    
            % local angle of attack
            AoA(i,j) = pitch_j - phi -AoA_zl;

            % CL and CD from airfoil data (compressibility corrected)
            [Cl, Cd] = airfoilCoeffs(AoA(i,j)*180/pi, polar, Mach(i,j));

            % record
            sigma(j) = sigma_j;
            pitch(j) = pitch_j;
            inflow_angle(j) = phi;
            prandlt_loss(j) = F;
 
            % incremental thrust
            % dCt_j(j) = sigma_j*Cl_alpha/2 * (pitch_j - phi) * r(j)^2 * dr; % without small angle approx
            
            % sectional force
            % Cd = 0.02; % airfoil
            if U_T <=0 % reverse flow
                dFz = 0;
                dFx = 0;
            else
                dFz2 = 0.5*rho_input*c_u(r(j))*(Cl_alpha0/sqrt(1-Mach(i,j)^2))*((pitch_j) * U_T^2 - U_P*U_T)*dr*R;
                dFz = 0.5*rho_input*U^2*c_u(r(j))*Cl*dr*R;
                dFx2 = 0.5*rho_input*c_u(r(j))*(Cl_alpha0/sqrt(1-Mach(i,j)^2))*((pitch_j) * U_T*U_P - U_P^2 + Cd/(Cl_alpha0/sqrt(1-Mach(i,j)^2))*U_T^2)*dr*R;
                dFx = 0.5*rho_input*U^2*c_u(r(j))*Cd*dr*R;
            end
            dFr = -beta(psi)*dFz;

            dT(i,j) = Nb*dFz;                             % Rotor thrust
            dH(i,j) = Nb*(dFx*sin(psi) + dFr*cos(psi));   % Rotor drag force
            dY(i,j) = Nb*(-dFx*cos(psi) + dFr*sin(psi));  % Rotor side force
            dQ(i,j) = Nb*dFx*r(j)*R;
            
            dMx(i,j) = Nb*dFz*r(j)*R*sin(psi);
            dMy(i,j) = Nb*dFz*r(j)*R*cos(psi);
            % sectional thrust coeff
            dCt_j(i,j) = dT(i,j)/(rho_input*Ae*Vtip^2);
            dCt_j2(i,j) = Nb*dFz2/(rho_input*Ae*Vtip^2);

            % sectional torque coeff
            dCq_j(i,j) = dQ(i,j)*Omega/(rho_input*Ae*Vtip^3);
            dCq_j2(i,j) = Nb*dFx2*r(j)*R*Omega/(rho_input*Ae*Vtip^3);

            % sectional rotor drag coeff
            dCH_j(i,j) = dH(i,j)/(rho_input*Ae*Vtip^2);

            % sectional rotor side force coeff
            dCY_j(i,j) = dY(i,j)/(rho_input*Ae*Vtip^2);

            % dCl_j(i,j) = Cl_alpha*(pitch_j - phi - AoA_zl);
            dCl_j(i,j) = Cl;
           
            
            % sectional power
            dCp_induced_j(i,j) = dCt_j(i,j) * lambda_induced(i,j); % induced power
            dCp_profile_j(i,j) = Nb*dFx*U_T/(rho_input*Ae*Vtip^3); % profile power
            dCp_j(i,j) = k_cr*dCp_induced_j(i,j) + dCp_profile_j(i,j);
 
        end
    
        % Sum sectional thrust to get overall Ct at the azimuth angle
        Ct_j(i) = sum(dCt_j(i,:));
        Cq_j(i) = sum(dCq_j(i,:));
        CH_j(i) = sum(dCH_j(i,:));
        CY_j(i) = sum(dCY_j(i,:));

        Mx_j(i) = sum(dMx(i,:));
        My_j(i) = sum(dMy(i,:));
        % Sum sectional power to get overall Cp at the azimuth angle
        Cp_j(i) = sum(dCp_j(i,:));
        Cp_induced_j(i) = sum(dCp_induced_j(i,:));
        Cp_profile_j(i) = sum(dCp_profile_j(i,:));

        % Ct_j2 = sum(dCt_j2)
        record_Ct(count) = Ct_j(i);

    end

    Ct_i(count) = sum(Ct_j*dpsi)/(2*pi);
    Cq_i = sum(Cq_j*dpsi)/(2*pi);
    CH_i(count) = sum(CH_j*dpsi)/(2*pi);
    CY_i(count) = sum(CY_j*dpsi)/(2*pi);

    Mx = sum(Mx_j*dpsi)/(2*pi);
    My = sum(My_j*dpsi)/(2*pi);

    % output
    T = Ct_i(count)*(rho_input*Ae*Vtip^2);
    H = CH_i(count)*(rho_input*Ae*Vtip^2);
    Y = CY_i(count)*(rho_input*Ae*Vtip^2);
    Q = Cq_i*R*(rho_input*Ae*Vtip^2);

    % Power coeff
    Cp_i(count) = sum(Cp_j*dpsi)/(2*pi);
    Cp_induced_i(count) = sum(Cp_induced_j*dpsi)/(2*pi);
    Cp_profile_i(count) = sum(Cp_profile_j*dpsi)/(2*pi);

    % Lift

%     theta_0_u_new = theta_0_u + ((2*(Ct_u_req - Ct_i(count))*pi*R/Cl_alpha0/Nb)) * 3*(c0_u*(3*TR_u/4+1/4))^-1;
% 
%     % Calculate residual
%     eps(count) = theta_0_u_new - theta_0_u;
%     eps2(count) = (Ct_u_req - Ct_i(count));
% 
%     % update 
%     theta_0_u = theta_0_u_new;
%     % lambda_total = lambda_total_j;
% 
%     record(count) = theta_0_u;
%     count = count+1;
% end

% Lift to drag ratio
L_to_D = Ct_u_req*mu_x/(Cp_i(end)); % rotor lift-to-drag
L_to_D2 = Ct_u_req*mu_x/(Cp_i(end)+0.5*S_f*Cd_f/Ae*mu_x^3); % helicopter lift-to-drag


% theta_0_u = theta_0_u;

% disp([' ----------- Upper rotor -----------']);
% disp(['Collective pitch angle initial estimate at Vx = ', num2str(Vx), 'm/s and Vc = ', 'm/s is: ', num2str(theta_0_u_initial*180/pi), ' degree.']);
% disp(['Converged value: ', num2str(180*theta_0_u_new/pi), ' degree.']);
% disp(['CT: ', num2str(Ct_i(end)),'; Required CT: ',num2str(Ct_u_req)]);
% disp(['CQ: ', num2str(Cq_i(end))]);
% disp(['CP: ', num2str(Cp_i(end)),'; CP_induced: ', num2str(Cp_induced_i(end)),'; CP_profile: ', num2str(Cp_profile_i(end))]);
% disp(['Rotor Lift-to-Drag: ', num2str(L_to_D)]);
% disp(['Helicopter Lift-to-Drag: ', num2str(L_to_D2)]);
% disp(['Number of iteration: ', num2str(count-1)]);
% disp([' -----------------------------------']);
% 
% % Setup polar grid
% azimuth = linspace(0, 2*pi, NN);
% r = linspace(r_0, 1, N);
% 
% % Adjust azimuth angle so 0 is at -Y, rotation is clockwise
% psi_rot = azimuth - pi/2;
% 
% % Generate meshgrid and transpose to match data
% [PSI_grid, R_grid] = meshgrid(psi_rot, r);   % [N, NN]
% PSI_grid = PSI_grid';                        % [NN, N]
% R_grid = R_grid';                            % [NN, N]
% 
% % Convert to Cartesian
% [X, Y] = pol2cart(PSI_grid, R_grid);
% 
% % Plotting
% scr = get(0, 'ScreenSize');
% figWidth = 800; figHeight = 700;
% figLeft = scr(3) - figWidth - 100;
% figBottom = (scr(4) - figHeight)/2;
% 
% figure('Position', [figLeft, figBottom, figWidth, figHeight]);
% 
% 
% % --- 1. AoA ---
% 
% % Create mask: 1 for valid (forward flow), NaN for reverse flow
% mask = double(U_t >= 0);     % forward flow = 1, reverse = 0
% mask(mask == 0) = NaN;       % convert 0 to NaN to exclude in plot
% 
% % Apply mask to your AoA or other data field
% AoA_masked = AoA .* mask;
% 
% subplot(3,2,1);
% surf(X, Y, AoA_masked*180/pi, 'EdgeColor', 'none');
% axis equal tight;
% title('Angle of Attack (AoA)', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo; caxis([-10 10]);   % <<< Limit color scale here
% view(2); box on;
% 
% % --- Overlay Reverse Flow Region (Solid Color) ---
% % Use mask where U_t < 0
% reverse_mask = U_t < 0;
% 
% % Create a patch layer: same shape as AoA, but filled only where reverse
% reverse_overlay = NaN(size(AoA));
% reverse_overlay(reverse_mask) = 1;  % All 1s in reverse flow
% hold on
% % Use surf again with solid color for overlay
% h_overlay = surf(X, Y, reverse_overlay, ...
%     'EdgeColor', 'none', ...
%     'FaceColor', [0.3 0.3 0.3], ...   % Gray color
%     'FaceAlpha', 0.8);                % Optional: 50% transparency
% 
% % --- 2. dCl_j ---
% subplot(3,2,2);
% surf(X, Y, dCq_j, 'EdgeColor', 'none');
% axis equal tight;
% title('Sectional dCq', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;

% % --- 3. lambda ---
% subplot(3,2,3);
% surf(X, Y, lambda_induced, 'EdgeColor', 'none');
% axis equal tight;
% title('Induced Inflow \lambda_i', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;
% 
% % --- 4. dCt (induced inflow) ---
% subplot(3,2,4);
% surf(X, Y, dCt_j, 'EdgeColor', 'none');
% axis equal tight;
% title('Sectional dCT', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;
% 
% % --- 5. U_t (In plane speed) ---
% subplot(3,2,5);
% surf(X, Y, U_t, 'EdgeColor', 'none');
% axis equal tight;
% title('In plane speed', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;
% 
% % --- 6. U_p (In plane speed) ---
% subplot(3,2,6);
% surf(X, Y, U_p, 'EdgeColor', 'none');
% axis equal tight;
% title('Out of plane speed', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;

% figure;
% subplot(1,2,1);
% surf(X, Y, Mach, 'EdgeColor', 'none');
% axis equal tight;
% title('Local Mach', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;
% 
% subplot(1,2,2);
% surf(X, Y, Reynolds, 'EdgeColor', 'none');
% axis equal tight;
% title('Local Reynolds No.', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;
% 
% figure
% % subplot(1,2,1)
% surf(X, Y, dCp_j, 'EdgeColor', 'none');
% axis equal tight;
% title('sectional dCp', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;

% subplot(1,2,2)
% surf(X, Y, dCt_j./dCH_j, 'EdgeColor', 'none');
% axis equal tight;
% title('Sectional L/D', 'FontWeight', 'bold');
% xlabel('X (r/R)'); ylabel('Y (r/R)');
% colorbar; colormap turbo;
% view(2); box on;

end