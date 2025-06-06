% Required parameters 
% subscripts of lat and lon indicate lateral and longitudinal respectively
% A: dynamics/state matrix
% B: input matrix
% C: sensor matrix
% D: direct term/feedforward matrix

% Input-output mapping
input_lon = 1;   % 1: theta_0 to theta, 2: theta_1s to theta
input_lat = 1;   % 1: theta_0, 2: theta_1c
output_lat = 1;  % 1: phi, 2: psi

% Open-loop plant
P_lon = ss(A_lon, B_lon, C_lon, D_lon);
P_lat = ss(A_lat, B_lat, C_lat, D_lat);

% Assign correct state names
P_lon.StateName = {'u-vel'; 
                   'w-vel'; 
                   'pitch rate';
                   'pitch angle'};

P_lat.StateName = {'v-vel'; 
                   'roll rate'; 
                   'yaw rate'; 
                   'roll angle';
                   'yaw angle'};

% Assign input variable names
P_lon.InputName = {'theta_0'; 
                   'theta_1s';
                   'delta_e'};

P_lat.InputName = {'Delta_theta';
                   'theta_1c';
                   'delta_r'};

% Assign output variable names
P_lon.OutputName = {'theta'};
P_lat.OutputName = {'phi'; 'psi'};

P_lon_tf = tf(P_lon); P_lon_tf = P_lon_tf(input_lon);
% Create transfer functions
P_lat_tf = tf(P_lat); P_lat_tf = P_lat_tf(input_lat, output_lat);

% Display open-loop poles
disp("Open-loop longitudinal poles")
damp(P_lon_tf)

disp("Open-loop lateral poles")
damp(P_lat_tf)


% Longitudinal pitch controller
% Pre-compensator and post-compensator design

% Low-pass pre-compensator
K1 = 4;  % Low pass filter gain
W1 = tf(K1, [1, 8]);