function geometry = Parameters(flightCondition)

% Rotor
geometry.GTOW = 3491;                          % gross take-off weight [kg]
geometry.Nb = 3;                               % number of blades
geometry.AR = 20;                              % aspect ratio
geometry.Vtip = 220;                           % tip speed [m/s]
geometry.R = 5.2992;                           % rotor radius [m]
geometry.r0 = 0.2;                             % non-dimensional root cutout ??
geometry.Cd0 = 0.007;                          % profile drag coefficient
geometry.Omega = geometry.Vtip/geometry.R;     % angular speed [rad/s]
geometry.z_R = [];                             % rotor separation [m]
geometry.tw_u = deg2rad(-11);                  % upper rotor twist [rad]
geometry.tw_l = deg2rad(-6);                   % lower rotor twist [rad]
geometry.e = 0.08;                             % hinge offset [1/R]
geometry.CLalpha = 6.2305;                     % lift curve slope [1/rad]
geometry.m = 8.88;                             % Mass distribution [kg/m]
geometry.Ib = (1/3)*geometry.m*geometry.R^3*(1-geometry.e)^3; 
geometry.v_beta = sqrt(1+1.5*geometry.e/(1-geometry.e)); % natural frequency
geometry.z_R = 0.52992; % rotor separation [m]

% Flexbeam
E1 = 139.2e9;
b1 = 0.015;
b2 = 0.13;
h1 = 0.03;
h2 = 0.016;
Ix1 = (1/12) * b1 * h1^3;
Ix2_local = (1/12) * b2 * h2^3;
d = (h1/2) - (h2/2);
A2 = b2 * h2;
Ix2 = Ix2_local + A2 * d^2;
Ix = Ix1 + Ix2;
geometry.K_beta =(3*E1*Ix)/(0.14*geometry.R);

% Derived Geometric Values
geometry.A = pi * geometry.R^2;                                       % rotor disc area [m^2]
geometry.Ae = geometry.A - pi * (geometry.r0 * geometry.R)^2;       % effective area [m^2]
geometry.cbar = geometry.R / geometry.AR;                            % mean chord [m]
geometry.sigma_e = geometry.Nb * geometry.cbar / (pi * geometry.R); % equivalent solidity

% Weight coefficient (non-dimensional)
geometry.Cw = geometry.GTOW * 9.81 / (1.225 * geometry.A * geometry.Vtip^2);

% Fuselage Drag
geometry.Cd_f = 0.12;  % fuselage drag coeff
geometry.S_f = 4;      % fuselage frontal area

% Horizontal Tail
geometry.HT.z = 0.682;             % vertical distance from ac to lower rotor hub [m]
geometry.HT.x = 5.5;               % horizontal distance from ac to lower rotor hub [m]
geometry.HT.c = 0.9;               % chord length [m]
geometry.HT.b = 5;                 % span [m]
geometry.HT.S = geometry.HT.c*geometry.HT.b;

% Vertical Tail
geometry.VT.z = 0.682; % vertical distance from ac to lower rotor hub [m]
geometry.VT.x = 5.5;   % horizontal distance from ac to lower rotor hub [m]
geometry.VT.c = 0.9;   % chord length [m]
geometry.VT.b = 1;     % span [m]
geometry.VT.S = geometry.VT.c*geometry.VT.b;

%% Moments of Inertia [kg m^2] & CG [m]

switch flightCondition

    case 1
    % 1 - full fuel, full payload
    geometry.flightCondition.m = 2.8676e+03;
    geometry.flightCondition.Ixx = 5.9503e+03;
    geometry.flightCondition.Iyy = 4.0466e+04;
    geometry.flightCondition.Izz = 3.4703e+04;
    geometry.flightCondition.Ixz = 1.2291e+04;
    geometry.flightCondition.xCG = 3.2496;
    geometry.flightCondition.yCG = -0.0215;
    geometry.flightCondition.zCG = 1.2354;
    
    case 2
    % 2 - zero fuel, full payload
    geometry.flightCondition.m = 2.1676e+03;
    geometry.flightCondition.Ixx = 5.7107e+03;
    geometry.flightCondition.Iyy = 3.2135e+04;
    geometry.flightCondition.Izz = 2.6611e+04;
    geometry.flightCondition.Ixz = 1.0898e+04;
    geometry.flightCondition.xCG = 3.2010;
    geometry.flightCondition.yCG = -0.0285;
    geometry.flightCondition.zCG = 1.4454;
    
    case 3
    % 3 - full fuel, zero payload
    geometry.flightCondition.m = 2.2614e+03;
    geometry.flightCondition.Ixx = 5.1058e+03;
    geometry.flightCondition.Iyy = 3.0526e+04;
    geometry.flightCondition.Izz = 2.5451e+04;
    geometry.flightCondition.Ixz = 9.8181e+03;
    geometry.flightCondition.xCG = 3.1349;
    geometry.flightCondition.yCG = -0.0153;
    geometry.flightCondition.zCG = 1.2687;
    
    case 4
    % 4 - 50% fuel, 50% payload
    geometry.flightCondition.m = 2.3376e+03;
    geometry.flightCondition.Ixx = 5.6409e+03;
    geometry.flightCondition.Iyy = 3.4489e+04;
    geometry.flightCondition.Izz = 2.8995e+04;
    geometry.flightCondition.Ixz = 1.1102e+04;
    geometry.flightCondition.xCG = 3.2599;
    geometry.flightCondition.yCG = -0.0264;
    geometry.flightCondition.zCG = 1.3532;

end

end