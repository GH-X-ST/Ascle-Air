function constants = getConstants()
%SETUPCONSTANTS Returns a struct with all constants used in rotor modeling

    % Parameters from initial sizing of rotor
    constants.GTOW = 3491;                        % kg
    constants.Nb = 3;                             % number of blades
    constants.AR = 20;                            % aspect ratio
    constants.Vtip = 220;                         % tip speed [m/s]
    constants.R = 5.2992;                         % rotor radius [m]
    constants.r0 = 0.2;                           % non-dimensional root cutout
    constants.Cd0 = 0.007;                        % profile drag coefficient
    constants.Omega = constants.Vtip/constants.R; % angular speed [rad/s]

    constants.theta_tw_u = deg2rad(-11);          % root to tip [rad]
    constants.theta_tw_l = deg2rad(-6);           % root to tip [rad]

    % Derived geometric values
    constants.A = pi * constants.R^2;                                       % rotor disc area [m^2]
    constants.Ae = constants.A - pi * (constants.r0 * constants.R)^2;       % effective area [m^2]
    constants.cbar = constants.R / constants.AR;                            % mean chord [m]
    constants.sigma_e = constants.Nb * constants.cbar / (pi * constants.R); % equivalent solidity

    % Operating conditions
    constants.rho = 1.225;    % air density at sea level [kg/m^3]
    constants.rho_cr = 1.058; % cruise condition density [kg/m^3]

    % Weight coefficient (non-dimensional)
    constants.Cw = constants.GTOW * 9.81 / (constants.rho * constants.Ae * constants.Vtip^2);

    % fuselage drag
    constants.Cd_f = 0.12; % fuselage drag coeff
    constants.S_f = 4;     % fuselage frontal area

    % BEMT descritisation
    constants.N = 100;
    constants.r = linspace(constants.r0,1,constants.N);

    % Flapping
    constants.e = 0.08;       % eq. flap hinge
    constants.theta_1c = 0.3; % deg;
    constants.theta_1s = -2;  % deg;  

end