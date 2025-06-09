function [Cl, Cd] = airfoilCoeffs(alpha_input, polar, Mach, degInput)
% AIRFOILCOEFFS  Returns Cl and Cd from AoA using XFOIL polar data
% alpha_input: scalar or array of angles of attack
% polar: struct with fields .alpha, .Cl, .Cd (from loadPolarData)
% degInput: true if alpha_input is in degrees (default: true)

    if nargin < 4
        degInput = true;
    end

    % Convert radians to degrees if needed
    if ~degInput
        alpha_input = rad2deg(alpha_input);
    end

    % Interpolate Cl and Cd
    Cl = interp1(polar.alpha, polar.Cl, alpha_input, 'linear', 'extrap');
    Cd = interp1(polar.alpha, polar.Cd, alpha_input, 'linear', 'extrap');

    % Compressibility correction
    Cl = Cl/sqrt(1-Mach^2);
    Cd = Cd/sqrt(1-Mach^2);
end