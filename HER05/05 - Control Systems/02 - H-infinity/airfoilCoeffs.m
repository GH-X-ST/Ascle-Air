function [Cl, Cd] = airfoilCoeffs(alpha_input, polar, Mach, degInput)
% AIRFOILCOEFFS  Returns Cl and Cd from AoA using XFOIL polar data
%   Any alpha_input outside the range polar.alpha will give Cl=Cd=0.
%
% alpha_input: scalar or array of angles of attack
% polar:       struct with fields .alpha, .Cl, .Cd (from loadPolarData)
% Mach:        Mach number for compressibility correction
% degInput:    true if alpha_input is in degrees (default: true)

    if nargin < 4
        degInput = true;
    end

    % Convert radians to degrees if needed
    if ~degInput
        alpha_input = rad2deg(alpha_input);
    end

    % Preallocate
    Cl = zeros(size(alpha_input));
    Cd = zeros(size(alpha_input));

    % Determine which inputs lie inside the polar’s alpha range
    amin = min(polar.alpha);
    amax = max(polar.alpha);
    inRange = alpha_input >= amin & alpha_input <= amax;

    % Only interpolate for in-range angles; out-of-range stay at zero
    if any(inRange)
        % linear interpolation (no extrapolation)
        Cl(inRange) = interp1(polar.alpha, polar.Cl,  ...
                              alpha_input(inRange), 'linear');
        Cd(inRange) = interp1(polar.alpha, polar.Cd,  ...
                              alpha_input(inRange), 'linear');
    end

    % Compressibility correction on the valid portion
    Cl(inRange) = Cl(inRange) / sqrt(1 - Mach^2);
    Cd(inRange) = Cd(inRange) / sqrt(1 - Mach^2);
end