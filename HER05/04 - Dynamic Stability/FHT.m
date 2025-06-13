function [X,Z] = FHT(u,w,q,rho,deltaE,geometry,l_HT,NACA0015Data)
% Forces and Moments of Horizontal Tailplane around ac
% Lift positive upwards
% Drag positive backwards
% positive u in forward flight
% positive w while descending
% deltaE [radian] positive when airfoil pitches down (nose-up pitching moment)

% Extract parameters
S = geometry.HT.S;                     % area [m^2]

% Calculate alpha
dw = w + q*l_HT;
phi = atan2d(dw,u);
alpha = phi - rad2deg(deltaE); % [deg]

if abs(alpha) > 19.75 
    L = 0;
    D = 0;
    % error('Horizontal tail angle of attack limit exceeded (maximum 19.75 deg): %.2f deg',alpha)
else
    % Find sectional coefficients
    [C_L,C_D] = NACA0015(alpha,NACA0015Data);

    % Calculate Total Forces
    const = .5 * rho * (u^2 + w^2) * S;
    L = C_L * const;
    D = C_D * const;
end

% Convert to body-fixed axes
X = L*sind(phi) - D*cosd(phi);
Z = -L*cosd(phi) - D*sind(phi);

end
