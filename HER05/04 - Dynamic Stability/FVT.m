function [X,Y] = FVT(u,v,r,rho,deltaR,geometry,l_VT,NACA0015Data)
% Forces and Moments of Vertical Tailplane around ac
% Lift positive rightwards
% Drag positive backwards
% positive u in forward flight
% positive v from right
% deltaR [radian] positive when airfoil points left (yaw right)

% Extract constants
c = geometry.VT.c;     % chord length [m]
b = geometry.VT.b;     % span [m]

% Planform area
S = b*c;

% Calculate alpha
dv = v + r*l_VT;
phi = atand(dv/u);
alpha = phi + rad2deg(deltaR); % [deg]

if abs(alpha) > 19.75
    L = 0;
    D = 0;
    % error(sprintf('Vertical tail angle of attack limit exceeded (maximum 19.75 deg): %.2f deg',alpha))
else

    % Find sectional coefficients
    [C_L,C_D] = NACA0015(alpha,NACA0015Data);

    % Calculate Total Forces and Moments
    const = .5 * rho * (u^2 + dv^2) * S;
    L = C_L * const;
    D = C_D * const;
end

% Convert to body-fixed axes
X = L*sind(phi) - D*cosd(phi);
Y = -L*cosd(phi) - D*sind(phi);

end