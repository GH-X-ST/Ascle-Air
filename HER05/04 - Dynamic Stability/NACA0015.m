function [CL,CD,CM] = NACA0015(alpha,NACA0015Data)
% Calculates sectoional drag, lift and pitching moment of empenage 
% components
% alpha must be between -19.75 deg and 19.75 deg

% Linear Interpolation
CL = fixed.interp1(NACA0015Data.alpha,NACA0015Data.Cl,alpha);
CD = fixed.interp1(NACA0015Data.alpha,NACA0015Data.Cl,alpha);
CM = fixed.interp1(NACA0015Data.alpha,NACA0015Data.Cm,alpha);

end