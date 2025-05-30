function [LVT,DVT] = FVT(u,v,w,rho,delR)
% Lift & Drag Force of Vertical Tailplane
% Lift positive rightwards
% Drag positive backwards
% positive u in forward flight
% positive v while drifting right
% negative w while climbing
% del R positive when airfoil points right (yaw left)

HTl = 0.9; % Horizontail Tail length (airfoil chord length)
HTt = HTl*0.15; % Horizontal Tail thickness
VTl = HTl; % Vertical Tail length (airfoil chord length) SAME AS HTl
VTw = 1; % Vertical Tail width (top to bottom length)

Area = VTl*VTw;

alpha = atand(v/u)+rad2deg(delR);

[CL,CD] = NACA0015(alpha);
LVT = CL*0.5*rho*u^2*Area;
DVT = CD*0.5*rho*u^2*Area;

end