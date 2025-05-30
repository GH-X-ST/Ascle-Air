function [LHT,DHT] = FHT(u,v,w,rho,delE)
% Lift & Drag Force of Horizontail Tailplane
% Lift positive upwards
% Drag positive backwards
% positive u in forward flight
% positive v while drifting right
% negative w while climbing
% del E positive when airfoil points up (pitch down)

TBD = 0.2; % Tail Boom Diametre
HTl = 0.9; % Horizontail Tail length (airfoil chord length)
HTw = 3; % Horizontal Tail width (VT to VT)
HTt = HTl*0.15; % Horizontal Tail thickness
VTl = HTl;
cutout = VTl*sind(20); %Loss in elevator width due to cutout of angle in degrees for VT

Area = HTl*(HTw-2*TBD-2*cutout);

alpha = atand(-w/u)+rad2deg(delE);

[CL,CD] = NACA0015(alpha);
LHT = CL*0.5*rho*u^2*Area;
DHT = CD*0.5*rho*u^2*Area;

end