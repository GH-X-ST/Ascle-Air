function [X,Y,Z,L,M,N] = RotorForces(states, ctrl, rho, geometry)
% 
% Outputs: T - thrust (+ve upwards)
%          H - Drag (+ve rearwards)
%          Y - Side force (+ve to the right)
%          Q - torque (+ve clockwise)
%          Mx - rolling moment (+ve right side down)
%          My - pitching moment (+ve nose up)

% Extract states
u = states(1);
v = states(2);
w = states(3);
p = states(4);
q = states(5);

% Extract Controls
theta_0 = ctrl(1);
Delta_theta = ctrl(2);
theta_1c = ctrl(3);
theta_1s = ctrl(4);

% Collective
theta_u = theta_0 + .5*Delta_theta;
theta_l = theta_0 - .5*Delta_theta;

options = optimoptions('fsolve','Display','off');

% Solve for upper rotor
test = @(C_T) BEMT(u,v,w,p,q,rho,C_T,theta_u,theta_1c,theta_1s,'upper',geometry) - C_T;
C_Tinit = geometry.Cw / 2;
C_Tu = fsolve(test,C_Tinit,options);
[~,Xu,Yu,Zu,Lu,Mu,Nu] = BEMT(u,v,w,p,q,rho,C_Tu,theta_u,theta_1c,theta_1s,'upper',geometry);

% Solve for lower rotor
test = @(C_T) BEMT(u,v,w,p,q,rho,C_T,theta_l,theta_1c,theta_1s,'lower',geometry) - C_T;
C_Tinit = geometry.Cw / 2;
C_Tl = fsolve(test,C_Tinit,options);
[~,Xl,Yl,Zl,Ll,Ml,Nl] = BEMT(u,v,w,p,q,rho,C_Tl,theta_l,theta_1c,theta_1s,'lower',geometry);

% Sum forces
X = Xu - Xl;
Y = Yu - Yl;
Z = Zu + Zl;
L = Lu - Ll;
M = Mu + Ml;
N = Nu - Nl;

end