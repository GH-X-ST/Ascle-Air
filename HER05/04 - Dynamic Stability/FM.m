function Forces = FM(state,ctrl,altitude,flightCondition)
% 
% Inputs:  state - [u v w p q r phi theta psi]
%          ctrl - [theta_0 Delta_theta theta_1c theta_1s delta_r delta_e]
%          altitude - altitude [ft]
%          flightCondition - 1 - full fuel, full payload
%                            2 - zero fuel, full payload
%                            3 - full fuel, zero payload
%                            4 - 50% fuel, 50% payload
% Outputs: Forces - [X Y Z L M N]
%          X - forward force [N]
%          Y - side force [N]
%          Z - vertical force [Nm]
%          L - rolling moment [Nm]
%          M - pitching moment [Nm]
%          N - yawing moment [Nm]

% Extract states
u = state(1);
v = state(2);
w = state(3);
p = state(4);
q = state(5);
r = state(6);

% Speed of sound and density
[~,~,~,rho] = atmosisa(convlength(altitude,'ft','m'));

% Constant parameters
geometry = Parameters(flightCondition);

%% Moment arms

theta_FP = atan2(w,u);
phi_F = atan2(v,w);

% CG location relative to lower rotor
x_CG = 3.2599 - geometry.flightCondition.xCG;
z_CG = 2.66847 - geometry.flightCondition.zCG;

% Tail arm
h_VT = z_CG - geometry.VT.z;
h_HT = z_CG - geometry.HT.z;
l_VT = x_CG + geometry.VT.x;
l_HT = x_CG + geometry.HT.x;
VT_l = sqrt(l_VT^2 + h_VT^2);
HT_l = sqrt(l_HT^2 + h_HT^2);


%% Rotors

[XR,YR,ZR,LR,MR,NR] = RotorForces(state, ctrl, rho, geometry);


%% Empenage

load NACA0015Data.mat NACA0015Data

% Horizontal tailplane
[X_HT,Z_HT] = FHT(u,w,q,rho,ctrl(6),geometry,HT_l,NACA0015Data);

% Vertical tailplane
[X_VT,Y_VT] = FVT(u,v,r,rho,ctrl(5),geometry,VT_l,NACA0015Data);


%% Fuselage

[L_F, D_F, Y_F, Mx_F, My_F] = Fuselage_Aero(u, v, w, rho);


%% Sum Forces & Moments

% X - Longitudinal Force
X = XR - D_F*cos(theta_FP) - L_F*sin(theta_FP) + X_HT + X_VT;

% Y - Lateral Force
Y = YR + Y_F + Y_VT;

% Z - Vertical Force
Z = ZR + D_F * sin(theta_FP) - L_F * cos(theta_FP) + Z_HT;

% L - Rolling moment
L = LR + Mx_F + Y_VT*h_VT;

% M - Pitching moment
M = MR + My_F + Z_HT*l_HT;

% N - Yawing moment
N = NR + Y_VT*l_VT;

Forces = [X Y Z L M N]';

end