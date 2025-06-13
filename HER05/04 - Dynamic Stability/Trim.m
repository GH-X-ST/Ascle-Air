function ctrl = Trim(states,altitude,flightCondition)
% Inputs:  x - state vector [u v w p q r phi theta psi]'
%          u - controls vector [theta_0 Delta_theta theta_1c theta_1s delta_r delta_e]'

% Extract states
u = states(1);
v = states(2);
w = states(3);
p = states(4);
q = states(5);
r = states(6);
phi = states(7);
theta = states(8);
psi = states(9);

% Extract mass and MOI
geometry = Parameters(flightCondition);
m = geometry.flightCondition.m;
I_xx = geometry.flightCondition.Ixx;
I_yy = geometry.flightCondition.Iyy;
I_zz = geometry.flightCondition.Izz;
I_xz = geometry.flightCondition.Ixz;
g = 9.81;

fun = @(ctrl) ...
    FM(states,[ctrl],altitude,flightCondition) ....
    -[      m*(q*w - r*v) + m*g*sin(theta);
            m*(r*u - p*w) - m*g*sin(theta)*cos(phi);
            m*(p*v - q*u) - m*g*cos(theta)*cos(phi);
            I_xz*(p*q) - (I_yy - I_zz)*q*r;
            -I_xz*(r^2 - p^2) - (I_zz - I_xx)*p*r;
            I_xz*(q*r) - (I_xx - I_yy)*p*q;              
            ];
x0 = [deg2rad(5) 0 0 0 0 0];

options = optimoptions('fsolve','Display','off');
% options = optimoptions('fsolve','Display','iter-detailed','FunctionTolerance',1e1,'StepTolerance',1e-3);
ctrl = fsolve(fun,x0,options);


end