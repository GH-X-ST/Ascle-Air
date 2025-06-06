% Some estimates for the state and control matricies
% Only valid in hover
% Magnitudes may be off

clear
clc
close all

% Final function will look like this with inputs of forward flight speed,
% climb rate etc.
function [A,B,A_lon,B_lon,A_lat,B_lat] = Linearise(x,u)
% Inputs:  x - state vector [u v w p q r phi theta psi]'
%          u - controls vector [theta_0 Delta_theta theta_1c theta_1s delta_r delta_e]'
% Outputs: A - state matrix (9x9)
%          B - Controls matrix (9x6)
%          A_lon - longitudinal states matrix (4x4) [u w q theta]'
%          B_lon - longitudinal controls matrix (4x3) [theta_0 theta_1s delta_e]'
%          A_lat - lateral states matrix (5x5) [v p r phi psi]'
%          B_lat - lateral controls matrix (5x3) [Delta_theta theta_1c delta_r]'


Xu = -0.01;
Xw = 0.01;
Xq = 3.608;
Zu = 0;
Zw = -0.03;
Zq = 0;
Mu = 0.01;
Mw = -0.005;
Mq = -3.75;
A_lon = [Xu Xw Xq -9.81;
         Zu Zw Zq 0    ;
         Mu Mw Mq 0    ;
         0  0  1  0    ];

Xtheta0 = 9.8;
Xtheta1s = -31.1;
Ztheta0 = -295.2;
Ztheta1s = 0;
Mtheta0 = 0;
Mtheta1s = 11;
B_lon = [Xtheta0 Xtheta1s 0;
         Ztheta0 Ztheta1s 0;
         Mtheta0 Mtheta1s 0;
         0       0        0];

Yv = -0.05;
Yp = -0.29;
Yr = 0;
Lv = -0.01;
Lp = -7.5;
Lr = 0.1;
Nv = -0.008;
Np = -1.5;
Nr = -0.25;
A_lat = [Yv Yp Yr 9.81 0;
         Lv Lp Lr 0    0;
         Nv Np Nr 0    0;
         0  1  0  0    0;
         0  0  1  0    0];

Ydeltatheta0 = 13.12;
Ytheta1c = -32.8;
Ldeltatheta0 = 5;
Ltheta1c = 25;
Ndeltatheta0 = -8;
Ntheta1c = -10;
B_lat = [Ydeltatheta0 Ytheta1c 0;
         Ldeltatheta0 Ltheta1c 0;
         Ndeltatheta0 Ntheta1c 0;
         0            0        0;
         0            0        0];

% These will be filled in on the final version
A = zeros(9,9);
B = zeros(9,6);

end

[~,~,A_lon,B_lon,A_lat,B_lat] = Linearise

evalslon = eig(A_lon);
figure
plot(real(evalslon),imag(evalslon),'x');

evalslat = eig(A_lat);
figure
plot(real(evalslat),imag(evalslat),'x')