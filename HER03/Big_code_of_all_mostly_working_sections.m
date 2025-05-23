%%%     LAST EDITED 22/05/2025
%%%     Mostly works as it is meant to

%%%  ASSUMPTIONS 
%%%  Blade is rigid
%%%  Blade is not moving
%%%  Maximum loading case is that of approaching blade in forward flight
%%%  Uniform beam with central CG - 8kg/m
%%%  Consider oncoming velocity as (omega*r)+V_inf
%%%  Twist linear from 10 deg to 0 deg
%%%  No sweep
%%%  No taper ratio
%%%  Upwards as positive
%%%  Both upper and lower blade are the same

clear
clc

%%% INPUTS
n = 100;                            % number of steps
V_inf = 70;                         % Maximum oncoming velocity, m/s
R = 0.8*(10.6/2);                   % Lower rotor radius, m
steps = R/n;                        % 
r = 0:steps:R;                      % Lower rotor radius distributions, m
rho = 1.058;                        % Density at max cruise altitude, kg/m^3
omega = 39.62;                      % rotational velocity at each point, ???
V_local = (omega.*r) + V_inf;       % local velocity at each point on the bottom blade, m/s
AR = 20;                            % aspect ratio
c_bar = R/AR;                       % Mean aerodynamic chord, m;
c = 0.265;                          % chord at each point along lower blade, m
twist = (-10/R).*r + 10;    % Twist distribution along blade, deg
AoA = 0 + twist;                    % Angle of attack at each location along blade;
mass = 8*9.81*R;                    % Approximate mass of lower blade

%%% Aerofoils
% RC(4)-10 - inboard
    % Cl = 0.1079*AoA+0.00768
% RC(6)-08 - outboard
    % Cl = 0.09223*AoA+0.0336
Cl_in = 0.1079*AoA(1)+0.00768;                  % Lift at root
Cl_out = 0.09223*AoA(n)+0.0336;                 % Lift at tip  
Cl_step = (Cl_out-Cl_in)/n;               % Step for linear twist therefore linear lift distribution
Cl = Cl_in:Cl_step:Cl_out;                   % lift coefficient distribution for the upper blade

%%% Finding the forces
lift = 0.5 * rho .* V_local.^2 .* (c.*steps) .* Cl; % lift along lower blade

weight = mass/(n+1);      % sectional weight of the lower blade, assuming 0 taper ratio

load_factor = 3.5;        % by CS-29, max is 3.5 min is -1.0

force_lower = load_factor*lift - weight;        % resultant forces

%%% Shear moment is found by summing all loads from that point until the
%%% tip of the blade
shear = [];
for i = 1: n+1
    d_shear = 0;
    for j = i:n+1
        d_shear = d_shear + force_lower(j);
    end
    shear(i) = d_shear;
end

%%% Bending moment is found from summing dM from the point on the blade to
%%% the tip, where dM is the shear at that point and the next 
d_bending = [];
for i = 1:n
    d_bending(i) = (shear(i) + shear(i+1))*steps/2;
end

bending = [];
for i = 1: n
    bend = 0;
    for j = i:n
        bend = bend + d_bending(j);
    end
    bending(i) = bend;
end
bending(n+1) = 0;


%%% MATERIAL PROPERTIES
%%% Mean values used, for ranges see Granta CES
% epoxy/HS carbon fibre UD prepreg
% PEEK/IM carbon fibre UD prepreg
% BMI/HS carbonn fibre UD prepreg
% phenolic/E-glass woven prepreg

E_epo = 140;        % Young's modulus, GPa
nu_epo = 0.33;      % Poisson ratio
G_epo = 5;          % Shear modulus, GPa
K_epo = 10;         % Bulk modulus, GPa
rho_epo = 1560;     % Density, kg/m^3
ten_y_epo = 1900;   % Tensile yield stress, MPa
shear_y_epo = 80;   % Shear yield stress, MPa, approx
com_y_epo = 1500;   % Compressive yield stress, MPa
cost_epo = 30;      % Price per kg, GBP

E_PEEK = 147;       % Young's modulus, GPa
nu_PEEK = 0.36;     % Poisson ratio
G_PEEK = 5.4;       % Shear modulus, GPa
K_PEEK = 15;        % Bulk modulus, GPa
rho_PEEK = 1560;    % Density, kg/m^3
ten_y_PEEK = 2420;  % Tensile yield stress, MPa
shear_y_PEEK = 95;  % Shear yield stress, MPa, approx
com_y_PEEK = 1100;  % Compressive yield stress, MPa
cost_PEEK = 85;     % Price per kg, GBP

E_BMI = 120;        % Young's modulus, GPa
nu_BMI = 0.32;      % Poisson ratio
G_BMI = 5.3;        % Shear modulus, GPa
K_BMI = 10;         % Bulk modulus, GPa
rho_BMI = 1590;     % Density, kg/m^3
ten_y_BMI = 1720;   % Tensile yield stress, MPa
shear_y_BMI = 95;   % Shear yield stress, MPa
com_y_BMI = 1250;   % Compressive yield stress, MPa
cost_BMI = 80;      % Price per kg, GBP

E_gla = 46;         % Young's modulus, GPa
nu_gla = 0.337;     % Poisson ratio
G_gla = 17;         % Shear modulus, GPa
% K_gla = 10;         % Bulk modulus, GPa
% rho_gla = 1590;     % Density, kg/m^3
% ten_y_gla = 1720;   % Tensile yield stress, MPa
shear_y_gla = 270;  % Shear yield stress, MPa, approx
% com_y_gla = 1250;   % Compressive yield stress, MPa
% cost_gla = 80;      % Price per kg, GBP


%%%     Wing box design from AVD
%%% Only takes skin into account
box_width = 0.35*c_bar;        % m, approx, can be varied
box_height = 0.09*c_bar;       % m, approx, can be varied
compress_load = bending./box_width./box_height;   % compressive load per unit length, N/m
skin_t_epo = (compress_load./3.62./(E_epo*10^9)*box_width.^2).^(1/3)*1000; % skin thickness, mm
crit_buckl_epo = compress_load./skin_t_epo;



%%% change the force into a perpendicular to the blade and a parallel to the blade components
%%% max bending moment and shear force at the root
% require sigma_c smaller than 1500 MPa = 1,500,000,000 N/m^2 = 1.5e9 N/m^2
M_y = bending(1)*cosd(twist(1));  % Acting parallel to y axis
M_x = bending(1)*sind(twist(1));  % Acting parallel to x axis

V_y = shear(1)*cosd(twist(1));    % Acting parallel to y axis
V_x = shear(1)*sind(twist(1));    % Acting parallel to x axis


%%%     Assuming blade as a rectangle, and taking components of force
t = 0.006;      % mm, skin thickness
t2 = 0.006;     % mm, spar thickness

h_bar = 0.0185;       % From the geometry, to get same cross-sectional area as the solid aerofoil

Ixx_rect = 1/12 * c_bar * h_bar^3 - 1/12 * (c_bar - 2*t) * (h_bar - 2*t)^3;
Ixx_spar_rect = 1/12 * (c_bar*0.35) * (c_bar*0.1)^3 - 1/12 * (c_bar*0.35 - 2*t2) * (c_bar*0.1 - 2*t2)^3;
Ixx_rect = Ixx_rect + Ixx_spar_rect;

sigma_t = -M_y * (h_bar/2) / Ixx_rect;  % maximum tensile stress in terms of t
sigma_c = -M_y * (h_bar/2) / Ixx_rect;  % maximum compressive stress in terms of t
% will be equal due to approximation of cross-section as a
% symmetric rectangle

% thickness 2 mm skin 4 mm spar gives 7.7e9 > 1.5e9
% thickness 5 mm gives sigma_c 4.6e9 > 1.5e9
% thickness 10 mm gives sigma_c 4.2e9 > 1.5e9


%%%     Assuming skin of aerofoil is modelled by a sixth/eigth degree
%%%     polynomial, can integrate to get Ixx
syms yu yl x y y_area y_in y_area_in yu_int yl_int
t = 0.004;      % m
t2 = 0.005;     % m
% allow neutral axis to be x-axis for now
yu = -5582*x^6 + 4864*x^5 - 1628*x^4 + 264*x^3 - 22.4*x^2 + 0.954*x + 0.002;
yl = 3.843*10^5*x^8 - 4.5*10^5*x^7 + 2.19*10^5*x^6 - 5.73*10^4*x^5 + ...
    8650*x^4 - 751*x^3 + 35.2*x^2 - 0.774*x - 0.0013;
yu_int = -5582*(((c_bar-2*t)/c_bar)*x)^6 + 4864*(((c_bar-2*t)/c_bar)*x)^5 - 1628*(((c_bar-2*t)/c_bar)*x)^4 + 264*(((c_bar-2*t)/c_bar)*x)^3 - 22.4*(((c_bar-2*t)/c_bar)*x)^2 + 0.954*(((c_bar-2*t)/c_bar)*x) + 0.002;
yl_int = 3.843*10^5*(((c_bar-2*t)/c_bar)*x)^8 - 4.5*10^5*(((c_bar-2*t)/c_bar)*x)^7 + 2.19*10^5*(((c_bar-2*t)/c_bar)*x)^6 - 5.73*10^4*(((c_bar-2*t)/c_bar)*x)^5 + 8650*(((c_bar-2*t)/c_bar)*x)^4 - 751*(((c_bar-2*t)/c_bar)*x)^3 + 35.2*(((c_bar-2*t)/c_bar)*x)^2 - 0.774*(((c_bar-2*t)/c_bar)*x) - 0.0013;

% https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/21cb09c67195d8b306de55214c2b2f87_spl10b.pdf
y = 1/3 * (yu^3 - yl^3);
y_in = 1/3 * (yu_int^3 - yl_int^3);

Ixx_int = double(int(y, 0, 0.265));
Ixx_int_in = double(int(y_in, 0, (c_bar-2*t)));

Ixx_int = Ixx_int - Ixx_int_in;

Ixx_spar = 1/12 * (c_bar*0.35) * (c_bar*0.1)^3 - 1/12 * (c_bar*0.35 - 2*t2) * (c_bar*0.1 - 2*t2)^3;
Ixx_int = Ixx_int + Ixx_spar;

sigma_t = -M_y * (h_bar/2) / Ixx_int;  % maximum tensile stress in terms of t
sigma_c = -M_y * (h_bar/2) / Ixx_int;  % maximum compressive stress in terms of t

% thickness = 1 mm gives sigma_c = 5.08e9 > 1.5e9
% thickness = 3 mm gives sigma_c = 2.1e9 > 1.5e9
% thickness = 4.5 mm gives sigma_c = 1.67e9 < 1.5e9
% thickness = 6 mm gives sigma_c = 1.67e9 < 1.5e9
% thickness = 4 spar = 5 gives sigma_c = 1.4e9

% boom-shear idealisation - turn the cross-sectional area a of the skin
% ~= thickness t * length l into n booms

% BROKEN CHANING N SHOULD NOT CHANGE THE ANSWER MUCH BUT IT DOES

%syms t A_tot A_boom Ixx r d
t = 0.003;
t2 = 0.004;
length = 0.540561;    % m
n = 10;     % must be even

% locate one boom at nose and one at tail
% i = 1 and i = n/2 are nose and tail respectively
% all others distributed evenly along c_bar
x_pos = [];
y_pos = [];
for i = 1:n
    if i == 1
        % located at nose
        x_pos(i) = 0;
        y_pos(i) = 0;
    elseif i < (n/2+1)
        % located on top surface
        x_pos(i) = (i-1)*(c_bar/5);
        y_pos(i) = -5582*x_pos(i)^6 + 4864*x_pos(i)^5 - 1628*x_pos(i)^4 + ...
            264*x_pos(i)^3 - 22.4*x_pos(i)^2 + 0.954*x_pos(i) + 0.002;
    elseif i == (n/2+1)
        % located at tail
        x_pos(i) = c_bar;
        y_pos(i) = 0;
    else
        % located on bottom surface
        x_pos(i) = c_bar - (i-6)*(c_bar/5);
        y_pos(i) = 3.843*10^5*x_pos(i)^8 - 4.5*10^5*x_pos(i)^7 + ...
            2.19*10^5*x_pos(i)^6 - 5.73*10^4*x_pos(i)^5 + 8650*x_pos(i)^4 -...
            751*x_pos(i)^3 + 35.2*x_pos(i)^2 - 0.774*x_pos(i) - 0.0013;
    end
end

% divide area up into n booms
A_tot = length * t;
A_boom = A_tot / n;
r = (A_boom/pi)^0.5;
d = 2*r;

% Iz = pi/64 d^4
Izz = (pi/64) * d^4;
Ixx_boom = 0;
for i = 1:n
    % Ix = Iz + A*y^2
    Ixx_boom = Ixx_boom + Izz + A_boom*(y_pos(i))^2;
end

% Add in spar
A_spar = 0.35*c_bar*0.1*c_bar - (0.35*c_bar-2*t2)*(0.1*c_bar-2*t2);

% model with 3 booms on top surface and 3 on lower (one in each corner and
% one in centre of each long surface)
x_spar = [0.15*c_bar, 0.325*c_bar, 0.5*c_bar, 0.5*c_bar, 0.325*c_bar, 0.15*c_bar];
y_spar = [0.065*c_bar, 0.065*c_bar, 0.065*c_bar, -0.035*c_bar, -0.035*c_bar, -0.035*c_bar];

A_spar_boom = A_spar/6;
r_spar = (A_spar_boom/pi)^0.5;
d_spar = 2*r_spar;

Izz_spar = (pi/64) * d_spar^4;
Ixx_spar = 0;
for i = 1:6
    Ixx_spar = Ixx_spar + Izz_spar + A_boom*(y_spar(i))^2;
end

Ixx_boom = Ixx_boom + Ixx_spar;

sigma_t = -M_y * (h_bar/2) / Ixx_boom;  % maximum tensile stress in terms of t
sigma_c = -M_y * (h_bar/2) / Ixx_boom;  % maximum compressive stress in terms of t

% thickness = 1 mm gives sigma_c = 7.3e9 > 1.5e9
% thickness = 3 mm gives sigma_c = 2.3e9 > 1.5e9
% thickness = 4.5 mm gives sigma_c = 1.45e9 < 1.5e9
% thickness = 3 spar = 5 gives sigma_c = 1.4e9
% thickness = 2.5 spar = 4 gives sigma_c = 1.7e9





