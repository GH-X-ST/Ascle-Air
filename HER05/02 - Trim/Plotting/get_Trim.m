function [theta_LR, theta_UR, theta_1s, theta_1c, phi, theta, delta_E, delta_R] = get_Trim(u, w, config)
% Introduction:
%   This function computes the hover, vertical climb and forward flight
%   trim, be aware that this function only work for pure forward flight
%   or pure vertical climb
%
% Author:
%   Hanchen Li (hl3422@ic.ac.uk)
%
% Inputs:
%   u         - forward flight velocity
%   w         - vertical velocity
%   config    - 'MTOW', 'MZFW', 'ZeroP', 'LinearPoint'
%
% Outputs:
%   theta_LR  - Lower rotor collective (rad)
%   theta_UR  - Upper rotor collective (rad)
%   theta_1s  - Longitudinal cyclic (rad)
%   theta_1c  - Lateral cyclic (rad)
%   phi       - Euler roll angle (rad)
%   theta     - Euler pitch angle (rad)
%   delta_E   - Elevator deflection (rad)
%   delta_R   - Rudder deflection (rad)

%% 0 Input
% Manually input data from optimise

% 0.1 Data Point

% Trim Point
v_F_t = 72.0222222; % 140 knots = 72.0222222 m/s

% Vertical Climb
% from 0 to 17.5566 m/s
v_C = [0, 4, 8, 12, 17.5566];

% Forward Flight
% from 0 to 84.7 m/s
v_F = [0, 25, 50, v_F_t, 84.7];

% 0.2 MTOW
% Hover
theta_LR_H_1 = 0;
theta_UR_H_1 = 0;
theta_1s_H_1 = 0;
theta_1c_H_1 = 0;
phi_H_1      = 0;
theta_H_1    = 0;
delta_E_H_1  = 0;
delta_R_H_1  = 0;

% Trim Point
theta_LR_t_1 = 0.5;
theta_UR_t_1 = 0.5;
theta_1s_t_1 = 0.5;
theta_1c_t_1 = 0.5;
phi_t_1      = 0.5;
theta_t_1    = 0.5;
delta_E_t_1  = 0.5;
delta_R_t_1  = 0.5;

% Vertical Climb
theta_LR_C_1 = [theta_LR_H_1, 0.3, 0.1, 0.2, 0.4];
theta_UR_C_1 = [theta_UR_H_1, 0.1, 0.2, 0.3, 0.4];

% Forward Flight
theta_LR_F_1 = [theta_LR_H_1, 0.3, 0.2, theta_LR_t_1, 0.6];
theta_UR_F_1 = [theta_UR_H_1, 0.1, 0.2, theta_UR_t_1, 0.5];
theta_1s_F_1 = [theta_1s_H_1, 0.1, 0.2, theta_1s_t_1, 0.5];
theta_1c_F_1 = [theta_1c_H_1, 0.1, 0.2, theta_1c_t_1, 0.5];
phi_F_1      = [phi_H_1, 0.1, 0.2, phi_t_1, 0.5];
theta_F_1    = [theta_H_1, 0.1, 0.2, theta_t_1, 0.5];
delta_E_F_1  = [delta_E_H_1, 0.1, 0.2, delta_E_t_1, 0.5];
delta_R_F_1  = [delta_R_H_1, 0.1, 0.2, delta_R_t_1, 0.5];

% 0.3 MZFW
% Hover
theta_LR_H_2 = 0;
theta_UR_H_2 = 0;
theta_1s_H_2 = 0;
theta_1c_H_2 = 0;
phi_H_2      = 0;
theta_H_2    = 0;
delta_E_H_2  = 0;
delta_R_H_2  = 0;

% Trim Point
theta_LR_t_2 = 0.6;
theta_UR_t_2 = 0.6;
theta_1s_t_2 = 0.6;
theta_1c_t_2 = 0.6;
phi_t_2      = 0.6;
theta_t_2    = 0.6;
delta_E_t_2  = 0.6;
delta_R_t_2  = 0.6;

% Vertical Climb
theta_LR_C_2 = [theta_LR_H_2, 0.2, 0.3, 0.4, 0.5];
theta_UR_C_2 = [theta_UR_H_2, 0.2, 0.3, 0.4, 0.5];

% Forward Flight
theta_LR_F_2 = [theta_LR_H_2, 0.2, 0.3, theta_LR_t_2, 0.7];
theta_UR_F_2 = [theta_UR_H_2, 0.2, 0.3, theta_UR_t_2, 0.7];
theta_1s_F_2 = [theta_1s_H_2, 0.2, 0.3, theta_1s_t_2, 0.7];
theta_1c_F_2 = [theta_1c_H_2, 0.2, 0.3, theta_1c_t_2, 0.7];
phi_F_2      = [phi_H_2, 0.2, 0.3, phi_t_2, 0.7];
theta_F_2    = [theta_H_2, 0.2, 0.3, theta_t_2, 0.7];
delta_E_F_2  = [delta_E_H_2, 0.2, 0.3, delta_E_t_2, 0.7];
delta_R_F_2  = [delta_R_H_2, 0.2, 0.3, delta_R_t_2, 0.7];

% 0.4 Maximum Fuel with Zero Payload
% Hover
theta_LR_H_3 = 0;
theta_UR_H_3 = 0;
theta_1s_H_3 = 0;
theta_1c_H_3 = 0;
phi_H_3      = 0;
theta_H_3    = 0;
delta_E_H_3  = 0;
delta_R_H_3  = 0;

% Trim Point
theta_LR_t_3 = 0.7;
theta_UR_t_3 = 0.7;
theta_1s_t_3 = 0.7;
theta_1c_t_3 = 0.7;
phi_t_3      = 0.7;
theta_t_3    = 0.7;
delta_E_t_3  = 0.6;
delta_R_t_3  = 0.6;

% Vertical Climb
theta_LR_C_3 = [theta_LR_H_3, 0.3, 0.4, 0.5, 0.6];
theta_UR_C_3 = [theta_UR_H_3, 0.3, 0.4, 0.5, 0.6];

% Forward Flight
theta_LR_F_3 = [theta_LR_H_3, 0.3, 0.4, theta_LR_t_3, 0.7];
theta_UR_F_3 = [theta_UR_H_3, 0.3, 0.4, theta_UR_t_3, 0.7];
theta_1s_F_3 = [theta_1s_H_3, 0.3, 0.4, theta_1s_t_3, 0.7];
theta_1c_F_3 = [theta_1c_H_3, 0.3, 0.4, theta_1c_t_3, 0.7];
phi_F_3      = [phi_H_3, 0.3, 0.4, phi_t_3, 0.7];
theta_F_3    = [theta_H_3, 0.3, 0.4, theta_t_3, 0.7];
delta_E_F_3  = [delta_E_H_3, 0.3, 0.4, delta_E_t_3, 0.7];
delta_R_F_3  = [delta_R_H_3, 0.3, 0.4, delta_R_t_3, 0.7];

% 0.5 Linearization Point
% Hover
theta_LR_H_4 = 0;
theta_UR_H_4 = 0;
theta_1s_H_4 = 0;
theta_1c_H_4 = 0;
phi_H_4      = 0;
theta_H_4    = 0;
delta_E_H_4  = 0;
delta_R_H_4  = 0;

% Trim Point
theta_LR_t_4 = 0.7;
theta_UR_t_4 = 0.7;
theta_1s_t_4 = 0.7;
theta_1c_t_4 = 0.7;
phi_t_4      = 0.7;
theta_t_4    = 0.7;
delta_E_t_4  = 0.6;
delta_R_t_4  = 0.6;

% Vertical Climb
theta_LR_C_4 = [theta_LR_H_4, 0.4, 0.5, 0.6, 0.7];
theta_UR_C_4 = [theta_UR_H_4, 0.4, 0.5, 0.6, 0.7];

% Forward Flight
theta_LR_F_4 = [theta_LR_H_4, 0.4, 0.5, theta_LR_t_4, 0.8];
theta_UR_F_4 = [theta_UR_H_4, 0.4, 0.5, theta_UR_t_4, 0.8];
theta_1s_F_4 = [theta_1s_H_4, 0.4, 0.5, theta_1s_t_4, 0.8];
theta_1c_F_4 = [theta_1c_H_4, 0.4, 0.5, theta_1c_t_4, 0.8];
phi_F_4      = [phi_H_4, 0.4, 0.5, phi_t_4, 0.8];
theta_F_4    = [theta_H_4, 0.4, 0.5, theta_t_4, 0.8];
delta_E_F_4  = [delta_E_H_4, 0.4, 0.5, delta_E_t_4, 0.8];
delta_R_F_4  = [delta_R_H_4, 0.4, 0.5, delta_R_t_4, 0.8];


%% 1 Linearisation
% 1.1 Hover and vertical climb
v_C_fine = linspace(min(v_C), max(v_C), 500);

% theta_LR
% MTOW
theta_LR_C_1_interp = interp1(v_C, theta_LR_C_1, v_C_fine, 'pchip');

% MZFW
theta_LR_C_2_interp = interp1(v_C, theta_LR_C_2, v_C_fine, 'pchip');

% Maximum Fuel with Zero Payload
theta_LR_C_3_interp = interp1(v_C, theta_LR_C_3, v_C_fine, 'pchip');

% Linearization Point
theta_LR_C_4_interp = interp1(v_C, theta_LR_C_4, v_C_fine, 'pchip');

% theta_UR
% MTOW
theta_UR_C_1_interp = interp1(v_C, theta_UR_C_1, v_C_fine, 'pchip');

% MZFW
theta_UR_C_2_interp = interp1(v_C, theta_UR_C_2, v_C_fine, 'pchip');

% Maximum Fuel with Zero Payload
theta_UR_C_3_interp = interp1(v_C, theta_UR_C_3, v_C_fine, 'pchip');

% Linearization Point
theta_UR_C_4_interp = interp1(v_C, theta_UR_C_4, v_C_fine, 'pchip');

% 1.2 Hover and forward Flight
v_F_fine = linspace(min(v_F), max(v_F), 500);

% theta_LR
% MTOW
theta_LR_F_1_interp = interp1(v_F, theta_LR_F_1, v_F_fine, 'pchip');

% MZFW
theta_LR_F_2_interp = interp1(v_F, theta_LR_F_2, v_F_fine, 'pchip');

% Maximum Fuel with Zero Payload
theta_LR_F_3_interp = interp1(v_F, theta_LR_F_3, v_F_fine, 'pchip');

% Linearization Point
theta_LR_F_4_interp = interp1(v_F, theta_LR_F_4, v_F_fine, 'pchip');

% theta_UR
% MTOW
theta_UR_F_1_interp = interp1(v_F, theta_UR_F_1, v_F_fine, 'pchip');

% MZFW
theta_UR_F_2_interp = interp1(v_F, theta_UR_F_2, v_F_fine, 'pchip');

% Maximum Fuel with Zero Payload
theta_UR_F_3_interp = interp1(v_F, theta_UR_F_3, v_F_fine, 'pchip');

% Linearization Point
theta_UR_F_4_interp = interp1(v_F, theta_UR_F_4, v_F_fine, 'pchip');

% theta_1s
% MTOW
theta_1s_F_1_interp = interp1(v_F, theta_1s_F_1, v_F_fine, 'pchip');

% MZFW
theta_1s_F_2_interp = interp1(v_F, theta_1s_F_2, v_F_fine, 'pchip');

% Maximum Fuel with Zero Payload
theta_1s_F_3_interp = interp1(v_F, theta_1s_F_3, v_F_fine, 'pchip');

% Linearization Point
theta_1s_F_4_interp = interp1(v_F, theta_1s_F_4, v_F_fine, 'pchip');

% theta_1c
% MTOW
theta_1c_F_1_interp = interp1(v_F, theta_1c_F_1, v_F_fine, 'pchip');

% MZFW
theta_1c_F_2_interp = interp1(v_F, theta_1c_F_2, v_F_fine, 'pchip');

% Maximum Fuel with Zero Payload
theta_1c_F_3_interp = interp1(v_F, theta_1c_F_3, v_F_fine, 'pchip');

% Linearization Point
theta_1c_F_4_interp = interp1(v_F, theta_1c_F_4, v_F_fine, 'pchip');

% phi
% MTOW
phi_F_1_interp = interp1(v_F, phi_F_1, v_F_fine, 'pchip');

% MZFW
phi_F_2_interp = interp1(v_F, phi_F_2, v_F_fine, 'pchip');

% Maximum Fuel with Zero Payload
phi_F_3_interp = interp1(v_F, phi_F_3, v_F_fine, 'pchip');

% Linearization Point
phi_F_4_interp = interp1(v_F, phi_F_4, v_F_fine, 'pchip');

% theta
% MTOW
theta_F_1_interp = interp1(v_F, theta_F_1, v_F_fine, 'pchip');

% MZFW
theta_F_2_interp = interp1(v_F, theta_F_2, v_F_fine, 'pchip');

% Maximum Fuel with Zero Payload
theta_F_3_interp = interp1(v_F, theta_F_3, v_F_fine, 'pchip');

% Linearization Point
theta_F_4_interp = interp1(v_F, theta_F_4, v_F_fine, 'pchip');

% delta_E
% MTOW
delta_E_F_1_interp = interp1(v_F, delta_E_F_1, v_F_fine, 'pchip');

% MZFW
delta_E_F_2_interp = interp1(v_F, delta_E_F_2, v_F_fine, 'pchip');

% Maximum Fuel with Zero Payload
delta_E_F_3_interp = interp1(v_F, delta_E_F_3, v_F_fine, 'pchip');

% Linearization Point
delta_E_F_4_interp = interp1(v_F, delta_E_F_4, v_F_fine, 'pchip');

% delta_R
% MTOW
delta_R_F_1_interp = interp1(v_F, delta_R_F_1, v_F_fine, 'pchip');

% MZFW
delta_R_F_2_interp = interp1(v_F, delta_R_F_2, v_F_fine, 'pchip');

% Maximum Fuel with Zero Payload
delta_R_F_3_interp = interp1(v_F, delta_R_F_3, v_F_fine, 'pchip');

% Linearization Point
delta_R_F_4_interp = interp1(v_F, delta_R_F_4, v_F_fine, 'pchip');

%% 2 Output

if config == 'MTOW'

    if u == 0

        if w <= 0

            [~, idx]   = min(abs(v_C_fine + w));
            
            theta_LR = theta_LR_C_1_interp(idx);
            theta_UR = theta_UR_C_1_interp(idx);
            theta_1s = 0;
            theta_1c = 0;
            phi      = 0;
            theta    = 0;
            delta_E  = 0;
            delta_R  = 0;

        else

            error('BadInput', ...
              ['Input must be climb,' ...
              ' got vertical descend %g.'], -w);

        end

    else

        if w == 0

            [~, idx]   = min(abs(v_F_fine - u));

            theta_LR = theta_LR_F_1_interp(idx);
            theta_UR = theta_UR_F_1_interp(idx);
            theta_1s = theta_1s_F_1_interp(idx);
            theta_1c = theta_1c_F_1_interp(idx);
            phi      = phi_F_1_interp(idx);
            theta    = theta_F_1_interp(idx);
            delta_E  = delta_E_F_1_interp(idx);
            delta_R  = delta_R_F_1_interp(idx);

        else

             error('BadInput', ...
              ['Input must be pure forward flight,' ...
              ' got vertical speed %g.'], w);

        end

    end

    
elseif config == 'MZFW'

    if u == 0

        if w <= 0

            [~, idx]   = min(abs(v_C_fine + w));
            
            theta_LR = theta_LR_C_2_interp(idx);
            theta_UR = theta_UR_C_2_interp(idx);
            theta_1s = 0;
            theta_1c = 0;
            phi      = 0;
            theta    = 0;
            delta_E  = 0;
            delta_R  = 0;

        else

            error('BadInput', ...
              ['Input must be climb,' ...
              ' got vertical descend %g.'], -w);

        end

    else

        if w == 0

            [~, idx]   = min(abs(v_F_fine - u));

            theta_LR = theta_LR_F_2_interp(idx);
            theta_UR = theta_UR_F_2_interp(idx);
            theta_1s = theta_1s_F_2_interp(idx);
            theta_1c = theta_1c_F_2_interp(idx);
            phi      = phi_F_2_interp(idx);
            theta    = theta_F_2_interp(idx);
            delta_E  = delta_E_F_2_interp(idx);
            delta_R  = delta_R_F_2_interp(idx);

        else

             error('BadInput', ...
              ['Input must be pure forward flight,' ...
              ' got vertical speed %g.'], w);

        end

    end
    
elseif config == 'ZeroP'
    
    if u == 0

        if w <= 0

            [~, idx]   = min(abs(v_C_fine + w));
            
            theta_LR = theta_LR_C_3_interp(idx);
            theta_UR = theta_UR_C_3_interp(idx);
            theta_1s = 0;
            theta_1c = 0;
            phi      = 0;
            theta    = 0;
            delta_E  = 0;
            delta_R  = 0;

        else

            error('BadInput', ...
              ['Input must be climb,' ...
              ' got vertical descend %g.'], -w);

        end

    else

        if w == 0

            [~, idx]   = min(abs(v_F_fine - u));

            theta_LR = theta_LR_F_3_interp(idx);
            theta_UR = theta_UR_F_3_interp(idx);
            theta_1s = theta_1s_F_3_interp(idx);
            theta_1c = theta_1c_F_3_interp(idx);
            phi      = phi_F_3_interp(idx);
            theta    = theta_F_3_interp(idx);
            delta_E  = delta_E_F_3_interp(idx);
            delta_R  = delta_R_F_3_interp(idx);

        else

             error('BadInput', ...
              ['Input must be pure forward flight,' ...
              ' got vertical speed %g.'], w);

        end

    end

elseif config == 'LinearPoint'
    
    if u == 0

        if w <= 0

            [~, idx]   = min(abs(v_C_fine + w));
            
            theta_LR = theta_LR_C_4_interp(idx);
            theta_UR = theta_UR_C_4_interp(idx);
            theta_1s = 0;
            theta_1c = 0;
            phi      = 0;
            theta    = 0;
            delta_E  = 0;
            delta_R  = 0;

        else

            error('BadInput', ...
              ['Input must be climb,' ...
              ' got vertical descend %g.'], -w);

        end

    else

        if w == 0

            [~, idx]   = min(abs(v_F_fine - u));

            theta_LR = theta_LR_F_4_interp(idx);
            theta_UR = theta_UR_F_4_interp(idx);
            theta_1s = theta_1s_F_4_interp(idx);
            theta_1c = theta_1c_F_4_interp(idx);
            phi      = phi_F_4_interp(idx);
            theta    = theta_F_4_interp(idx);
            delta_E  = delta_E_F_4_interp(idx);
            delta_R  = delta_R_F_4_interp(idx);

        else

             error('BadInput', ...
              ['Input must be pure forward flight,' ...
              ' got vertical speed %g.'], w);

        end

    end

else
    
    error('BadInput', 'Wrong config input');

end

end