function [A,B,A_lon,B_lon,A_lat,B_lat] = Linearise(state,altitude,flightCondition,ctrl)
% Inputs:  state - state vector [u v w p q r phi theta psi]
%          altitude - ISA altitude [ft]
%          flightCondition - 1 - full fuel, full payload
%                            2 - zero fuel, full payload
%                            3 - full fuel, zero payload
%                            4 - 50% fuel, 50% payload
%          ctrl (optional) - controls vector [theta_0 Delta_theta theta_1c theta_1s delta_r delta_e]
% Outputs: A - state matrix (9x9)
%          B - Controls matrix (9x6)
%          A_lon - longitudinal state matrix (4x4) [u w q theta]'
%          B_lon - longitudinal controls matrix (4x3) [theta_0 theta_1s delta_e]'
%          A_lat - lateral state matrix (5x5) [v p r phi psi]'
%          B_lat - lateral controls matrix (5x3) [Delta_theta theta_1c delta_r]'

arguments
        state (1,9) double
        altitude (1,1) double
        flightCondition
        ctrl (1,6) double = [1e10 0 0 0 0 0]
end

g = 9.81;

% Extract state and controls
u = state(1);
v = state(2);
w = state(3);
p = state(4);
q = state(5);
r = state(6);
phi = state(7);
theta = state(8);
psi = state(9);

% Calculate trimmed controls if needed
if ctrl(1) == 1e10
    ctrl = Trim(state,altitude,flightCondition);
end

% Calculate forces at base condition
u_b = FM(state,ctrl,altitude,flightCondition);

% Perturbations
for ii = 1:6
    stateperturbations(ii) = max(1e-3,1e-3*state(ii));
    ctrlperturbations(ii) = max(1e-3,1e-3*ctrl(ii));
end
stateperturbation = diag(stateperturbations);
ctrlperturbation = diag(ctrlperturbations);

% Calculate forces at perturbed condition
for ii = 1:6
    statep(ii,:) = [state(1:6)+stateperturbation(ii,:) state(7:9)];
    u_pA(ii,:) = FM(statep(ii,:),ctrl,altitude,flightCondition);
    ctrlp(ii,:) = ctrl + ctrlperturbation(ii,:);
    u_pB(ii,:) = FM(state,ctrlp(ii,:),altitude,flightCondition);
end

% Calculate derivatives
for ii = 1:6
    A(:,ii) = (u_pA(:,ii) - u_b(ii)) / stateperturbations(ii);
    B(:,ii) = (u_pB(:,ii) - u_b(ii)) / ctrlperturbations(ii);
end

% Coupling of L and N
geometry = Parameters(flightCondition);
Ixx = geometry.flightCondition.Ixx;
Izz = geometry.flightCondition.Izz;
Ixz = geometry.flightCondition.Ixz;
Atemp = (A(4,:)-Ixz*A(6,:)./Izz)*(Ixx+Ixz^2/Izz)/Ixx;
Atemp(2,:) = (A(6,:)-Ixz*A(4,:)./Ixx)*(Izz+Ixz^2/Ixx)/Izz;
A(4,:) = Atemp(1,:);
A(6,:) = Atemp(2,:);

% Kinematic equations
A1 = [  0, -g*cos(theta), 0;
        g*cos(phi)*cos(theta), -g*sin(phi)*sin(theta), 0;
        -g*cos(theta)*sin(phi), -g*cos(phi)*sin(theta), 0;
        0, 0, 0;
        0, 0, 0;
        0, 0, 0];
A = [A A1;
    0 0 0 1 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 1 0 0 0];

% Fill out B with zeros
B = [B; zeros(3,6)];

% Longitudinal decomposition
lonA_idx = [1 3 5 8];
lonB_idx = [1 4 6];
A_lon = A(lonA_idx,lonA_idx);
B_lon = B(lonA_idx,lonB_idx);

% Lateral decomposition
latA_idx = [2 4 6 7 9];
latB_idx = [2 3 5];
A_lat = A(latA_idx,latA_idx);
B_lat = B(latA_idx,latB_idx);

end