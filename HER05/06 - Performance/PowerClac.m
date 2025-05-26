function output = PowerClac(Weight,h,T,alp ,V_forward,input)


omega = input(1);
R = input(2);
A = input(3);
rho_0 = input(4);
kappa_ii = input(5);
kappa_int = input(6);
f = input(7);
Paux_total = input(8);
sigema = input(9);
Cd0 = input(10);
kappa_02 = input(11);
kappa_04 = input(12);
gamma = input(13);
R_0 = input(14);
t_c = input(15);
sweep = input(16);
M_dd = input(17);
kappa_trans = input(18);
ISA_factor = input(19);
T_trans = input(20);
mu = V_forward ./ (omega * R);

if ISA_factor == 1
    T1_T2 = ((15 - 0.001981 * h) + T_trans) / ((35 - 1.981 * h * 10^-3) + T_trans);
    T = T + 20;
else
    T1_T2 = 1;
end

% Pi clac
CT = Weight/(2*rho_0*(1-6.876*10^-6*h)^4.265*T1_T2*A*omega^2*R^2);
lambda_initial = sqrt(CT/2);
lambda_n = lambda_initial;
lambda = zeros(size(V_forward));
for j = 1:length(V_forward)
    for i = 1:100
        lambda_n1 = mu(j)*tan(alp(j)) + CT / (2*sqrt(mu(j)^2 + lambda_n^2));
        lambda_n = lambda_n1;
    end
    lambda(j) = lambda_n1;
end
kappa_i = kappa_ii .* cosh(7.5 .* mu.^2);
Pi = (kappa_int * Weight/2 * omega * R * 2) .* lambda .* kappa_i;


Cpp = (0.5 * f / A) .* mu.^3;
Pp = Cpp .* (rho_0*(1-6.876*10^-6*h)^4.265*T1_T2 * A * omega^3 * R^3);



Pc = 0.3 * Weight *ones(size(V_forward));


Paux = Paux_total *ones(size(V_forward));


% Po calc
Cpo_pure = sigema * Cd0 .* (1 + kappa_02.*mu.^2 + kappa_04.*mu.^4) ./ 8;
h_M = h;
M_19 = (omega * R + V_forward) ./ sqrt(gamma * R_0 * T);
M_hat = (M_19.^2 - 1) ./ (1.79 .* M_19.^(4/3) .* t_c^(2/3));
M_hover = (omega * R) / sqrt(gamma * R_0 * T);
delta_Cpo = zeros(size(V_forward));
for i = 1:length(V_forward)
    if M_19(i) >= M_dd
        delta_Cpo1 = sigema * (0.3 * (1+mu(i))^(5/2) * t_c^(5/2) * (M_hat(i)+1)^2);
    else
        delta_Cpo1 = 0;
    end
    % Assumption: tip relief only occurs at RC6-08 part (85%---100%)
    M_r = M_hover*(1+mu(i))*cos(sweep);
    if M_r >= M_dd
        psi1 = asin(mu(i)^-1 * (M_dd/M_hover - 1));
        psi2 = 180 - psi1;
        f = @(psi,r) (r + mu(i)*sin(psi))^3 * r * (0.165*M_hover*(r+mu(i)*sin(psi))*cos(sweep)-0.0895);
        delta_Cpo2 = sigema * integral2(f,psi1,psi2,@(psi) M_dd/M_hover - mu(i)*sin(psi),1) / (4 * pi);
    else
        delta_Cpo2 = 0;
    end
    delta_Cpo(i) = delta_Cpo2 + delta_Cpo1;
end
Po = (Cpo_pure+delta_Cpo) .* (2 * A * omega^3 * R^3 * rho_0*(1-6.876*10^-6*h_M)^4.265 *T1_T2);

PPP = (Pi + Po + Pp + Paux) .* kappa_trans;


output = PPP;
end