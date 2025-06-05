% Inputs (example values for HER25)
MTOW = 3500;       % kg
MaxSpeed = 259.28;    % km/h
RotorDiam = 10.44;    % m




% Coefficients from literature (example only — adjust with calibration)


a = 1e-4;           % scaling factor
b1 = 1.1;           % MTOW exponent
b2 = 0.3;           % Max speed exponent
b3 = 0.2;           % Rotor diameter exponent






% Cost estimation
cer_cost = a * (MTOW^b1) * (MaxSpeed^b2) * (RotorDiam^b3);
fprintf('Estimated CER Cost: $%.2f million\n', cer_cost);




% ===================================================================
% COMPETITOR-BASED COSTING MODEL (Market Benchmarking)
% ===================================================================




% Competitor data (MTOW in kg, MaxSpeed in km/h, Cost in $ million)
comp_data = [2980, 259, 5.5;   % Airbus H135
           3400, 287, 6.2;   % Bell 429
           3175, 285, 6.5;   % Leonardo AW109
           3000, 270, 5.9
           3800, 259, 9.7];  % Additional sample




% HER25 design inputs
HER25_MTOW = 3500;
HER25_Speed = 259.28;




% Extract inputs
X_comp = comp_data(:, 1:2);    % MTOW and Speed
Y_comp = comp_data(:, 3);      % Cost




% Add intercept term
X_comp_aug = [ones(size(X_comp,1),1), X_comp];




% Linear regression model
beta_comp = (X_comp_aug' * X_comp_aug) \ (X_comp_aug' * Y_comp);




% Predict HER25 cost
HER25_pred_cost = [1, HER25_MTOW, HER25_Speed] * beta_comp;




fprintf('HER25 Estimated Market Cost (based on competitors): $%.2f million\n', HER25_pred_cost);




% Plot
figure;
scatter(X_comp(:,1), Y_comp, 80, 'filled');
hold on;
scatter(HER25_MTOW, HER25_pred_cost, 100, 'r', 'filled');
xlabel('MTOW (kg)'); ylabel('Cost ($ million)');
title('HER25 Cost vs Competitors (Market-Based Estimation)');
legend('Competitors', 'HER25 Estimate');
grid on;








% ===================================================================
% NASA/ODU PROCESS-BASED COST MODEL (First-order process velocity model)
% ===================================================================
% Placeholder values provided – update with actual parameters




% Inputs (Example placeholders)
V0 = 5;               % Steady-state process velocity (parts/hr)
tau = 0.5;            % Time constant for ramp-up (hr)
t_process = 2;        % Actual process time per part (hr)
cost_modulus = 1.2;   % Relative complexity factor
labor_rate = 50;      % $/hr
machine_rate = 150;   % $/hr
num_parts = 50;




% Dynamic velocity model cost per part
effective_velocity = V0 * (1 - exp(-t_process / tau));
time_per_part = 1 / effective_velocity;
unit_cost = cost_modulus * (labor_rate + machine_rate) * time_per_part;
total_cost = unit_cost * num_parts;




fprintf('Process-Based Total Cost: $%.2f\n', total_cost);








% ===================================================================
% ACTIVITY-BASED COSTING (ABC) MODEL
% ===================================================================




% Inputs: Each activity and its consumption by the product
activities = {'Machining', 'Layup', 'Assembly'};
time_hours = [12, 8, 10];       % time per activity (hrs)
cost_rates = [120, 100, 90];    % cost per hour for each activity ($/hr)




% Total cost calculation
abc_cost = sum(time_hours .* cost_rates);




fprintf('ABC Estimated Cost: $%.2f\n', abc_cost);





% ===================================================================
% BREAKEVEN ANALYSIS (with dynamic sales forecast)
% ===================================================================

% Inputs
fixed_cost = 180;           % Development + tooling ($ million)
unit_cost = 7.0;            % Production cost per unit ($ million)
sell_price = 8.1;           % Selling price per unit ($ million)

% Sales model parameters (logistic growth)
K = 50;                     % Max sales capacity (units/year)
P0 = 10;                    % Initial year sales
r = 0.35;                   % Growth rate
t_max = 15;                 % Simulate over 15 years

% Simulate sales over time using logistic growth
sales_forecast = zeros(1, t_max);
for t = 1:t_max
    sales_forecast(t) = K / (1 + ((K - P0) / P0) * exp(-r * t));
end

% Calculate cumulative units and profit
cumulative_units = 0;
cumulative_margin = 0;
margin_per_unit = sell_price - unit_cost;
breakeven_units = fixed_cost / margin_per_unit;
breakeven_year = 0;

for t = 1:t_max
    yearly_units = sales_forecast(t);
    cumulative_units = cumulative_units + yearly_units;
    cumulative_margin = cumulative_margin + yearly_units * margin_per_unit;
    if cumulative_margin >= fixed_cost && breakeven_year == 0
        breakeven_year = t;
    end
end

fprintf('Breakeven Units (Logistic Model): %.0f\n', breakeven_units);
fprintf('Breakeven Achieved in Year: %d\n', breakeven_year);

% Plot 1: Breakeven Sales Forecast Over Time
figure;
plot(1:t_max, sales_forecast, '-o', 'LineWidth', 2);
xlabel('Year');
ylabel('Units Sold');
title('HER25 Projected Sales per Year (Logistic Growth)');
grid on;

% ===================================================================
% LEARNING CURVE COST ESTIMATION
% ===================================================================

C1 = 7;                    % Cost of first unit ($ million)
n = 136;                   % Unit number to evaluate
learning_rate = 0.96;

b = log2(learning_rate);
Cn = C1 * n^b;

fprintf('Cost of unit %d with learning rate %.2f: $%.2f million\n', n, learning_rate, Cn);

% Plot 2: Learning Curve
units = 1:136;
learning_curve = C1 * units.^b;
figure;
plot(units, learning_curve, 'r-', 'LineWidth', 2);
xlabel('Unit Number');
ylabel('Unit Cost ($ million)');
title('Learning Curve: Cost per Unit vs Unit Number');
grid on;

% ===================================================================
% LIFE CYCLE COSTING (LCC)
% ===================================================================

units_total = 136;
C_dev = 180;              
C_prod = unit_cost * units_total;   
C_ops = 2.0 * units_total;      
C_maint = 1.0 * units_total;    
C_disposal = 0.2 * units_total; 

LCC = C_dev + C_prod + C_ops + C_maint + C_disposal;
fprintf('Estimated Life Cycle Cost (LCC): $%.2f million\n', LCC);

% Plot 3: Life Cycle Cost Breakdown
figure;
bar([C_dev, C_prod, C_ops, C_maint, C_disposal], 'FaceColor', [0.2 0.6 0.5]);
set(gca, 'xticklabel', {'Dev', 'Prod', 'Ops', 'Maint', 'Disposal'});
ylabel('Cost ($ million)');
title('Life Cycle Cost Breakdown');
grid on;

% ===================================================================
% DEVELOPMENT COST PHASING WITH INFLATION
% ===================================================================

base_cost = 180;       
inflation = 0.03;      
years = 8;

phase_dist = [0.20, 0.18, 0.12, 0.10, 0.10, 0.08, 0.12, 0.10]; 

year_costs = zeros(1, years);
for i = 1:years
   year_costs(i) = base_cost * phase_dist(i) * (1 + inflation)^(i-1);
end

total_inflated = sum(year_costs);
fprintf('Total Inflated Cost over %d years: $%.2f million\n', years, total_inflated);

