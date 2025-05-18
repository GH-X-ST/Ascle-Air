%% Pie Chart of Funding Sources
categories = {'Donations','Legacies','Corporate','NHS/Grants'};
values = [80, 15, 3, 2]; % Example percentages for UK air ambulance charities
figure;
p = pie(values);
title('Air Ambulance Charity Funding Sources');

% Remove text labels from pie slices
for i = 1:length(p)
    if isprop(p(i), 'String')
        p(i).String = '';
    end
end

% Add legend instead of overlapping labels
legend(categories, 'Location', 'southoutside', 'Orientation', 'horizontal');

legend(["Donations", "Legacies", "Corporate", "NHS/Grants"], "Location", "none", "Orientation", "horizontal", "Position", [0.0068 0.0376 0.9750, 0.0559])

%% Radar Chart: HER25 vs Competitors
metrics = {'Noise','Emissions','Maintenance','Range','Modularity','Price'};
HER25  = [9 9 8 7 9 9];
AW169  = [7 6 7 9 8 6];
H145   = [8 6 8 6 7 5];
Bell429 = [6 6 7 7 7 8];

% Close the loop
HER25(end+1)  = HER25(1);
AW169(end+1)  = AW169(1);
H145(end+1)   = H145(1);
Bell429(end+1)= Bell429(1);
metrics{end+1} = metrics{1}; % Add first label again to close loop

theta = linspace(0, 2*pi, length(HER25));
figure;
polaraxes;
hold on;
polarplot(theta, HER25, '-r', 'LineWidth', 2);
polarplot(theta, AW169, '-y', 'LineWidth', 2);
polarplot(theta, H145, '-m', 'LineWidth', 2);
polarplot(theta, Bell429, '-b', 'LineWidth', 2);
thetaticks(0:60:300);
thetaticklabels(metrics);
rticks([2 4 6 8]);
title('HER25 vs Competitors (Qualitative)');
legend('HER25','AW169','H145','Bell 429','Location','southoutside', 'Orientation', 'horizontal');
hold off;

legend(["HER25", "AW169", "H145", "Bell 429"], "Location", "none", "Orientation", "horizontal", "Position", [0.0982 0.0148 0.8120, 0.0559])

%% Line Chart of UK HEMS Missions
years = 2015:2024;
missions = [15000, 15500, 16000, 17000, 18000, 16000, 17000, 18000, 19500, 20000];
figure;
plot(years, missions, '-o','LineWidth',2);
xlabel('Year'); ylabel('Annual HEMS Missions');
title('UK Air Ambulance Missions (Simulated)');
grid on;

%% Break-even Bar Chart
variants = categorical({'Lite','Base','Plus'});
breakEvenUnits = [20, 30, 35]; % Example values
figure;
bar(variants, breakEvenUnits, 0.6);
ylabel('Units to Break-even');
title('Break-even Units for HER25 Variants');




% Radar chart for HER25 vs competitors (Qualitative Comparison)
attributes = {'Speed','Range','Payload','Safety','Cost','Versatility'};
% Qualitative scores (1-10) for HER25, AW109, H145, Bell429
data = [8, 7, 8, 9, 6, 8;   % HER25 – strong safety (twin+coax), high versatility, moderate cost
        9, 6, 6, 8, 5, 7;   % AW109 – very fast, but smaller payload, higher cost
        7, 8, 7, 9, 6, 8;   % H145 – good range/payload, very safe (modern twin), but more expensive
        8, 7, 6, 8, 5, 6];  % Bell 429 – decent speed, mid range, cost similar to AW109
% Close the radar data loop by repeating first value at end for plotting
P = size(attributes,2);
theta = linspace(0, 2*pi, P+1);
data_closed = [data, data(:,1)];  % append first column at end for each row
% Plot
figure;
polaraxes; hold on;
competitors = {'HER25','AW109','Airbus H145','Bell 429'};
colors = lines(size(data,1));
for i = 1:size(data,1)
    polarplot(theta, data_closed(i,:), '-o', 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', competitors{i});
end
thetalim([0 360]); % set 0 at right
% Set radial axis limits and labels
rlim([0 10]); rticks([0 5 10]); % scores from 0 to 10
thetaticks(theta(1:end-1)*180/pi); thetaticklabels(attributes);
legend('Location','southoutside');
title('HER25 vs Competitors – Attribute Radar Chart');



% Regional demand potential by segment (relative index values)
regions = {'Europe','Asia-Pac','Americas'};
segments = {'EMS','Law Enf','VIP'};
% Hypothetical demand indices (could represent fleet growth or investment potential)
demand = [90, 70, 50;   % EMS: Europe high (established), APAC growing, Americas large but mature
          50, 80, 60;   % Law Enforcement: Europe moderate, APAC high need, Americas also high
          40, 60, 80];  % VIP/Corporate: Europe moderate, APAC growing wealthy class, Americas highest
figure;
bar(demand');  % transpose so that each region becomes a group on x-axis
set(gca, 'XTickLabel', regions);
legend(segments, 'Location','northoutside','Orientation','horizontal');
ylabel('Relative Demand Index');
title('Regional Demand Potential by Segment');



% Price vs Range scatter for HER25 and competitors
models = {'HER25','AW109','Airbus H145','Bell 429','Ka-226','Hill HX50','Jaunt eVTOL'};
range_km = [650, 830, 650, 760, 600, 1296, 160];  % representative ranges (some ferry range)
price_usd_m = [7, 6.5, 9, 7.5, 5.5, 1.0, 2.5];     % price in USD millions (approx estimates)
figure;
scatter(range_km, price_usd_m, 100, 'filled'); grid on; hold on;
text(range_km+5, price_usd_m, models, 'FontSize',8);  % label points (offset slightly)
xlabel('Range (km)'); ylabel('Price (USD millions)');
title('Price vs Range: HER25 and Competing Aircraft');
