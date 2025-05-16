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
