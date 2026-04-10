nn = 10;

%% Run
[alpha_ml, beta_ml, alpha_aq, beta_aq, alpha_qm, beta_qm] = PlotNPA(3, nn);

%% Plot
figure; hold on;
plot(beta_ml,  alpha_ml,  'b-', 'LineWidth', 2, 'DisplayName', 'ML');
plot(beta_aq,  alpha_aq,  'r-', 'LineWidth', 2, 'DisplayName', 'AQ');
plot(beta_qm,  alpha_qm,  'p-', 'LineWidth', 2, 'DisplayName', 'QM');
xlabel('q_2', 'FontSize', 13);
ylabel('q_1', 'FontSize', 13);
title('Comparison of ML, AQ, QM regions', 'FontSize', 14);
legend('show', 'Location', 'northeast');
grid on;
hold off;


%% Save Results to File
% 1. Ensure all vectors are column vectors and combine them
% We use the ' (transpose) to flip them from rows to columns
data_matrix = [beta_ml',  alpha_ml', ...
               beta_aq',  alpha_aq', ...
               beta_qm',  alpha_qm'];

% 2. Create a Table with descriptive headers
colNames = {'beta_ML', 'alpha_ML', 'beta_AQ', 'alpha_AQ', ...
            'beta_QM', 'alpha_QM'};
T = array2table(data_matrix, 'VariableNames', colNames);

% 3. Write to a file
filename = 'slice_data.xlsx';
writetable(T, filename);