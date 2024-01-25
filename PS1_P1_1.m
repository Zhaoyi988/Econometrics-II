% 1.1 A Monte Carlo exercise
%% 1 
M = 5000;
T = 280;
Phi = 0.9;
y_mean = 3;
std = 3;

mu = y_mean * (1 - Phi);
gammasqr = (std^2) * (1 - Phi^2);
y = zeros(M, T);

rng(42) 

for i = 1:M
        epsilon = 0 + gammasqr * randn(1, T);
    for j = 1:T-1
        y(i,1) = y_mean + epsilon(1);
        y(i,j + 1) = mu + Phi * y(i, j) + epsilon(j + 1);
    end
end

%% 2
mu_hat = zeros(M, 1);
phi_hat = zeros(M, 1);

for i = 1:M
    X = [ones(T-1, 1) y(i, 1:T-1)'];
    y_t = y(i, 2:T)';
    beta_hat = X\y_t;
    mu_hat(i) = beta_hat(1);
    phi_hat(i) = beta_hat(2);
 
end

figure;
subplot(2, 1, 1);
histogram(mu_hat, 'FaceColor', [31/255, 153/255, 166/255], 'Normalization', 'probability', 'EdgeColor', 'w');
hold on;
line([mu mu], [0 0.25], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-');
title('\mu Estimates');
xlabel('\mu');
legend('OLS Estimates', 'True Value');

subplot(2, 1, 2);
histogram(phi_hat, 'FaceColor', [31/255, 153/255, 166/255], 'Normalization', 'probability', 'EdgeColor', 'w');
hold on;
line([Phi Phi], [0 0.25], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '-');
title('\Phi Estimates');
xlabel('\Phi');
legend('OLS Estimates', 'True Value');

saveas(gcf, '1.1_2.png');



%% 3


M = 5000;   
phi_values = [0.90, 0.95, 0.97, 0.99];  
T_values = [40, 80, 120, 280];     
rng(42);  

unconditional_mean = 3;
unconditional_std = 3;

results_table = table();

for phi_idx = 1:length(phi_values) % Iterate all phi_values
        Phi = phi_values(phi_idx);
    
    for T_idx = 1:length(T_values) % Iterate all T values
        T = T_values(T_idx);
      
        mu = y_mean * (1 - Phi);
        gammasqr = (std^2) * (1 - Phi^2);
        y = zeros(M, T);
        
        for i = 1:M
                epsilon = 0 + gammasqr * randn(1, T);
            for j = 1:T-1
                y(i,1) = y_mean + epsilon(1);
                y(i,j + 1) = mu + Phi * y(i, j) + epsilon(j + 1);
            end
        end % Create a matrix of y for each combination of Phi and T

        phi_hat = zeros(M, 1);
        
        for i = 1:M
            X = [ones(T-1, 1) y(i, 1:T-1)'];
            y_t = y(i, 2:T)';
            beta_hat = X\y_t;
            phi_hat(i) = beta_hat(2);
         
        end % For y calculate the mean estimated phi_hat
            k = mean(phi_hat);
            result_entry = table(Phi, T, k);
            results_table = [results_table; result_entry];
    end
end
disp(results_table);

latex_file = 'C:\Users\Circle\Dropbox\applications\Overleaf\Assignment1\Table.tex';
fid = fopen(latex_file, 'w');

fprintf(fid, '\\begin{table}\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{ccc}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'True $\\phi$ & Sample Size $T$ & Mean OLS Estimate of $\\phi$ \\\\ \n');
fprintf(fid, '\\hline\n');

for row_idx = 1:height(results_table)
    fprintf(fid, '%.2f & %d & %.4f \\\\ \n', results_table.Phi(row_idx), results_table.T(row_idx), results_table.k(row_idx));
end


fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Mean OLS Estimate of $\\phi$ for Different True $\\phi$ Values and Sample Sizes}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);

% 1.2 Forecasting Under Structural Breaks
%% 1
rng(8);  
T = 280;
phi = 0.9;
sigma_y = 3;
sigma_e = sqrt((sigma_y^2) * (1 - phi^2));
t_series = datetime(1948, 3, 31):calmonths(3):datetime(2017, 12, 31);

mu_series = zeros(T, 1); % Generate mu series
for t = 1:T
    if t_series(t) < datetime(1973, 1, 1)
        mu_series(t, 1) = 5 * (1 - phi);
    elseif t_series(t) < datetime(1996, 1, 1)
        mu_series(t, 1) = 0 * (1 - phi);
    elseif t_series(t) < datetime(2005, 1, 1)
        mu_series(t, 1) = 4.5 * (1 - phi);
    elseif t_series(t) < datetime(2018, 1, 1)
        mu_series(t, 1) = 0.5 * (1 - phi);
    end
end

%% 2
y = zeros(T, 1); % generate y series
epsilon = sigma_e * randn(T);
y(1) = mu_series(1)/(1 - phi) + epsilon(1);

for t = 2:T
    y(t) = mu_series(t) + phi * y(t - 1) + epsilon(t);
end

figure;
hold on;
legend;
plot(t_series, y, 'LineWidth', 1.5, 'DisplayName', 'y');
plot(t_series, mu_series / (1 - phi), 'DisplayName','Unconditional Mean');
xlabel('Year');
legend('show');
saveas(gcf, '1.2_1.png');

%% 3
end_est = sum(t_series < datetime(1990, 1, 1)); % The t index of the end of estimation period 1989Q4
start_pred = end_est + 1; % The t index of the start of prediction period 1990Q1
pred_period = 12;

for t = end_est:(T-1)

    % Estimate with all past data points
    X = y(1:(t-1), 1);
    Y = y(2:t, 1);
    model1 = fitlm(X, Y);

    % Predict
    ypred1 = predAR1(model1, y(t, 1), pred_period);
    if t == end_est
        Ypred1 = ypred1;
    else
        Ypred1 = [Ypred1 ypred1];
    end
end

%% 4
for t = end_est:(T-1)

    % Estimate with the past 40 data points
    X = y((t-40):(t-1), 1);
    Y = y((t-39):t, 1);
    model2 = fitlm(X, Y);

    % Predict
    ypred2 = predAR1(model2, y(t, 1), pred_period);
    if t == end_est
        Ypred2 = ypred2;
    else
        Ypred2 = [Ypred2 ypred2];
    end
end

%% 5
for t = end_est:(T-1)
   
    % Predict using the most recent value
    ypred3 = ones(pred_period, 1) * y(t, 1);
    if t == end_est
        Ypred3 = ypred3;
    else
        Ypred3 = [Ypred3 ypred3];
    end
end



