% Emmanouil Savvakis, 9093
% Ioannis Gkouzoumas, 9057

clc
clear

% The starting and ending point of First Wave of Covid_19 for every country
% was produced using: https://www.worldometers.info. Some countries had
% either NaN or negative values, which were restored.

% Reading data from the given files.
Input1 = readtable('Covid19Confirmed.xlsx', 'basic', true);
Input2 = readtable('Covid19Deaths.xlsx', 'basic', true);

% 'A' country and 10 other European countries.
% A. Cases and Deaths of Austria. (10/3/2020 - 20/5/2020)
Austria = table2array([Input1(8, 73:144); Input2(8, 73:144)]);

% Fixing max deaths value (it's not on 22/4/2020, it's on 8/4).
% By changing this cell, max value falls on 8/4.
Austria(2, 44) = 19;

% 1. Cases and Deaths of Finland. (12/3/2020 - 3/6/2020)
Finland = table2array([Input1(47, 75:158); Input2(47, 75:158)]);

% 2. Cases and Deaths of Spain. (3/3/2020 - 2/5/2020)
Spain = table2array([Input1(130, 66:126); Input2(130, 66:126)]);

% Fixing negative point. (24/5/2020)
Spain(1, 47) = round(mean([Spain(1, 44:46), Spain(1, 48:50)]));

% Fixing max cases value (20/4/2020) and value of the cell previously
% listed as max (26/4).
Spain(1, 18) = 10656;
Spain(1, 24) = 6353;

% 3. Cases and Deaths of Germany. (12/3/2020 - 16/5/2020)
Germany = table2array([Input1(52, 75:140); Input2(52, 75:140)]);

% Fixing max deaths value (8/4/2020).
Germany(2, 28) = 333;

% 4. Cases and Deaths of Portugal. (18/3/2020 - 10/5/2020)
Portugal = table2array([Input1(113, 81:134); Input2(113, 81:134)]);

% Fixing negative point. (2/5/2020)
Portugal(1, 46) = round(mean([Portugal(1, 43:45), Portugal(1, 47:49)]));

% Fixing max cases value (10/4/2020).
Portugal(1, 24) = 1726;

% Fixing max deaths value (3/4/2020).
% Also fixing 24/4 and 3/5.
Portugal(2, 17) = 37;
Portugal(2, 38) = 34;
Portugal(2, 47) = 20;

% 5. Cases and Deaths of France. (16/3/2020 - 28/4/2020)
France = table2array([Input1(48, 79:124); Input2(48, 79:124)]);

% Fixing max deaths value (15/4/2020).
% Also fixing 3/4/2020
France(2, 31) = 1437;
France(2, 19) = 1119;

% Multiple Linear Regression for every country, using as independent variables
% the shifted data: x(t), x(t-1), ...., x(t-20). Also, find the optimal
% selection among them to fit the data best (Selection is made from stepwiselm function).
% These Models and the Linear Regression Model of Exercise 5 are compared
% using Mean Squared Error.
% Col.1: Model from Ex.5,   Col.2: Multiple Lin. Reg. Model,    Col.3:
% Stepwise Model.
mse_A = zeros(1, 3);
mse_1 = zeros(1, 3);
mse_2 = zeros(1, 3);
mse_3 = zeros(1, 3);
mse_4 = zeros(1, 3);
mse_5 = zeros(1, 3);

% Offset that creates the best Linear Regression Model for every country,
% as shown from Exercise 5.
offset = [13 17 7 13 5 7];

% Confidence level and critical value (diagnostic plot).
alpha = 0.05;
zcrit = norminv(1-alpha/2);

%% A. Austria.
% Independent variables for Austria's Multiple Linear Regression Model.
% x(t), x(t-1), ...., x(t-20).
x = zeros( size(Austria, 2), 21);
y = Austria(2, :)';

% Residuals of the 3 Models.
e_A = zeros(length(y), 3);

% Adjusted R^2 of the 3 Models.
adjR2_A = zeros(1, 3);

% MSE of Linear Regressiong Model from Ex.5, Multiple Linear Regression
% Model with 21 independent var. and stepwise Model.
for t = 0:20
        
    x(:, t+1) = offset_data_fun(Austria(1, :), t);
    
end

% Optimal Regression Model (Ex. 5).
X = [ones(length(x), 1) x(:, offset(1))];
beta_A = regress(y, X);

y_pred_A = X * beta_A;
mse_A(1) = immse(y_pred_A, y);

% Residuals.
e_A(:, 1) = y_pred_A - y;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
adjR2_A(1) =  1 - ((n-1)/(n-2)) * (sum(e_A(:, 1).^2)) / (sum((y - mu_y).^2));

% Regression Model including all 21 independent variables.
X = [ones(length(x), 1) x];
beta_A_ml = regress(y, X);

% Predictive y from Multiple Lin. Regression.
y_pred_ml_A = X * beta_A_ml;
mse_A(2) = immse(y_pred_ml_A, y);

% Residuals.
e_A(:, 2) = y_pred_ml_A - y;

% Adjusted R^2.
adjR2_A(2) =  1 - ((n-1)/(n-2)) * (sum(e_A(:, 2).^2)) / (sum((y - mu_y).^2));

% Stepwise Model (Using the most significant independent variables.)
[beta_A_step, ~, ~, mdl_A, stats_A] = stepwisefit(x, y, 'Display', 'off');
b0 = stats_A.intercept;

% Keep only significant variables.
i = find(mdl_A == 1);
x_step = x(:, i);

y_pred_step_A = [ones(length(x_step), 1) x_step] * ([b0; beta_A_step(i)]);
mse_A(3) = immse(y_pred_step_A, y);

% Residuals.
e_A(:, 3) = y_pred_step_A - y;

% Adjusted R^2.
adjR2_A(3) =  1 - ((n-1)/(n-2)) * (sum(e_A(:, 3).^2)) / (sum((y - mu_y).^2));

% Standarised Error.
std_e_A = std(e_A);
stand_e_A = e_A ./ std_e_A;


% Diagnostic Plots for each Model.
figure('NumberTitle', 'off', 'Name', 'A. Austria');
clf

% Ex.5 Model.
subplot(1, 3, 1); 
scatter(y_pred_A, stand_e_A(:, 1))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', ['Linear Model, Offset: ', num2str(offset(1))]})
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Multiple Linear Regression Model.
subplot(1, 3, 2);
scatter(y_pred_ml_A, stand_e_A(:, 2))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Multiple var. (21) Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Stepwise Regression Model.
subplot(1, 3, 3);
scatter(y_pred_step_A, stand_e_A(:, 3))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Stepwise Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')


%% 1. Finland.
% Independent variables for Finland's Multiple Linear Regression Model.
% x(t), x(t-1), ...., x(t-20).
x = zeros( size(Finland, 2), 21);
y = Finland(2, :)';

% Residuals of the 3 Models.
e_1 = zeros(length(y), 3);

% Adjusted R^2 of the 3 Models.
adjR2_1 = zeros(1, 3);

% MSE of Linear Regressiong Model from Ex.5, Multiple Linear Regression
% Model with 21 independent var. and stepwise Model.
for t = 0:20
        
    x(:, t+1) = offset_data_fun(Finland(1, :), t);
    
end

% Optimal Regression Model (Ex. 5).
X = [ones(length(x), 1) x(:, offset(2))];
beta_1 = regress(y, X);

y_pred_1 = X * beta_1;
mse_1(1) = immse(y_pred_1, y);

% Residuals.
e_1(:, 1) = y_pred_1 - y;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
adjR2_1(1) =  1 - ((n-1)/(n-2)) * (sum(e_1(:, 1).^2)) / (sum((y - mu_y).^2));

% Regression Model including all 21 independent variables.
X = [ones(length(x), 1) x];
beta_1_ml = regress(y, X);

% Predictive y from Multiple Lin. Regression.
y_pred_ml_1 = X * beta_1_ml;
mse_1(2) = immse(y_pred_ml_1, y);

% Residuals.
e_1(:, 2) = y_pred_ml_1 - y;

% Adjusted R^2.
adjR2_1(2) =  1 - ((n-1)/(n-2)) * (sum(e_1(:, 2).^2)) / (sum((y - mu_y).^2));

% Stepwise Model (Using the most significant independent variables.)
[beta_1_step, ~, ~, mdl_1, stats_1] = stepwisefit(x, y, 'Display', 'off');
b0 = stats_1.intercept;

% Keep only significant variables.
i = find(mdl_1 == 1);
x_step = x(:, i);

y_pred_step_1 = [ones(length(x_step), 1) x_step] * ([b0; beta_1_step(i)]);
mse_1(3) = immse(y_pred_step_1, y);

% Residuals.
e_1(:, 3) = y_pred_step_1 - y;

% Adjusted R^2.
adjR2_1(3) =  1 - ((n-1)/(n-2)) * (sum(e_1(:, 3).^2)) / (sum((y - mu_y).^2));

% Standarised Error.
std_e_1 = std(e_1);
stand_e_1 = e_1 ./ std_e_1;

% Diagnostic Plots for each Model.
figure('NumberTitle', 'off', 'Name', '1. Finland');
clf

% Ex.5 Model.
subplot(1, 3, 1); 
scatter(y_pred_1, stand_e_1(:, 1))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', ['Linear Model, Offset: ', num2str(offset(2))]})
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Multiple Linear Regression Model.
subplot(1, 3, 2);
scatter(y_pred_ml_1, stand_e_1(:, 2))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Multiple var. (21) Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Stepwise Regression Model.
subplot(1, 3, 3);
scatter(y_pred_step_1, stand_e_1(:, 3))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Stepwise Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

%% 2. Spain
% Independent variables for Spain's Multiple Linear Regression Model.
% x(t), x(t-1), ...., x(t-20).
x = zeros( size(Spain, 2), 21);
y = Spain(2, :)';

% Residuals of the 3 Models.
e_2 = zeros(length(y), 3);

% Adjusted R^2 of the 3 Models.
adjR2_2 = zeros(1, 3);

% MSE of Linear Regressiong Model from Ex.5, Multiple Linear Regression
% Model with 21 independent var. and stepwise Model.
for t = 0:20
        
    x(:, t+1) = offset_data_fun(Spain(1, :), t);
    
end

% Optimal Regression Model (Ex. 5).
X = [ones(length(x), 1) x(:, offset(3))];
beta_2 = regress(y, X);

y_pred_2 = X * beta_2;
mse_2(1, 1) = immse(y_pred_2, y);

% Residuals.
e_2(:, 1) = y_pred_2 - y;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
adjR2_2(1) =  1 - ((n-1)/(n-2)) * (sum(e_2(:, 1).^2)) / (sum((y - mu_y).^2));

% Regression Model including all 21 independent variables.
X = [ones(length(x), 1) x];
beta_2_ml = regress(y, X);

% Predictive y from Multiple Lin. Regression.
y_pred_ml_2 = X * beta_2_ml;
mse_2(1, 2) = immse(y_pred_ml_2, y);

% Residuals.
e_2(:, 2) = y_pred_ml_2 - y;

% Adjusted R^2.
adjR2_2(2) =  1 - ((n-1)/(n-2)) * (sum(e_2(:, 2).^2)) / (sum((y - mu_y).^2));

% Stepwise Model (Using the most significant independent variables.)
[beta_2_step, ~, ~, mdl_2, stats_2] = stepwisefit(x, y, 'Display', 'off');
b0 = stats_2.intercept;

% Keep only significant variables.
i = find(mdl_2 == 1);
x_step = x(:, i);

y_pred_step_2 = [ones(length(x_step), 1) x_step] * ([b0; beta_2_step(i)]);
mse_2(1, 3) = immse(y_pred_step_2, y);

% Residuals.
e_2(:, 3) = y_pred_step_2 - y;

% Adjusted R^2.
adjR2_2(3) =  1 - ((n-1)/(n-2)) * (sum(e_2(:, 3).^2)) / (sum((y - mu_y).^2));

% Standarised Error.
std_e_2 = std(e_2);
stand_e_2 = e_2 ./ std_e_2;

% Diagnostic Plots for each Model.
figure('NumberTitle', 'off', 'Name', '2. Spain');
clf

% Ex.5 Model.
subplot(1, 3, 1); 
scatter(y_pred_2, stand_e_2(:, 1))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', ['Linear Model, Offset: ', num2str(offset(3))]})
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Multiple Linear Regression Model.
subplot(1, 3, 2);
scatter(y_pred_ml_2, stand_e_2(:, 2))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Multiple var. (21) Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Stepwise Regression Model.
subplot(1, 3, 3);
scatter(y_pred_step_2, stand_e_2(:, 3))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Stepwise Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')


%% 3. Germany
% Independent variables for Germany's Multiple Linear Regression Model.
% x(t), x(t-1), ...., x(t-20).
x = zeros( size(Germany, 2), 21);
y = Germany(2, :)';

% Residuals of the 3 Models.
e_3 = zeros(length(y), 3);

% Adjusted R^2 of the 3 Models.
adjR2_3 = zeros(1, 3);

% MSE of Linear Regressiong Model from Ex.5, Multiple Linear Regression
% Model with 21 independent var. and stepwise Model.
for t = 0:20
        
    x(:, t+1) = offset_data_fun(Germany(1, :), t);
    
end

% Optimal Regression Model (Ex. 5).
X = [ones(length(x), 1) x(:, offset(4))];
beta_3 = regress(y, X);

y_pred_3 = X * beta_3;
mse_3(1) = immse(y_pred_3, y);

% Residuals.
e_3(:, 1) = y_pred_3 - y;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
adjR2_3(1) =  1 - ((n-1)/(n-2)) * (sum(e_3(:, 1).^2)) / (sum((y - mu_y).^2));

% Regression Model including all 21 independent variables.
X = [ones(length(x), 1) x];
beta_3_ml = regress(y, X);

% Predictive y from Multiple Lin. Regression.
y_pred_ml_3 = X * beta_3_ml;
mse_3(2) = immse(y_pred_ml_3, y);

% Residuals.
e_3(:, 2) = y_pred_ml_3 - y;

% Adjusted R^2.
adjR2_3(2) =  1 - ((n-1)/(n-2)) * (sum(e_3(:, 2).^2)) / (sum((y - mu_y).^2));

% Stepwise Model (Using the most significant independent variables.)
[beta_3_step, ~, ~, mdl_3, stats_3] = stepwisefit(x, y, 'Display', 'off');
b0 = stats_3.intercept;

% Keep only significant variables.
i = find(mdl_3 == 1);
x_step = x(:, i);

y_pred_step_3 = [ones(length(x_step), 1) x_step] * ([b0; beta_3_step(i)]);
mse_3(3) = immse(y_pred_step_3, y);

% Residuals.
e_3(:, 3) = y_pred_step_3 - y;

% Adjusted R^2.
adjR2_3(3) =  1 - ((n-1)/(n-2)) * (sum(e_3(:, 3).^2)) / (sum((y - mu_y).^2));

% Standarised Error.
std_e_3 = std(e_3);
stand_e_3 = e_3 ./ std_e_3;

% Diagnostic Plots for each Model.
figure('NumberTitle', 'off', 'Name', '3. Germany');
clf

% Ex.5 Model.
subplot(1, 3, 1); 
scatter(y_pred_3, stand_e_3(:, 1))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', ['Linear Model, Offset: ', num2str(offset(4))]})
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Multiple Linear Regression Model.
subplot(1, 3, 2);
scatter(y_pred_ml_3, stand_e_3(:, 2))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Multiple var. (21) Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Stepwise Regression Model.
subplot(1, 3, 3);
scatter(y_pred_step_3, stand_e_3(:, 3))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Stepwise Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')


%% 4. Portugal.
% Independent variables for Portugal's Multiple Linear Regression Model.
% x(t), x(t-1), ...., x(t-20).
x = zeros( size(Portugal, 2), 21);
y = Portugal(2, :)';

% Residuals of the 3 Models.
e_4 = zeros(length(y), 3);

% Adjusted R^2 of the 3 Models.
adjR2_4 = zeros(1, 3);

% MSE of Linear Regressiong Model from Ex.5, Multiple Linear Regression
% Model with 21 independent var. and stepwise Model.
for t = 0:20
        
    x(:, t+1) = offset_data_fun(Portugal(1, :), t);
    
end

% Optimal Regression Model (Ex. 5).
X = [ones(length(x), 1) x(:, offset(5))];
beta_4 = regress(y, X);

y_pred_4 = X * beta_4;
mse_4(1) = immse(y_pred_4, y);

% Residuals.
e_4(:, 1) = y_pred_4 - y;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
adjR2_4(1) =  1 - ((n-1)/(n-2)) * (sum(e_4(:, 1).^2)) / (sum((y - mu_y).^2));

% Regression Model including all 21 independent variables.
X = [ones(length(x), 1) x];
beta_4_ml = regress(y, X);

% Predictive y from Multiple Lin. Regression.
y_pred_ml_4 = X * beta_4_ml;
mse_4(2) = immse(y_pred_ml_4, y);

% Residuals.
e_4(:, 2) = y_pred_ml_4 - y;

% Adjusted R^2.
adjR2_4(2) =  1 - ((n-1)/(n-2)) * (sum(e_4(:, 2).^2)) / (sum((y - mu_y).^2));

% Stepwise Model (Using the most significant independent variables.)
[beta_4_step, ~, ~, mdl_4, stats_4] = stepwisefit(x, y, 'Display', 'off');
b0 = stats_4.intercept;

% Keep only significant variables.
i = find(mdl_4 == 1);
x_step = x(:, i);

y_pred_step_4 = [ones(length(x_step), 1) x_step] * ([b0; beta_4_step(i)]);
mse_4(3) = immse(y_pred_step_4, y);

% Residuals.
e_4(:, 3) = y_pred_step_4 - y;

% Adjusted R^2.
adjR2_4(3) =  1 - ((n-1)/(n-2)) * (sum(e_4(:, 3).^2)) / (sum((y - mu_y).^2));

% Standarised Error.
std_e_4 = std(e_4);
stand_e_4 = e_4 ./ std_e_4;

% Diagnostic Plots for each Model.
figure('NumberTitle', 'off', 'Name', '4. Portugal');
clf

% Ex.5 Model.
subplot(1, 3, 1); 
scatter(y_pred_4, stand_e_4(:, 1))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', ['Linear Model, Offset: ', num2str(offset(5))]})
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Multiple Linear Regression Model.
subplot(1, 3, 2);
scatter(y_pred_ml_4, stand_e_4(:, 2))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Multiple var. (21) Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Stepwise Regression Model.
subplot(1, 3, 3);
scatter(y_pred_step_4, stand_e_4(:, 3))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Stepwise Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

%% 5. France.
% Independent variables for France's Multiple Linear Regression Model.
% x(t), x(t-1), ...., x(t-20).
x = zeros( size(France, 2), 21);
y = France(2, :)';

% Residuals of the 3 Models.
e_5 = zeros(length(y), 3);

% Adjusted R^2 of the 3 Models.
adjR2_5 = zeros(1, 3);

% MSE of Linear Regressiong Model from Ex.5, Multiple Linear Regression
% Model with 21 independent var. and stepwise Model.
for t = 0:20
        
    x(:, t+1) = offset_data_fun(France(1, :), t);
    
end

% Optimal Regression Model (Ex. 5).
X = [ones(length(x), 1) x(:, offset(6))];
beta_5 = regress(y, X);

y_pred_5 = X * beta_5;
mse_5(1) = immse(y_pred_5, y);

% Residuals.
e_5(:, 1) = y_pred_5 - y;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
adjR2_5(1) =  1 - ((n-1)/(n-2)) * (sum(e_5(:, 1).^2)) / (sum((y - mu_y).^2));

% Regression Model including all 21 independent variables.
X = [ones(length(x), 1) x];
beta_5_ml = regress(y, X);

% Predictive y from Multiple Lin. Regression.
y_pred_ml_5 = X * beta_5_ml;
mse_5(2) = immse(y_pred_ml_5, y);

% Residuals.
e_5(:, 2) = y_pred_ml_5 - y;

% Adjusted R^2.
adjR2_5(2) =  1 - ((n-1)/(n-2)) * (sum(e_5(:, 2).^2)) / (sum((y - mu_y).^2));

% Stepwise Model (Using the most significant independent variables.)
[beta_5_step, ~, ~, mdl_5, stats_5] = stepwisefit(x, y, 'Display', 'off');
b0 = stats_5.intercept;

% Keep only significant variables.
i = find(mdl_5 == 1);
x_step = x(:, i);

y_pred_step_5 = [ones(length(x_step), 1) x_step] * ([b0; beta_5_step(i)]);
mse_5(3) = immse(y_pred_step_5, y);

% Residuals.
e_5(:, 3) = y_pred_step_5 - y;

% Adjusted R^2.
adjR2_5(3) =  1 - ((n-1)/(n-2)) * (sum(e_5(:, 3).^2)) / (sum((y - mu_y).^2));

% Standarised Error.
std_e_5 = std(e_5);
stand_e_5 = e_5 ./ std_e_5;

% Diagnostic Plots for each Model.
figure('NumberTitle', 'off', 'Name', '5. France');
clf

% Ex.5 Model.
subplot(1, 3, 1); 
scatter(y_pred_5, stand_e_5(:, 1))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', ['Linear Model, Offset: ', num2str(offset(6))]})
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Multiple Linear Regression Model.
subplot(1, 3, 2);
scatter(y_pred_ml_5, stand_e_5(:, 2))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Multiple var. (21) Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

% Stepwise Regression Model.
subplot(1, 3, 3);
scatter(y_pred_step_5, stand_e_5(:, 3))
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1])
hold on
plot([ax(1) ax(2)], zcrit*[1, 1])
title({'Diagnostic Plot', 'Stepwise Model'} )
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

%% Prints:
fprintf('Adjusted R^2 for every Model of every country: (order: Linear Model, Multiple Lin. Model, Stepwise Model)\n')
fprintf('A. Austria:\t\t%1.4f\t%1.4f\t%1.4f\n', adjR2_A)
fprintf('1. Finland:\t\t%1.4f\t%1.4f\t%1.4f\n', adjR2_1)
fprintf('2. Spain:\t\t%1.4f\t%1.4f\t%1.4f\n', adjR2_2)
fprintf('3. Germany:\t\t%1.4f\t%1.4f\t%1.4f\n', adjR2_3)
fprintf('4. Portugal:\t%1.4f\t%1.4f\t%1.4f\n', adjR2_4)
fprintf('5. France:\t\t%1.4f\t%1.4f\t%1.4f\n', adjR2_5)
fprintf('\nMean Squared Error for every Model of every country: (order: Linear Model, Multiple Lin. Model, Stepwise Model)\n')
fprintf('A. Austria:\t\t%1.3f\t\t%1.3f\t\t%1.3f\n', mse_A)
fprintf('1. Finland:\t\t%1.3f\t\t%1.3f\t\t%1.3f\n', mse_1)
fprintf('2. Spain:\t\t%1.3f\t%1.3f\t%1.3f\n', mse_2)
fprintf('3. Germany:\t\t%1.3f\t%1.3f\t\t%1.3f\n', mse_3)
fprintf('4. Portugal:\t%1.3f\t\t%1.3f\t\t%1.3f\n', mse_4)
fprintf('5. France:\t\t%1.3f\t%1.3f\t%1.3f\n', mse_5)

%% Comments:
% As shown in the console, the model that fits best to the data of every country
% is the Multiple Linear Regression Model that includes all 21 independent variables
% (obviously Adjusted R^2 is coherent with MSE).
% This Model seems to be good for most of the countries. Only Finland
% seem to have insufficient performance for this Model, while the rest of the
% countries got a high Adjusted R^2 value.
% Also, Stepwise Model adapts better to the data than Simple Linear
% Regression Model.
% As for MSE the reason why Spain, Germany, France have such big values is
% that their Daily Deaths are much higher than other countries, and so
% their error values get bigger values.
% In conclusion, Multiple Linear Regression Model of 21 independent var.
% can make way better predictions that Simple Linear Model for every
% country.
% Diagnostic Plots of every country still show points outside diagnostic boundaries.
% Multiple Linear Regression and Stepwise Model improve Diagnostic Test for
% some countries.

%% This function generates the data set that comes from a certain value of
% t. t is the time offset given and x is column vector.
function x = offset_data_fun(data_x, t)

% Shift x vector to the right by offset.
x = zeros(length(data_x), 1);
x(t+1 : end) = data_x(1 : end-t);

end
