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

%% Linear Regression Model (y = b1*x + b0), where independent variable is Daily
% Cases (x) and dependent var. is Daily Deaths (y).
% Create a Model for every set of x(n - t)*y(n) (or equally x(n)*y(n + t)),
% where t is an offset in [0, 20], n in [[1, data's columns]], for every country.
% Compare the Models of every country, find out which offset creates the Model that
% fits best to the data. 
% Also, for every offset, a  graph of standarised error is shown.

% beta matrices have the regression coefficients of each country for every
% offset in [0, 20]. row 1 <-- b0, row2 <-- b1
beta_A = zeros(2, 21);
beta_1 = zeros(2, 21);
beta_2 = zeros(2, 21);
beta_3 = zeros(2, 21);
beta_4 = zeros(2, 21);
beta_5 = zeros(2, 21);

% mse contains Mean Squared Error for every Lin. Regression Model.
% Each row corresponds to a country.
mse = zeros(6, 21);

% Confidence level and critical value (diagnostic plot).
alpha = 0.05;
zcrit = norminv(1-alpha/2);

for t = 0:20
    
% A. Regression Models for Austria and MSE for every model.   
    x = offset_data_fun(Austria(1, :), t);
    y = Austria(2, :)';
    
    X = [ones(length(x), 1) x];
    beta_A(:, t+1) = regress(y, X);
    
    y_pred = X * beta_A(:, t+1);
    mse(1, t+1) = immse(y_pred, y);
    
% 1. Regression Models for Finland and MSE for every model.
    x = offset_data_fun(Finland(1, :), t);
    y = Finland(2, :)';
    
    X = [ones(length(x), 1) x];
    beta_1(:, t+1) = regress(y, X);
    
    y_pred = X * beta_1(:, t+1);
    mse(2, t+1) = immse(y_pred, y);
    
% 2. Regression Models for Spain and MSE for every model.
    x = offset_data_fun(Spain(1, :), t);
    y = Spain(2, :)';
    
    X = [ones(length(x), 1) x];
    beta_2(:, t+1) = regress(y, X);
    
    y_pred = X * beta_2(:, t+1);
    mse(3, t+1) = immse(y_pred, y);

% 3. Regression Models for Germany and MSE for every model.
    x = offset_data_fun(Germany(1, :), t);
    y = Germany(2, :)';
    
    X = [ones(length(x), 1) x];
    beta_3(:, t+1) = regress(y, X);
    
    y_pred = X * beta_3(:, t+1);
    mse(4, t+1) = immse(y_pred, y);
    
% 4. Regression Models for Portugal and MSE for every model.
    x = offset_data_fun(Portugal(1, :), t);
    y = Portugal(2, :)';
    
    X = [ones(length(x), 1) x];
    beta_4(:, t+1) = regress(y, X);
    
    y_pred = X * beta_4(:, t+1);
    mse(5, t+1) = immse(y_pred, y);
    
% 5. Regression Models for France and MSE for every model.    
    x = offset_data_fun(France(1, :), t);
    y = France(2, :)';
     
    X = [ones(length(x), 1) x];
    beta_5(:, t+1) = regress(y, X);
    
    y_pred = X * beta_5(:, t+1);
    mse(6, t+1) = immse(y_pred, y);
    
end

% Find which offset creates the best Linear Regression Model for every
% country (Offset = Index - 1).
index = zeros(6, 1);

% Adjusted R^2.
AdjR2 = zeros(6, 1);

for i = 1:size(mse, 1)

    [~, index(i)] = min(mse(i, :));
    
end

% Find Standarised Error for these Models.
%% A. Austria. 
x = offset_data_fun(Austria(1, :), index(1) - 1);
y = Austria(2, :)';

X = [ones(length(x), 1) x];
y_pred_A = X * beta_A(:, index(1));

e_A = y_pred_A - y;
std_e_A = std(e_A);
stand_e_A = e_A ./ std_e_A;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
AdjR2(1) =  1 - ((n-1)/(n-2)) * (sum(e_A.^2)) / (sum((y - mu_y).^2));

% Regression plot.
figure(1)
clf
subplot(2, 1, 1);
scatter(x, y)
hold on
ax = axis;
x_ax = ax(1): 0.5 : ax(2);
X_plot = [ones(length(x_ax), 1) x_ax'];
plt = plot(X_plot * beta_A(:, index(1)), 'LineWidth', 1.2);
hold off
title({'Austria','Scatter Plot of Daily Cases - Deaths'})
xlabel('Daily Cases')
ylabel('Daily Deaths')
legend(plt, 'Linear Regression curve')

% Diagnostic Plot.
subplot(2, 1, 2); 
scatter(y_pred_A, stand_e_A)
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1], 'r--')
hold on
plot([ax(1) ax(2)], zcrit*[1, 1], 'r--')
title('Diagnostic Plot')
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

%% 1. Finland.
x = offset_data_fun(Finland(1, :), index(2)-1);
y = Finland(2, :)';

X = [ones(length(x), 1) x];
y_pred_1 = X * beta_1(:, index(2));

e_1 = y_pred_1 - y;
std_e_1 = std(e_1);
stand_e_1 = e_1 ./ std_e_1;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
AdjR2(2) =  1 - ((n-1)/(n-2)) * (sum(e_1.^2)) / (sum((y - mu_y).^2));


% Regression plot.
figure(2)
clf
subplot(2, 1, 1);
scatter(x, y)
hold on
ax = axis;
x_ax = ax(1): 0.5 : ax(2);
X_plot = [ones(length(x_ax), 1) x_ax'];
plt = plot(X_plot * beta_1(:, index(2)), 'LineWidth', 1.2);
hold off
title({'Finland','Scatter Plot of Daily Cases - Deaths'})
xlabel('Daily Cases')
ylabel('Daily Deaths')
legend(plt, 'Linear Regression curve')

% Diagnostic Plot.
subplot(2, 1, 2); 
scatter(y_pred_1, stand_e_1)
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1], 'r--')
hold on
plot([ax(1) ax(2)], zcrit*[1, 1], 'r--')
title('Diagnostic Plot')
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')


%% 2. Spain.
x = offset_data_fun(Spain(1, :), index(3)-1);
y = Spain(2, :)';

X = [ones(length(x), 1) x];
y_pred_2 = X * beta_2(:, index(3));

e_2 = y_pred_2 - y;
std_e_2 = std(e_2);
stand_e_2 = e_2 ./ std_e_2;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
AdjR2(3) =  1 - ((n-1)/(n-2)) * (sum(e_2.^2)) / (sum((y - mu_y).^2));

% Regression plot.
figure(3)
clf
subplot(2, 1, 1);
scatter(x, y)
hold on
ax = axis;
x_ax = ax(1): 0.5 : ax(2);
X_plot = [ones(length(x_ax), 1) x_ax'];
plt = plot(X_plot * beta_2(:, index(3)), 'LineWidth', 1.2);
hold off
title({'Spain','Scatter Plot of Daily Cases - Deaths'})
xlabel('Daily Cases')
ylabel('Daily Deaths')
legend(plt, 'Linear Regression curve')

% Diagnostic Plot.
subplot(2, 1, 2); 
scatter(y_pred_2, stand_e_2)
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1], 'r--')
hold on
plot([ax(1) ax(2)], zcrit*[1, 1], 'r--')
title('Diagnostic Plot')
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

%% 3. Germany. 
x = offset_data_fun(Germany(1, :), index(4)-1);
y = Germany(2, :)';

X = [ones(length(x), 1) x];
y_pred_3 = X * beta_3(:, index(4));

e_3 = y_pred_3 - y;
std_e_3 = std(e_3);
stand_e_3 = e_3 ./ std_e_3;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
AdjR2(4) =  1 - ((n-1)/(n-2)) * (sum(e_3.^2)) / (sum((y - mu_y).^2));

% Regression plot.
figure(4)
clf
subplot(2, 1, 1);
scatter(x, y)
hold on
ax = axis;
x_ax = ax(1): 0.5 : ax(2);
X_plot = [ones(length(x_ax), 1) x_ax'];
plt = plot(X_plot * beta_3(:, index(4)), 'LineWidth', 1.2);
hold off
title({'Germany','Scatter Plot of Daily Cases - Deaths'})
xlabel('Daily Cases')
ylabel('Daily Deaths')
legend(plt, 'Linear Regression curve')

% Diagnostic Plot.
subplot(2, 1, 2); 
scatter(y_pred_3, stand_e_3)
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1], 'r--')
hold on
plot([ax(1) ax(2)], zcrit*[1, 1], 'r--')
title('Diagnostic Plot')
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')


%% 4. Portugal. 
x = offset_data_fun(Portugal(1, :), index(5)-1);
y = Portugal(2, :)';

X = [ones(length(x), 1) x];
y_pred_4 = X * beta_4(:, index(5));

e_4 = y_pred_4 - y;
std_e_4 = std(e_4);
stand_e_4 = e_4 ./ std_e_4;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
AdjR2(5) =  1 - ((n-1)/(n-2)) * (sum(e_4.^2)) / (sum((y - mu_y).^2));

% Regression plot.
figure(5)
clf
subplot(2, 1, 1);
scatter(x, y)
hold on
ax = axis;
x_ax = ax(1): 0.5 : ax(2);
X_plot = [ones(length(x_ax), 1) x_ax'];
plt = plot(X_plot * beta_4(:, index(5)), 'LineWidth', 1.2);
hold off
title({'Germany','Scatter Plot of Daily Cases - Deaths'})
xlabel('Daily Cases')
ylabel('Daily Deaths')
legend(plt, 'Linear Regression curve')

% Diagnostic Plot.
subplot(2, 1, 2); 
scatter(y_pred_4, stand_e_4)
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1], 'r--')
hold on
plot([ax(1) ax(2)], zcrit*[1, 1], 'r--')
title('Diagnostic Plot')
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

%% 5. France. 
x = offset_data_fun(France(1, :), index(6)-1);
y = France(2, :)';

X = [ones(length(x), 1) x];
y_pred_5 = X * beta_5(:, index(6));

e_5 = y_pred_5 - y;
std_e_5 = std(e_5);
stand_e_5 = e_5 ./ std_e_5;

% Adjusted R^2.
mu_y = mean(y);
n = length(y);
AdjR2(6) =  1 - ((n-1)/(n-2)) * (sum(e_5.^2)) / (sum((y - mu_y).^2));

% Regression plot.
figure(6)
clf
subplot(2, 1, 1);
scatter(x, y)
hold on
ax = axis;
x_ax = ax(1): 0.5 : ax(2);
X_plot = [ones(length(x_ax), 1) x_ax'];
plt = plot(X_plot * beta_5(:, index(6)), 'LineWidth', 1.2);
hold off
title({'Germany','Scatter Plot of Daily Cases - Deaths'})
xlabel('Daily Cases')
ylabel('Daily Deaths')
legend(plt, 'Linear Regression curve')

% Diagnostic Plot.
subplot(2, 1, 2);
scatter(y_pred_5, stand_e_5)
hold on
ax = axis;
plot([ax(1) ax(2)], -zcrit*[1, 1], 'r--')
hold on
plot([ax(1) ax(2)], zcrit*[1, 1], 'r--')
hold on
title('Diagnostic Plot')
xlabel('Predicted Daily Deaths')
ylabel('Standarised Error')

%% Print:
fprintf('Best Linear Regression Models.\n')
fprintf('A. Austria:\t\tOffset: %d,\t\tAdjR^2: %1.4f\n', index(1)-1, AdjR2(1))
fprintf('1. Finland:\t\tOffset: %d,\t\tAdjR^2: %1.4f\n', index(2)-1, AdjR2(2))
fprintf('2. Spain:\t\tOffset: %d,\t\tAdjR^2: %1.4f\n', index(3)-1, AdjR2(3))
fprintf('3. Germany:\t\tOffset: %d,\t\tAdjR^2: %1.4f\n', index(4)-1, AdjR2(4))
fprintf('4. Portugal:\tOffset: %d,\t\tAdjR^2: %1.4f\n', index(5)-1, AdjR2(5))
fprintf('5. France:\t\tOffset: %d,\t\tAdjR^2: %1.4f\n', index(6)-1, AdjR2(6))

%% Comments:
% Adjusted R^2 values show that Linear Regression Model can't make precise
% predictions of Daily Deaths, depending on Daily Cases.
% For some countries (Spain, Germany) the Model can bring moderate
% results, but for the rest it's not useful at all.
% Diagnostic Plots confirm previous remarks. For every country there are
% at least 2 points that are outside diagnostic boundaries [-1.96, 1.96].

% Assumptions that can be made are:
% Austria: Non-Normal distributed data.
% Finland: Extreme value of error and independent variable.
% Spain: Non-Normal distributed data.
% Germany: Non-Normal distributed data.
% Portugal: Non-Normal distributed data.
% France: Extreme value of error and independent variable.

% The optimal offsets for every country (in this Exercise) seems to satisfy
% Ex.4 offsets. That means that only 50% of the offsets are in the
% confidence interval of Ex.3, and none of them was 14.

%% This function generates the data set that comes from a certain value of
% t. t is the time offset given and x is column vector.
function x = offset_data_fun(data_x, t)

% Shift x vector to the right by offset.
x = zeros(length(data_x), 1);
x(t+1 : end) = data_x(1 : end-t);

end

