% Emmanouil Savvakis, 9093
% Ioannis Gkouzoumas, 9057

clc
clear

% The starting and ending point of First and Second Wave of Covid_19 for every country
% was produced using: https://www.worldometers.info. Some countries had
% either NaN or negative values, which were restored.

% Reading data from the given files.
Input1 = readtable('Covid19Confirmed.xlsx', 'basic', true);
Input2 = readtable('Covid19Deaths.xlsx', 'basic', true);

% First Wave of Covid-19.

% 'A' country and 10 other European countries.
% A. Cases and Deaths of Austria. (10/3/2020 - 20/5/2020)
Austria_w1 = table2array([Input1(8, 73:144); Input2(8, 73:144)]);

% Fixing max deaths value (it's not on 22/4/2020, it's on 8/4).
% By changing this cell, max value falls on 8/4.
Austria_w1(2, 44) = 19;

% 1. Cases and Deaths of Finland. (12/3/2020 - 3/6/2020)
Finland_w1 = table2array([Input1(47, 75:158); Input2(47, 75:158)]);

% 2. Cases and Deaths of Spain. (3/3/2020 - 2/5/2020)
Spain_w1 = table2array([Input1(130, 66:126); Input2(130, 66:126)]);

% Fixing negative point. (24/5/2020)
Spain_w1(1, 47) = round(mean([Spain_w1(1, 44:46), Spain_w1(1, 48:50)]));

% Fixing max cases value (20/4/2020) and value of the cell previously
% listed as max (26/4).
Spain_w1(1, 18) = 10656;
Spain_w1(1, 24) = 6353;

% 3. Cases and Deaths of Germany. (12/3/2020 - 16/5/2020)
Germany_w1 = table2array([Input1(52, 75:140); Input2(52, 75:140)]);

% Fixing max deaths value (8/4/2020).
Germany_w1(2, 28) = 333;

% 4. Cases and Deaths of Portugal. (18/3/2020 - 10/5/2020)
Portugal_w1 = table2array([Input1(113, 81:134); Input2(113, 81:134)]);

% Fixing negative point. (2/5/2020)
Portugal_w1(1, 46) = round(mean([Portugal_w1(1, 43:45), Portugal_w1(1, 47:49)]));

% Fixing max cases value (10/4/2020).
Portugal_w1(1, 24) = 1726;

% Fixing max deaths value (3/4/2020).
% Also fixing 24/4 and 3/5.
Portugal_w1(2, 17) = 37;
Portugal_w1(2, 38) = 34;
Portugal_w1(2, 47) = 20;

% 5. Cases and Deaths of France. (16/3/2020 - 28/4/2020)
France_w1 = table2array([Input1(48, 79:124); Input2(48, 79:124)]);

% Fixing max deaths value (15/4/2020).
% Also fixing 3/4/2020
France_w1(2, 31) = 1437;
France_w1(2, 19) = 1119;

% Second Wave of Covid-19.

% 'A' country and 10 other European countries.
% A. Cases and Deaths of Austria. (23/10/2020 - 13/12/2020)
Austria_w2 = table2array([Input1(8, 300:351); Input2(8, 300:351)]);

% 1. Cases and Deaths of Finland. (9/9/2020 - 20/10/2020)
Finland_w2 = table2array([Input1(47, 256:297); Input2(47, 256:297)]);

% 2. Cases and Deaths of Spain. (9/7/2020 - 27/9/2020)
Spain_w2 = table2array([Input1(130, 194:274); Input2(130, 194:274)]);

% Fixing zero values (Weekends).
for i = 0:11
    
% Daily Cases.
    part = round(Spain_w2(1, 4+i*7) / 3);
    Spain_w2(1, 4+i*7) = part;
    Spain_w2(1, (4+i*7) - 1) = part;
    Spain_w2(1, (4+i*7) - 2) = part;
    
% Daily Deaths.
    part = round(Spain_w2(2, 4+i*7) / 3);
    Spain_w2(2, 4+i*7) = part;
    Spain_w2(2, (4+i*7) - 1) = part;
    Spain_w2(2, (4+i*7) - 2) = part;
          
end

% Fixing negative value (11/8/2020).
Spain_w2(2, 34) = 0;

% 3. Cases and Deaths of Germany. (5/10/2020 - 29/11/2020)
Germany_w2 = table2array([Input1(52, 282:337); Input2(52, 282:337)]);

% 4. Cases and Deaths of Portugal. (5/10/2020 - 1/12/2020)
Portugal_w2 = table2array([Input1(113, 282:339); Input2(113, 282:339)]);

% 5. Cases and Deaths of France. (1/9/2020 - 22/11/2020)
France_w2 = table2array([Input1(48, 248:330); Input2(48, 248:330)]);


%% Apply the Model with the best performance in Ex.6 for every country on the
% data of First and Second Wave.


%% A. Austria.
% Model: Multiple Linear Regression using 21 independent variables.

% FIRST WAVE.
% Normalize data.
Austria_w1(1, :) = Austria_w1(1, :) ./ sum(abs(Austria_w1(1, :)));
Austria_w1(2, :) = Austria_w1(2, :) ./ sum(abs(Austria_w1(2, :)));

x_train_A = zeros( size(Austria_w1, 2), 21);
y_train_A = Austria_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_A = zeros(1, 2);

for t = 0:20
        
    x_train_A(:, t+1) = offset_data_fun(Austria_w1(1, :), t);
    
end

X = [ones(length(x_train_A), 1) x_train_A];
beta_A = regress(y_train_A, X);

y_pred_train_A = X * beta_A;

% Prediction Error.
e_A_1 = y_pred_train_A - y_train_A;

% Adjusted R^2.
n = length(y_train_A);
mu_y_train_A = mean(y_train_A);

adjR2_A(1) = 1 - ((n-1)/(n-2)) * (sum(e_A_1.^2)) / (sum((y_train_A - mu_y_train_A).^2));

% SECOND WAVE.
% Normalize data.
Austria_w2(1, :) = Austria_w2(1, :) ./ sum(abs(Austria_w2(1, :)));
Austria_w2(2, :) = Austria_w2(2, :) ./ sum(abs(Austria_w2(2, :)));

x_test_A = zeros( size(Austria_w2, 2), 21);
y_test_A = Austria_w2(2, :)';

for t = 0:20
        
    x_test_A(:, t+1) = offset_data_fun(Austria_w2(1, :), t);
    
end

X = [ones(length(x_test_A), 1) x_test_A];
y_pred_test_A = X * beta_A;

% Prediction Error.
e_A_2 = y_pred_test_A - y_test_A;

% Adjusted R^2.
n = length(y_test_A);
mu_y_test_A = mean(y_test_A);

adjR2_A(2) = 1 - ((n-1)/(n-2)) * (sum(e_A_2.^2)) / (sum((y_test_A - mu_y_test_A).^2));


%% 1. Finland.
% Model: Multiple Linear Regression using 21 independent variables.

% FIRST WAVE.
% Normalize data.
Finland_w1(1, :) = Finland_w1(1, :) ./ sum(abs(Finland_w1(1, :)));
Finland_w1(2, :) = Finland_w1(2, :) ./ sum(abs(Finland_w1(2, :)));

x_train_1 = zeros( size(Finland_w1, 2), 21);
y_train_1 = Finland_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_1 = zeros(1, 2);

for t = 0:20
        
    x_train_1(:, t+1) = offset_data_fun(Finland_w1(1, :), t);
    
end

X = [ones(length(x_train_1), 1) x_train_1];
beta_1 = regress(y_train_1, X);

y_pred_train_1 = X * beta_1;

% Residuals.
e_1_1 = y_pred_train_1 - y_train_1;

% Adjusted R^2.
n = length(y_train_1);
mu_y_train_1 = mean(y_train_1);

adjR2_1(1) = 1 - ((n-1)/(n-2)) * (sum(e_1_1.^2)) / (sum((y_train_1 - mu_y_train_1).^2));

% SECOND WAVE.
% Normalize data.
Finland_w2(1, :) = Finland_w2(1, :) ./ sum(abs(Finland_w2(1, :)));
Finland_w2(2, :) = Finland_w2(2, :) ./ sum(abs(Finland_w2(2, :)));

x_test_1 = zeros( size(Finland_w2, 2), 21);
y_test_1 = Finland_w2(2, :)';

for t = 0:20
        
    x_test_1(:, t+1) = offset_data_fun(Finland_w2(1, :), t);
    
end

X = [ones(length(x_test_1), 1) x_test_1];
y_pred_test_1 = X * beta_1;

% Residuals.
e_1_2 = y_pred_test_1 - y_test_1;

% Adjusted R^2.
n = length(y_test_1);
mu_y_test_1 = mean(y_test_1);

adjR2_1(2) = 1 - ((n-1)/(n-2)) * (sum(e_1_2.^2)) / (sum((y_test_1 - mu_y_test_1).^2));


%% 2. Spain.
% Model: Multiple Linear Regression using 21 independent variables.

% FIRST WAVE.
% Normalize data.
Spain_w1(1, :) = Spain_w1(1, :) ./ sum(abs(Spain_w1(1, :)));
Spain_w1(2, :) = Spain_w1(2, :) ./ sum(abs(Spain_w1(2, :)));

x_train_2 = zeros( size(Spain_w1, 2), 21);
y_train_2 = Spain_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_2 = zeros(1, 2);

for t = 0:20
        
    x_train_2(:, t+1) = offset_data_fun(Spain_w1(1, :), t);
    
end

X = [ones(length(x_train_2), 1) x_train_2];
beta_2 = regress(y_train_2, X);

y_pred_train_2 = X * beta_2;

% Residuals.
e_2_1 = y_pred_train_2 - y_train_2;

% Adjusted R^2.
n = length(y_train_2);
mu_y_train_2 = mean(y_train_2);

adjR2_2(1) = 1 - ((n-1)/(n-2)) * (sum(e_2_1.^2)) / (sum((y_train_2 - mu_y_train_2).^2));

% SECOND WAVE.
% Normalize data.
Spain_w2(1, :) = Spain_w2(1, :) ./ sum(abs(Spain_w2(1, :)));
Spain_w2(2, :) = Spain_w2(2, :) ./ sum(abs(Spain_w2(2, :)));

x_test_2 = zeros( size(Spain_w2, 2), 21);
y_test_2 = Spain_w2(2, :)';

for t = 0:20
        
    x_test_2(:, t+1) = offset_data_fun(Spain_w2(1, :), t);
    
end

X = [ones(length(x_test_2), 1) x_test_2];
y_pred_test_2 = X * beta_2;

% Residuals.
e_2_2 = y_pred_test_2 - y_test_2;

% Adjusted R^2.
n = length(y_test_2);
mu_y_test_2 = mean(y_test_2);

adjR2_2(2) = 1 - ((n-1)/(n-2)) * (sum(e_2_2.^2)) / (sum((y_test_2 - mu_y_test_2).^2));


%% 3. Germany.
% Model: Multiple Linear Regression using 21 independent variables.

% FIRST WAVE.
% Normalize data.
Germany_w1(1, :) = Germany_w1(1, :) ./ sum(abs(Germany_w1(1, :)));
Germany_w1(2, :) = Germany_w1(2, :) ./ sum(abs(Germany_w1(2, :)));

x_train_3 = zeros( size(Germany_w1, 2), 21);
y_train_3 = Germany_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_3 = zeros(1, 2);

for t = 0:20
        
    x_train_3(:, t+1) = offset_data_fun(Germany_w1(1, :), t);
    
end

X = [ones(length(x_train_3), 1) x_train_3];
beta_3 = regress(y_train_3, X);

y_pred_train_3 = X * beta_3;

% Residuals.
e_3_1 = y_pred_train_3 - y_train_3;

% Adjusted R^2.
n = length(y_train_3);
mu_y_train_3 = mean(y_train_3);

adjR2_3(1) = 1 - ((n-1)/(n-2)) * (sum(e_3_1.^2)) / (sum((y_train_3 - mu_y_train_3).^2));

% SECOND WAVE.
% Normalize data.
Germany_w2(1, :) = Germany_w2(1, :) ./ sum(abs(Germany_w2(1, :)));
Germany_w2(2, :) = Germany_w2(2, :) ./ sum(abs(Germany_w2(2, :)));

x_test_3 = zeros( size(Germany_w2, 2), 21);
y_test_3 = Germany_w2(2, :)';

for t = 0:20
        
    x_test_3(:, t+1) = offset_data_fun(Germany_w2(1, :), t);
    
end

X = [ones(length(x_test_3), 1) x_test_3];
y_pred_test_3 = X * beta_3;

% Residuals.
e_3_2 = y_pred_test_3 - y_test_3;

% Adjusted R^2.
n = length(y_test_3);
mu_y_test_3 = mean(y_test_3);

adjR2_3(2) = 1 - ((n-1)/(n-2)) * (sum(e_3_2.^2)) / (sum((y_test_3 - mu_y_test_3).^2));


%% 4. Portugal.
% Model: Multiple Linear Regression using 21 independent variables.

% FIRST WAVE.
% Normalize data.
Portugal_w1(1, :) = Portugal_w1(1, :) ./ sum(abs(Portugal_w1(1, :)));
Portugal_w1(2, :) = Portugal_w1(2, :) ./ sum(abs(Portugal_w1(2, :)));

x_train_4 = zeros( size(Portugal_w1, 2), 21);
y_train_4 = Portugal_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_4 = zeros(1, 2);

for t = 0:20
        
    x_train_4(:, t+1) = offset_data_fun(Portugal_w1(1, :), t);
    
end

X = [ones(length(x_train_4), 1) x_train_4];
beta_4 = regress(y_train_4, X);

y_pred_train_4 = X * beta_4;

% Prediction Error.
e_4_1 = y_pred_train_4 - y_train_4;

% Adjusted R^2.
n = length(y_train_4);
mu_y_train_4 = mean(y_train_4);

adjR2_4(1) = 1 - ((n-1)/(n-2)) * (sum(e_4_1.^2)) / (sum((y_train_4 - mu_y_train_4).^2));

% SECOND WAVE.
% Normalize data.
Portugal_w2(1, :) = Portugal_w2(1, :) ./ sum(abs(Portugal_w2(1, :)));
Portugal_w2(2, :) = Portugal_w2(2, :) ./ sum(abs(Portugal_w2(2, :)));

x_test_4 = zeros( size(Portugal_w2, 2), 21);
y_test_4 = Portugal_w2(2, :)';

for t = 0:20
        
    x_test_4(:, t+1) = offset_data_fun(Portugal_w2(1, :), t);
    
end

X = [ones(length(x_test_4), 1) x_test_4];
y_pred_test_4 = X * beta_4;

% Residuals.
e_4_2 = y_pred_test_4 - y_test_4;

% Adjusted R^2.
n = length(y_test_4);
mu_y_test_4 = mean(y_test_4);

adjR2_4(2) = 1 - ((n-1)/(n-2)) * (sum(e_4_2.^2)) / (sum((y_test_4 - mu_y_test_4).^2));


%% 5. France.
% Model: Multiple Linear Regression using 21 independent variables.

% FIRST WAVE.
% Normalize data.
France_w1(1, :) = France_w1(1, :) ./ sum(abs(France_w1(1, :)));
France_w1(2, :) = France_w1(2, :) ./ sum(abs(France_w1(2, :)));

x_train_5 = zeros( size(France_w1, 2), 21);
y_train_5 = France_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_5 = zeros(1, 2);

for t = 0:20
        
    x_train_5(:, t+1) = offset_data_fun(France_w1(1, :), t);
    
end

X = [ones(length(x_train_5), 1) x_train_5];
beta_5 = regress(y_train_5, X);

y_pred_train_5 = X * beta_5;

% Residuals.
e_5_1 = y_pred_train_5 - y_train_5;

% Adjusted R^2.
n = length(y_train_5);
mu_y_train_5 = mean(y_train_5);

adjR2_5(1) = 1 - ((n-1)/(n-2)) * (sum(e_5_1.^2)) / (sum((y_train_5 - mu_y_train_5).^2));

% SECOND WAVE.
% Normalize data.
France_w2(1, :) = France_w2(1, :) ./ sum(abs(France_w2(1, :)));
France_w2(2, :) = France_w2(2, :) ./ sum(abs(France_w2(2, :)));

x_test_5 = zeros( size(France_w2, 2), 21);
y_test_5 = France_w2(2, :)';

for t = 0:20
        
    x_test_5(:, t+1) = offset_data_fun(France_w2(1, :), t);
    
end

X = [ones(length(x_test_5), 1) x_test_5];
y_pred_test_5 = X * beta_5;

% Residuals.
e_5_2 = y_pred_test_5 - y_test_5;

% Adjusted R^2.
n = length(y_test_5);
mu_y_test_5 = mean(y_test_5);

adjR2_5(2) = 1 - ((n-1)/(n-2)) * (sum(e_5_2.^2)) / (sum((y_test_5 - mu_y_test_5).^2));


%% Print the results.
fprintf('Adjusted R^2 of Mul. Lin. Rgression Model (21 independent var.)\n')
fprintf('\nA. Austria:\nAdjusted R^2: Training Set (First Wave): %1.4f,\t Testing Set (Second Wave): %1.4f\n', adjR2_A)
fprintf('\n1. Finland:\nAdjusted R^2: Training Set (First Wave): %1.4f,\t Testing Set (Second Wave): %1.4f\n', adjR2_1)
fprintf('\n2. Spain:\nAdjusted R^2: Training Set (First Wave): %1.4f,\t Testing Set (Second Wave): %1.4f\n', adjR2_2)
fprintf('\n3. Germany:\nAdjusted R^2: Training Set (First Wave): %1.4f,\t Testing Set (Second Wave): %1.4f\n', adjR2_3)
fprintf('\n4. Portugal:\nAdjusted R^2: Training Set (First Wave): %1.4f,\t Testing Set (Second Wave): %1.4f\n', adjR2_4)
fprintf('\n5. France:\nAdjusted R^2: Training Set (First Wave): %1.4f,\t Testing Set (Second Wave): %1.4f\n', adjR2_5)

%% Comments:
% Some countries have data of greater order of magnitude in the Second Wave than those
% in the First Wave. To face this problem, data are normalised before any processing
% and therefore, results can be compared with each other.
% As shown in the console, training sets create Regression Models that
% completely insufficient. The only exception is Portugal, whose Model can
% make even better predictions in the Second Wave than it did in the First.
% Finland has a negative Adjusted R^2 value, which is translated to 0.

%% This function generates the data set that comes from a certain value of
% t. t is the time offset given and x is column vector.
function x = offset_data_fun(data_x, t)

% Shift x vector to the right by offset.
x = zeros(length(data_x), 1);
x(t+1 : end) = data_x(1 : end-t);

end
