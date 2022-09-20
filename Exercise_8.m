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


%% Use a method to reduce the number of independent variables of the Multiple
% Linear Regression Model (21 ind. variables). In this case, PLR (Principal
% Component Regression) is applied.
% Source: https://www.mathworks.com/matlabcentral/answers/406353-how-can-i-test-pcr-model
n = 21;

%% A. Austria.

% FIRST WAVE.
% Normalize data.
Austria_w1(1, :) = Austria_w1(1, :) ./ sum(abs(Austria_w1(1, :)));
Austria_w1(2, :) = Austria_w1(2, :) ./ sum(abs(Austria_w1(2, :)));

x_train_A = zeros( size(Austria_w1, 2), n);
y_train_A = Austria_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_A = zeros(1, 2);

for t = 0:n-1
       
    x_train_A(:, t+1) = offset_data_fun(Austria_w1(1, :), t);
    
end

% Apply PCA to training data.
[PCA_loadings_A, PCA_scores_A, PCA_eig_A] = pca(x_train_A, 'Economy', false);

% Scree plot
figure(1)
clf
plot(1:n, PCA_eig_A, 'ko-')
title({'A. Austria', 'Training data Scree Plot'})
xlabel('Index')
ylabel('Eigenvalues')

% The first 3 eigenvalues are selected.
eig_id = 3;

% Regression coefficients, using the 3 significant eigenvectors.
betaPCR_A = regress(y_train_A - mean(y_train_A), PCA_scores_A(:, 1:eig_id));

% Transform Regressiong coeff. back to initial vector space.
betaPCR_A = PCA_loadings_A(:, 1:eig_id) * betaPCR_A;

% Add intercept (b0 coefficient).
betaPCR_A = [mean(y_train_A) - mean(x_train_A) * betaPCR_A; betaPCR_A];

% Predicted values.
y_pred_train_A = [ones(length(x_train_A), 1) x_train_A] * betaPCR_A;

% Residuals.
e_A_1 = y_pred_train_A - y_train_A;

% Adjusted R^2.
m = length(y_train_A);
adjR2_A(1) = 1 - ((m-1)/(m-2)) * (sum(e_A_1.^2)) / (sum((y_train_A - mean(y_train_A)).^2));

% SECOND WAVE.
% Normalize data.
Austria_w2(1, :) = Austria_w2(1, :) ./ sum(abs(Austria_w2(1, :)));
Austria_w2(2, :) = Austria_w2(2, :) ./ sum(abs(Austria_w2(2, :)));

x_test_A = zeros( size(Austria_w2, 2), n);
y_test_A = Austria_w2(2, :)';

for t = 0:n-1
        
    x_test_A(:, t+1) = offset_data_fun(Austria_w2(1, :), t);

end

% Use regression coefficients of training data to predictions
% Predicted values.
y_pred_test_A = [ones(length(x_test_A), 1) x_test_A] * betaPCR_A;

% Residuals.
e_A_2 = y_pred_test_A - y_test_A;

% Adjusted R^2.
m = length(y_test_A);
adjR2_A(2) = 1 - ((m-1)/(m-2)) * (sum(e_A_2.^2)) / (sum((y_test_A - mean(y_test_A)).^2));


%% 1. Finland.

% FIRST WAVE.
% Normalize data.
Finland_w1(1, :) = Finland_w1(1, :) ./ sum(abs(Finland_w1(1, :)));
Finland_w1(2, :) = Finland_w1(2, :) ./ sum(abs(Finland_w1(2, :)));

x_train_1 = zeros( size(Finland_w1, 2), n);
y_train_1 = Finland_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_1 = zeros(1, 2);

for t = 0:n-1
       
    x_train_1(:, t+1) = offset_data_fun(Finland_w1(1, :), t);
    
end

% Apply PCA to training data.
[PCA_loadings_1, PCA_scores_1, PCA_eig_1] = pca(x_train_1, 'Economy', false);

% Scree plot
figure(2)
clf
plot(1:n, PCA_eig_1, 'ko-')
title({'1. Finland', 'Training data Scree Plot'})
xlabel('Index')
ylabel('Eigenvalues')

% The first 3 eigenvalues are selected.
eig_id = 3;

% Regression coefficients, using the 3 significant eigenvectors.
betaPCR_1 = regress(y_train_1 - mean(y_train_1), PCA_scores_1(:, 1:eig_id));

% Transform Regressiong coeff. back to initial vector space.
betaPCR_1 = PCA_loadings_1(:, 1:eig_id) * betaPCR_1;

% Add intercept (b0 coefficient).
betaPCR_1 = [mean(y_train_1) - mean(x_train_1) * betaPCR_1; betaPCR_1];

% Predicted values.
y_pred_train_1 = [ones(length(x_train_1), 1) x_train_1] * betaPCR_1;

% Residuals.
e_1_1 = y_pred_train_1 - y_train_1;

% Adjusted R^2.
m = length(y_train_1);
adjR2_1(1) = 1 - ((m-1)/(m-2)) * (sum(e_1_1.^2)) / (sum((y_train_1 - mean(y_train_1)).^2));

% SECOND WAVE.
% Normalize data.
Finland_w2(1, :) = Finland_w2(1, :) ./ sum(abs(Finland_w2(1, :)));
Finland_w2(2, :) = Finland_w2(2, :) ./ sum(abs(Finland_w2(2, :)));

x_test_1 = zeros( size(Finland_w2, 2), n);
y_test_1 = Finland_w2(2, :)';

for t = 0:n-1
        
    x_test_1(:, t+1) = offset_data_fun(Finland_w2(1, :), t);

end

% Use regression coefficients of training data to predictions
% Predicted values.
y_pred_test_1 = [ones(length(x_test_1), 1) x_test_1] * betaPCR_1;

% Residuals.
e_1_2 = y_pred_test_1 - y_test_1;

% Adjusted R^2.
m = length(y_test_1);
adjR2_1(2) = 1 - ((m-1)/(m-2)) * (sum(e_1_2.^2)) / (sum((y_test_1 - mean(y_test_1)).^2));


%% 2. Spain.

% FIRST WAVE.
% Normalize data.
Spain_w1(1, :) = Spain_w1(1, :) ./ sum(abs(Spain_w1(1, :)));
Spain_w1(2, :) = Spain_w1(2, :) ./ sum(abs(Spain_w1(2, :)));

x_train_2 = zeros( size(Spain_w1, 2), n);
y_train_2 = Spain_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_2 = zeros(1, 2);

for t = 0:n-1
       
    x_train_2(:, t+1) = offset_data_fun(Spain_w1(1, :), t);
    
end

% Apply PCA to training data.
[PCA_loadings_2, PCA_scores_2, PCA_eig_2] = pca(x_train_2, 'Economy', false);

% Scree plot
figure(3)
clf
plot(1:n, PCA_eig_2, 'ko-')
title({'2. Spain', 'Training data Scree Plot'})
xlabel('Index')
ylabel('Eigenvalues')

% The first 3 eigenvalues are selected.
eig_id = 3;

% Regression coefficients, using the 3 significant eigenvectors.
betaPCR_2 = regress(y_train_2 - mean(y_train_2), PCA_scores_2(:, 1:eig_id));

% Transform Regressiong coeff. back to initial vector space.
betaPCR_2 = PCA_loadings_2(:, 1:eig_id) * betaPCR_2;

% Add intercept (b0 coefficient).
betaPCR_2 = [mean(y_train_2) - mean(x_train_2) * betaPCR_2; betaPCR_2];

% Predicted values.
y_pred_train_2 = [ones(length(x_train_2), 1) x_train_2] * betaPCR_2;

% Residuals.
e_2_1 = y_pred_train_2 - y_train_2;

% Adjusted R^2.
m = length(y_train_2);
adjR2_2(1) = 1 - ((m-1)/(m-2)) * (sum(e_2_1.^2)) / (sum((y_train_2 - mean(y_train_2)).^2));

% SECOND WAVE.
% Normalize data.
Spain_w2(1, :) = Spain_w2(1, :) ./ sum(abs(Spain_w2(1, :)));
Spain_w2(2, :) = Spain_w2(2, :) ./ sum(abs(Spain_w2(2, :)));

x_test_2 = zeros( size(Spain_w2, 2), n);
y_test_2 = Spain_w2(2, :)';

for t = 0:n-1
        
    x_test_2(:, t+1) = offset_data_fun(Spain_w2(1, :), t);

end

% Use regression coefficients of training data to predictions
% Predicted values.
y_pred_test_2 = [ones(length(x_test_2), 1) x_test_2] * betaPCR_2;

% Residuals.
e_2_2 = y_pred_test_2 - y_test_2;

% Adjusted R^2.
m = length(y_test_2);
adjR2_2(2) = 1 - ((m-1)/(m-2)) * (sum(e_2_2.^2)) / (sum((y_test_2 - mean(y_test_2)).^2));


%% 3. Germany.

% FIRST WAVE.
% Normalize data.
Germany_w1(1, :) = Germany_w1(1, :) ./ sum(abs(Germany_w1(1, :)));
Germany_w1(2, :) = Germany_w1(2, :) ./ sum(abs(Germany_w1(2, :)));

x_train_3 = zeros( size(Germany_w1, 2), n);
y_train_3 = Germany_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_3 = zeros(1, 2);

for t = 0:n-1
       
    x_train_3(:, t+1) = offset_data_fun(Germany_w1(1, :), t);
    
end

% Apply PCA to training data.
[PCA_loadings_3, PCA_scores_3, PCA_eig_3] = pca(x_train_3, 'Economy', false);

% Scree plot
figure(4)
clf
plot(1:n, PCA_eig_3, 'ko-')
title({'3. Germany', 'Training data Scree Plot'})
xlabel('Index')
ylabel('Eigenvalues')

% The first 3 eigenvalues are selected.
eig_id = 3;

% Regression coefficients, using the 3 significant eigenvectors.
betaPCR_3 = regress(y_train_3 - mean(y_train_3), PCA_scores_3(:, 1:eig_id));

% Transform Regressiong coeff. back to initial vector space.
betaPCR_3 = PCA_loadings_3(:, 1:eig_id) * betaPCR_3;

% Add intercept (b0 coefficient).
betaPCR_3 = [mean(y_train_3) - mean(x_train_3) * betaPCR_3; betaPCR_3];

% Predicted values.
y_pred_train_3 = [ones(length(x_train_3), 1) x_train_3] * betaPCR_3;

% Residuals.
e_3_1 = y_pred_train_3 - y_train_3;

% Adjusted R^2.
m = length(y_train_3);
adjR2_3(1) = 1 - ((m-1)/(m-2)) * (sum(e_3_1.^2)) / (sum((y_train_3 - mean(y_train_3)).^2));

% SECOND WAVE.
% Normalize data.
Germany_w2(1, :) = Germany_w2(1, :) ./ sum(abs(Germany_w2(1, :)));
Germany_w2(2, :) = Germany_w2(2, :) ./ sum(abs(Germany_w2(2, :)));

x_test_3 = zeros( size(Germany_w2, 2), n);
y_test_3 = Germany_w2(2, :)';

for t = 0:n-1
        
    x_test_3(:, t+1) = offset_data_fun(Germany_w2(1, :), t);

end

% Use regression coefficients of training data to predictions
% Predicted values.
y_pred_test_3 = [ones(length(x_test_3), 1) x_test_3] * betaPCR_3;

% Residuals.
e_3_2 = y_pred_test_3 - y_test_3;

% Adjusted R^2.
m = length(y_test_3);
adjR2_3(2) = 1 - ((m-1)/(m-2)) * (sum(e_3_2.^2)) / (sum((y_test_3 - mean(y_test_3)).^2));


%% 4. Portugal.

% FIRST WAVE.
% Normalize data.
Portugal_w1(1, :) = Portugal_w1(1, :) ./ sum(abs(Portugal_w1(1, :)));
Portugal_w1(2, :) = Portugal_w1(2, :) ./ sum(abs(Portugal_w1(2, :)));

x_train_4 = zeros( size(Portugal_w1, 2), n);
y_train_4 = Portugal_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_4 = zeros(1, 2);

for t = 0:n-1
       
    x_train_4(:, t+1) = offset_data_fun(Portugal_w1(1, :), t);
    
end

% Apply PCA to training data.
[PCA_loadings_4, PCA_scores_4, PCA_eig_4] = pca(x_train_4, 'Economy', false);

% Scree plot
figure(5)
clf
plot(1:n, PCA_eig_4, 'ko-')
title({'4. Portugal', 'Training data Scree Plot'})
xlabel('Index')
ylabel('Eigenvalues')

% The first 3 eigenvalues are selected.
eig_id = 3;

% Regression coefficients, using the 3 significant eigenvectors.
betaPCR_4 = regress(y_train_4 - mean(y_train_4), PCA_scores_4(:, 1:eig_id));

% Transform Regressiong coeff. back to initial vector space.
betaPCR_4 = PCA_loadings_4(:, 1:eig_id) * betaPCR_4;

% Add intercept (b0 coefficient).
betaPCR_4 = [mean(y_train_4) - mean(x_train_4) * betaPCR_4; betaPCR_4];

% Predicted values.
y_pred_train_4 = [ones(length(x_train_4), 1) x_train_4] * betaPCR_4;

% Residuals.
e_4_1 = y_pred_train_4 - y_train_4;

% Adjusted R^2.
m = length(y_train_4);
adjR2_4(1) = 1 - ((m-1)/(m-2)) * (sum(e_4_1.^2)) / (sum((y_train_4 - mean(y_train_4)).^2));

% SECOND WAVE.
% Normalize data.
Germany_w2(1, :) = Germany_w2(1, :) ./ sum(abs(Germany_w2(1, :)));
Germany_w2(2, :) = Germany_w2(2, :) ./ sum(abs(Germany_w2(2, :)));

x_test_4 = zeros( size(Germany_w2, 2), n);
y_test_4 = Germany_w2(2, :)';

for t = 0:n-1
        
    x_test_4(:, t+1) = offset_data_fun(Germany_w2(1, :), t);

end

% Use regression coefficients of training data to predictions
% Predicted values.
y_pred_test_4 = [ones(length(x_test_4), 1) x_test_4] * betaPCR_4;

% Residuals.
e_4_2 = y_pred_test_4 - y_test_4;

% Adjusted R^2.
m = length(y_test_4);
adjR2_4(2) = 1 - ((m-1)/(m-2)) * (sum(e_4_2.^2)) / (sum((y_test_4 - mean(y_test_4)).^2));


%% 5. France.

% FIRST WAVE.
% Normalize data.
France_w1(1, :) = France_w1(1, :) ./ sum(abs(France_w1(1, :)));
France_w1(2, :) = France_w1(2, :) ./ sum(abs(France_w1(2, :)));

x_train_5 = zeros( size(France_w1, 2), n);
y_train_5 = France_w1(2, :)';

% Adjusted R^2 of train and test data.
adjR2_5 = zeros(1, 2);

for t = 0:n-1
       
    x_train_5(:, t+1) = offset_data_fun(France_w1(1, :), t);
    
end

% Apply PCA to training data.
[PCA_loadings_5, PCA_scores_5, PCA_eig_5] = pca(x_train_5, 'Economy', false);

% Scree plot
figure(6)
clf
plot(1:n, PCA_eig_5, 'ko-')
title({'5. France', 'Training data Scree Plot'})
xlabel('Index')
ylabel('Eigenvalues')

% The first 3 eigenvalues are selected.
eig_id = 3;

% Regression coefficients, using the 3 significant eigenvectors.
betaPCR_5 = regress(y_train_5 - mean(y_train_5), PCA_scores_5(:, 1:eig_id));

% Transform Regressiong coeff. back to initial vector space.
betaPCR_5 = PCA_loadings_5(:, 1:eig_id) * betaPCR_5;

% Add intercept (b0 coefficient).
betaPCR_5 = [mean(y_train_5) - mean(x_train_5) * betaPCR_5; betaPCR_5];

% Predicted values.
y_pred_train_5 = [ones(length(x_train_5), 1) x_train_5] * betaPCR_5;

% Residuals.
e_5_1 = y_pred_train_5 - y_train_5;

% Adjusted R^2.
m = length(y_train_5);
adjR2_5(1) = 1 - ((m-1)/(m-2)) * (sum(e_5_1.^2)) / (sum((y_train_5 - mean(y_train_5)).^2));

% SECOND WAVE.
% Normalize data.
France_w2(1, :) = France_w2(1, :) ./ sum(abs(France_w2(1, :)));
France_w2(2, :) = France_w2(2, :) ./ sum(abs(France_w2(2, :)));

x_test_5 = zeros( size(France_w2, 2), n);
y_test_5 = France_w2(2, :)';

for t = 0:n-1
        
    x_test_5(:, t+1) = offset_data_fun(France_w2(1, :), t);

end

% Use regression coefficients of training data to predictions.
% Predicted values.
y_pred_test_5 = [ones(length(x_test_5), 1) x_test_5] * betaPCR_5;

% Residuals.
e_5_2 = y_pred_test_5 - y_test_5;

% Adjusted R^2.
m = length(y_test_5);
adjR2_5(2) = 1 - ((m-1)/(m-2)) * (sum(e_5_2.^2)) / (sum((y_test_5 - mean(y_test_5)).^2));

% Sum of absolute residuals of Multiple Linear Regression Model for every
% country from Ex.7.
% row 1 --> First Wave, row 2 --> Second Wave, columns = [A, 1, 2, 3, 4, 5]
e_mul = [0.2962, 0.6273, 0.107, 0.1684, 0.1562, 0.1585;
         0.4099, 1.5312, 0.56, 0.3424, 0.1724, 0.6382];

% Sum of absolute residuals of PCR Model.
e_pcr = [sum(abs(e_A_1)), sum(abs(e_1_1)), sum(abs(e_2_1)), sum(abs(e_3_1)), sum(abs(e_4_1)), sum(abs(e_5_1));
        sum(abs(e_A_2)), sum(abs(e_1_2)), sum(abs(e_2_2)), sum(abs(e_3_2)), sum(abs(e_4_2)), sum(abs(e_5_2))];

%% Prints:
fprintf('Adjusted R^2 of PCR Model.\n')
fprintf('A. Austria:\nAdjusted R^2: Training Set: %1.4f,\t Testing Set: %1.4f', adjR2_A)
fprintf('\n1. Finland:\nAdjusted R^2: Training Set: %1.4f,\t Testing Set: %1.4f', adjR2_1)
fprintf('\n2. Spain:\nAdjusted R^2: Training Set: %1.4f,\t Testing Set: %1.4f', adjR2_2)
fprintf('\n3. Germany:\nAdjusted R^2: Training Set: %1.4f,\t Testing Set: %1.4f', adjR2_3)
fprintf('\n4. Portugal:\nAdjusted R^2: Training Set: %1.4f,\t Testing Set: %1.4f', adjR2_4)
fprintf('\n5. France:\nAdjusted R^2: Training Set: %1.4f,\t Testing Set: %1.4f\n\n', adjR2_5)

fprintf('\nResiduals of Multiple Linear Regression Model / PCR Model\n')
fprintf('A. Austria:\nTraining Set: %1.4f | %1.4f,\t Testing Set: %1.4f | %1.4f', e_mul(1, 1), e_pcr(1, 1), e_mul(2, 1), e_pcr(2, 1))
fprintf('\n1. Finland:\nTraining Set: %1.4f | %1.4f,\t Testing Set: %1.4f | %1.4f', e_mul(1, 2), e_pcr(1, 2), e_mul(2, 2), e_pcr(2, 2))
fprintf('\n2. Spain:\nTraining Set: %1.4f | %1.4f,\t Testing Set: %1.4f | %1.4f', e_mul(1, 3), e_pcr(1, 3), e_mul(2, 3), e_pcr(2, 3))
fprintf('\n3. Germany:\nTraining Set: %1.4f | %1.4f,\t Testing Set: %1.4f | %1.4f', e_mul(1, 4), e_pcr(1, 4), e_mul(2, 4), e_pcr(2, 4))
fprintf('\n4. Portugal:\nTraining Set: %1.4f | %1.4f,\t Testing Set: %1.4f | %1.4f', e_mul(1, 5), e_pcr(1, 5), e_mul(2, 5), e_pcr(2, 5))
fprintf('\n5. France:\nTraining Set: %1.4f | %1.4f,\t Testing Set: %1.4f | %1.4f\n', e_mul(1, 6), e_pcr(1, 6), e_mul(2, 6), e_pcr(2, 6))

%% Comments:
% Some countries have data of greater order of magnitude in the Second Wave than those
% in the First Wave. To face this problem, data are normalised before any processing
% and therefore, results can be compared with each other.
% As shown in the console, PCR Model makes great First Wave predictions for
% Spain and Germany and has moderate results for Portugal, France and
% Austria. Finland's Model doesn't fit to the data satisfactorily.
% Applying the same Model for each country seems to make insufficient
% predictions for every case. Germany's Model has the best performance, but
% it's still nott good enough, while all the other countries way worse. 
% Just like Ex.7, Finland has a negative Adjusted R^2, which corresponds to zero.
% Comparing PRC Model to Multiple Linear Regression Model of 21 independent
% variables, the second one can predict First Wave Daily Deaths more
% effectively. As for Second Wave, in some cases PCR Model works better (Austria,
% France), almost the same (Finland, Spain, Germany) and worse (Portugal) than Mul.
% Lin. Regression Model.
% As far as Residuals (sum of absolute residuals) of Multiple Linear Regression 
% Model and PCR Model are concerned, the first one fits training set better
% than the second in all cases but one (Finland). It has to be mentioned
% though that the performance is close for most of the countries. About the
% testing data, PCR Model makes better predictions for Austria, Finland,
% Germany and France. In conclusion, Mult. Linear Regression Model
% describes the training data better than PCR Model, but PCR predicts more
% accurately, in general, test set values.

%% This function generates the data set that comes from a certain value of
% t. t is the time offset given and x is column vector.
function x = offset_data_fun(data_x, t)

% Shift x vector to the right by offset.
x = zeros(length(data_x), 1);
x(t+1 : end) = data_x(1 : end-t);

end