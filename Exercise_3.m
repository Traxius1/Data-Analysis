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

% 3. Cases and Deaths of the United Kingdom. (18/3/2020 - 16/6/2020)
Un_Kingdom = table2array([Input1(147, 81:181); Input2(147, 81:181)]);

% Fixing max cases value (10/4/2020).
Un_Kingdom(1, 24) = 4858;

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

% 6. Cases and Deaths of Germany. (12/3/2020 - 16/5/2020)
Germany = table2array([Input1(52, 75:140); Input2(52, 75:140)]);

% Fixing max deaths value (8/4/2020).
Germany(2, 28) = 333;

% 7. Cases and Deaths of Switzerland. (13/3/2020 - 30/4/2020)
Switzerland = table2array([Input1(134, 76:124); Input2(134, 76:124)]);

% Fixing max cases value (20/3/2020) and the cell previously listed
% as max (27/3).
Switzerland(1, 8) = 1393;
Switzerland(1, 15) = 1117;

% Fixing max deaths value (8/4/2020) and the cell previously listed
% as max (31/3). Also fixing 15/4
Switzerland(2, 27) = 63;
Switzerland(2, 19) = 58;
Switzerland(2, 34) = 46;

% 8. Cases and Deaths of Norway. (4/3/2020 - 16/5/2020)
Norway = table2array([Input1(103, 67:140); Input2(103, 67:140)]);

% Fixing max deaths value (20/4/2020).
Norway(2, 48) = 16;

% 9. Cases and Deaths of Italy. (3/3/2020 - 19/5/2020)
Italy = table2array([Input1(67, 66:143); Input2(67, 66:143)]);

% 10. Cases and Deaths of Ireland. (18/3/2020 - 18/5/2020)
Ireland = table2array([Input1(65, 81:142); Input2(65, 81:142)]);

% Fixing max deaths value (10/4/2020).
Ireland(1, 24) = 1508;

% Fixing max deaths value (24/4/2020) and the cell previously listed
% as max (25/4).
Ireland(2, 38) = 214;
Ireland(2, 39) = 48;

% First row contains the predicted date of the maximum Daily Cases for every country.
% Second row contains the predicted date of the maximum Daily Deaths for every country.
% Third row contains the number of days between the maximum Daily Cases and max. Daily Deaths.
dt = zeros(3, 11);

% First row contains the distribution of Daily Cases data for every country.
% Second row contains the distribution of Daily Deaths data for every country.
dist = strings([2, 11]);

% Call the 'fun' function for every country (for both Daily Cases and Deaths).
[dist(1, 1), dt(1, 1)] = fun(Austria(1, :));
[dist(2, 1), dt(2, 1)] = fun(Austria(2, :));
 
[dist(1, 2), dt(1, 2)] = fun(Finland(1, :));
[dist(2, 2), dt(2, 2)] = fun(Finland(2, :));

[dist(1, 3), dt(1, 3)] = fun(Spain(1, :));
[dist(2, 3), dt(2, 3)] = fun(Spain(2, :));

[dist(1, 4), dt(1, 4)] = fun(Un_Kingdom(1, :));
[dist(2, 4), dt(2, 4)] = fun(Un_Kingdom(2, :));

[dist(1, 5), dt(1, 5)] = fun(Portugal(1, :));
[dist(2, 5), dt(2, 5)] = fun(Portugal(2, :));

[dist(1, 6), dt(1, 6)] = fun(France(1, :));
[dist(2, 6), dt(2, 6)] = fun(France(2, :));

[dist(1, 7), dt(1, 7)] = fun(Germany(1, :));
[dist(2, 7), dt(2, 7)] = fun(Germany(2, :));

[dist(1, 8), dt(1, 8)] = fun(Switzerland(1, :));
[dist(2, 8), dt(2, 8)] = fun(Switzerland(2, :));

[dist(1, 9), dt(1, 9)] = fun(Norway(1, :));
[dist(2, 9), dt(2, 9)] = fun(Norway(2, :));

[dist(1, 10), dt(1, 10)] = fun(Italy(1, :));
[dist(2, 10), dt(2, 10)] = fun(Italy(2, :));

[dist(1, 11), dt(1, 11)] = fun(Ireland(1, :));
[dist(2, 11), dt(2, 11)] = fun(Ireland(2, :));

% Difference between max. Daily Cases and max. Daily Deaths.
dt(3, :) = diff(dt(1:2, :));

% It's assumed that Daily Deaths max must be later than Daily Cases max.
% The United Kingdom dt is negative, which is not accepted.
% So, it's excluded from the c.i. calculation.

% New dt set.
new_dt = dt(3, :);
new_dt(4) = [];

% Sample is small (10 observations) and it doesn't come from Normal Distribution
% as shown in the figure (histogram - boxplot).
figure('NumberTitle', 'off', 'Name', 'Days difference between Daily Cases and Daily Deaths max. for 11 countries.');
clf
subplot(1, 2, 1)
histogram(new_dt)
title('Histogram')
xlabel('Days Difference')
ylabel('Counts')

subplot(1, 2, 2)
boxplot(new_dt)
title('Boxplot')
xlabel('Column Number')
ylabel('Days Difference')

% It can't be assumed that the sample comes from Norm. Distr.

% Parametric 95% confidence interval for the mean time between the date of
% Daily Cases and Daily Deaths maximization.
[~, ~, ci, ~] = ttest(new_dt);
fprintf('Parametric 95%% Confidence Intervals: [%1.3f, %1.3f]\n', ci(1:2))

% Bootstrap 95% confidence intervals for the mean time between the date of
% Daily Cases and Daily Deaths maximization.
% Number of bootstrap samples is 1000.
n_boot = 10000;
a = 0.05;
b_sample = bootstrp(n_boot, @mean, new_dt');
sorted_b_sample = sort(b_sample);

% Bootstrap confidence interval's boundaries.
low = fix((n_boot + 1)*(a/2));
up = n_boot - low;
ci_b = [sorted_b_sample(low, :); sorted_b_sample(up, :)];
fprintf('Bootstrap 95%% Confidence Intervals: [%1.3f, %1.3f]\n', ci_b(1:2))

if ci_b(1) < 14 && 14 < ci_b(2)
   fprintf('It can be assumed with %d%% confidence that time offset between max Daily Cases and Deaths is 14 days.\n', (1-a)*100) 
else
   fprintf('It can not be assumed with %d%% confidence that time offset between max Daily Cases and Deaths is 14 days.\n', (1-a)*100) 
end
    
%% Comments:
% As said before, parametric confidence interval can not be applied to
% this case, as shown in the figure, because sample is small (< 30) and it
% come from Normal Distribution.
% Bootstrap confidence interval include 14, so it can be said with
% confidence of 95% that the mean time between the maximization of Daily
% Cases - Deaths is 14 days.
% Sample was small (10 values), so confidence interval isn't really
% precise.


%% This function returns the distribution that fits the given data best.
% Also returns the day that the distribution predicts the max valie.
function [distribution, t] = fun(data)

% mse contains the Mean Squared Error for each distribution.
mse = zeros(1, 5);
days = length(data);

% Distributions used to fit the data:
% 1 -> Normal distribution.
% 2 -> Exponential distribution.
% 3 -> Gamma distribution.
% 4 -> Weibull distribution.
% 5 -> Poisson distribution.

% Normalize data. Every cell is divided with the sum of absolute values.
data_n = data ./ sum(abs(data));

% Fitting days' data, using function's input data as frequency, to parametric distribution.
% Check Goodness-of-fit using Mean Squared Error.

% Normal Distribution.
pd1 = fitdist((1:days)' , 'Normal', 'Frequency', data);
pdf1 = pdf(pd1, 1:days);
mse(1, 1) = immse(pdf1, data_n);

% Exponential Distribution.
pd2 = fitdist((1:days)' , 'Exponential', 'Frequency', data);
pdf2 = pdf(pd2, 1:days);
mse(1, 2) = immse(pdf2, data_n);

% Gamma Distribution.
pd3 = fitdist((1:days)' , 'Gamma', 'Frequency', data);
pdf3 = pdf(pd3, 1:days);
mse(1, 3) = immse(pdf3, data_n);

% Weibull Distribution.
pd4 = fitdist((1:days)' , 'Weibull', 'Frequency', data);
pdf4 = pdf(pd4, 1:days);
mse(1, 4) = immse(pdf4, data_n);

% Poisson Distribution.
pd5 = fitdist((1:days)' , 'Poisson', 'Frequency', data);
pdf5 = pdf(pd5, 1:days);
mse(1, 5) = immse(pdf5, data_n);

% Sort Mean Squared Error in ascending order.
[~, index] = sort(mse);
dist = {'Normal', 'Exponential', 'Gamma', 'Weibull', 'Poisson'};

% The distribution that fits the data best is returned from the function.
distribution = dist{index(1)};

% Matrix with all pdfs calculated before.
pdf_all = [pdf1; pdf2; pdf3; pdf4; pdf5];

% Date that the distribution model predicts the max value of the data occurs.
[~, t] = max( pdf_all(index(1), :) );

end