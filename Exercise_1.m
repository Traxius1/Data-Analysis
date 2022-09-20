clc
clear

% --- Introduction ---
% The country assigned to our team is Austria.

https://www.worldometers.info/coronavirus/country/austria/?fbclid=IwAR3EB2LWK_viMZEFs1o8mssKhRnXMY7NMu7g3ACIiPjTGeao4ntQGjzBXj4
% shows Austria's number of cases relative to time. We thought that we
% should define the starting of the First Wave where the curve rises
% and the ending at the point that the curve maintains at a constant value.
% According to the graph, the final decision about the length of the First 
% Covid-19 Wave is assumed to occur 10/3/2020 - 20/5/2020.

% Reading data from the given files.
Input1 = readtable('Covid19Confirmed.xlsx', 'basic', true);
Input2 = readtable('Covid19Deaths.xlsx', 'basic', true);

% Collecting the data needed.
cases = table2array(Input1(8, 73:144));
deaths = table2array(Input2(8, 73:144));
days = length(cases);

% Fixing max deaths value (it's not on 22/4/2020, it's on 8/4).
% By changing this cell, max value falls on 8/4.
deaths(44) = 19;

% mse contains the Mean Squared Error for each distribution.
% row 1 <- Daily Cases.   row 2 <- Daily Deaths.
mse = zeros(2, 5);

% Distributions used to fit the data:
% 1 -> Normal distribution.
% 2 -> Exponential distribution.
% 3 -> Gamma distribution.
% 4 -> Weibull distribution.
% 5 -> Poisson distribution.

%% Cases' Data.

% Normalize data. Every cell is divided with the sum of absolute values.
cases_n = cases ./ sum(abs(cases));
deaths_n = deaths ./ sum(abs(deaths));

% Fitting days' data, using daily cases as frequency, to parametric distribution.
% Check Goodness-of-fit using Mean Squared Error.

% Normal Distribution.
pd1_c = fitdist((1:days)' , 'Normal', 'Frequency', cases);
pdf1_c = pdf(pd1_c, 1:days);
mse(1, 1) = immse(pdf1_c, cases_n);

% Exponential Distribution.
pd2_c = fitdist((1:days)' , 'Exponential', 'Frequency', cases);
pdf2_c = pdf(pd2_c, 1:days);
mse(1, 2) = immse(pdf2_c, cases_n);

% Gamma Distribution.
pd3_c = fitdist((1:days)' , 'Gamma', 'Frequency', cases);
pdf3_c = pdf(pd3_c, 1:days);
mse(1, 3) = immse(pdf3_c, cases_n);

% Weibull Distribution.
pd4_c = fitdist((1:days)' , 'Weibull', 'Frequency', cases);
pdf4_c = pdf(pd4_c, 1:days);
mse(1, 4) = immse(pdf4_c, cases_n);

% Poisson Distribution.
pd5_c = fitdist((1:days)' , 'Poisson', 'Frequency', cases);
pdf5_c = pdf(pd5_c, 1:days);
mse(1, 5) = immse(pdf5_c, cases_n);

%% Deaths' Data.

% Fitting days' data, using daily cases as frequency, to parametric distribution.
% Check Goodness-of-fit using Mean Squared Error.

% Normal Distribution.
pd1_d = fitdist((1:days)' , 'Normal', 'Frequency', deaths);
pdf1_d = pdf(pd1_d, 1:days);
mse(2, 1) = immse(pdf1_d, deaths_n);

% Exponential Distribution.
pd2_d = fitdist((1:days)' , 'Exponential', 'Frequency', deaths);
pdf2_d = pdf(pd2_d, 1:days);
mse(2, 2) = immse(pdf2_d, deaths_n);

% Gamma Distribution.
pd3_d = fitdist((1:days)' , 'Gamma', 'Frequency', deaths);
pdf3_d = pdf(pd3_d, 1:days);
mse(2, 3) = immse(pdf3_d, deaths_n);

% Weibull Distribution.
pd4_d = fitdist((1:days)' , 'Weibull', 'Frequency', deaths);
pdf4_d = pdf(pd4_d, 1:days);
mse(2, 4) = immse(pdf4_d, deaths_n);

% Poisson Distribution.
pd5_d = fitdist((1:days)' , 'Poisson', 'Frequency', deaths);
pdf5_d = pdf(pd5_d, 1:days);
mse(2, 5) = immse(pdf5_d, deaths_n);


%% Print the result.
distributions = {'Normal', 'Exponential', 'Gamma', 'Weibull', 'Poisson'};
[~, index_c] = sort(mse(1, :));
[~, index_d] = sort(mse(2, :));
fprintf('Distribution that fits the Daily Cases data best is the %s distribution.\n', distributions{index_c(1)});
fprintf('Distribution that fits the Daily Deaths data best is the %s distribution.\n', distributions{index_d(1)});

% Daily Cases figure.
figure;
clf
bar(cases_n)
hold on
plot(1:days, pdf3_c, 'LineWidth', 2)
title('Fitting of Gamma Distribution to Daily Cases distribution.')
xlabel('Time (days)')
ylabel('Daily Cases')

% Other Distribution fittings.
figure('NumberTitle', 'off', 'Name', 'Other Distribution fittings (Daily Cases).');
clf
subplot(2, 2, 1)
bar(cases_n)
hold on
plot(1:days, pdf1_c, 'LineWidth', 1.5)
title('Normal Distribution')
xlabel('Time (days)')
ylabel('Daily Cases')

subplot(2, 2, 2)
bar(cases_n)
hold on
plot(1:days, pdf2_c, 'LineWidth', 1.5)
title('Exponential Distribution')
xlabel('Time (days)')
ylabel('Daily Cases')

subplot(2, 2, 3)
bar(cases_n)
hold on
plot(1:days, pdf4_c, 'LineWidth', 1.5)
title('Weibull Distribution')
xlabel('Time (days)')
ylabel('Daily Cases')

subplot(2, 2, 4)
bar(cases_n)
hold on
plot(1:days, pdf5_c, 'LineWidth', 1.5)
title('Poisson Distribution')
xlabel('Time (days)')
ylabel('Daily Cases')


% Daily Deaths figure.
figure;
clf
bar(deaths_n)
hold on
plot(1:days, pdf3_d, 'LineWidth', 2)
title('Fitting of Gamma Distribution to Daily Deaths distribution.')
xlabel('Time (days)')
ylabel('Daily Deaths')

% Other Distribution fittings.
figure('NumberTitle', 'off', 'Name', 'Other Distribution fittings (Daily Deaths).');
clf
subplot(2, 2, 1)
bar(deaths_n)
hold on
plot(1:days, pdf1_d, 'LineWidth', 1.5)
title('Normal Distribution')
xlabel('Time (days)')
ylabel('Daily Deaths')

subplot(2, 2, 2)
bar(deaths_n)
hold on
plot(1:days, pdf2_d, 'LineWidth', 1.5)
title('Exponential Distribution')
xlabel('Time (days)')
ylabel('Daily Deaths')

subplot(2, 2, 3)
bar(deaths_n)
hold on
plot(1:days, pdf4_d, 'LineWidth', 1.5)
title('Weibull Distribution')
xlabel('Time (days)')
ylabel('Daily Deaths')

subplot(2, 2, 4)
bar(deaths_n)
hold on
plot(1:days, pdf5_d, 'LineWidth', 1.5)
title('Poisson Distribution')
xlabel('Time (days)')
ylabel('Daily Deaths')

%% Comments:
% As shown in console, Daily Cases and Deaths are can be predicted from the
% same distribution, which is Gamma.
