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

% Mean Square Error for every country. First row contains MSE for daily
% cases and second row contains MSE for daily deaths.
mse = zeros(2, 10);

% Days vector contains the number of columns of every country's data.
days = zeros(1, 10);

% Countries.
id = [47, 130, 147, 113, 48, 52, 134, 103, 67, 65];
countries = table2array(Input1(id, 1));

% Distribution selected from Exersice 1.
dist = 'Gamma';

% Index for simple transition between vectors' columns.
i = 1;

%% 1. Goodness-of-fit: Finland.
x = Finland(1, :);
y = Finland(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd1_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf1_c = pdf(pd1_c, 1:days(i));
mse(1, i) = immse(pdf1_c, cases_n);

pd1_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf1_d = pdf(pd1_d, 1:days(i));
mse(2, i) = immse(pdf1_d, deaths_n);

%% 2. Goodness-of-fit: Spain.
i = i + 1;
x = Spain(1, :);
y = Spain(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd2_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf2_c = pdf(pd2_c, 1:days(i));
mse(1, i) = immse(pdf2_c, cases_n);

pd2_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf2_d = pdf(pd2_d, 1:days(i));
mse(2, i) = immse(pdf2_d, deaths_n);

%% 3. Goodness-of-fit: the United Kingdom.
i = i + 1;
x = Un_Kingdom(1, :);
y = Un_Kingdom(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd3_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf3_c = pdf(pd3_c, 1:days(i));
mse(1, i) = immse(pdf3_c, cases_n);

pd3_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf3_d = pdf(pd3_d, 1:days(i));
mse(2, i) = immse(pdf3_d, deaths_n);

%% 4. Goodness-of-fit: Portugal.
i = i + 1;
x = Portugal(1, :);
y = Portugal(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd4_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf4_c = pdf(pd4_c, 1:days(i));
mse(1, i) = immse(pdf4_c, cases_n);

pd4_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf4_d = pdf(pd4_d, 1:days(i));
mse(2, i) = immse(pdf4_d, deaths_n);

%% 5. Goodness-of-fit: France.
i = i + 1;
x = France(1, :);
y = France(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd5_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf5_c = pdf(pd5_c, 1:days(i));
mse(1, i) = immse(pdf5_c, cases_n);

pd5_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf5_d = pdf(pd5_d, 1:days(i));
mse(2, i) = immse(pdf5_d, deaths_n);

%% 6. Goodness-of-fit: Germany.
i = i + 1;
x = Germany(1, :);
y = Germany(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd6_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf6_c = pdf(pd6_c, 1:days(i));
mse(1, i) = immse(pdf6_c, cases_n);

pd6_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf6_d = pdf(pd6_d, 1:days(i));
mse(2, i) = immse(pdf6_d, deaths_n);

%% 7. Goodness-of-fit: Switzerland.
i = i + 1;
x = Switzerland(1, :);
y = Switzerland(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd7_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf7_c = pdf(pd7_c, 1:days(i));
mse(1, i) = immse(pdf7_c, cases_n);

pd7_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf7_d = pdf(pd7_d, 1:days(i));
mse(2, i) = immse(pdf7_d, deaths_n);

%% 8. Goodness-of-fit: Norway.
i = i + 1;
x = Norway(1, :);
y = Norway(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd8_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf8_c = pdf(pd8_c, 1:days(i));
mse(1, i) = immse(pdf8_c, cases_n);

pd8_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf8_d = pdf(pd8_d, 1:days(i));
mse(2, i) = immse(pdf8_d, deaths_n);

%% 9. Goodness-of-fit: Italy.
i = i + 1;
x = Italy(1, :);
y = Italy(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd9_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf9_c = pdf(pd9_c, 1:days(i));
mse(1, i) = immse(pdf9_c, cases_n);

pd9_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf9_d = pdf(pd9_d, 1:days(i));
mse(2, i) = immse(pdf9_d, deaths_n);

%% 10. Goodness-of-fit: Ireland.
i = i + 1;
x = Ireland(1, :);
y = Ireland(2, :);

days(i) = length(x);
cases_n = x ./ (sum(abs(x)));
deaths_n = y ./ (sum(abs(y)));

pd10_c = fitdist((1:days(i))' , dist, 'Frequency', x);
pdf10_c = pdf(pd10_c, 1:days(i));
mse(1, i) = immse(pdf10_c, cases_n);

pd10_d = fitdist((1:days(i))' , dist, 'Frequency', y);
pdf10_d = pdf(pd10_d, 1:days(i));
mse(2, i) = immse(pdf10_d, deaths_n);

%% Sort Daily Cases Mean Squared Error in ascending order.
[~, index_c] = sort(mse(1, :));
fprintf('Fitness of %s distribution on Daily Cases (Countries shown from best to worst order):\n', dist)
disp(countries(index_c))

% Sort Daily Deaths Mean Squared Error in ascending order.
[~, index_d] = sort(mse(2, :));
fprintf('\nFitness of %s distribution on Daily Deaths (Countries shown from best to worst order):\n', dist)
disp(countries(index_d))


%% Comments:
% Daily Cases:
% Gamma distribution seems to fit the same or even better to the 10 countries selected
% than country A Austria, as shown from MSE matrix (first row).

% Daily Deaths:
% Gamma distribution seems to fit better to almost half of the 10 countries selected
% than country A Austria, while the rest got a higher MSE (second row).


