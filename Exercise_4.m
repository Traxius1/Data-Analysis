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

% Calculate the maximum correlation between Daily Cases and Deaths.
% If x <- Daily Cases, y <- Daily Deaths, find the offset t that maximizes
% Pearson coefficient r(x(n)*y(n+t)) and t in [0, 20], n in [1, data's columns].
% It assumed that Daily Cases flow must be past Daily Deaths.

% In this case, matrix r is a matrix with columns for every 20 values of
% Pearson coeff. and every row contains the results for a different
% country.
r = zeros(6, 21);
n = 21;

for t = 0:n-1
    
% A. Pearson coefficient for Austria.
    [x, y] = offset_data_fun(Austria(1, :) ,Austria(2, :), t);
    corr_m = corrcoef(x, y);
    r(1, t+1) = corr_m(1, 2);

% 1. Pearson coefficient for Finland.
    [x, y] = offset_data_fun(Finland(1, :), Finland(2, :), t);
    corr_m = corrcoef(x, y);
    r(2, t+1) = corr_m(1, 2);

% 2. Pearson coefficient for Spain.
    [x, y] = offset_data_fun(Spain(1, :) ,Spain(2, :), t);
    corr_m = corrcoef(x, y);
    r(3, t+1) = corr_m(1, 2);

% 3. Pearson coefficient for Germany.
    [x, y] = offset_data_fun(Germany(1, :), Germany(2, :), t);
    corr_m = corrcoef(x, y);
    r(4, t+1) = corr_m(1, 2);

% 4. Pearson coefficient for Portugal.
    [x, y] = offset_data_fun(Portugal(1, :), Portugal(2, :), t);
    corr_m = corrcoef(x, y);
    r(5, t+1) = corr_m(1, 2);

% 5. Pearson coefficient for France.
    [x, y] = offset_data_fun(France(1, :), France(2, :), t);
    corr_m = corrcoef(x, y);
    r(6, t+1) = corr_m(1, 2);
    
end

% r_max containts the max Pearson coefficient for every country.
% The offset that maximizes Pearson coeff. is (index-1).
% x and y are row vectors.
r_max = zeros(1, 6);
offset = zeros(1, 6);

for i = 1:size(r, 1)
    
    [r_max(i), index] = max(r(i, :));
    offset(i) = index - 1;
    
end

%% Plots
% A. Austria.
figure
clf
plot(0:n-1, r(1, :), 'o-')
hold on
ax = axis;
plt(1) = plot(7.7*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plot(14.4*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plt(2) = plot(offset(1)*[1, 1], [ax(3), r_max(1)], 'k--');
title({'Austria'; 'Pearson coeff. values for offset in [0, 20]'})
ylabel('Pearson Coefficient')
xlabel('Offset')
legend(plt, 'Offset confidence interval', 'Max Pearson coeff')

% 1. Finland.
figure
clf
plot(0:n-1, r(2, :), 'o-')
hold on
ax = axis;
plt(1) = plot(7.7*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plot(14.4*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plt(2) = plot(offset(2)*[1, 1], [ax(3), r_max(2)], 'k--');
title({'Finland'; 'Pearson coeff. values for offset in [0, 20]'})
ylabel('Pearson Coefficient')
xlabel('Offset')
legend(plt, 'Offset confidence interval', 'Max Pearson coeff')

% 2. Spain.
figure
clf
plot(0:n-1, r(3, :), 'o-')
hold on
ax = axis;
plt(1) = plot(7.7*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plot(14.4*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plt(2) = plot(offset(3)*[1, 1], [ax(3), r_max(3)], 'k--');
title({'Spain'; 'Pearson coeff. values for offset in [0, 20]'})
ylabel('Pearson Coefficient')
xlabel('Offset')
legend(plt, 'Offset confidence interval', 'Max Pearson coeff')

% 3. Germany.
figure
clf
plot(0:n-1, r(4, :), 'o-')
hold on
ax = axis;
plt(1) = plot(7.7*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plot(14.4*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plt(2) = plot(offset(4)*[1, 1], [ax(3), r_max(4)], 'k--');
title({'Spain'; 'Pearson coeff. values for offset in [0, 20]'})
ylabel('Pearson Coefficient')
xlabel('Offset')
legend(plt, 'Offset confidence interval', 'Max Pearson coeff')

% 4. Portugal.
figure
clf
plot(0:n-1, r(5, :), 'o-')
hold on
ax = axis;
plt(1) = plot(7.7*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plot(14.4*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plt(2) = plot(offset(5)*[1, 1], [ax(3), r_max(5)], 'k--');
title({'Portugal'; 'Pearson coeff. values for offset in [0, 20]'})
ylabel('Pearson Coefficient')
xlabel('Offset')
legend(plt, 'Offset confidence interval', 'Max Pearson coeff')

% 5. France.
figure
clf
plot(0:n-1, r(6, :), 'o-')
hold on
ax = axis;
plt(1) = plot(7.7*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plot(14.4*[1, 1], [ax(3), ax(4)], 'r--');
hold on
plt(2) = plot(offset(6)*[1, 1], [ax(3), r_max(4)], 'k--');
title({'France'; 'Pearson coeff. values for offset in [0, 20]'})
ylabel('Pearson Coefficient')
xlabel('Offset')
legend(plt, 'Offset confidence interval', 'Max Pearson coeff')

% Print offset used for max Pearson Coeff.
fprintf('A. Austria\tOffset: %d, \tPearson Coeff.: %1.3f\n', offset(1), r_max(1))
fprintf('1. Finland\tOffset: %d, \tPearson Coeff.: %1.3f\n', offset(2), r_max(2))
fprintf('2. Spain\tOffset: %d, \t\tPearson Coeff.: %1.3f\n', offset(3), r_max(3))
fprintf('3. Germany\tOffset: %d, \tPearson Coeff.: %1.3f\n', offset(4), r_max(4))
fprintf('4. Portugal\tOffset: %d, \t\tPearson Coeff.: %1.3f\n', offset(5), r_max(5))
fprintf('5. France\tOffset: %d, \t\tPearson Coeff.: %1.3f\n', offset(6), r_max(6))


%% Comments.
% The only countries which show a kind of linear relation of Daily Cases -
% Deaths are Spain and Germany, using Offset of 5 and 13 respectively.
% The other countries can't have Daily Cases - Deaths linearly related.

% Offsets calculated in this exercise don't really support Exercise 3
% conclusion. Only 3/6 are in the confidence interval, and there wasn't an
% offset equal to 14.
% This happens because confidence interval calculated came from a really small
% sample (10 values), so it's not really precise.


%% This function generates the data set that comes from a certain value of
% t. t is the time offset between Daily Cases and Deaths.
function [x, y] = offset_data_fun(data_x, data_y, t)

% Shift y vector to the right by offset.
y = data_y(t+1 : end);
x = data_x(1:length(y));

end

