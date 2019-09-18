%CONSUMP_SAMPLES Sample for consumption and bodyweight.
%
% Define empirical distributions for body weight and lettuce consumption to
% sample from for risk estimations. Saves the results in a file called
% `consumpsamp2.mat`.

clear

n = 100000;
% Values from Kahn and Stralka 2009, Table 6 (last row, 'All Ages')
percentiles = [0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99];
values = [8, 15, 22, 52, 67, 81, 95, 104, 122];
% Together, the percentile and weight values can be used to define a
% cumulative distribution function (CDF). Sampling from a CDF is carried
% out by first generating a large number of uniformly distributed random
% numbers in [0,1]. The body weight samples are obtained by interpolating
% the inverse CDF function (percentiles on X axis, body weight on Y axis)
% at these samples.
xx = rand(n, 1);
% Sample from the (inverse) cumulative distribution function.
% Shape-preserving piecewise cubic interpolation used.
bwsample = interp1(percentiles, values, xx, 'phcip'); % (kg bodywt)

% Get the per capita consumption samples
% This is from CSFII survey, page 133 (Table B-10a. Consumer-only intake of
% lettuce (g/(kg bodywt*day) as consumed)
% Using the same argument as for body weight.
percentiles = [0, 0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1];
values = [0, 0.026, 0.045, 0.063, 0.138, 0.379, 0.722, 1.213, 1.617, 2.798, 7.395];
perconssample = interp1(percentiles, values, xx, 'phcip'); % (kg bodywt)
perconssample = perconssample(randperm(n, n));

% Multiply both to get consumption distribution
consumpsamp = bwsample .* perconssample;
%% Save results
save('Data/consumpsamp2.mat', 'consumpsamp');
