%=====================================================================

% Purpose: To generate the sequence of observations Z

% Arguments
% Input(s):
% M - number of samples of recieved signal
% d - distance between target and each sensor
% mu - threshold to determine hypothesis, i.e.,
% target present or absent
% Input(s):
% result - vector of observations where each element
% of the vector is 1 (if energy > mu) or 0
%
% Assumptions: vector 'd', which is an input vector should
% have no value less than 0
%
% Function declaration:
% function result = Z_k(M,d,mu)
%
%=====================================================================
function result = Z_k(M,d,mu)
len_d = length(d);
noise_mean = zeros(1,M);
signal_mean = zeros(1,M);
mean_H1 = noise_mean + signal_mean;
noise_var = 10^(-4);
for i = 1:len_d;
if d(i) > 0
signal_var = 1/d(i)^2;
var_H1 = noise_var + signal_var;
rho = mean_H1 + randn(1,M) * sqrt(var_H1);
energy = abs(rho * rho');
if energy > mu
result(i) = 1;
else
result(i) = 0;
end
elseif d(i) == 0
result(i) = 1;
else
error('Distance from sensor to target is less than zero!!!');
end
end
