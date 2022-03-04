%% Optimal Estimation - Homework 2 - Problem 4

clear
clc

%% Part A

numd = 0.25 * [1 -0.8];
dend = [1 -1.9 0.95];

u = randn(1000,1);
numSamps = length(u);

y = dlsim(numd, dend, u);

y_ = zeros(numSamps,1);

A = -0.25;
B = -A * 0.8;
C = -1.90; 
D = 0.95;

% for i = 3:numSamps
% 
%     y_(i) = A * u(i-1) + B * u(i-2) - C * y_(i-1) - D * y_(i-2);
% 
% end

sigma = 0.01;
noise = sigma * randn(1000,1);
Y_ = y_ + noise;
Y = y + noise;