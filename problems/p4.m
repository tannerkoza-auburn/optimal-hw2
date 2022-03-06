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

sigma = 0.01;
noise = sigma * randn(1000,1);
Y_ = y_ + noise;
Y = y + noise;

H = [u(2:end-1) -u(1:end-2) y(2:end-1) -y(1:end-2)];

est = (H' * H)^-1 * H' * y(3:end);