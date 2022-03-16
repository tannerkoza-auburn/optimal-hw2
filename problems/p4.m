%% Optimal Estimation - Homework 2 - Problem 4

clear
clc

%% Part A & B

numd = 0.25 * [1 -0.8];
dend = [1 -1.9 0.95];

dt = 1/100;

TF = tf(numd,dend,1);

u = randn(1000,1);
numSamps = length(u);

y = dlsim(numd, dend, u);

sigma = 0.01;
noise = sigma * randn(1000,1);
Y = y + noise;

H = [u(2:end-1) -u(1:end-2) Y(2:end-1) -Y(1:end-2)];

est = (H' * H)^-1 * H' * Y(3:end);

TFS = tf([est(1) -est(2)], [1 -est(3) est(4)],1);

figure
bode(TF)
hold on
bode(TFS)
legend('Simulated','SYS. ID')

% SNR 
snr = std(Y)/sigma;

figure
plot(y)
hold on
plot(Y)
title('System Output: Clean & Noisy')
xlabel('Samples')
ylabel('Magnitude')
legend('Clean','Noisy')

%% Part C

numSims = 10;
est = zeros(4,numSims);

for i = 1:numSims

    u = randn(1000,1);
    numSamps = length(u);

    y = dlsim(numd, dend, u);

    sigma = 0.01;
    noise = sigma * randn(1000,1);
    Y = y + noise;

    H = [u(2:end-1) -u(1:end-2) Y(2:end-1) -Y(1:end-2)];

    est(:,i) = (H' * H)^-1 * H' * Y(3:end);

end

P = sigma * (H' * H)^-1;

mean_est = mean(est,2);
std_est = std(est,0,2);

%% Part D

numd = 0.25 * [1 -0.8];
dend = [1 -1.9 0.95];

dt = 1/100;

TF = tf(numd,dend,1);

u = randn(1000,1);
numSamps = length(u);

y = dlsim(numd, dend, u);

sigma = 0.5;
noise = sigma * randn(1000,1);
Y = y + noise;

H = [u(2:end-1) -u(1:end-2) Y(2:end-1) -Y(1:end-2)];

est = (H' * H)^-1 * H' * Y(3:end);

TFS = tf([est(1) -est(2)], [1 -est(3) est(4)],1);

figure
bode(TF)
hold on
bode(TFS)
legend('Simulated','SYS. ID')

% SNR 
snr = std(Y)/sigma;

figure
plot(y)
hold on
plot(Y)
title('System Output: Clean & Noisy')
xlabel('Samples')
ylabel('Magnitude')
legend('Clean','Noisy')

numSims = 10;
est = zeros(4,numSims);

for i = 1:numSims

    u = randn(1000,1);
    numSamps = length(u);

    y = dlsim(numd, dend, u);

    sigma = 0.5;
    noise = sigma * randn(1000,1);
    Y = y + noise;

    H = [u(2:end-1) -u(1:end-2) Y(2:end-1) -Y(1:end-2)];

    est(:,i) = (H' * H)^-1 * H' * Y(3:end);

end

P = sigma * (H' * H)^-1;

mean_est = mean(est,2);
std_est = std(est,0,2);