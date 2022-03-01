%% Optimal Estimation - Homework 2 - Tanner Koza

clear
clc
close all

%% Problem 3 - Part A & B

% Time Initialization
dt = 0.1;
t_end = 30;
t = 0:dt:t_end;
numSamp = length(t);

% Monte Carlo Initialization
numSims = 10000;

% Noise & Frequency Initialization
sigma = 0.3; % deg/s
freq = 2;
omega = freq * (2 * pi) ; % rads/s

% Arbitrary Coefficient Initialization
a = 3;
b = 10;

% Least Squares Initialization
estSamps = 10; % # of samples used in estimate
R = sigma * eye(estSamps);

% Preallocation
r = zeros(numSamp,1); 
g = zeros(numSamp,1);
est = zeros(numSims,2);


for i = 1:numSims

    n = sigma * randn(numSamp,1);

    for k = 1:numSamp
       
        r(k) = 100 * sin(omega * t(k));
    
        g(k) = a * r(k) + b + n(k); % degs/s
    
    end

H = [r(1:estSamps) ones(estSamps,1)];
est(i,:) = (H' * R^-1 * H)^-1 * H' * R^-1 * g(1:estSamps);

P = (H' * R^-1 * H)^-1 ;
end

mean_est = mean(est);
std_est = std(est); % Monte Carlo Standard Deviation

std_a = sqrt(P(1,1)); % Theoretical Standard Deviation
std_b = sqrt(P(2,2));

figure
plot(1:numSims, est(:,1), '*')
hold on
yline(mean_est(1),'k','LineWidth',3)
yline(mean_est(1) + std_est(1),'r','LineWidth',3)
yline(mean_est(1) - std_est(1),'r','LineWidth',3)
title('Monte Carlo a Estimate (10 Samples)')
xlabel('Simulation #')
ylabel('Coefficient Estimate')

figure
plot(1:numSims, est(:,2), '*')
hold on
yline(mean_est(2),'k','LineWidth',3)
yline(mean_est(2) + std_est(2),'r','LineWidth',3)
yline(mean_est(2) - std_est(2),'r','LineWidth',3)
title('Monte Carlo b Estimate (10 Samples)')
xlabel('Simulation #')
ylabel('Coefficient Estimate')

clearvars

%% Problem 3 - Part C

% Time Initialization
dt = 0.1;
t_end = 300;
t = 0:dt:t_end;
numSamp = length(t);

% Monte Carlo Initialization
numSims = 1000;

% Noise & Frequency Initialization
sigma = 0.3; % deg/s
freq = 2;
omega = freq * (2 * pi) ; % rads/s

% Arbitrary Coefficient Initialization
a = 3;
b = 10;

% Least Squares Initialization
estSamps = 1000; % # of samples used in estimate

% Preallocation
r = zeros(numSamp,1); 
g = zeros(numSamp,1);
est = zeros(numSims,2);


for i = 1:numSims

    n = sigma * randn(numSamp,1);

    for k = 1:numSamp
       
        r(k) = 100 * sind(omega * t(k));
    
        g(k) = a * r(k) + b + n(k);
    
    end

H = [r(1:estSamps) ones(estSamps,1)];
est(i,:) = (H' * H)^-1 * H' * g(1:estSamps);

end

mean_est = mean(est);
std_est = std(est);

figure
plot(1:numSims, est(:,1), '*')
hold on
yline(mean_est(1),'k','LineWidth',3)
yline(mean_est(1) + std_est(1),'r','LineWidth',3)
yline(mean_est(1) - std_est(1),'r','LineWidth',3)
title('Monte Carlo a Estimate (1000 Samples)')
xlabel('Simulation #')
ylabel('Coefficient Estimate')

figure
plot(1:numSims, est(:,2), '*')
hold on
yline(mean_est(2),'k','LineWidth',3)
yline(mean_est(2) + std_est(2),'r','LineWidth',3)
yline(mean_est(2) - std_est(2),'r','LineWidth',3)
title('Monte Carlo b Estimate (1000 Samples)')
xlabel('Simulation #')
ylabel('Coefficient Estimate')

clearvars

%% Problem 3 - Part D

%TODO: MODIFY FOR RLS ALGORITHIM

% Time Initialization
dt = 0.1;
t_end = 30;
t = 0:dt:t_end;
numSamp = length(t);

% Monte Carlo Initialization
numSims = 1000;

% Noise & Frequency Initialization
sigma = 0.3; % deg/s
omega = 450; % degs/s

% Arbitrary Coefficient Initialization
a = 3;
b = 10;

% Least Squares Initialization
estSamps = 10; % # of samples used in estimate

% Preallocation
r = zeros(numSamp,1); 
g = zeros(numSamp,1);
est = zeros(numSims,2);


for i = 1:numSims

    n = sigma * randn(numSamp,1);

    for k = 1:numSamp
       
        r(k) = 100 * sind(omega * t(k));
    
        g(k) = a * r(k) + b + n(k);
    
    end

H = [r(1:estSamps) ones(estSamps,1)];
est(i,:) = (H' * H)^-1 * H' * g(1:estSamps);

end

mean_est = mean(est);
std_est = std(est);