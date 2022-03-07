%% Optimal Estimation - Homework 2 - Problem 3

clear
clc

%% Parts A & B

% Time Initialization
dt = 0.1;
t_end = 30;
t = 0:dt:t_end;
numSamps = length(t);

% Monte Carlo Initialization
numSims = 1000;

% Noise & Frequency Initialization
sigma = 0.3; % deg/s
var = sigma^2;
freq = 2;
omega = freq * (2 * pi) ; % rads/s

% Arbitrary Coefficient Initialization
a = 3;
b = 10;

% Least Squares Initialization
estSamps = 10; % # of samples used in estimate
R = var * eye(estSamps);

% Preallocation
r = zeros(numSamps,1); 
g = zeros(numSamps,1);
est = zeros(numSims,2);


for i = 1:numSims

    n = sigma * randn(numSamps,1);

    for k = 1:numSamps
       
        r(k) = 100 * sin(omega * t(k));
    
        g(k) = a * r(k) + b + n(k); % degs/s
    
    end

H = [r(1:estSamps) ones(estSamps,1)];
est(i,:) = (H' * H)^-1 * H' * g(1:estSamps);

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

%% Part C

% Time Initialization
dt = 0.1;
t_end = 300;
t = 0:dt:t_end;
numSamps = length(t);

% Monte Carlo Initialization
numSims = 1000;

% Noise & Frequency Initialization
sigma = 0.3; % deg/s
var = sigma^2;
freq = 2;
omega = freq * (2 * pi) ; % rads/s

% Arbitrary Coefficient Initialization
a = 3;
b = 10;

% Least Squares Initialization
estSamps = 1000; % # of samples used in estimate
R = var * eye(estSamps);

% Preallocation
r = zeros(numSamps,1); 
g = zeros(numSamps,1);
est = zeros(numSims,2);


for i = 1:numSims

    n = sigma * randn(numSamps,1);

    for k = 1:numSamps
       
        r(k) = 100 * sin(omega * t(k));
    
        g(k) = a * r(k) + b + n(k); % degs/s
    
    end

H = [r(1:estSamps) ones(estSamps,1)];
est(i,:) = (H' * H)^-1 * H' * g(1:estSamps);

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

%% Part D

% Time Initialization
dt = 0.1;
t_end = 30; % s
t = 0:dt:t_end;
numSamps = length(t);

% Monte Carlo Initialization
numSims = 1000;

% Noise & Frequency Initialization
sigma = 0.3; % deg/s
var = sigma^2;
freq = 2;
omega = freq * (2 * pi) ; % rads/s

% Arbitrary Coefficient Initialization
a = 3;
b = 10;

% Recursive Least Squares Initialization
est = zeros(2,1);
R = var;
P = eye(2);

% Preallocation
est_a = zeros(numSamps,numSims);
est_b = zeros(numSamps,numSims);
std_a = zeros(numSamps,numSims);
std_b = zeros(numSamps,numSims);


for i = 1:numSims

    n = sigma * randn(numSamps,1);
    r = 100 * sin(omega .* t);
    g = a * r' + b + n;

    for j = 1:numSamps

        H = [r(j) 1];

        K = P * H' * (H * P * H' + R)^ -1;

        est =  est + K * (g(j) - H * est);

        est_a(j,i) = est(1);
        est_b(j,i) = est(2);

        P = (P^-1 + H' * R^-1 * H)^ -1;

        std_a(j,i) = sqrt(P(1,1));
        std_b(j,i) = sqrt(P(2,2)); 

    end

end

mean_est_a = mean(est_a,2);
mean_est_b = mean(est_b,2);

th_std_a = mean(std_a,2);
th_std_b = mean(std_b,2);

mc_std_a = std(est_a,0,2);
mc_std_b = std(est_b,0,2);

figure
plot(t,mean_est_a)
hold on
plot(t,mean_est_b)
title('RLS Coefficient Estimate')
xlabel('Time (s)')
ylabel('Coefficient Estimate')
legend('a Estimate','b Estimate')
ylim([-inf 11])

figure
ax(1) = subplot(2,1,1);
plot(t,th_std_a, 'Parent', ax(1))
hold on
plot(t,th_std_b, 'Parent', ax(1))
title('Standard Deviation of Estimate Error')
xlabel('Time (s)')
ylabel('Standard Deviation')
legend('a Error Std.','b Error Std.')

ax(2) = subplot(2,1,2);
plot(t,abs(a - mean_est_a),'Parent', ax(2))
hold on
plot(t,abs(b - mean_est_b), 'Parent', ax(2))
title('Actual Error of Estimate')
xlabel('Time (s)')
ylabel('Error (Residual)')
legend('a Error','b Error')

linkaxes(ax,'y');
