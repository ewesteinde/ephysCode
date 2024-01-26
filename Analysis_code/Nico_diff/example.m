% Examplary usages of RDIFF, function implementing various methods
% for regularised numerical differentiation.
% 
% Jakub Wagner, February 2, 2020
% Institute of Radioelectronics and Multimedia Technology
% Warsaw University of Technology

% =========================================================================
% Example #1: Differentiation of a single sequence of data

% Test function and its derivative
f = @(t) 1 + 2/3*exp(-((t-3)/3).^6).^12 + (t.^3)/30;
f1 = @(t) -16*exp(-((t-3)/3).^6).^12 .* ((t-3)/3).^5 + (t.^2)/10;

% Generation of synthetic data
sigma = 0.0001;                 % variance of errors corrupting the data
N = 51;                         % number of data points
t = linspace(0, 3, N)';         % abscissa
e = sqrt(sigma)*randn(N,1);     % measurement errors
x = f(t) + e;                   % sequence of data to be differentiated

% Differentiation
x1 = rdiff(x, t);

% Results
figure(1); clf; hold on
plot(t, f(t), 'Color', [0 0 0])
plot(t, f1(t), 'Color', [0 0.4 0.7])
plot(t, x, 'o:', 'MarkerSize', 3, 'Color', [0 0 0])
plot(t, x1, 'o:', 'MarkerSize', 3, 'Color', [0 0.4 0.7])
hold off; grid on; box on; ylim([-0.5 3]); xlabel('\itt')
legend('f({\itt})', 'f^{(1)}({\itt})', 'Data', 'Estimates of derivative', 'Location', 'NorthWest')

% =========================================================================
% Example #2: Differentiation of several sequences of data simultaneously

% Generation of synthetic data
R = 10;                         % number of sequences
E = sqrt(sigma)*randn(N,R);     % measurement errors
X = repmat(f(t), [1 R]) + E;    % sequences of data to be differentiated

% Differentiation
X1 = rdiff(X, t);

% Results
figure(2); clf; hold on
plot(t, f1(t), 'Color', [0 0 0])
plot(t, X1, 'o:', 'MarkerSize', 3, 'Color', [0 0.4 0.7])
hold off; grid on; box on; ylim([-0.5, 1.5])
legend('f^{(1)}({\itt})', 'Estimates of derivative', 'Location', 'NorthWest')

% =========================================================================
% Example #3: Comparison of performance of various methods

% Compared methods
methods = {'totalvar', 'tikhonov', 'central', 'landweber', 'nbcos'};
NM = numel(methods);

% Reference values of the derivative
X1Ref = repmat(f1(t), [1 R]);

% Signal-to-noise ratio in the estimates of the derivative
snr1 = zeros(R,NM);

for nm = 1:NM   % for each method
    
    % Differentiation
    X1 = rdiff(X, t, methods{nm}, sigma);
    
    % Signal-to-noise ratio
    snr1(:,nm) = 10*log10(sum(X1Ref.^2) ./ sum((X1 - X1Ref).^2));
end

% Results
figure(3); clf
plot(1:NM, snr1', '*', 'MarkerSize', 6, 'Color', [0 0.4 0.7])
grid on; xlim([0.5 NM+0.5]); ylabel('Signal-to-noise ratio');
ax = gca; ax.XTick = 1:NM; ax.XTickLabel = methods;

% =========================================================================
% Example #4: Estimate derivative at different abscissa

% Abscissa at which derivative is to be estimated
t1 = linspace(0,3,14);

% Differentiation
x1 = rdiff(x, t, 'nbcos', sigma, 'gcv', t1);

% Results
figure(4); clf; hold on
plot(t, f(t), 'Color', [0 0 0])
plot(t, f1(t), 'Color', [0 0.4 0.7])
plot(t, x, 'o:', 'MarkerSize', 3, 'Color', [0 0 0])
plot(t1, x1, 'o:', 'MarkerSize', 3, 'Color', [0 0.4 0.7])
hold off; grid on; box on; ylim([-0.5 3]); xlabel('\itt')
legend('f({\itt})', 'f^{(1)}({\itt})', 'Data', 'Estimates of derivative', 'Location', 'NorthWest')