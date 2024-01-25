function x1 = rdiff_landweber(x, t, sigma, strategy, x1max)
% RDIFF_LANDWEBER   Regularised numerical differentiation using method
%                   based on Landweber's iterative algorithm.
%    X1 = RDIFF_LANDWEBER(X, T) returns an estimate of the first derivative
%    of X; if T is a scalar, it is used as spacing between abscissa; if T
%    is a vector, its elements are used as the abscissa. If X is a matrix,
%    RDIFF_LANDWEBER differentiates each of its columns independently.
%    
%    X1 = RDIFF_LANDWEBER(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF_LANDWEBER(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_LANDWEBER(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the number of iterations selected among the following
%    ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'ncp'     based on normalised cumulative periodogram
%    By default, the 'dp' strategy is used if SIGMA is provided, and the
%    'ncp' strategy otherwise.
%    
%    X1 = RDIFF_LANDWEBER(X, T, METHOD, SIGMA, STRATEGY, X1MAX) uses X1MAX
%    as an a priori upper bound for absolute value of the first derivative.
%    
%    Based on:
%    P. C. Hansen, "Discrete Inverse Problems: Insight and Algorithms",
%    Society for Industrial and Applied Mathematics, 2010 (Chapter 6).
% 
%    Jakub Wagner, February 2, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

% Check size of x
if isrow(x)
    x = x';
    xrow = true;
end
if isscalar(x)
    error('x must be a vector or a 2-D matrix.')
end
[N, R] = size(x);

% Check size of t
if isscalar(t), t = t*(0:N-1)'; end
if ~isvector(t) || length(t) ~= N
    error('for size(x) = [%d, %d], t should be a scalar or a vector of length %d; size(t) = [%d, %d]', size(x,1), size(x,2), N, size(t,1), size(t,2))
end
if isrow(t); t = t'; end

% Select default strategy
if ~exist('strategy', 'var') || isempty(strategy)
    if exist('sigma', 'var') && ~isempty(sigma)
        strategy = 'dp';
    else
        strategy = 'ncp';
    end
end

% Check availability and size of sigma
if strcmpi(strategy, 'dp')
    if ~exist('sigma', 'var') || isempty(sigma)
        warning('Cannot use ''%s'' strategy when sigma is not specified; using ''ncp'' strategy instead.', strategy)
        strategy = 'ncp';
    else
        if numel(sigma) == 1
            sigma = sigma * eye(N);
        elseif any(size(sigma) ~= N)
            error('Wrong size of sigma: [%d, %d] (should be scalar or %dx%d matrix)', size(sigma, 1), size(sigma, 2), N, N)
        end
    end
end

% Get quadrature matrix
Q = qmatrix(t);

% Modify matrix Q to add equation x'(1) = x'(2)
Q(1, 1:2) = [1, -1];

% Maximum number of iterations
I = 1500;

% Constants used during optimisation of number of iterations
switch lower(strategy)
    case 'dp'
        uDP = 2;
        trS = trace(sigma);
    case 'ncp'
        L = floor(N/2);
        l = repmat((1:L)' / L, [1 R]);
    otherwise
        error('Unknown strategy: %s', strategy)
end

% Criterion for optimising number of iterations
crit = inf(I+1, R);

% Estimates of the derivative
x1 = zeros(N,R);

% First singular value needed to determine parameter omega; svds is faster
% than svd for N larger than ~100.
if N > 100
    sv1 = svds(Q,1);
else
    [~, S] = svd(Q);
    sv1 = S(1,1);
end
omega = 0.99 * 2 / sv1^2;

% Optional constraint on the values of the derivative
if exist('x1max', 'var') && ~isempty(x1max)
    P = @(x) x .* (abs(x) <= x1max) + sign(x) .* x1max .* (abs(x) > x1max);
else
    P = @(x) x;
end

for r = 1:R     % for each sequence of data to be differentiated
    
    % Estimates of the derivative in subsequent iterations
    X1 = zeros(N, I+1);
    
    % Make sure first data point equals zero
    xp = x(:,r) - x(1,r);
    
    for i = 1:I     % in each iteration of the Landweber algorithm
        
        % Estimate the derivative
        X1(:,i+1) = P(X1(:,i) + omega * Q' * (xp - Q*X1(:,i)));
        
        % Approximate the data
        xa = Q * X1(:, i+1);
        
        % Compute criteria for optimisation of the number of iterations
        switch lower(strategy)
            case 'dp'
                crit(i+1,r) = abs(sum((xa - xp).^2, 1) - uDP*trS);
            case 'ncp'
                residual = xa - xp;
                rf = fft(residual);
                rp = abs(rf(1:L)).^2;
                rncp = cumsum(rp) / sum(rp);
                crit(i+1,r) = norm(rncp - l);
        end
    end
    
    % Select optimum number of iterations
    ISel = find(crit(:,r) == min(crit(:,r)), 1, 'first');
    x1(:,r) = X1(:,ISel);
end

% Return row vector if x was one
if exist('xrow', 'var') && xrow
    x1 = x1';
end