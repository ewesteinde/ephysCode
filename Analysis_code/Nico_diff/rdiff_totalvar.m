function x1 = rdiff_totalvar(x, t, sigma, strategy)
% RDIFF_TOTALVAR   Regularised numerical differentiation using method based
%                  on constraining the 1-norm of the first derivative.
%    X1 = RDIFF_TOTALVAR(X, T) returns an estimate of the first derivative
%    of X; if T is a scalar, it is used as spacing between abscissa; if T
%    is a vector, its elements are used as the abscissa. If X is a matrix,
%    RDIFF_TOTALVAR differentiates each of its columns independently.
%    
%    X1 = RDIFF_TOTALVAR(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF_TOTALVAR(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_TOTALVAR(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the constraint on the 1-norm of the first derivative
%    selected among the following ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'ncp'     based on normalised cumulative periodogram
%    By default, the 'ncp' strategy is used.
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
if ~exist('t', 'var') || isempty(t), t = (1:N)'; end
if isscalar(t), t = t*(0:N-1)'; end
if ~isvector(t) || length(t) ~= N
    error('for size(x) = [%d, %d], t should be a scalar or a vector of length %d; size(t) = [%d, %d]', size(x,1), size(x,2), N, size(t,1), size(t,2))
end
if isrow(t); t = t'; end

% Select default strategy
if ~exist('strategy', 'var') || isempty(strategy)
    strategy = 'ncp';
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

% Get progressive-difference matrix
D = pdmatrix(t, 1);

% Considered values of regularisation parameter
as = logspace(-5, -1, 12);
Na = numel(as);

% Constants
e = 1e-8;
IMax = 100;
iStopThr = 5e-2;
switch lower(strategy)
    case 'dp'
        uDP = 2;
        trS = trace(sigma);
    case 'ncp'
        Nl = floor(N/2);
        l = (1:Nl)' / Nl;
    otherwise
        error('Unknown strategy: %s', strategy)
end

% Make sure first data point equals 0
xp = x - repmat(x(1,:), [N 1]);

% Estimates of the derivative
x1 = zeros(N, R);

for r = 1:R     % for each sequence of data to be differentiated
    xpr = xp(:,r);
    
    % Criterion for optimisation of regularisation parameter
    crit = zeros(Na, 1);
    
    % Estimates of the derivative for different values
    % of the regularisation parameter
    X1a = zeros(N, Na);
    
    % Criterion for stopping the iteration
    iStopFun = @(x1, a) norm(Q*x1 - xpr, 2) + a*norm(x1, 1);
    
    for na = 1:Na   % for each considered value of the reg. par.
        a = as(na);
        
        % Estimates of the derivative in subsequent iterations
        X1i = zeros(N, IMax+1);
        
        % Iterative implementation of constraint on total variation
        % based on:
        % - R. Chartrand, "Numerical differentiation of noisy, nonsmooth
        %   data", ISRN Applied Mathematics, vol. 2011, pp. 1-11, 2011
        % - C. R. Vogel and M. E. Oman, "Iterative methods for total
        %   variation denoising", SIAM Journal on Scientific Computing,
        %   vol. 17, no. 1, pp. 227-238, 1996
        for i = 1:IMax
            E = diag(1./sqrt((D*X1i(:,i)).^2 + e));
            L = D'*E*D;
            H = Q'*Q + a*L;
            g = Q'*(Q*X1i(:,i) - xpr) + a*L*X1i(:,i);
            s = H \ -g;
            X1i(:,i+1) = X1i(:,i) + s;
            if (iStopFun(X1i(:,i), a) - iStopFun(X1i(:,i+1),a)) / norm(s) <= iStopThr
                break;
            end
        end
        X1a(:,na) = X1i(:,i+1);
        
        % Vector of approximation residuals
        residual = Q*X1a(:,na) - xpr;
        
        % Compute criteria for optimisation of regularisation parameter
        switch lower(strategy)
            case 'dp'
                crit(na) = abs(sum(residual.^2, 1) - uDP*trS);
            case 'ncp'
                rf = fft(residual);
                rp = abs(rf(1:Nl)).^2;
                rncp = cumsum(rp) / sum(rp);
                crit(na) = sqrt(sum((rncp - l).^2));
        end
    end
    
    % Select optimum value of regularisation parameter
    naSel = find(crit == min(crit), 1, 'last');
    x1(:,r) = X1a(:, naSel);
end

% Return row vector if x was one
if exist('xrow', 'var') && xrow
    x1 = x1';
end