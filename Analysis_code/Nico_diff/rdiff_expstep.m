function x1 = rdiff_expstep(x, t, sigma, strategy)
% RDIFF_EXPSTEP   Regularised numerical differentiation using method
%                 based on modified central-difference formula with
%                 extended differentiation step; in contrary to standard
%                 central-difference formula, here subsequences of data
%                 are approximated using functions of following form:
%                    f(t) = exp(p0 + p1*t + p2*t^2) + p3
%                 instead of 2nd degree algebraic polynomials.
%    X1 = RDIFF_EXPSTEP(X, T) returns an estimate of the first derivative
%    of the N-dimensional vector X; if T is a scalar, it is used as spacing
%    between abscissa; if T is an N-dimensional vector, its elements are
%    used as the abscissa. This method is based on the assumption that the
%    abscissa are equidistant. It is required that N >= 5. If X is a
%    matrix, RDIFF_EXPSTEP differentiates each column of X independently.
%    
%    X1 = RDIFF_EXPSTEP(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF_EXPSTEP(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_EXPSTEP(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the differentiation step selected among the following
%    ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'ncp'     based on normalised cumulative periodogram
%    By default, the 'dp' strategy is used if SIGMA is provided, and the
%    'ncp' strategy otherwise.
%    
%    Jakub Wagner, February 2, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

% Make sure data are positive
if any(x <= 0)
    x = x - min(x) + 1;
end

% Check size of x
if isrow(x)
    x = x';
    xrow = true;
end
[N, R] = size(x);
if N < 5
    error('At least 5 data points are necessary for ''ce3ds'' method (given %d)', N)
end

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
        if numel(sigma) ~= 1 && any(size(sigma) ~= N)
            error('Wrong size of sigma: [%d, %d] (should be a scalar or a %dx%d matrix)', size(sigma, 1), size(sigma, 2), N, N)
        end
    end
end

% Initialise vector/matrix for estimates of the derivative
x1 = zeros(N, R);

% Considered values of differentiation step
K = min([120, floor((N-1)/2)]);
ks = 2:K;
Nk = numel(ks);

% Sampling period
dt = (t(N) - t(1)) / (N-1);

% Check sigma
if strcmpi(strategy, 'dp')
    if isempty(sigma)
        error('%s strategy requires Sigma', strategy)
    else
        if numel(sigma) ~= 1 && any(size(sigma) ~= N)
            error('Wrong size of Sigma: [%d, %d] (should be [1, 1] or [%d, %d])', size(sigma, 1), size(sigma, 2), N, N)
        end
    end
end

% Criterion values for optimisation of differentiation step
crit = inf(Nk, 1);

% Scaling factor for discrepancy principle
uDP = 2;

for r = 1:R     % for each differentiated sequence of data
    for n = 1:N     % for each time instant
        
        % Estimates of the derivative corresponding to different values
        % of the differentiation step
        x1n = zeros(Nk,1);
        
        for nk = 1:Nk   % for each considered value of diff. step
            k = ks(nk);
            
            % Approximate the data and estimate the derivative
            if n < k+1
                in = 1:n+k;         % indices of approximated data
                xn = x(in, r);      % approximated subsequence of data
                                    % parameters of approximating function
                p0 = log(xn(n));
                if xn(n+k)/xn(n) > exp(1)^0.5
                    p1 = log(xn(n+k)/xn(n))/k/dt;
                    p2 = 0;
                else
                    p1(1) = (1+sqrt(1-2*log(xn(n+k)/xn(n))))/k/dt;
                    p1(2) = (1-sqrt(1-2*log(xn(n+k)/xn(n))))/k/dt;
                    [~, ip1] = min(abs(p1));
                    p1 = p1(ip1);
                    p2 = -p1^2 / 2;
                end
                x1n(nk) = xn(n) * p1;   % estimate of derivative
                                        % result of approximation
                xa = exp(p0 + p1*(t(in)-t(n)) + p2*(t(in)-t(n)).^2);
                
            elseif n <= N-k
                in = n-k:n+k;       % indices of approximated data
                xn = x(in, r);      % approximated subsequence of data
                                    % parameters of approximating function
                p0 = log(xn(k+1));
                p1 = 1/2/k/dt * log(xn(2*k+1)/xn(1));
                p2 = 1/2/k^2/dt^2 * log(xn(2*k+1)*xn(1)/xn(k+1)^2);
                x1n(nk) = xn(k+1) * p1; % estimate of derivative
                                        % result of approximation
                xa = exp(p0 + p1*(t(in)-t(n)) + p2*(t(in)-t(n)).^2);
                
            elseif n > N-k
                in = n-k:N;         % indices of approximated data
                xn = x(in, r);      % approximated subsequence of data
                                    % parameters of approximating function
                p0 = log(xn(k+1));
                if xn(k+1)/xn(1) < exp(1)^-0.5
                    p1 = log(xn(k+1)/xn(1))/k/dt;
                    p2 = 0;
                else
                    p1(1) = (-1+sqrt(1+2*log(xn(k+1)/xn(1))))/k/dt;
                    p1(2) = (-1-sqrt(1+2*log(xn(k+1)/xn(1))))/k/dt;
                    [~, ip1] = min(abs(p1));
                    p1 = p1(ip1);
                    p2 = -p1^2 / 2;
                end
                x1n(nk) = xn(k+1) * p1; % estimate of derivative
                                        % result of approximation
                xa = exp(p0 + p1*(t(in)-t(n)) + p2*(t(in)-t(n)).^2);
            end
            
            % Length of analysed sequence for current value of k
            Ni = length(in);
            
            % Compute criteria for optimisation of k
            switch lower(strategy)
                case 'dp'
                    if numel(sigma) == 1
                        crit(nk) = abs(sum((xa - xn).^2) - uDP*Ni*sigma) / Ni;
                    else
                        crit(nk) = abs(sum((xa - xn).^2) - uDP*trace(sigma(in, in))) / Ni;
                    end
                case 'ncp'
                    L = floor(Ni/2);
                    l = (1:L)'/L;
                    rf = fft(xa - xn);
                    rp = abs(rf(1:L)).^2;
                    rncp = cumsum(rp)/sum(rp);
                    crit(nk) = sum((rncp-l).^2)/L;
                otherwise
                    error('Unknown strategy: %s', strategy)
            end
        end
        
        % Select estimate of derivative corresponding to minimum of criterion
        nkSel = find(crit == min(crit), 1, 'last');
        x1(n, r) = x1n(nkSel);
    end
end

% Return row vector if x was one
if exist('xrow', 'var') && xrow
    x1 = x1';
end