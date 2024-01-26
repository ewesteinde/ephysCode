function x1 = rdiff_sinstep(x, t, sigma, strategy)
% RDIFF_SINSTEP   Regularised numerical differentiation using method based
%                 on modified central-difference formula with extended
%                 differentiation step; in contrary to standard
%                 central-difference formula, here subsequences of data
%                 are approximated using functions of following form:
%                    f(t) = p0 + p1*sin(t) + p2*cos(t);
%                 instead of 2nd degree algebraic polynomials.
%    X1 = RDIFF_SINSTEP(X, T) returns an estimate of the first derivative
%    of the N-dimensional vector X; if T is a scalar, it is used as spacing
%    between abscissa; if T is an N-dimensional vector, its elements are
%    used as the abscissa. This method is based on the assumption that the
%    abscissa are equidistant. It is required that N >= 5. If X is a
%    matrix, RDIFF_SINSTEP differentiates each column of X independently.
%    
%    X1 = RDIFF_SINSTEP(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF_SINSTEP(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_SINSTEP(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the differentiation step selected among the following
%    ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'gcv'     based on generalised cross validation
%       'ncp'     based on normalised cumulative periodogram
%       'sure'    based on Stein's unbiased risk estimator (requires SIGMA)
%    By default, the 'gcv' strategy is used.
%    
%    Jakub Wagner, February 2, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

% Check size of x
if isrow(x)
    x = x';
    xrow = true;
end
[N, R] = size(x);
if N < 5
    error('At least 5 data points are necessary for ''cs3ds'' method (given %d)', N)
end

% Check size of t
if isscalar(t), t = t*(0:N-1)'; end
if ~isvector(t) || length(t) ~= N
    error('for size(x) = [%d, %d], t should be a scalar or a vector of length %d; size(t) = [%d, %d]', size(x,1), size(x,2), N, size(t,1), size(t,2))
end
if isrow(t); t = t'; end

% Select default strategy
if ~exist('strategy', 'var') || isempty(strategy)
    strategy = 'gcv';
end

% Check availability and size of sigma
if any(strcmpi(strategy, {'dp', 'sure'}))
    if ~exist('sigma', 'var') || isempty(sigma)
        warning('Cannot use ''%s'' strategy when sigma is not specified; using ''gcv'' strategy instead.', strategy)
        strategy = 'gcv';
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
ks = unique(round(linspace(1, K, 30)));
Nk = numel(ks);

% Sampling period
dt = (t(N) - t(1)) / (N-1);

% Generate "hat" matrices
H0 = cell(Nk,1);
H0a = cell(Nk,Nk);
H0b = cell(Nk,Nk);
for nk = 1:Nk
    k = ks(nk);
    ct = cos((-k:k)*dt);
    st = sin((-k:k)*dt);
    H0{nk} = zeros(2*k+1, 2*k+1);
    H0{nk}(:,1) = (ct-1)/2/(cos(k*dt)-1) - st/2/sin(k*dt);
    H0{nk}(:,k+1) = (cos(k*dt)-ct)/(cos(k*dt)-1);
    H0{nk}(:,2*k+1) = (ct-1)/2/(cos(k*dt)-1) + st/2/sin(k*dt);
    for n = 1:k
        st = sin((-n+1:k)*dt);
        H0a{nk,n} = zeros(n+k,n+k);
        H0a{nk,n}(:,n) = 1-st/sin(k*dt);
        H0a{nk,n}(:,n+k) = st/sin(k*dt);
        u = N-k+n;
        st = sin((-k:N-u)'*dt);
        Ni = k+N-u+1;
        H0b{nk,n} = zeros(Ni,Ni);
        H0b{nk,n}(:,1) = -st/sin(k*dt);
        H0b{nk,n}(:,k+1) = 1+st/sin(k*dt);
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
        
        for nk = 1:Nk     % for each considered value of diff. step
            k = ks(nk);
            
            % Approximate the data and estimate the derivative
            if n < k+1
                in = 1:n+k;             % indices of approximated data
                xn = x(in, r);          % approximated subsequence of data
                H = H0a{nk, n};         % "hat" matrix
                xa = H * xn;            % result of approximation
                                        % estimate of derivative
                x1n(nk) = (xn(n+k)-xn(n))/sin(k*dt);
            elseif n <= N-k
                in = n-k:n+k;           % indices of approximated data
                xn = x(in, r);          % approximated subsequence of data
                H = H0{nk};             % "hat" matrix
                xa = H * xn;            % result of approximation
                                        % estimate of derivative
                x1n(nk) = (xn(end) - xn(1)) / 2 / sin(k*dt);
            elseif n > N-k
                in = n-k:N;             % indices of approximated data
                xn = x(in, r);          % approximated subsequence of data
                H = H0b{nk, k-N+n};     % "hat" matrix
                xa = H * xn;            % result of approximation
                                        % estimate of derivative
                x1n(nk) = (xn(k+1)-xn(1)) / sin(k*dt);
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
                case 'gcv'
                    crit(nk) = sum((xa - xn).^2) / trace(eye(Ni) - H)^2;
                case 'ncp'
                    L = floor(Ni/2);
                    l = (1:L)'/L;
                    rf = fft(xa - xn);
                    rp = abs(rf(1:L)).^2;
                    rncp = cumsum(rp)/sum(rp);
                    crit(nk) = sum((rncp-l).^2)/L;
                case 'sure'
                    if numel(sigma) == 1
                        crit(nk) = (sum(xa.^2) + sum(xn.^2) - Ni*sigma - 2*xn'*xa + 2*sigma*trace(H)) / Ni;
                    else
                        Sigma12 = sqrt(sigma(in, in));
                        crit(nk) = (sum(xa.^2) + sum(xn.^2) - trace(sigma(in, in)) - 2*xn'*xa + 2*trace(Sigma12 * H * Sigma12)) / Ni;
                    end
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