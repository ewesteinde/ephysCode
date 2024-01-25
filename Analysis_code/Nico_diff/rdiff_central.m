function x1 = rdiff_central(x, t, sigma, strategy)
% RDIFF_CENTRAL   Regularised numerical differentiation using method
%                 based on central-difference formula with extended
%                 differentiation step.
%    X1 = RDIFF_CENTRAL(X, T) returns an estimate of the first derivative
%    of the N-dimensional vector X; if T is a scalar, it is used as spacing
%    between abscissa; if T is an N-dimensional vector, its elements are
%    used as the abscissa. This method is based on the assumption that the
%    abscissa are equidistant. It is required that N >= 5. If X is a
%    matrix, RDIFF_CENTRAL differentiates each column of X independently.
%    
%    X1 = RDIFF_CENTRAL(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF_CENTRAL(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_CENTRAL(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the differentiation step selected among the following
%    ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'gcv'     based on generalised cross validation
%       'ncp'     based on normalised cumulative periodogram
%       'sure'    based on Stein's unbiased risk estimator (requires SIGMA)
%       'lp'      based on: S. Lu and S. Pereverzev, "Numerical
%                 differentiation from a viewpoint of regularization
%                 theory", Mathematics of Computation, vol. 75, no. 256,
%                 pp. 1853-1870, 2006 (requires SIGMA)
%    By default, the 'lp' strategy is used if SIGMA is provided, and the
%    'gcv' strategy otherwise.
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
    error('At least 5 data points are necessary for ''central'' method (given %d)', N)
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
        strategy = 'lp';
    else
        strategy = 'gcv';
    end
end

% Check availability and size of sigma
if any(strcmpi(strategy, {'dp', 'sure', 'lp'}))
    if ~exist('sigma', 'var') || isempty(sigma)
        warning('Cannot use ''%s'' strategy when sigma is not specified; using ''gcv'' strategy instead.', strategy)
        strategy = 'gcv';
    else
        if numel(sigma) ~= 1 && any(size(sigma) ~= N)
            error('Wrong size of sigma: [%d, %d] (should be scalar or %dx%d matrix)', size(sigma, 1), size(sigma, 2), N, N)
        end
    end
end

% Initialise vector/matrix for estimates of the derivative
x1 = zeros(N, R);

if strcmpi(strategy, 'lp')
    % Strategy for optimisation of differentiation step based on the 2006
    % paper by Lu & Pereverzev (see above for complete reference).
    
    % Maximum differentiation step considered
    K = max([2, min([120, floor((N-1)/3)])]);
    
    % Considered values of differentiation step
    kBase = 1.5;
    ikMax = ceil(log(K)/log(kBase)-1);
    ks = ceil(kBase.^(0:ikMax));
    Nk = numel(ks);
    
    if numel(sigma) == 1
        delta = 3*sqrt(sigma);
    end
    
    for r = 1:R     % for each differentiated sequence of data
        for n = 1:N     % for each time instant
            
            x1n = zeros(Nk, 1); % Estimates of derivative at n obtained using different values of k
            dts = zeros(Nk, 1); % Differentiation steps corresponding to different values of k
            
            for nk = 1:Nk   % for each considered value of differentiation step
                
                k = ks(nk); % differentiation step
                
                % Indices of data points analysed for current value of k
                if n < k+1
                    in = [n, n+k];
                elseif n <= N-k
                    in = [n-k, n+k];
                elseif n > N-k
                    in = [n-k, n];
                end
                
                % Estimate derivative
                x1n(nk) = (x(in(2), r) - x(in(1), r)) / (t(in(2)) - t(in(1)));
                
                % Check condition for further extending diff. step
                if numel(sigma) > 1
                    delta = 3*sqrt(max(diag(sigma(in(1):in(2), in(1):in(2)))));
                end
                dts(nk) = (t(in(2)) - t(in(1)));
                crit = abs(x1n(nk) - x1n(1:nk-1));
                thr = 2*delta ./ (dts(1:nk-1));
                
                % If threshold exceeded, keep previous value of dt and stop
                if any(crit > thr)
                    nkOpt = nk-1;
                    break
                else
                    % If not, consider larger step
                    nkOpt = nk;
                end
            end
            
            % Select estimate of derivative corresponding to optimum step
            x1(n,r) = x1n(nkOpt);
        end
    end
else
% Other strategies for optimisation of differentiation step
    
    % Maximum differentiation step considered
    K = min([120, floor((N-1)/2)]);
    % (i.e. at most 2*K+1 data points approximated by a single polynomial)
    
    % Generate "hat" matrices
    cdmatrices = generatecdmatrices(K);
    
    % Criterion values for optimisation of differentiation step
    crit = inf(K, 1);
    
    % Scaling factor for discrepancy principle
    uDP = 2;
    
    for r = 1:R     % for each differentiated sequence of data
        
        % Matrix with estimates of the derivative corresponding to different
        % time instants and different values of the differentiation step
        x1r = zeros(N,K);
        
        for n = 1:N     % for each time instant
            for k = 1:K     % for each considered value of diff. step
                
                % Approximate the data and estimate the derivative
                if n < k+1
                                            % estimate of derivative
                    x1r(n,k) = (x(n+k, r) - x(n, r)) / (t(n+k) - t(n));
                    in = 1:n+k;             % indices of approximated data
                    xn = x(in, r);          % approximated subsequence
                    H = zeros(n+k, n+k);    % "hat" matrix
                    H(:,n+k) = (-n+1:k)'/k;
                    H(:,n) = 1 - H(:,n+k);
                    xa = H * xn;            % result of approximation
                    
                elseif n <= N-k
                                            % estimate of derivative
                    x1r(n,k) = (x(n+k, r) - x(n-k, r)) / (t(n+k) - t(n-k));
                    in = n-k:n+k;           % indices of approximated data
                    xn = x(in, r);          % approximated subsequence
                    H = cdmatrices{k};      % "hat" matrix
                    xa = H * xn;            % result of approximation
                elseif n > N-k
                                            % estimate of derivative
                    x1r(n,k) = (x(n, r) - x(n-k, r)) / (t(n) - t(n-k));
                    in = n-k:N;             % indices of approximated data
                    xn = x(in, r);          % approximated subsequence
                    H = zeros(N-n+k+1, N-n+k+1); % "hat" matrix
                    H(:,1) = (k:-1:n-N)'/k;
                    H(:,k+1) = 1 - H(:,1);
                    xa = H * xn;            % result of approximation
                end
                
                % Length of analysed subsequence for current value of k
                Ni = length(in);
                
                % Compute criteria for optimisation of k
                switch lower(strategy)
                    case 'dp'
                        if numel(sigma) == 1
                            crit(k) = abs(sum((xa - xn).^2) - uDP*Ni*sigma) / Ni;
                        else
                            crit(k) = abs(sum((xa - xn).^2) - uDP*trace(sigma(in, in))) / Ni;
                        end
                    case 'gcv'
                        crit(k) = sum((xa - xn).^2) / trace(eye(Ni) - H)^2;
                    case 'ncp'
                        L = floor(Ni/2);
                        l = (1:L)'/L;
                        rf = fft(xa - xn);
                        rp = abs(rf(1:L)).^2;
                        rncp = cumsum(rp)/sum(rp);
                        crit(k) = sum((rncp-l).^2)/L;
                    case 'sure'
                        if numel(sigma) == 1
                            crit(k) = (sum(xa.^2) + sum(xn.^2) - Ni*sigma - 2*xn'*xa + 2*sigma*trace(H)) / Ni;
                        else
                            Sigma12 = sqrt(sigma(in, in));
                            crit(k) = (sum(xa.^2) + sum(xn.^2) - trace(sigma(in, in)) - 2*xn'*xa + 2*trace(Sigma12 * H * Sigma12)) / Ni;
                        end
                    otherwise
                        error('Unknown strategy: %s', strategy)
                end
            end
            
            % Select estimate of derivative corresponding to minimum of criterion
            nkSel = find(crit == min(crit), 1, 'last');
            x1(n, r) = x1r(n,nkSel);
        end
    end
end

% Return row vector if x was one
if exist('xrow', 'var') && xrow
    x1 = x1';
end

% End main function
% =========================================================================

% =========================================================================
% Generation of matrices corresponding to 2nd-order-polynomial
% approximation based on 3 points
function cdm = generatecdmatrices(KMax)

cdm = cell(KMax, 1);

for K = 1:KMax
    
    N = 2*K+1;
    cdm{K} = zeros(N,N);
    
    k = (-K:K)';
    ksq = k.^2;
    
    cdm{K}(:,1) = -k / (2*K) + ksq / (2*K^2);
    cdm{K}(:,K+1) = 1 - ksq / K^2;
    cdm{K}(:,N) = k / (2*K) + ksq / (2*K^2);
    
end