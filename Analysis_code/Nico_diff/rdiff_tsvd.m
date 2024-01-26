function x1 = rdiff_tsvd(x, t, sigma, strategy)
% RDIFF_TSVD   Regularised numerical differentiation method based on
%              truncated singular value decomposition.
%    X1 = RDIFF_TSVD(X, T) returns an estimate of the first derivative
%    of X; if T is a scalar, it is used as spacing between abscissa; if T
%    is a vector, its elements are used as the abscissa. If X is a matrix,
%    RDIFF_TSVD differentiates each of its columns independently.
%    
%    X1 = RDIFF_TSVD(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA as
%    an estimate of the variance of the errors corrupting the elements of X
%    under the assumption that those errors are independent and identically
%    distributed. If SIGMA is an NxN matrix, it is used as the covariance
%    matrix of those errors.
%    
%    X1 = RDIFF_TSVD(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_TSVD(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the regularisation parameter selected among the
%    following ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'gcv'     based on generalised cross validation
%       'ncp'     based on normalised cumulative periodogram
%       'sure'    based on Stein's unbiased risk estimator (requires SIGMA)
%    By default, the 'dp' strategy is used if SIGMA is provided, and the
%    'ncp' strategy otherwise.
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
    if exist('sigma', 'var') && ~isempty(sigma)
        strategy = 'dp';
    else
        strategy = 'ncp';
    end
end

% Check availability and size of sigma
if any(strcmpi(strategy, {'dp', 'sure'}))
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

% Considered values of regularisation parameter
nThs = (1:N-1)';
NnTh = numel(nThs);

% Constants used during optimisation of regularisation parameters
switch lower(strategy)
    case 'dp'
        uDP = 2;
        trS = trace(sigma);
    case 'gcv'
        
    case 'ncp'
        L = floor(N/2);
        l = repmat((1:L)' / L, [1 R]);
    case 'sure'
        Sigma12 = sqrt(sigma);
    otherwise
        error('Unknown strategy: %s', strategy)
end

% Criterion values for optimisation of regularisation parameter
crit = inf(NnTh, R);

% Estimates of the derivative
X1 = zeros(N, R, NnTh);

% Make sure first data point equals 0
xp = x - repmat(x(1,:), [N 1]);

% Singular value decomposition
[U, S, V] = svd(Q);
sv = diag(S);

% Optimise regularisation parameter
for nnTh = 1:NnTh   % for each considered value of regularisation parameter
    nTh = nThs(nnTh);
    
    % Differentiate and approximate data
    S1i = diag([1 ./ sv(1:nTh); zeros(N - nTh, 1)]);
    X1(:,:,nnTh) = V*S1i*U'*xp;
    xa = Q * X1(:,:,nnTh);
    
    % Compute criteria for optimisation of regularisation parameter
    switch lower(strategy)
        case 'dp'
            crit(nnTh,:) = abs(sum((xa - xp).^2, 1) - uDP*trS);
        case 'gcv'
            crit(nnTh,:) = sum((xa - xp).^2, 1) / trace(eye(N) - Q*V*S1i*U')^2;
        case 'ncp'
            residual = xa - xp;
            rf = fft(residual,[],1);
            rp = abs(rf(1:L,:)).^2;
            rncp = cumsum(rp,1) ./ repmat(sum(rp,1), [L 1]);
            crit(nnTh,:) = sqrt(sum((rncp - l).^2, 1));
        case 'sure'
            if R > 1
                crit(nnTh,:) = sum(xa.^2, 1) - 2*diag(xp'*xa)' + 2*trace(Sigma12 * (Q*V*S1i*U') * Sigma12);
            else
                crit(nnTh,:) = sum(xa.^2, 1) - 2*xp'*xa + 2*trace(Sigma12 * (Q*V*S1i*U') * Sigma12);
            end
    end
end

% Select optimum value of regularisation parameter
x1 = zeros(N,R);
for r = 1:R
    nnThSel = find(crit(:,r) == min(crit(:,r)), 1, 'first');
    x1(:,r) = X1(:, r, nnThSel);
end

% Return row vector if x was one
if exist('xrow', 'var') && xrow
    x1 = x1';
end