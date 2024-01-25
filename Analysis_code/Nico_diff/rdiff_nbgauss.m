function x1 = rdiff_nbgauss(x, t, sigma, strategy, t1)
% RDIFF_NBGAUSS   Regularised numerical differentiation using method based
%                 on approximation of data using linear combinations
%                 of Gaussian basis functions.
%    X1 = RDIFF_NBGAUSS(X, T) returns an estimate of the first derivative
%    of X; if T is a scalar, it is used as spacing between abscissa; if T
%    is a vector, its elements are used as the abscissa. If X is a matrix,
%    RDIFF_NBGAUSS differentiates each of its columns independently.
%    
%    X1 = RDIFF_NBGAUSS(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF_NBGAUSS(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_NBGAUSS(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the number of basis functions, and the parameter
%    controlling their scale, selected among the following ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'gcv'     based on generalised cross validation
%       'ncp'     based on normalised cumulative periodogram
%       'sure'    based on Stein's unbiased risk estimator (requires SIGMA)
%    By default, the 'sure' strategy is used if SIGMA is provided, and the
%    'gcv' strategy otherwise.
% 
%    X1 = RDIFF_NBGAUSS(X, T, SIGMA, STRATEGY, T1), estimates
%    the derivative at abscissa specified by T1.
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
        strategy = 'sure';
    else
        strategy = 'gcv';
    end
end

% Check availability and size of sigma
if any(strcmpi(strategy, {'dp', 'sure'}))
    if ~exist('sigma', 'var') || isempty(sigma)
        warning('Cannot use ''%s'' strategy when sigma is not specified; using ''gcv'' strategy instead.', strategy)
        strategy = 'gcv';
    else
        if numel(sigma) ~= 1 && any(size(sigma) ~= N)
            error('Wrong size of sigma: [%d, %d] (should be scalar or %dx%d matrix)', size(sigma, 1), size(sigma, 2), N, N)
        end
    end
end

% Check availability and size of t1
if ~exist('t1', 'var') || isempty(t1)
    t1 = t;
end
if ~isvector(t)
    error('t1 should be a vector; size(t1) = [%d, %d]', size(t1,1), size(t1,2))
end
if isrow(t1); t1 = t1'; end
M = length(t1);

% Basis function
phi = @(t,s) 1 / sqrt(2*pi) / s * exp(- (t.^2 ./ (2*s.^2)));
phi1 = @(t,s) -t / sqrt(2*pi) / s^3 .* exp(- (t.^2 ./ (2*s.^2)));

% Considered values of regularisation parameters
cs = linspace(0.5, 2.5, 5)';
Ks = unique(round(linspace(3,round(N/2),N)))';
NK = numel(Ks);
Nc = numel(cs);

% Constants used during optimisation of regularisation parameters
switch lower(strategy)
    case 'dp'
        uDP = 2;
        if isscalar(sigma), trS = N*sigma;
        else, trS = trace(sigma); end
    case 'gcv'
        
    case 'ncp'
        L = floor(N/2);
        l = repmat((1:L)' / L, [1 R]);
    case 'sure'
        sigma12 = sqrt(sigma);
    otherwise
        error('Unknown strategy: %s', strategy)
end

% Criterion values for optimisation of regularisation parameters
crit = inf(NK, R, Nc);

% Estimates of the derivative
X1 = zeros(M, R, NK, Nc);

% Optimise regularisation parameters
for nK = 1:NK   % for each considered number of basis functions
    K = Ks(nK);
    
    % Centres of basis functions
    tp = linspace(t(1), t(N), K);
    
    for nc = 1:Nc   % for each considered value of scaling parameter
        c = cs(nc);
        s = c * (t(end) - t(1)) / (K-1);
        
        % Generate matrix of basis functions
        Phi = zeros(N, K);
        for k = 1:K
            Phi(:, k) = phi(t - tp(k), s);
        end
        
        if cond(Phi) < 1e10
            
            % Approximate data
            p = Phi \ x;
            xa = Phi * p;
            
            % Generate matrix of derivatives of basis functions
            Phi1 = zeros(M, K);
            for k = 1:K
                Phi1(:, k) = phi1(t1 - tp(k), s);
            end
            
            % Differentiate
            X1(:,:,nK,nc) = Phi1 * p;
            
            % Compute criteria for optimisation of regularisation parameters
            switch lower(strategy)
                case 'dp'
                    crit(nK,:,nc) = abs(sum((xa - x).^2, 1) - uDP*trS);
                case 'gcv'
                    crit(nK,:,nc) = sum((xa - x).^2, 1) / trace(eye(N) - Phi*(Phi\eye(N)))^2;
                case 'ncp'
                    residual = xa - x;
                    rf = fft(residual,[],1);
                    rp = abs(rf(1:L,:)).^2;
                    rncp = cumsum(rp,1) ./ repmat(sum(rp,1), [L 1]);
                    crit(nK,:,nc) = sqrt(sum((rncp - l).^2, 1));
                case 'sure'
                    if R > 1
                        if numel(sigma) == 1
                            crit(nK,:,nc) = sum(xa.^2, 1) - 2*diag(x'*xa)' + 2*sigma*trace(Phi * (Phi \ eye(N)));
                        else
                            crit(nK,:,nc) = sum(xa.^2, 1) - 2*diag(x'*xa)' + 2*trace(sigma12 * (Phi * (Phi \ eye(N))) * sigma12);
                        end
                    else
                        if numel(sigma) == 1
                            crit(nK,:,nc) = sum(xa.^2, 1) - 2*x'*xa + 2*sigma*trace(Phi * (Phi \ eye(N)));
                        else
                            crit(nK,:,nc) = sum(xa.^2, 1) - 2*x'*xa + 2*trace(sigma12 * (Phi * (Phi \ eye(N))) * sigma12);
                        end
                    end
            end
        end
    end
end

% Select optimum values of regularisation parameters
x1 = zeros(M,R);
for r = 1:R
    critr = squeeze(crit(:,r,:));
    
    KsCrit = repmat(Ks, [1 Nc]);
    KsCrit(critr > min(critr(:))) = Inf;
    
    csCrit = repmat(cs', [NK 1]);
    csCrit(KsCrit > min(KsCrit(:))) = -Inf;
    
    [~, pSel] = max(csCrit(:));
    [nKSel, ncSel] = ind2sub([NK, Nc], pSel);
    x1(:,r) = X1(:,r,nKSel, ncSel);
end

% Return row vector if x was one
if exist('xrow', 'var') && xrow
    x1 = x1';
end