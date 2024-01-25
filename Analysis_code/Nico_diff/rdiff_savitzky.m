function x1 = rdiff_savitzky(x, t, sigma, strategy)
% RDIFF_SAVITZKY   Regularised numerical differentiation using method based
%                  on Savitzky-Golay filter.
%    X1 = RDIFF_SAVITZKY(X, T) returns an estimate of the first derivative
%    of the N-dimensional vector X; if T is a scalar, it is used as spacing
%    between abscissa; if T is an N-dimensional vector, its elements are
%    used as the abscissa. This method is based on the assumption that the
%    abscissa are equidistant. It is required that N >= 5. If X is a
%    matrix, RDIFF_SAVITZKY differentiates each column of X independently.
%    
%    X1 = RDIFF_SAVITZKY(X, T, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF_SAVITZKY(X, T, SIGMA, STRATEGY) or
%    X1 = RDIFF_SAVITZKY(X, T, [], STRATEGY) uses the strategy for the
%    optimisation of the degrees of the algebraic polynomials, implicitly
%    approximating subsequences of data, selected among the following ones:
%       'dp'      based on discrepancy principle (requires SIGMA)
%       'ncp'     based on normalised cumulative periodogram
%       'sure'    based on Stein's unbiased risk estimator (requires SIGMA)
%    By default, the 'sure' strategy is used if SIGMA is provided, and the
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
[N, R] = size(x);
if N < 5
    error('At least 5 data points are necessary for ''savitzky'' method (given %d)', N)
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
        strategy = 'sure';
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
        if numel(sigma) ~= 1 && any(size(sigma) ~= N)
            error('Wrong size of sigma: [%d, %d] (should be scalar or %dx%d matrix)', size(sigma, 1), size(sigma, 2), N, N)
        end
    end
end

% Filter width
W = min([floor((N-1)/2), 45]);
% (number of data points approximated by a single polynomial equals 2*W+1)

% Generate matrices with coefficients of Savizky-Golay filters
persistent sgmatrices
if isempty(sgmatrices) || W > size(sgmatrices.H0, 1)
    sgmatrices = generatesgmatrices(W);
end

% Initialise vector/matrix for estimates of the derivative
x1 = zeros(N, R);

% Number of points approximated by single polynomial
Ni = 2*W+1;

% Maximum degree of approximating polynomial
NK = min([50 Ni-1]);

% Criterion values for optimisation of degree of approximating polynomial
crit = inf(NK, 1);

% Scaling factor for discrepancy principle
uDP = 2;

for n = 1:N     % for each time instant
    
    % Get indices of data points to be approximated
    if n < W+1
        in = 1:Ni;
        n1 = n;
    elseif n <= N-W
        in = n-W:n+W;
        n1 = W+1;
    elseif n >= N-W
        in = N-Ni+1:N;
        n1 = Ni-(N-n);
    end
    dt = (t(in(end)) - t(in(1)))/2/W;
    
    for r = 1:R     % for each differentiated sequence of data
        
        % Subset of data to be approximated and differentiated
        xn = x(in, r);
        
        % Estimates of the derivative for different degrees
        % of approximating polynomial
        x1n = zeros(NK, 1);
        
        for K = 1:NK
            
            % Approximate data
            H = sgmatrices.H0{W,K};
            xa = H * xn;
            
            % Estimate derivative
            x1n(K) = sgmatrices.H1{W,K}(n1, :) * xn / dt;
            
            % Compute criteria for optimisation of degree of approximating
            % polynomial K
            switch lower(strategy)
                case 'dp'
                    if numel(sigma) == 1
                        crit(K) = abs(sum((xa - xn).^2) - uDP*N*sigma) / Ni;
                    else
                        crit(K) = abs(sum((xa - xn).^2) - uDP*trace(sigma(in, in))) / Ni;
                    end
                case 'ncp'
                    L = floor(Ni/2);
                    l = (1:L)'/L;
                    rf = fft(xa - xn);
                    rp = abs(rf(1:L)).^2;
                    rncp = cumsum(rp)/sum(rp);
                    crit(K) = sum((rncp-l).^2)/L;
                case 'sure'
                    if numel(sigma) == 1
                        crit(K) = (sum(xa.^2) + sum(xn.^2) - Ni*sigma - 2*xn'*xa + 2*sigma*trace(H)) / Ni;
                    else
                        Sigma12 = sqrt(sigma(in, in));
                        crit(K) = (sum(xa.^2) + sum(xn.^2) - trace(sigma(in, in)) - 2*xn'*xa + 2*trace(Sigma12 * H * Sigma12)) / Ni;
                    end
                otherwise
                    error('Unknown strategy: %s', strategy)
            end
        end
        
        % Select estimate of derivative corresponding to minimum of criterion
        [~, KSel] = min(crit);
        x1(n, r) = x1n(KSel);
    end
end

% Return row vector if x was one
if exist('xrow', 'var') && xrow
    x1 = x1';
end

% End main function
% =========================================================================

% =========================================================================
% Generation of matrices with coefficients of Savitzky-Golay filters
function sgm = generatesgmatrices(W)
% Based on: P. A. Gorry, "General least-squares smoothing and
% differentiation by the convolution (Savitzky-Golay) method," Analytical
% Chemistry, vol. 62, no. 6, pp. 570-573, 1990.

sgm.H0 = cell(W, 2*W+1);
sgm.H1 = cell(W, 2*W+1);

for w = 1:W
    
    % Compute Gram polynomials
    P0 = zeros(w+1, 2*w+1);
    P1 = zeros(w+1, 2*w+1);
    P0(2,:) = 1;
    i = -w:w;
    K = 2*w;
    for k = 1:K
        P0(k+2,:) = (4*k-2)/k/(2*w-k+1) * i.*P0(k+1,:) - (k-1)*(2*w+k)/k/(2*w-k+1) * P0(k,:);
        P1(k+2,:) = (4*k-2)/k/(2*w-k+1) * (i.*P1(k+1,:) + P0(k+1,:)) - (k-1)*(2*w+k)/k/(2*w-k+1) * P1(k,:);
    end
    P0 = P0(2:end,:);
    P1 = P1(2:end,:);
    
    % Compute coefficients of Savitzky-Golay filters
    sgm.H0{w,1} = 1/(2*w+1) * repmat(P0(1,1:2*w+1), [2*w+1, 1]) .* repmat(P0(1,1:2*w+1)', [1, 2*w+1]);
    sgm.H1{w,1} = 1/(2*w+1) * repmat(P0(1,1:2*w+1), [2*w+1, 1]) .* repmat(P1(1,1:2*w+1)', [1, 2*w+1]);
    for k = 1:K
        sgm.H0{w,k+1} = sgm.H0{w,k} + (2*k+1)*gf(2*w,k)/gf(2*w+k+1,k+1) * repmat(P0(k+1,1:2*w+1), [2*w+1, 1]) .* repmat(P0(k+1,1:2*w+1)', [1, 2*w+1]);
        sgm.H1{w,k+1} = sgm.H1{w,k} + (2*k+1)*gf(2*w,k)/gf(2*w+k+1,k+1) * repmat(P0(k+1,1:2*w+1), [2*w+1, 1]) .* repmat(P1(k+1,1:2*w+1)', [1, 2*w+1]);
    end
end

sgm.H0 = sgm.H0(:,2:end);
sgm.H1 = sgm.H1(:,2:end);

% =========================================================================
% Generalised factorial (used for computing coefficients of S.-G. filters)
function r = gf(a,b)

r = 1;
for f = (a-b+1):a
    r = r*f;
end