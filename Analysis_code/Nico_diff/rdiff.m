function x1 = rdiff(x, t, method, sigma, strategy, t1, x1max)
% RDIFF   Regularised numerical differentiation.
%    X1 = RDIFF(X, T) returns an estimate of the first derivative of X;
%    if T is a scalar, it is used as spacing between abscissa; if T is
%    a vector, its elements are used as the abscissa. If X is a matrix,
%    RDIFF differentiates each of its columns independently.
%    
%    X1 = RDIFF(X, T, METHOD) uses method selected among following ones:
%       'totalvar'  based on total-variation regularisation
%       'tikhonov'  based on Tikhonov regularisation
%       'savitzky'  based on Savitzky-Golay filter
%       'kalman'    based on Kalman filter
%       'central'   based on central-difference formula with extended
%                   differentiation step
%       'sinstep'   based on modified central-difference formula
%       'expstep'   based on modified central-difference formula
%       'nbgauss'   based on smoothing approximation using Gaussian basis
%                   functions
%       'nbcos'     based on smoothing approximation using even powers
%                   of cosine functions
%       'tsvd'      based on truncated singular value decomposition
%       'landweber' based on Landweber iteration
%    By default, the 'totalvar' method is used.
%    
%    X1 = RDIFF(X, T, METHOD, SIGMA), where SIGMA is a scalar, uses SIGMA
%    as an estimate of the variance of the errors corrupting the elements
%    of X under the assumption that those errors are independent and
%    identically distributed. If SIGMA is an NxN matrix, it is used as the
%    covariance matrix of those errors.
%    
%    X1 = RDIFF(X, T, METHOD, SIGMA, STRATEGY) uses a selected strategy for
%    the optimisation of the regularisation parameter(s); see documentation
%    of M-files corresponding to particular methods ('rdiff_[METHOD].m')
%    for lists of available strategies.
%    
%    X1 = RDIFF(X, T, METHOD, SIGMA, STRATEGY, T1), if METHOD is 'nbgauss'
%    or 'nbcos', estimates the derivative at abscissa specified by T1.
%    For methods other than 'nbgauss' and 'nbcos', T1 is ignored.
%    
%    X1 = RDIFF(X, T, METHOD, SIGMA, STRATEGY, [], X1MAX), if METHOD is
%    'landweber', uses X1MAX as an a priori upper bound for absolute value
%    of the first derivative. For methods other than 'landweber', X1MAX is
%    ignored.
% 
%    Jakub Wagner, February 3, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

% Check input arguments
if ~exist('method', 'var') || isempty(method)
    if exist('t1', 'var') && ~isempty(t1)
        method = 'nbgauss';
    elseif exist('x1max', 'var') && ~isempty(x1max)
        method = 'landweber';
    else
        method = 'totalvar';
    end
end
if ~exist('sigma', 'var')
    sigma = [];
end
if ~exist('strategy', 'var')
    strategy = '';
end
if ~exist('t1', 'var')
    t1 = [];
end
if ~any(strcmpi(method, {'nbgauss', 'nbcos'})) && ~isempty(t1)
    warning('Cannot use parameter ''t1'' with method ''%s''.', method)
end
if ~exist('x1max', 'var')
    x1max = [];
elseif ~strcmpi(method, 'landweber') && ~isempty(x1max)
    warning('Parameter ''x1max'' can only be used with method ''landweber''.')
end

% Differentiate
switch lower(method)
        
    case 'savitzky'
        x1 = rdiff_savitzky(x, t, sigma, strategy);
        
    case 'central'
        x1 = rdiff_central(x, t, sigma, strategy);
        
    case 'sinstep'
        x1 = rdiff_sinstep(x, t, sigma, strategy);
        
    case 'expstep'
        x1 = rdiff_expstep(x, t, sigma, strategy);
        
    case 'nbgauss'
        x1 = rdiff_nbgauss(x, t, sigma, strategy, t1);
        
    case 'nbcos'
        x1 = rdiff_nbcos(x, t, sigma, strategy, t1);
        
    case 'tikhonov'
        x1 = rdiff_tikhonov(x, t, sigma, strategy);
        
    case 'totalvar'
        x1 = rdiff_totalvar(x, t, sigma, strategy);
        
    case 'tsvd'
        x1 = rdiff_tsvd(x, t, sigma, strategy);
        
    case 'landweber'
        x1 = rdiff_landweber(x, t, sigma, strategy, x1max);
        
    case 'kalman'
        x1 = rdiff_kalman(x, t, sigma, strategy);
        
    otherwise
        error('Unknown method: %s', method)
end