function D = pdmatrix(t, k)
% PDMATRIX   Matrix of coefficients of progressive-difference formula.
%    D = PDMATRIX(N) generates an NxN matrix whose elements are the
%    coefficients of the progressive-difference formula.
%    The first derivative of an N-dimensional vector X can be estimated as:
%    X1 = PDMATRIX(N) * X / DT
%    where DT is the spacing of the data contained in X.
%    
%    D = PDMATRIX(T), where T is an N-dimensional vector, generates such a
%    matrix for differentiation of vectors which correspond to abscissa
%    specified by T, i.e.:
%    X1 = PDMATRIX(T) * X
%    The elements of T do not need to be equidistant.
%    
%    D = PDMATRIX(N, K) or D = PDMATRIX(T, K) generates a matrix
%    representative of the progressive-difference formula applied K times.
%    The K-th derivative of an N-dimensional vector X can be estimated as:
%    XK = PDMATRIX(N) * X / (DT^K)
%    or:
%    XK = PDMATRIX(T) * X
%    
%    Based on: J. J. Stickel, "Data smoothing and numerical differentiation
%    by a regularization method", Computers & Chemical Engineering,
%    vol. 34, no. 4, pp. 467-475, 2010.
%    
%    Jakub Wagner, February 2, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

if ~exist('k', 'var'), k = 1; end   % default order of differentiation = 1
if isscalar(t), t = (1:t)'; end     % number of points specified instead
                                    % of vector of abscissa
N = length(t);

if k == 0, D = eye(N); return; end % "order = 0" -> identity matrix

% Matrix of inverses of spacings between abscissa
Vd = zeros(N-k,N-k);
for n = (ceil(k/2)+1) : (N-floor(k/2))
    i = n - ceil(k/2);
    Vd(i,i) = 1/(t(n+floor(k/2)) - t(n-ceil(k/2)));
end

% Matrix of "unscaled" progressive-difference formula coefficients
Dh1 = zeros(N-k,N-k+1);
for n = 1:N-k
    Dh1(n,n) = -1;
    Dh1(n,n+1) = 1;
end

% Recursive implementation for order higher than 1
if k == 1
    D = Vd*Dh1;
else
    D = k*Vd*(Dh1*pdmatrix(t,k-1));
end