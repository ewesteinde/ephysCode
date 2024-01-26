function Q = qmatrix(t)
% QMATRIX   Matrix of coefficients of 2nd-order Newton-Cotes quadrature.
%    Q = QMATRIX(N), where N is a scalar, returns an NxN matrix whose
%    elements are the coefficients of the second-order Newton-Cotes
%    quadrature corresponding to equidistant abscissa.
% 
%    Q = QMATRIX(T), where T is an N-dimensional vector, returns an NxN
%    matrix whose elements are the coefficients of the second-order
%    Newton-Cotes quadrature corresponding to the abscissa specified by T.
%    The elements of T do not need to be equidistant.
%    
%    The integral of an N-dimensional vector X can be estimated as:
%    I = DT * QMATRIX(N) * X
%    where DT is the spacing between the abscissa or:
%    I = QMATRIX(T) * X
%    where it is assumed that X represents the values of the integrated
%    function evaluated at the abscissa T.
%    
%    Jakub Wagner, February 2, 2020
%    Institute of Radioelectronics and Multimedia Technology
%    Warsaw University of Technology

if isscalar(t), t = (1:t)'; end     % number of points specified instead
                                    % of vector of abscissa
N = length(t);
Q = zeros(N,N);
for n = 2:N
    Q(n:N, [n-1, n]) = Q(n:N, [n-1, n]) + (t(n)-t(n-1))/2;
end