function [ m ] = CoxianDistMoment( rates, k, alpha )
%PHMoment Summary of this function goes here
%   rates: the 2 x n array: first row is deterioration rate; 
%          second row is failure rate.
%   k: the 1 x p array indicates the k_th moment to generate
%   m: the 1 x p array generate moments
    n = size(rates, 2);
    p = size(k, 2);
    m = zeros(1,p);
    S = diag(-rates(1,:)-rates(2,:))+diag(rates(1,1:end-1),1);  % sub-ordinate matrix
    if ~exist('alpha', 'var')
        alpha = [1, zeros(1,n-1)];
    end

    for i = 1:p
       m(i) = (-1)^k(i)*factorial(k(i))*alpha*S^(-k(i))*ones(n,1);
    end


end

