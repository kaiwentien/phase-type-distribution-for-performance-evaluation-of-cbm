function [lam, alpha, rates] = CoxianDist(mean, b)
a = 1 - b;
A = [1, cumprod(a(1:end-1))];
coef_1 = sum((1:length(b)).*A.*b);
lam = coef_1/mean;
alpha = [1, zeros(1, length(b)-1)];
rates = [a*lam; b*lam];
end