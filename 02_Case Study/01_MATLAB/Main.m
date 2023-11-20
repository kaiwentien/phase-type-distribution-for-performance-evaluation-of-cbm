clear;clc;
% Heller machine with three different tempers
b1 = [0.1611 0.1289 0.0967 0.0644 0.0322 0.1933 0.2256 0.2578 0.2900 1.0000];            % increase failure rate
b2 = [0.1240 0.1240 0.1240 0.1240 0.1240 0.1240 0.1240 0.1240 0.1240 1.0000];            % tube-shaped 1
b3 = [0.1616 0.1077 0.0539 0.2154 0.2693 0.3231 0.3770 0.4308 0.4847 1.0000];            % constant failure rate
%%
mttf = 385; % average time to failure of three machines
%b = [];
% source of variation
m_T0 = 0.2;
v_T0 = 0;
m_R = 9*1;
v_R = 0.25*9;
m_M = 3;
v_M = 9;
%%
[~, alpha, rate] = CoxianDist(mttf, b);
rate = [rate;ones(1,10)];

%%
% get the rates
[lam1, alpha1, rate1] = CoxianDist(mttf, b1);
[lam2, alpha2, rate2] = CoxianDist(mttf, b2);
[lam3, alpha3, rate3] = CoxianDist(mttf, b3);
rate1 = [rate1;ones(1,10)];
rate2 = [rate2;ones(1,10)];
rate3 = [rate3;ones(1,10)];

%% CBM data
iter = 10;
output=zeros(iter, 9);

for i = 1:iter
    output(i,1) = i;
    [output(i,1), output(i,2), output(i,3), output(i,4), output(i,5), output(i,6), output(i,7),output(i,8), output(i,9)] ...
    = CBMPolicy(rate3, alpha3, m_T0, v_T0, m_R, v_R, m_M, v_M, i);
end

%% PM data -offline

iter = 15;
output = zeros(iter, 9);
for i = 1:iter
    output(i,1) = i;
    [output(i,1), output(i,2), output(i,3), output(i,4), output(i,5), output(i,6), output(i,7),output(i,8), output(i,9)] ...
    = TBMPolicy(rate, alpha, m_T0, v_T0, m_R, v_R, 0, 0, i*24);
end
%% PM data -online

iter = 30;
output = zeros(iter, 9);
for i = 1:iter
    output(i,1) = i;
    [output(i,1), output(i,2), output(i,3), output(i,4), output(i,5), output(i,6), output(i,7),output(i,8), output(i,9)] ...
    = TBMPolicy(rate, alpha, m_T0, v_T0, m_R, v_R, m_M, v_M, i*24);
end
%% probability plot






