%% RTF
clc;clear;
data = 'example';
seed = 100;
reps = 1000;
maintType = 1;
dv = Inf;
cv = 1.5;
output = [];
[maintCycles, reason4stop] = SimuMaintCycles2( data, reps, maintType, dv, seed, cv);
%avgMaintCycles = mean(maintCycles);
%varMaintCycles = var(maintCycles);
%cv2MaintCycles = varMaintCycles./avgMaintCycles.^2;
%output = [avgMaintCycles(1), cv2MaintCycles(1), avgMaintCycles(2), cv2MaintCycles(2)];
output = [output,maintCycles];
clear data reps maintType dv seed cv
%cor_coef = corrcoef(maintCycles(:,1),maintCycles(:,2))
%% PM
clc;clear;
data = 'example';
seed = 100;
reps = 1000;
maintType = 1;
tm = 24;
cv = 1.5;
output = [];
for i = 1:30
    dv = tm*i;
    [maintCycles, reason4stop] = SimuMaintCycles2( data, reps, maintType, dv, seed, cv);
    avgMaintCycles = mean(maintCycles);
    varMaintCycles = var(maintCycles);
    cv2MaintCycles = varMaintCycles./avgMaintCycles.^2;
    output = [output; avgMaintCycles(1), cv2MaintCycles(1), avgMaintCycles(2), cv2MaintCycles(2)];
    %output = [output,maintCycles];
end
clear data reps maintType dv seed cv i
%% CBM
clc;clear;
data = 'example';
seed = 100;
reps = 1000;
maintType = 2;
theta = 1;
cv = 1.5;
output = [];
for i = 1:10
    dv = theta*i;
    [maintCycles, reason4stop] = SimuMaintCycles2( data, reps, maintType, dv, seed, cv);
    %avgMaintCycles = mean(maintCycles);
    %varMaintCycles = var(maintCycles);
    %cv2MaintCycles = varMaintCycles./avgMaintCycles.^2;
    %output = [output; avgMaintCycles(1), cv2MaintCycles(1), avgMaintCycles(2), cv2MaintCycles(2)];
    output = [output,maintCycles];
end
clear data reps maintType dv seed cv i
%%
order = [];
for i = 1:10
    order = [order;randperm(2)]
end
