%%
clc;clear;
load('example.mat');
rng(81);
tm = Inf;
theta = 6;
reps = 10000;
MT = gamrnd(m_M^2/v_M, v_M/m_M,reps,1);
RT = gamrnd(m_R^2/v_R, v_R/m_R,reps,1);
time_in_states = zeros(reps, n);
ttf = zeros(reps, 1);
ttr = zeros(reps, 1);
reason4stop = zeros(reps,1);
start_state = 1;
maint_state = 0;
start_state_count = zeros(1,n);
ipf = 0;
%% sampleing exp cv =1
DT = exprnd(repmat(1./rate(1,:),reps,1),[reps,size(rate,2)]);
FT = exprnd(repmat(1./rate(2,:),reps,1),[reps,size(rate,2)]);
%% sampling gamma cv = 1.25
cv = 1.25;
alpha_1 = repmat(1/cv^2, reps, size(rate,2));
beta_1 = repmat(cv^2*1./rate(1,:), reps, 1);
alpha_2 = repmat(1/cv^2, reps, size(rate,2));
beta_2 = repmat(cv^2*1./rate(2,:), reps, 1);
DT = gamrnd(alpha_1,beta_1,[reps,size(rate,2)]);
FT = gamrnd(alpha_2,beta_2,[reps,size(rate,2)]);
clear alpha_1 alpha_2 beta_1 beta_2
%% TBM
for i = 1:reps
    is_failed = false;
    is_maint = false;
    for j = start_state:n
        if ~is_failed && ~ is_maint
            if DT(i,j) > FT(i,j)
                time_in_states(i,j) = FT(i,j);
                ttf(i,1) = ttf(i,1) + FT(i,j);
                if ttf(i,1) >= tm
                    is_maint = true;
                    maint_state = j;
                    reason4stop(i,1) = 2;
                    time_in_states(i,j) = tm - (ttf(i,1) - FT(i,j));
                    ttf(i,1) = tm;
                else
                    is_failed = true;
                    maint_state = n;
                    reason4stop(i,1) = 1;
                end
            else
                time_in_states(i,j) = DT(i,j);
                ttf(i,1) = ttf(i,1) + DT(i,j);
                if ttf(i,1) >= tm             
                    is_maint = true;
                    maint_state = j;
                    reason4stop(i,1) = 2;
                    time_in_states(i,j) = tm - (ttf(i,1) - DT(i,j));
                    ttf(i,1) = tm;
                else
                end
            end
        else
        end
        ttr(i,1) = (reason4stop(i,1)==1)*RT(i,1)+(reason4stop(i,1)==2)*0;
%         ttr(i,1) = (reason4stop(i,1)==1)*RT(i,1)+(reason4stop(i,1)==2)*MT(i,1); 
    end
    start_state = binornd(maint_state-1,ipf)+1;
    start_state_count(start_state) = start_state_count(start_state) + 1;
end
clear i j is_failed is_maint ans start_state maint_state
%% CBM
for i = 1:reps
    is_failed = false;
    is_maint = false;
    if start_state<=theta
        for j = start_state:n
            if ~is_failed && ~ is_maint
                if DT(i,j) > FT(i,j)
                    time_in_states(i,j) = FT(i,j);
                    ttf(i,1) = ttf(i,1) + FT(i,j);
                    is_failed = true;
                    maint_state = n;
                    reason4stop(i,1) = 1;
                else
                    time_in_states(i,j) = DT(i,j);
                    ttf(i,1) = ttf(i,1) + DT(i,j);
                    if j>=theta            
                        is_maint = true;
                        maint_state = j;
                        reason4stop(i,1) = 2;
                    else
                    end
                end
            else
            end
            ttr(i,1) = is_failed*RT(i,1)+is_maint*MT(i,1);    
        end
    else
        ttf(i,1)=0;
        ttr(i,1)=MT(i,1);
        maint_state = start_state;
    end
    start_state = binornd(maint_state-1,ipf)+1;
    start_state_count(start_state) = start_state_count(start_state) + 1;
end
clear i j is_failed is_maint ans start_state maint_state
%%
hf = histogram(ttf);
cf = [cumsum([0,hf.Values]/reps)',hf.BinEdges'];
str_cf = 'random.continuous(';
for i = 1:size(cf,1)
    str_cf = strcat(str_cf,string(cf(i,2)), ',', string(cf(i,1)),',');
end
str_cf = strcat(str_cf,')');
%%
hr = histogram(ttr);
cr = [cumsum([0,hr.Values]/reps)',hr.BinEdges'];
str_cr = 'random.continuous(';
for i = 1:size(cr,1)
    str_cr = strcat(str_cr,string(cr(i,2)), ',', string(cr(i,1)),',');
end
str_cr = strcat(str_cr,')');
