function [maintCycles,reason4Stop] = SimuMaintCycles( data, reps, maintType, dv, seed, cv)
    %SimuMaintCycles :Simulate maintenance cycles using Coxian distribution
    % using [sojourn times, failure probabiliy]
    %   data.mat: the data set should have rate, m_M, v_M, m_R, v_R
    %   reps: number of replication running in this simulation
    %   maintType: 1-TBM, 2-CBM
    %   dv: the maint decisions: TBM: maint interval, CBM: maint threshold
    %   cv (optional): the distribution of time in each stage is set to be a
    %   gamma distribution determined by rate and cv. if the input of cv is
    %   empty, then cv = 1 (exponential)
    %   maintCycles: col1: uptimes, col2:downtimes
    %   reason4Stop: 1: corrective maint, 2: planned maint
    if ~exist('cv','var')
        cv = 1;
    end
    load(data);                                 % load data
    n = size(rate,2);                           % number of stages of HI
    rng(seed);                                  % random seed
    tm = Inf;                                   % maintenance interval
    theta = 6;                                  % maintenance threshold
    uptimes = zeros(reps, 1);                   % uptimes
    downtimes = zeros(reps, 1);                 % downtimes
    reason4Stop = zeros(reps,1);                % 1: corrective, 2: preventive
    MT = gamrnd(m_M^2/v_M, v_M/m_M,reps,1);     % sampling of PT time
    RT = gamrnd(m_R^2/v_R, v_R/m_R,reps,1);     % sampling of CT time
    % create deterioation process
    m_S = 1./sum(rate(1:2,:));                  % sojourn times in states
    fail_p = rate(2,:)./sum(rate(1:2,:));       % failure prob in each states
    alpha_ST = repmat(1/cv^2, reps, size(m_S,2));
    beta_ST = repmat(cv^2*m_S, reps, 1);
    ST = gamrnd(alpha_ST,beta_ST,[reps,size(m_S,2)]);
    % find the closest one in each row (the first fail state)
    [~, fail_state] = max(rand(reps,size(m_S,2)) < repmat(fail_p,reps,1),[],2);
    % calculate failure times
    ttf = zeros(reps, 1);
    for i = 1:reps
        for j = 1: fail_state(i,1)
            ttf(i,1) = ttf(i,1) + ST(i,j);
        end
    end
    for i = 1: reps
        if maintType == 1
            tm = dv;
            if ttf(i,1) <= tm
                % do corrective maintenance
                reason4Stop(i,1) = 1;
                uptimes(i,1) = ttf(i,1);
                downtimes(i,1) = RT(i,1);
            else
                % do planned maintenance
                reason4Stop(i,1) = 2;
                uptimes(i,1) = tm;
                downtimes(i,1) = 0; % offline maint
                %downtimes(i,1) = MT(i,1); % online maint
            end
        elseif maintType == 2
            theta = dv;
            if fail_state(i,1) <= theta
                % do corrective maintenance
                reason4Stop(i,1) = 1;
                uptimes(i,1) = ttf(i,1);
                downtimes(i,1) = RT(i,1);                
            else
                % do planned maintenance
                reason4Stop(i,1) = 2;
                uptimes(i,1) = sum(ST(i,1:theta));
                downtimes(i,1) = MT(i,1); 
            end
        else
        end
    end
    maintCycles = [uptimes, downtimes];
end

