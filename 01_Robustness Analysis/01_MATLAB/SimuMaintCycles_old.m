function [maintCycles,reason4Stop] = SimuMaintCycles( data, reps, maintType, dv, seed, cv)
    %SimuMaintCycles :Simulate maintenance cycles using Coxian distribution
    % using [deterioration times, failure times]
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
    load(data);                             % load data
    n = size(rate,2);                       % number of stages of HI
    rng(seed);                              % random seed
    tm = Inf;                               % maintenance interval
    theta = 6;                              % maintenance threshold
    start_state = 1;                        % starting state
    maint_state = 0;                        % maintenance state
    time_in_states = zeros(reps, n);        % time in each HI states
    uptimes = zeros(reps, 1);                   % uptimes
    downtimes = zeros(reps, 1);                   % downtimes
    reason4Stop = zeros(reps,1);            % 1: corrective, 2: preventive
    start_state_count = zeros(1,n);         % count the states to start
    ipf = 0;                                % imperfect factor 0-perfect
    MT = gamrnd(m_M^2/v_M, v_M/m_M,reps,1); % sampling of PT time
    RT = gamrnd(m_R^2/v_R, v_R/m_R,reps,1); % sampling of CT time
    a_1 = repmat(1/cv^2, reps, size(rate,2));   % scale of gamma - DT  
    b_1 = repmat(cv^2*1./rate(1,:), reps, 1);   % rate of gamma - DT
    a_2 = repmat(1/cv^2, reps, size(rate,2));   % scale of gamma - FT
    b_2 = repmat(cv^2*1./rate(2,:), reps, 1);   % rate of gamma - FT
    DT = gamrnd(a_1,b_1,[reps,size(rate,2)]);   % sampling of DT
    FT = gamrnd(a_2,b_2,[reps,size(rate,2)]);   % sampling of FT
    %
    if maintType == 1
        tm = dv;
        for i = 1:reps
            is_failed = false;
            is_maint = false;
            for j = start_state:n
                if ~is_failed && ~ is_maint
                    if DT(i,j) > FT(i,j)
                        time_in_states(i,j) = FT(i,j);
                        uptimes(i,1) = uptimes(i,1) + FT(i,j);
                        if uptimes(i,1) >= tm
                            is_maint = true;
                            maint_state = j;
                            reason4Stop(i,1) = 2;
                            time_in_states(i,j) = tm - (uptimes(i,1) - FT(i,j));
                            uptimes(i,1) = tm;
                        else
                            is_failed = true;
                            maint_state = n;
                            reason4Stop(i,1) = 1;
                        end
                    else
                        time_in_states(i,j) = DT(i,j);
                        uptimes(i,1) = uptimes(i,1) + DT(i,j);
                        if uptimes(i,1) >= tm             
                            is_maint = true;
                            maint_state = j;
                            reason4Stop(i,1) = 2;
                            time_in_states(i,j) = tm - (uptimes(i,1) - DT(i,j));
                            uptimes(i,1) = tm;
                        else
                        end
                    end
                else
                end
                % offline maintenance
                downtimes(i,1) = (reason4Stop(i,1)==1)*RT(i,1)+(reason4Stop(i,1)==2)*0;
            end
            start_state = binornd(maint_state-1,ipf)+1;
            start_state_count(start_state) = start_state_count(start_state) + 1;
        end
    elseif maintType == 2
        theta = dv;
        for i = 1:reps
            is_failed = false;
            is_maint = false;
            if start_state<=theta
                for j = start_state:n
                    if ~is_failed && ~ is_maint
                        if DT(i,j) > FT(i,j)
                            time_in_states(i,j) = FT(i,j);
                            uptimes(i,1) = uptimes(i,1) + FT(i,j);
                            is_failed = true;
                            maint_state = n;
                            reason4Stop(i,1) = 1;
                        else
                            time_in_states(i,j) = DT(i,j);
                            uptimes(i,1) = uptimes(i,1) + DT(i,j);
                            if j>=theta            
                                is_maint = true;
                                maint_state = j;
                                reason4Stop(i,1) = 2;
                            else
                            end
                        end
                    else
                    end
                    downtimes(i,1) = is_failed*RT(i,1)+is_maint*MT(i,1);    
                end
            else
                uptimes(i,1)=0;
                downtimes(i,1)=MT(i,1);
                maint_state = start_state;
            end
            start_state = binornd(maint_state-1,ipf)+1;
            start_state_count(start_state) = start_state_count(start_state) + 1;
        end
    else
        uptimes = 0;
        downtimes = 0;
    end
    maintCycles = [uptimes, downtimes];
end

