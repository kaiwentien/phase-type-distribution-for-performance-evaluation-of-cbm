function [m_T, c2_T, m_W, c2_W, pm, A, m_Te, c2_Te, m_GR] = TBMPolicy(rates, alpha, m_T0, v_T0, m_R, v_R, m_M, v_M, tm)
    % Output
    %   m_T: mean of uptime in one maintenance cycle
    %   c2_T: CV2 of uptime in one maintenance cycle
    %   m_W: mean to downtime in one maintenance cycle
    %   c2_W: CV2 of downtime in one maintenance cycle
    %   pm: probability of maintenance among all stops
    %   A: availability of machine
    %   m_Te: mean of effective process time
    %   c2_Te: CV2 of effective process time
    %   m_GR: mean of good rate
    %
    % Input
    %   rates: matrix includes 3 rows: 
    %   1: deterioration rates, 2: failure rates, 3: good part rates
    %   alpha: prob of initial state
    %   m_T0: mean of null process time
    %   v_T0: variance of null process time
    %   m_R: mean of repair time
    %   v_R: variance of repair time
    %   m_M: mean of maintenance time
    %   v_M: variance of maintenance time
    %   tm: maintenance interval
    
    %% Step 1. compute the moments of uptime and downtime in one machine life cycle.
    % construct phase type distribution
    n = size(rates,2);                                          % number of states 
    S = diag(-rates(1,:)-rates(2,:))+diag(rates(1,1:end-1),1);  % phase-type dist matrix
    GR = rates(3,:);                                            % good rate in each state
    % compute trucated expectation value of TTF.
    trc_F = - tm * alpha * expm(S*tm) * ones(n,1)...
            + alpha * expm(S*tm) * mpower(S,-1) * ones(n,1)...
            - alpha * mpower(S,-1)*ones(n,1);
%     trc_F2 = - power(tm,2) * alpha * expm(S*tm)*ones(n,1)...
%              + 2 * tm * alpha * expm(S*tm) * mpower(S,-1) * ones(n,1)...
%              - 2 * alpha * expm(S*tm) * mpower(S,-2) * ones(n,1)...
%              + 2 * alpha *mpower(S,-2)*ones(n,1);
    
    % uptime moments
    m_T = alpha * expm(S*tm)*mpower(S,-1)*ones(n,1)...
             - alpha * mpower(S,-1) * ones(n,1);                % mean of uptimes
    m_T2 = 2 * tm * alpha * expm(S*tm) * mpower(S,-1)* ones(n,1)...
            - 2 * alpha * expm(S*tm) * mpower(S,-2) * ones(n,1)...
            + 2 * alpha *mpower(S,-2)*ones(n,1);                % mean square of uptimes
    v_T = m_T2-m_T^2;                                           % variance of uptimes
    c2_T = v_T/m_T^2;                                           % CV2 of uptimes
    
    % downtime moments
    pm = alpha * expm(S*tm) * ones(n,1);                        % prob of maint when machine is stopped
    m_W = pm * m_M + (1-pm)* m_R;                               % mean of downtimes
    m_W2 = pm * (v_M+m_M^2) + (1-pm)* (v_R+m_R^2);              % mean square of downtimes
    v_W = m_W2-m_W^2;                                           % variance of downtimes
    c2_W = v_W/m_W^2;                                           % CV2 of downtimes
    % interaction moments
    m_TW = tm * m_M * pm + m_R * trc_F;                         % interaction of uptime and downtime
    % availability
    A = m_T/(m_T + m_W);                                        % availability
    
    %% Step 2. compute the moments of effcient process time (approximated)
    [m_Te, c2_Te] = EPTApprox(m_T, v_T, m_W, v_W, m_TW, m_T0, v_T0);
    
    %% Step 3. compute the good part rate
    S0 = -S*ones(n,1);
    Q_a = S+S0.*[ones(n,1),zeros(n,n-1)];
    p_m = expm(Q_a*tm);
    t_pei_a = (alpha*p_m)';
    m_GR = GR*t_pei_a;
end