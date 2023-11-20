function [m_T, c2_T, m_W, c2_W, pm, A, m_Te, c2_Te, m_GR] = RTFPolicy( rates, alpha, m_T0, v_T0, m_R, v_R)
    % Output
    %   m_T: mean of uptime in one maintenance cycle
    %   c2_T: CV2 of uptime in one maintenance cycle
    %   m_W: mean to downtime in one maintenance cycle
    %   c2_W: CV2 of downtime in one maintenance cycle
    %   pm: probability of maintenance among all stops
    %   A: availability of machine
    %   m_Te: mean of effective process time
    %   c2_Te: CV2 of effective process time
    %   m_GR: mean of good part rate
    %
    % Input
    %   rates: matrix includes 3 rows: 
    %   1: deterioration rates, 2: failure rates, 3: good part rates
    %   alpha: prob of initial state
    %   m_T0: mean of null process time
    %   v_T0: variance of null process time
    %   m_R: mean of repair time
    %   v_R: variance of repair time
    %

    
    %% Step 1. compute the moments of uptime and downtime in one machine life cycle.
    % construct phase type distribution
    n = size(rates,2);                                          % number of states 
    S = diag(-rates(1,:)-rates(2,:))+diag(rates(1,1:end-1),1);  % phase-type dist matrix
    GR = rates(3,:);                                            % good rate in each state
    % uptime moments
    m_T = - alpha * mpower(S,-1) * ones(n,1);                   % mean of uptimes
    m_T2 = 2 * alpha * mpower(S,-2) * ones(n,1);                % mean square of uptimes
    v_T = m_T2-m_T^2;                                           % variance of uptimes
    c2_T = v_T/m_T^2;                                           % CV2 of uptimes
    % downtime moments
    pm = 0;                                                     % prob of maint when machine is stopped
    m_W = m_R;                                                  % mean of downtimes
    m_W2 = v_R+m_R^2;                                           % mean square of downtimes
    v_W = m_W2-m_W^2;                                           % variance of downtimes
    c2_W = v_W/m_W^2;                                           % CV2 of downtimes
    % interaction moments
    m_TW = m_T * m_W;                                           % interaction of uptime and downtime
    % availability
    A = m_T/(m_T+m_W);                                          % availability

    
    %% Step 2. compute the moments of effcient process time (approximated)
    [m_Te, c2_Te] = EPTApprox(m_T, v_T, m_W, v_W, m_TW, m_T0, v_T0);
    
    %% Step 3. compute the good part rate
    S0 = -S*ones(n,1);
    Q_a = S+S0.*[ones(n,1),zeros(n,n-1)];
    QE_a = [Q_a';ones(1,size(Q_a,1))];
    B_a = [zeros(size(Q_a,2),1);1];
    t_pei_a = QE_a\B_a;
    m_GR = GR*t_pei_a;
end