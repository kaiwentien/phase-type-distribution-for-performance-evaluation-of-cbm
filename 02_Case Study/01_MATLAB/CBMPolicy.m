function [m_T, c2_T, m_W, c2_W, pm, A, m_Te, c2_Te, m_GR] = CBMPolicy(rates, alpha, m_T0, v_T0, m_R, v_R, m_M, v_M, th)
    % Output
    %   m_T: mean of uptime in one maintenance cycle
    %   c2_T: CV2 of uptime in one maintenance cycle
    %   m_W: mean to downtime in one maintenance cycle
    %   c2_W: CV2 of downtime in one maintenance cycle
    %   pm: probability of maintenance among all stops
    %   A: availability of machine
    %   m_Te: mean of effective process time
    %   c2_Te: CV2 of effective process time
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
    %   th: state threshold (once encounter th+1 states, stop machine)
    %   Note: this is only used for calculate Coxian deterioration model
    %% Step 1. compute the moments of uptime and downtime in one machine life cycle.
    % construct phase type distribution
    n = size(rates,2);                                          % number of states 
    S = diag(-rates(1,:)-rates(2,:))+diag(rates(1,1:end-1),1);  % phase-type dist matrix
    GR = rates(3,:);                                            % good rate in each state
    % Compute new matrix with respect to the threshlod state
    n_cbm = th;
    S_cbm = S(1:th,1:th);                                       % new phasetype matrix when CBM
    alpha_cbm = alpha(1,1:th);                                  % new initial distribution when CBM
    cum_alpha = cumsum(alpha,2,'reverse');                      % cumsum alpha
    alpha_cbm(th) = cum_alpha(th);                              % make sure sum of alpha is 1
    % uptime moments
    m_T = - alpha_cbm * mpower(S_cbm,-1) * ones(n_cbm,1);       % mean of uptimes
    m_T2 = 2* alpha_cbm * mpower(S_cbm,-2) * ones(n_cbm,1);     % mean square of uptimes
    v_T = m_T2 - m_T^2;                                         % variance of uptimes
    c2_T = v_T/m_T^2;                                           % CV2 of uptimes
    % downtime moments
    %pm = prod(rates(1,1:th)./(rates(1,1:th)+rates(2,1:th)));    % prob of maint when machine is stopped
    tran_prob = rates(1,1:th)./(rates(1,1:th)+rates(2,1:th));
    pm = 0;
    for i = 1:th
        pm = pm + alpha_cbm(i)*prod(tran_prob(i:th));
    end
    m_R2 = v_R + m_R^2;                                         % mean square of repair times
    m_M2 = v_M + m_M^2;                                         % mean square of maint times
    m_W = pm * m_M + (1-pm)* m_R;                               % mean of downtimes
    m_W2 = pm * m_M2 + (1-pm)* m_R2;                            % mean square of downtimes
    v_W = m_W2 - m_W^2;                                         % variance of downtimes
    c2_W = v_W/m_W^2;                                           % CV2 of downtimes
    % interaction of moments
    m_TW = m_T * m_W;                                           % interaction of uptime and downtime
    % availability
    A = m_T/(m_T + m_W);                                        %availability
    
    %% Step 2. compute the moments of effcient process time (approximated)
    [m_Te, c2_Te] = EPTApprox(m_T, v_T, m_W, v_W, m_TW, m_T0, v_T0);
    
    %% Step 3. compute the good part rate
%     S0 = -S_cbm*ones(th,1);
%     Q_a = S_cbm+S0.*[ones(th,1),zeros(th,th-1)];
%     QE_a = [Q_a';ones(1,size(Q_a,1))];
%     B_a = [zeros(size(Q_a,2),1);1];
%     t_pei_a = QE_a\B_a;
%     m_GR = GR(1:th)*t_pei_a;
    m_GR = 1;
end