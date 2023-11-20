function [m_Te, c2_Te] = EPTApprox(m_T, v_T, m_W, v_W, m_TW, m_T0, v_T0, index)
    %% check number of input variables
    if nargin < 7
       error('Not enough input arguments')
    elseif nargin==7
       % use the best approximation
       index = 1;
    elseif nargin > 8
       error('Too many input arguments')
    end  
    % we have three types of approximation now...
    switch index
        
        case 1  % Approx in Factory Physics: assume the uptime distribution is exponential.
            A = m_T/(m_W+m_T);
            m_Te = m_T0/A;
            v_Te = v_T0/A^2 + (m_W^2+v_W)*(1-A)*m_T0/(A*m_W);
            c2_Te = v_Te/m_Te^2;
            
        case 2  % Approx in Brown&Solomen's paper: assume the process time go to infinity.
            m_T2 = v_T + m_T^2;
            m_W2 = v_W + m_W^2;
            cof_a = power(m_T,-1)*m_W;                                           % coefficient of mean
            cof_c = m_T^(-3)*m_T2*m_W^2 - 2*m_T^(-2)*m_TW*m_W + m_T^(-1)*m_W2;   % coefficient of variance
            m_Te = (1 + cof_a) * m_T0;
            v_Te = (1 + cof_a) * v_T0 + cof_c * m_T0;
            c2_Te = v_Te/m_Te^2;
            
        case 3 % self-define approximation, may fit when maintenance and repair is an rare event, using residual life distribution
            p_rul = 1 - alpha*expm(S*m_T0)*S^-1*ones(n,1)/(alpha*S^-1*ones(n,1))    % compute prob that the machine will fail in next m_T0
            m_Te = m_T0*(1-p_rul) + (m_T0 + m_R)*p_rul
            m_Te2 = m_T0^2*(1-p_rul) + (m_T0 + m_R)^2*p_rul
            v_Te = m_Te2-m_Te^2
            c2_Te = v_Te/m_Te^2;
        
        otherwise
            m_Te = 0;
            c2_Te = 0;
            disp('no policies determined...');
    end
            

    %% Approx in Tien's: assume the maintenance and repair is a rare event.
end