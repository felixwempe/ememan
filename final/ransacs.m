function [ m_1, m_2 ] = ransacs( m_1, m_2, acc )
%RANSAC algorithm for the estimation of the essential matrix.
%   
% Initialisation for U,V
InGuess = eye(3)*[0 -1 0; 1 0 0; 0 0 0];
[U_init, ~, V_init] = svd(InGuess);
E_0 = [1 0 0; 0 1 0; 0 0 0];

bool = false(length(m_1),1);

% Iterate the Ransac algorithm.

for i=1:50
    % Get 10 random points.
    test = ceil(rand(10,1)*length(m_1));
    m_1_test = m_1(test,:);
    m_2_test = m_2(test,:);
    % Get the estimated essential matrix.
    [U_est, V_est, ~] = Helmke(U_init, V_init, m_1_test, m_2_test, 5e-06, 'helmke');
    E_est = U_est * E_0 * V_est';
    
    % Compute how appropriate the estimation is.
    test_res = m_1* E_est * m_2';
    test_res = diag(test_res);
    test_bol = abs(test_res)<acc;
    if sum(test_bol)>100 %if many points match the estimation add to consensus set.
        bool = bool|test_bol;
    end
    % Iterate.
    
end
% Return the good points.
m_1 = m_1(bool,:);
m_2 = m_2(bool,:);

end

