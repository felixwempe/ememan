function [U_kk, V_kk, iter] = Helmke(U_init, V_init, m_1, m_2, accuracy, ...
        algorithm, U_prev, V_prev)
%The function computed the Helmke algorithm for different types of Energy
%functions.
%   Input parameters:
%       U_init:     Initialization for U.
%       V_init:     Initialization for V.
%       M1:         matched points from image 1.
%       M2:         matched points from image 2.
%       accuracy:   Accuracy that has to be reached to converge.
%       algorithm:  Determines the Energy function that is used. Available
%                   input parameters for algorithm are: 'helmke', 'huber',
%                   'smooth'.
%   Output parameters:
%       U:      Estimated U of essential matrix E= U*E_0*V'
%       V:      Estimated V.
%       iter:   iterations needed to converge.
%       time:   time used for the algorithm.
lambda = 1e-01;

iter = 500;
delta = 1e-03;
% Matrices we need throughout the algorithm.
max_it = 20;

E_0 = [1 0 0; 0 1 0; 0 0 0];

Q_x = [0 0 0; 0 0 -1; 0 1 0];
Q_y = [0 0 1; 0 0 0; -1 0 0];
Q_z = [0 -1 0; 1 0 0; 0 0 0];

Q_1 = [ (1/sqrt(2))*Q_x(:), (1/sqrt(2))*Q_y(:), (1/2)*Q_z(:), zeros(9,2)];
Q_2 = [ zeros(9,2), -(1/2)*Q_z(:), (1/sqrt(2))*Q_x(:), (1/sqrt(2))*Q_y(:)];

U_k = U_init;
V_k = V_init;
% Get the gradient of the initial guess.
if strcmp(algorithm, 'helmke')
    M = calc_M(m_1, m_2);
    [grad, J] = get_gradient_helmke(U_k, V_k, M, Q_1, Q_2);
elseif strcmp(algorithm, 'smooth')
    M = calc_M(m_1, m_2);
    E_est = U_prev * E_0 * V_prev';
    [grad, J] = get_gradient_smooth(U_k, V_k, E_est, M, Q_1, Q_2, lambda);
else
    error('The input for algorithm is not valid.')
end

% Iterate the algorithm.
for i=1:max_it
    % Get the Hessian.
    if strcmp(algorithm, 'helmke')
        [Hessian, H_hat] = get_hessian_helmke(U_k, V_k, J, M, Q_1, Q_2);
    elseif strcmp(algorithm, 'smooth')
        [Hessian, H_hat] = get_hessian_smoothed(U_k, V_k, J, E_est, M, Q_1, Q_2, lambda);
    end
    
    % Get the opimal direction.
    if min(eig(Hessian))> accuracy
        x_opt = Hessian'\(-grad);
    else
        x_opt = H_hat'\(-grad);
    end
    
    % Project back on the essential manifold.
    U_kk = U_k*expm([0 -x_opt(3)/sqrt(2) x_opt(2);...
        x_opt(3)/sqrt(2) 0 -x_opt(1); -x_opt(2) x_opt(1) 0]);
    V_kk = V_k*expm([0 x_opt(3)/sqrt(2) x_opt(5);...
        -x_opt(3)/sqrt(2) 0 -x_opt(4); -x_opt(5) x_opt(4) 0]);
    
    % Calculate the new gradient.
    if strcmp(algorithm, 'helmke')
        [grad, J] = get_gradient_helmke(U_kk, V_kk, M, Q_1, Q_2);
    elseif strcmp(algorithm, 'smooth')
        [grad, J] = get_gradient_smooth(U_kk, V_kk, E_est, M, Q_1, Q_2, lambda);
    end
    
    % Check for convergence.
    if norm(grad) < accuracy
        iter = i;
        break
    end
%     out_txt = sprintf('At iteration: %d the norm is: %d', i, norm(grad));
%     disp(out_txt)
    % If not iterate.
    U_k = U_kk;
    V_k = V_kk;
end


