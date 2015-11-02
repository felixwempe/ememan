function [ U_est, V_est, numb_it ] = helmke_alg_smoothed( U_init, V_init,...
    U_prev, V_prev, acc, m_1, m_2, max_it )
%% Smoothed Helmke algorithm., U_prev_prev, V_prev_prev
% Estimates the essential matrix based on the algorithm provided in 
% Helmkes paper. Using additionally the information of the step before to
% trajectory.
%   INPUT: Initial guess for the essential matrix, Essential matrix of the 
%   2 steps before, accuracy, matched points of frame1 and frame2. And max 
%   number of iteration.
if nargin <10
    max_it = 20;
end
% Fixed matrices that are used throughout the algorithm.
E_0 = [1 0 0; 0 1 0; 0 0 0];

Q_x = [0 0 0; 0 0 -1; 0 1 0];
Q_y = [0 0 1; 0 0 0; -1 0 0];
Q_z = [0 -1 0; 1 0 0; 0 0 0];

Q_1 = [ (1/sqrt(2))*Q_x(:), (1/sqrt(2))*Q_y(:), (1/2)*Q_z(:), zeros(9,2)];
Q_2 = [ zeros(9,2), -(1/2)*Q_z(:), (1/sqrt(2))*Q_x(:), (1/sqrt(2))*Q_y(:)];

M = calc_M(m_1, m_2);

%% Calculate the Estimated Essentiel Matrix.
% Omeg_1 = logm(U_prev_prev' * U_prev);
% Omeg_2 = logm(V_prev_prev * V_prev');
% 
% % Get the direction from Ess_prev_prev to Ess_prev.
% x_est = sqrt(2) * [Omeg_1(3,2); Omeg_1(1,3); (sqrt(2)/2)*(Omeg_1(2,1) + Omeg_2(2,1)); ...
%     Omeg_2(2,3); Omeg_2(3,1)];
% 
% % Use the direction and progress it from Ess_prev to get an estimation for
% % the Essential matrix.
% U_est = U_prev*expm([0 -x_est(3)/sqrt(2) x_est(2);...
%     x_est(3)/sqrt(2) 0 -x_est(1); -x_est(2) x_est(1) 0]);
% V_est = V_prev*expm([0 x_est(3)/sqrt(2) x_est(5);...
%     -x_est(3)/sqrt(2) 0 -x_est(4); -x_est(5) x_est(4) 0]);
% 
% Ess_est = U_est * E_0 * V_est';

Ess_est = U_prev * E_0 * V_prev';

%% First projection.
U_k = U_init;
V_k = V_init;

% Project the intial guess on the essential manifold.
E_k = U_k * E_0 * V_k';

% Calculate the gradient of the initial guess. Since the cost function
% changed to f_(E) = f(E)+D(E,Ess_est), where D is a second order distance
% function, we have to change the gradient (compared to the original Helmke
% algorithm).
J = kron(V_k, U_k) * ( kron(E_0, eye(3)) * Q_1 - kron(eye(3), E_0) * Q_2);
grad = J' * M * E_k(:) + 2*J'*(E_k(:)-Ess_est(:));

%% Iterate the algorithm
for i=1:max_it
    %% Calculate the Newton step.
    % First calculate the Hessian.
    H_hat = J'* M * J;
    
    D_vec = kron(V_k', U_k')* M * E_k(:);
    D = reshape(D_vec, [3,3]);
    H_til = [Q_1', Q_2'] * [-(kron(D*E_0, eye(3))), kron(D, E_0);...
        kron(D', E_0), -kron(E_0*D,eye(3))] * [Q_1; Q_2];
    
    % Also we have to change the Hessian. Here the notation is very
    % similiar, since also the derivatives are similiar.
    H_D_hat = 2* (J'*J);
    
    D_D_vec = kron(V_k', U_k')* (E_k(:)-Ess_est(:));
    D_D = reshape(D_D_vec, [3,3]);
    H_D_til = [Q_1', Q_2'] * [-(kron(D_D*E_0, eye(3))), kron(D_D, E_0);...
        kron(D_D', E_0), -kron(E_0*D_D,eye(3))] * [Q_1; Q_2];

    H = H_hat + H_til + H_D_hat + H_D_til;
    
    % Check if Hessian close to singular.
    if min(eig(H))> acc
        x_opt = H'\(-grad);
    else
        x_opt = H_hat'\(-grad);
    end
   
    %% Project back on the essential manifold.
    U_kk = U_k*expm([0 -x_opt(3)/sqrt(2) x_opt(2);...
        x_opt(3)/sqrt(2) 0 -x_opt(1); -x_opt(2) x_opt(1) 0]);
    V_kk = V_k*expm([0 x_opt(3)/sqrt(2) x_opt(5);...
        -x_opt(3)/sqrt(2) 0 -x_opt(4); -x_opt(5) x_opt(4) 0]);
    
    E_kk = U_kk * E_0 * V_kk';
    
   
    %% Check if convergence criteria is fulfilled.
    % Calculate the new gradient.
    J = kron(V_kk, U_kk) * (kron(E_0, eye(3))*Q_1 - kron(eye(3), E_0) * Q_2);
    grad = J' * M * E_kk(:) + 2*J'*(E_kk(:)-Ess_est(:));
    if norm(grad) < acc
        break
    end
    % If not iterate.
    U_k = U_kk;
    V_k = V_kk;
    E_k = E_kk;
    numb_it =i;
end
U_est = U_kk;
V_est = V_kk;
end

