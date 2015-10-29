function [ U_est, V_est, numb_it ] = helmke_alg( U_init, V_init, acc, m_1, m_2, max_it )
%Estimates the essential matrix based on the algorithm provided in Helmkes
%paper.
%   INPUT: Initial guess for the essential matrix, accuracy, matched points
%   of frame1 and frame2. And max number of iteration.
if nargin <6
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

%% Project the initial Guess onto the manifold.
U_k = U_init;
V_k = V_init;

% Project the intial guess on the essential manifold.
E = U_k * E_0 * V_k';

% Calculate the gradient of the initial guess.
J = kron(V_k, U_k) * ( kron(E_0, eye(3)) * Q_1 - kron(eye(3), E_0) * Q_2);
grad = J' * M * E(:);

%% Iterate the algorithm
for i=1:max_it
    %% Calculate the Newton step
    % First calculate the Hessian.
    H_hat = J'* M * J;
    
    D_vec = kron(V_k', U_k')* M * E(:);
    D = reshape(D_vec, [3,3]);
    H_til = [Q_1', Q_2'] * [-(kron(D*E_0, eye(3))), kron(D, E_0); ...
        kron(D', E_0), -kron(E_0*D,eye(3))] * [Q_1; Q_2];
    
    H = H_hat + H_til;
    
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
    grad = J' * M * E_kk(:);
    if norm(grad) < acc
        break
    end
    % If not iterate.
    U_k = U_kk;
    V_k = V_kk;
    E = E_kk;
    numb_it = i;
end
U_est = U_kk;
V_est = V_kk;
end

