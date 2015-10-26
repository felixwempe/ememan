function [ EssMatGuess, numb_iter_used ] = Helmke_Huber_alg( InGuess, acc, m_1, m_2, delta, max_it )
%Estimates the essential matrix based on the algorithm provided in Helmkes
%paper.
%   INPUT: Initial guess for the essential matrix, accuracy, matched points
%   of frame1 and frame2. And max number of iteration.
% For the Huber loss function the function to be minimized consists of a
% sum, where the data points are split on both parts of the sum. therefore
% we have to split the points at each iteration.
if nargin <6
    max_it = 20;
end
if nargin <5
    delta = 1;
end
% Fixed matrices that are used throughout the algorithm.
E_0 = [1 0 0; 0 1 0; 0 0 0];

Q_x = [0 0 0; 0 0 -1; 0 1 0];
Q_y = [0 0 1; 0 0 0; -1 0 0];
Q_z = [0 -1 0; 1 0 0; 0 0 0];

Q_1 = [ (1/sqrt(2))*Q_x(:), (1/sqrt(2))*Q_y(:), (1/2)*Q_z(:), zeros(9,2)];
Q_2 = [ zeros(9,2), -(1/2)*Q_z(:), (1/sqrt(2))*Q_x(:), (1/sqrt(2))*Q_y(:)];

numb_points = length(m_1);
%% Project the initial guess onto the manifold.
% Calculate the single value decomposition of the initial guess.
[U_k, ~, V_k] = svd(InGuess);

% Project the intial guess on the essential manifold.
E = U_k * E_0 * V_k';

% Here we split the data points the first time.

[M_1, M_2, numb_good] = calc_Huber_M(E, m_1, m_2, delta);

% Calculate the gradient of the initial guess. Therefore we need to know
% the sign of the "bad" points.
Sign = sign(M_2 * E(:));
J = kron(V_k, U_k) * ( kron(E_0, eye(3)) * Q_1 - kron(eye(3), E_0) * Q_2);
grad = (numb_good/numb_points)* J' * M_1 * E(:) + (delta/numb_points) * ...
    J' * M_2' * Sign;

%% Iterate the algorithm
for i=1:max_it
    %% Calculate the Newton step
    % First calculate the Hessian.
    H_hat = J'* M_1 * J;
    
    D_vec_1 = kron(V_k', U_k')* M_1 * E(:);
    D_1 = reshape(D_vec_1, [3,3]);
    H_til_1 = [Q_1', Q_2'] * [-(kron(D_1*E_0, eye(3))), kron(D_1, E_0); kron(D_1', E_0), -kron(E_0*D_1,eye(3))] * [Q_1; Q_2];
    
    D_vec_2 = kron(V_k', U_k')* M_2' * Sign;
    D_2 = reshape(D_vec_2, [3,3]);
    H_til_2 = [Q_1', Q_2'] * [-(kron(D_2*E_0, eye(3))), kron(D_2, E_0); kron(D_2', E_0), -kron(E_0*D_2,eye(3))] * [Q_1; Q_2];

    
    H = H_hat + H_til_1 + H_til_2;
    
    % Check if Hessian close to singular.
    if min(eig(H))> acc
        x_opt = H'\(-grad);
%         disp('Newton')
    else
        x_opt = H_hat'\(-grad);
%         disp('Gauss-Newton')
    end
   
    %% Project back on the essential manifold.
    U_kk = U_k*expm([0 -x_opt(3)/sqrt(2) x_opt(2); x_opt(3)/sqrt(2) 0 -x_opt(1); -x_opt(2) x_opt(1) 0]);
    V_kk = V_k*expm([0 x_opt(3)/sqrt(2) x_opt(5); -x_opt(3)/sqrt(2) 0 -x_opt(4); -x_opt(5) x_opt(4) 0]);
    E_kk = U_kk * E_0 * V_kk';
    
   
    %% Check if convergence criteria is fulfilled.
    % Calculate the new gradient.
    [M_1, M_2, numb_good] = calc_Huber_M(E_kk, m_1, m_2, delta);
    Sign = sign(M_2 * E_kk(:));
    
    J = kron(V_kk, U_kk) * (kron(E_0, eye(3))*Q_1 - kron(eye(3), E_0) * Q_2);
    grad = (numb_good/numb_points)* J' * M_1 * E_kk(:) + (delta/numb_points) * ...
    J' * M_2' * Sign;

    if norm(grad) < acc
        break
    end
    % If not terminated iterate.
    U_k = U_kk;
    V_k = V_kk;
    E = E_kk;
    numb_iter_used = i;
end
% fprintf('Number of iterations used to converge: %s \n', i)
EssMatGuess = E_kk;
end

