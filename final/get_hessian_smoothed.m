function [Hessian, H_hat] = get_hessian_smoothed(U_k, V_k, J, E_est, M, Q_1, Q_2)
%This function computes the hessian of the smoothed energy function.
E_0 = [1 0 0; 0 1 0; 0 0 0];

E_k = U_k * E_0 * V_k';

H_hat = J'* M * J;

D_vec = kron(V_k', U_k')* M * E_k(:);
D = reshape(D_vec, [3,3]);
H_til = [Q_1', Q_2'] * [-(kron(D*E_0, eye(3))), kron(D, E_0);...
    kron(D', E_0), -kron(E_0*D,eye(3))] * [Q_1; Q_2];

H_D_hat = 2* (J'*J);

D_D_vec = kron(V_k', U_k')* (E_k(:)-E_est(:));
D_D = reshape(D_D_vec, [3,3]);
H_D_til = [Q_1', Q_2'] * [-(kron(D_D*E_0, eye(3))), kron(D_D, E_0);...
    kron(D_D', E_0), -kron(E_0*D_D,eye(3))] * [Q_1; Q_2];

Hessian = H_hat + H_til + H_D_hat + H_D_til;

H_hat = H_hat + H_D_hat;

end