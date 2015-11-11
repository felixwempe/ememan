function  [Hessian, H_hat] = get_hessian_helmke(U_k, V_k, J, M, Q_1, Q_2)
% Computes the Hessian.
E_0 = [1 0 0; 0 1 0; 0 0 0];

H_hat = J'* M * J;
E = U_k * E_0 * V_k';

D_vec = kron(V_k', U_k')* M * E(:);
D = reshape(D_vec, [3,3]);
H_til = [Q_1', Q_2'] * [-(kron(D*E_0, eye(3))), kron(D, E_0); ...
    kron(D', E_0), -kron(E_0*D,eye(3))] * [Q_1; Q_2];

Hessian = H_hat + H_til;

end