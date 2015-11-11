function [Hessian, H_hat] = get_hessian_huber(U_k, V_k, E_k, J, M_1, M_2, ...
            Sign, Q_1, Q_2, numb_good, numb_all, delta)
%This function computes the Hessian for the huber energy function.
E_0 = [1 0 0; 0 1 0; 0 0 0];

H_hat = J'* M_1 * J;

D_vec_1 = kron(V_k', U_k')* M_1 * E_k(:);
D_1 = reshape(D_vec_1, [3,3]);
H_til_1 = [Q_1', Q_2'] * [-(kron(D_1*E_0, eye(3))), kron(D_1, E_0);...
    kron(D_1', E_0), -kron(E_0*D_1,eye(3))] * [Q_1; Q_2];

D_vec_2 = kron(V_k', U_k')* M_2' * Sign;
D_2 = reshape(D_vec_2, [3,3]);
H_til_2 = [Q_1', Q_2'] * [-(kron(D_2*E_0, eye(3))), kron(D_2, E_0);...
    kron(D_2', E_0), -kron(E_0*D_2,eye(3))] * [Q_1; Q_2];


Hessian = (numb_good/numb_all)*(H_hat+H_til_1)+(delta/numb_all)*H_til_2;
        
end