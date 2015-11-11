function [grad, J] = get_gradient_helmke(U_k, V_k, M, Q_1, Q_2)
%This function computes the gradient for the helmke algorithm.
E_0 = [1 0 0; 0 1 0; 0 0 0];
E = U_k * E_0 * V_k';

J = kron(V_k, U_k) * ( kron(E_0, eye(3)) * Q_1 - kron(eye(3), E_0) * Q_2);
grad = J' * M * E(:);

end