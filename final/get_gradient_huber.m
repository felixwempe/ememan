function [grad, J, Sign] = get_gradient_huber(U_k, V_k, E_k, M_1, M_2, ...
    numb_good, numb_all, Q_1, Q_2, delta)
% This function computes the gradient for the huber energy function.
E_0 = [1 0 0; 0 1 0; 0 0 0];

Sign = sign(M_2 * E_k(:));
J = kron(V_k, U_k) * ( kron(E_0, eye(3)) * Q_1 - kron(eye(3), E_0) * Q_2);

grad = (numb_good/numb_all)* J' * M_1 * E_k(:) + (delta/numb_all) * ...
    J' * M_2' * Sign;

end