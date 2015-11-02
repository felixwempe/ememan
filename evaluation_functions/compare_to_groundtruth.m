function [ dist_R, dist_t ] = compare_to_groundtruth( U_est, V_est, P_groundtruth )
%The function compares for a given essential matrix the rotation and
%translation to the groundtruth data (P=(R|t)).
%   Input values are the estimated rotation matrices U and V, s.t.
%   E=U*E_0*V' and the groundtruth data.
%   Output is the distance between the two rotation matrices measured as
%   ||log(R_1'*R_2)||_F

% Some constant matrices used.

E_0 = [1 0 0; 0 1 0; 0 0 0];
W = [0 -1 0; 1 0 0; 0 0 1];

% First calculate the rotation matrix of the estimated essential matrix.
cross_t = U_est * W * E_0 * U_est';
t_est = [cross_t(3,2); cross_t(1,3); cross_t(2,1)];
R_est = U_est * W' * V_est';

% Then compare it to the groundtruth rotation matrix.
R_gt = P_groundtruth(1:3,1:3);
t_gt = P_groundtruth(1:3,4);

dist_R = norm(logm(R_est'*R_gt),'fro');

dist_t = acos((t_est.*t_gt)'*ones(3,1)/(norm(t_est)*norm(t_gt)));
% Furthermore compare the translation vector.

end

