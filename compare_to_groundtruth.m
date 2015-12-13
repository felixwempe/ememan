function [ dist_R, dist_t, R_est, t_est ] = compare_to_groundtruth( U_est, V_est, P, P_prev)
%The function compares for a given essential matrix the rotation and
%translation to the groundtruth data (P=(R|t)).
%   Input values are the estimated rotation matrices U and V, s.t.
%   E=U*E_0*V', the previous translation and rotation and the groundtruth
%   data.
%   Output is the distance between the two rotation matrices measured as
%   ||log(R_1'*R_2)||_F

% Some constant matrices used.

E_0 = [1 0 0; 0 1 0; 0 0 0];
W = [0 -1 0; 1 0 0; 0 0 1];

% First calculate the rotation matrix of the estimated essential matrix.
cross_t = U_est * W * E_0 * U_est';
t_est = [cross_t(3,2); cross_t(1,3); cross_t(2,1)];
R_est_1 = U_est * W' * V_est';
R_est_2 = U_est * W * V_est';

% if nargin==6
%     % Find the right translation vector (comparing to the previous motion).
%     if norm(t_est-t_prev)>norm(t_est+t_prev)
%         t_est = -t_est;
%     end
%     
%     
%     % Find the right rotation matrix (comparing to the previous motion).
%     dist_1 = norm(logm(R_est_1'*R_prev),'fro');
%     dist_2 = norm(logm(R_est_2'*R_prev),'fro');
%     
%     if dist_1<dist_2
%         R_est = R_est_1;
%     else
%         R_est = R_est_2;
%     end
% else
%     R_est = R_est_1;
% end

if t_est(3)<0
    t_est = -t_est;
end

ang_R_1 = norm(logm(R_est_1),'fro');
ang_R_2 = norm(logm(R_est_2),'fro');

if ang_R_1 < ang_R_2
    R_est = R_est_1;
else
    R_est = R_est_2;
end

% Then compare it to the groundtruth rotation matrix.
R_gt_prev = P_prev(1:3,1:3);
t_gt_prev = P_prev(1:3,4);
R_gt_act = P(1:3,1:3);
t_gt_act = P(1:3,4);

R_gt = R_gt_prev' * R_gt_act;
t_gt = R_gt_prev' * (-t_gt_act + t_gt_prev);

dist_R = norm(logm(R_est'*R_gt),'fro');
% if dist_R>1
%     R_est_1 = U_est * W * V_est';
%     dist_R = norm(logm(R_est_1'*R_gt),'fro');
% 
% end
% Furthermore compare the translation vector.
dist_t = acos(t_est'*t_gt/(norm(t_est)*norm(t_gt)));
% if dist_t>1
%     t_est = -t_est;
%     dist_t = pi - acos((t_est.*t_gt)'*ones(3,1)/(norm(t_est)*norm(t_gt)));
% 
% end
end

