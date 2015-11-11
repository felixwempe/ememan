function [] = plot_rigid_motion(R, t, t_gt)
%This function plots the rigid motion of the estimated rotations and
%translations.
P = [eye(3),zeros(3,1);0 0 0 1];
vec_est = P(1:3,4);
vec_gt = P(1:3,4);
for i=1:length(R)
    P = P* [R(i).R ,t(i).t; 0 0 0 1];
    vec_est = [vec_est, P(1:3,4)*(norm(t_gt(i).t)/norm(P(1:3,4)))];
    vec_gt = [ vec_gt, t_gt(i).t];
end
figure
hold on
plot(vec_est(1,:), vec_est(2,:), 'r')
plot(vec_gt(1,:), vec_gt(2,:), 'b');
hold off
end