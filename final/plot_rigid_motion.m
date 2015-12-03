function [] = plot_rigid_motion(Rotations, translations, translations_groundtruth)
%This function plots the rigid motion of the estimated rotations and
%translations.
P = [eye(3),zeros(3,1);0 0 0 1];
vec_est = P(1:3,4);
vec_gt = P(1:3,4);
for i=1:length(Rotations)
    P = P * [Rotations(i).R ,translations(i).t; 0 0 0 1];
    vec_est = [vec_est, P(1:3,4)*(norm(translations_groundtruth(i).t)/norm(P(1:3,4)))];
    vec_gt = [ vec_gt, translations_groundtruth(i).t];
end
figure
hold on
plot(vec_est(1,:), vec_est(3,:), 'r')
plot(vec_gt(1,:), vec_gt(3,:), 'b')
hold off
end