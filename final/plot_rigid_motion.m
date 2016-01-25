function [ret_est, ret_gt] = plot_rigid_motion(Rotations, translations, P_gt, numb_frames, seq, lambda)
%This function plots the rigid motion of the estimated rotations and
%translations.
P_est = [eye(3),zeros(3,1);0 0 0 1];
vec_est = P_est(1:3,4);
vec_gt = P_est(1:3,4);
for i=1:numb_frames
    if i>1
        translation_gt_rel = P_gt(i-1).P(1:3,1:3)'*(P_gt(i).P(1:3,4)-P_gt(i-1).P(1:3,4));
    else
        translation_gt_rel = P_gt(i).P(1:3,4);
    end
    P_est = P_est * [Rotations(i).R ,translations(i).t*(norm(translation_gt_rel)/norm(translations(i).t)); 0 0 0 1];
    vec_est = [vec_est, P_est(1:3,4)];
    vec_gt = [ vec_gt, P_gt(i).P(1:3,4)];
end
figure
hold on
plot(vec_est(1,:), vec_est(3,:), 'r', 'LineWidth', 3)
plot(vec_gt(1,:), vec_gt(3,:), 'b', 'LineWidth', 3)
% title('Rigid motion')
h_legend = legend('estimated', 'groundtruth');
set(h_legend,'FontSize',20);

ret_est = vec_est([1,3],:);
ret_gt = vec_gt([1,3],:);
hold off
% Save coordinates to files.
fileid = fopen(['plots/Seq', num2str(seq) ,'ransac', lambda, 'est.dat'], 'w');
fprintf(fileid, '%f \t %f\n', ret_est);
fclose(fileid);

fileid = fopen(['plots/Seq', num2str(seq) ,'ransac', lambda, 'gt.dat'], 'w');
fprintf(fileid, '%f \t %f\n', ret_gt);
fclose(fileid);
end