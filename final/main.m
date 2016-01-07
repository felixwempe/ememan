%% This script evaluates the helmke algorithm.
clear
%% Set the fixed values.
P0 = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00 0.000000000000e+00];

Calib_mat = P0(1:3, 1:3);

InGuess = eye(3)*[0 -1 0; 1 0 0; 0 0 0];
[U_init, ~, V_init] = svd(InGuess);


accuracy = 0.0000001;
algorithm = 'smooth';
numb_frames =4500;


%% Read the data. 
% M1 and M2 are structs containing in each cell the points from image1 and
% image2 for which the essentiell matrix shall be estimated.

load('results/matched_points_00_1.mat');
load('results/matched_points_00_2.mat');
M1=M1(2:end);
M2=M2(2:end);
%% Read the groundtruth data.
% Into a struct P containing rotation and translation groundtruth
% informations of each image pair.

load('results/gt_data_00.mat');
%% Run the algorithm on all frames.

[Estimate, Rotation, Translation] = evaluate_alg( M1, M2, P, U_init, V_init,...
    accuracy, algorithm, numb_frames);

%% Plot the results.

plot_rigid_motion(Rotation, Translation, P, numb_frames);

%% Save the results to directory.
% if ~isdir('save_dir')
%     mkdir(save_dir);
% end
% 
% save([save_dir, '/Final_res_' sequence, '_', algorithm , '.mat'], 'Estimate')
% save([save_dir, '/Final_rot_' sequence, '_', algorithm , '.mat'], 'Rotation')
% save([save_dir, '/Final_trans_' sequence, '_', algorithm , '.mat'], 'Translation')
