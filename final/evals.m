%% Evaluate everything.

% Fixed parameters.
P0 = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00 0.000000000000e+00];

Calib_mat = P0(1:3, 1:3);

InGuess = eye(3)*[0 -1 0; 1 0 0; 0 0 0];
[U_init, ~, V_init] = svd(InGuess);

algorithms = {'helmke';'huber';'smooth'};
rans = [true, false];


accuracy = 0.0000001;

%Iteration over sequence, algorithm, ransac.
for i=1:3
    % Get matched points and groundtruth data for sequence i-1.
    load(['results/gt_data_0', num2str(i-1), '.mat']);
    load(['results/matched_points_0', num2str(i-1), '_1.mat']);
    load(['results/matched_points_0', num2str(i-1), '_2.mat']);
    M1=M1(2:end);
    M2=M2(2:end);
    numb_frames = length(M1);
    for j=1:3
        % Get algorithm
        algorithm = algorithms(j,:);
        for k=1:2
            bol_rans = rans(k);
            [Estimate, Rotation, Translation] = evaluate_alg( M1, M2, P, ...
                U_init, V_init, accuracy, cellstr(algorithm), numb_frames, bol_rans);
            % Save the results.
            save(['eval/results_seq0', num2str(i-1),'_', char(algorithm), '_rans_'...
                , num2str(k), '.mat'], 'Estimate', 'Rotation', 'Translation');
            clear Estimate Rotation Translation
        end
    end
    clear M1 M2 P
end