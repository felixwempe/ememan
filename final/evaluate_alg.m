function [ Estimate, Rotation, Translation ] = evaluate_alg( M1, M2, P, ...
    U_init, V_init, accuracy, algorithm, numb_frames)
%This function evaluates different types of the algorithm presented in the
%paper of Helmke. 
%   The input parameters are:
%   data_dir:   data directory for which the algorithm shall be evaluated.
%   gt_dir:     groundtruth data directory belonging to the given data.
%   sequence:   Number of sequence of kitti data.
%   U_init:     Initialization of U.
%   V_init:     Initialization of V.
%   Calib_mat:  Calibration matrix of the camera belonging to the data.
%   accuracy:   Accuracy the algorithm should attain
%   algorithm:  Which type of Helmke's algorithm is used.
%   save_dir:   Directory to save the data to.
%   numb_frames:How many data frames are used (leave empty for all)

% Estimate saves the results of each iteration.
Estimate(numb_frames) = struct;
Rotation(numb_frames) = struct;
Translation(numb_frames) = struct;
for i=1:numb_frames
    if (strcmp(algorithm, 'smooth')) && (i>1)
        [U, V, iter, time] = Helmke(U_init, V_init, M1(i).m, M2(i).m, accuracy, ...
            algorithm, U_prev, V_prev);
    else
        [U,V, iter, time] = Helmke(U_init, V_init, M1(i).m, M2(i).m, accuracy, ...
            algorithm);
    end
    if i>1
        [dist_R, dist_t, R_est, t_est] = compare_to_groundtruth(U, V, P(i-1).P, P(i).P);
    else
        [dist_R, dist_t, R_est, t_est] = compare_to_groundtruth(U,V, ...
            [eye(3), zeros(3,1);zeros(1,3), 1], P(i).P);
    end
    
    %Save results.
    Rotation(i).R = R_est;
    Translation(i).t = t_est;
    Estimate(i).U = U;
    Estimate(i).V = V;
    Estimate(i).iter = iter;
    Estimate(i).time = time;
    Estimate(i).dist_R = dist_R;
    Estimate(i).dist_t = dist_t;
end

end

