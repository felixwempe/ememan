function [ Estimate, Rotation, Translation ] = evaluate_alg( M1, M2, P, ...
    U_init, V_init, accuracy, algorithm, numb_frames, bol_ransac)
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
%   numb_frames:How many data frames are used 
%   bol_ransac: Boolean to descide if ransac is used or not.
bol = exist('bol_ransac');
acc = 1e-04;
if ~bol
    bol_ransac = false;
end

% Estimate saves the results of each iteration.
Estimate(numb_frames) = struct;
Rotation(numb_frames) = struct;
Translation(numb_frames) = struct;
for i=1:numb_frames
    tic
    if bol_ransac
        for j=1:50
            [m_1, m_2] = ransacs(M1(i).m, M2(i).m, acc*10^(j-1));
            if length(m_1)>50
                break
            end
        end
    else
        m_1 = M1(i).m;
        m_2 = M2(i).m;
    end
    if (strcmp(algorithm, 'smooth')) && (i>1)
        [U, V, iter] = Helmke(U_init, V_init, m_1, m_2, accuracy, ...
            algorithm, U_prev, V_prev);
        U_prev =U;
        V_prev =V;
    elseif strcmp(algorithm, 'smooth')
        [U,V, iter] = Helmke(U_init, V_init, m_1, m_2, accuracy, ...
            'helmke');%first iteration of smooth algorithm.
        U_prev =U;
        V_prev =V;
    else
        [U,V, iter] = Helmke(U_init, V_init, m_1, m_2, accuracy, ...
            algorithm);
    end
    time = toc;
    
    [dist_R, dist_t, R_est, t_est] = compare_to_groundtruth(U, V,...
        P(i).P, P(i+1).P);
    
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

