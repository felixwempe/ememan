%% Skript to evaluate the speed, accuracy, etc of Helmke paper on a set of images.

% Set the values: Accuracy, maximal Iterations per step, Calibration
% matrix of the frames.
accuracy = 0.000001;


P0 = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00 0.000000000000e+00];


Calib = P0(1:3, 1:3);
Calib_inv = inv(Calib);

% Get the initialisation for the essential matrix.
% (Here it is a trivial initialisation, assuming the car always moves
% forward) E = R*[t], for the rotation matrix R=Id and a t=[1 0 0])
InGuess = eye(3)*[0 -1 0; 1 0 0; 0 0 0];
[U_init, ~, V_init] = svd(InGuess);

%% First we have to read in the whole set of images. 
% (Maybe only read in a subset and then on the flow more images, or use a 
% loop to read images only when they are used and forget about the old 
% ones)

% Get all filenames of a specified directory.
direct = '/home/felix/Universität/Praktikum_Andreas_Neufeld/dataset/sequences/00/image_0/';
filenames = dir(direct);
get_groundtruth;    % Reads the groundtruth data to a struct called P

% Read in the first image. And compute the Sift features.
frame1 = single(imread([direct, filenames(3).name]));
[point1, descr1] = vl_sift(frame1);

Estimate(length(filenames)-2) = struct;
% Start a loop over all images.
for i=4:length(filenames)-5
    %% Read in the current image.
    frame2 = single(imread([direct, filenames(i).name]));

    %% Compute the Sift features and match them.
    [point2, descr2] = vl_sift(frame2);
    [matches, score] = vl_ubcmatch(descr1, descr2);
    % Select 200 best matches.
    [~, idx] = sort(score);
    good_matches = matches(:,idx(1:200))';
    
    m_1 = point1(1:2,good_matches(:,1))';
    m_2 = point2(1:2,good_matches(:,2))';
    
    % Use the calibration matrix to normalize the image points.
    m_1_norm = (Calib_inv * [m_1, ones(length(m_1),1)]')';
    m_2_norm = (Calib_inv * [m_2, ones(length(m_2),1)]')';

    % Use the Helmke algorithm to compute a good essential matrix. Save the
    % Essential matrix of each step as well as the computation time.
    tic
    [U_huber,V_huber, numb_it_huber] = Helmke_Huber_alg(U_init,V_init, accuracy, m_1_norm, m_2_norm);
    time_huber = toc;
    
    Estimate(i-3).Huber_U = U_huber;
    Estimate(i-3).Huber_V = V_huber;
    Estimate(i-3).Huber_time = time_huber;
    Estimate(i-3).Huber_iter = numb_it_huber;
    
    tic
    [U_helmke, V_helmke, numb_it_helmke] = helmke_alg(U_init, V_init, accuracy, m_1_norm, m_2_norm);
    time_helmke = toc;
    
    Estimate(i-3).Helmke_U = U_helmke;
    Estimate(i-3).Helmke_V = V_helmke;
    Estimate(i-3).Helmke_time = time_helmke;
    Estimate(i-3).Helmke_iter = numb_it_helmke;
    
    if i>4
       tic
       [U_smooth, V_smoot, numb_it_smooth] = helmke_alg_smoothed(U_init, V_init, ...
           Estimate(i-4).Smooth_U, Estimate(i-4).Smooth_V, accuracy, m_1_norm, m_2_norm);
       time_smooth = toc;
       
       Estimate(i-3).Smooth_U = U_smooth;
       Estimate(i-3).Smooth_V = V_smooth;
       Estimate(i-3).Smooth_time = time_smooth;
       Estimate(i-3).Smooth_iter = numb_it_smooth;
    else
       U_smooth = U_helmke;
       V_smooth = V_helmke;
       Estimate(i-3).Smooth_U = U_helmke;
       Estimate(i-3).Smooth_V = V_helmke;
       Estimate(i-3).Smooth_time = time_helmke;
       Estimate(i-3).Smooth_iter = numb_it_helmke; 
    end
    
    % Compare the estimated essential manifold to the groundtruth.
    [dist_R_huber, dist_t_huber] = compare_to_groundtruth(U_huber, V_huber, P(i-2).P);
    Estimate(i-3).Huber_dist_R = dist_R_huber;
    Estimate(i-3).Huber_dist_t = dist_t_huber;
    
    [dist_R_helmke, dist_t_helmke] = compare_to_groundtruth(U_helmke, V_helmke, P(i-2).P);
    Estimate(i-3).Helmke_dist_R = dist_R_helmke;
    Estimate(i-3).Helmke_dist_t = dist_t_helmke;
    
    [dist_R_smooth, dist_t_smooth] = compare_to_groundtruth(U_smooth, V_smooth, P(i-2).P);
    Estimate(i-3).Smooth_dist_R = dist_R_smooth;
    Estimate(i-3).Smooth_dist_t = dist_t_smooth;
    
    
    % Go to next frame.
    frame1 = frame2;
    point1 = point2;
    descr1 = descr2;
end
% Save the results.
save('results/Final_res_00', 'Estimate');