%% Skript to evaluate the speed, accuracy, etc of Helmke paper on a set of images.

% Set the values: Accuracy, maximal Iterations per step, Calibration
% matrix of the frames.
accuracy = 0.000001;


P0 = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00 0.000000000000e+00];

Calib = P0(1:3, 1:3);
Calib_inv = inv(Calib);

%% First we have to read in the whole set of images. 
% (Maybe only read in a subset and then on the flow more images, or use a 
% loop to read images only when they are used and forget about the old 
% ones)

% Get all filenames of a specified directory.
direct = '/home/felix/Universität/Praktikum_Andreas_Neufeld/dataset/sequences/00/image_0/';
filenames = dir(direct);

% Read in the first image. And compute the Sift features.
frame1 = single(imread([direct, filenames(3).name]));
[point1, descr1] = vl_sift(frame1);

Estimate(length(filenames)-2) = struct;
% Start a loop over all images.
for i=4:20 %length(filenames)
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
    
    % Get the initialisation for the essential matrix.
    % (Here it is a trivial initialisation, assuming the car always moves
    % forward) E = R*[t], for the rotation matrix R=Id and a t=[1 0 0])
    InGuess = eye(3)*[0 -1 0; 1 0 0; 0 0 0];
    
    % Use the calibration matrix to normalize the image points.
    m_1_norm = (Calib_inv * [m_1, ones(length(m_1),1)]')';
    m_2_norm = (Calib_inv * [m_2, ones(length(m_2),1)]')';

    % Use the Helmke algorithm to compute a good essential matrix. Save the
    % Essential matrix of each step as well as the computation time.
    tic
    if i <6
        [E, numb_it] = helmke_alg(InGuess, accuracy, m_1_norm, m_2_norm);
    else
        [E, numb_it] = helmke_alg_smoothed(InGuess, Estimate(i-4).essmat,...
            Estimate(i-5).essmat, accuracy, m_1_norm, m_2_norm);
    end
    time = toc;
    Estimate(i-3).essmat = E;
    Estimate(i-3).time = time;
    Estimate(i-3).iter = numb_it;
    
    % Go to next frame.
    frame1 = frame2;
    point1 = point2;
    descr1 = descr2;
end
% Save the results.
%save('results/estimations_Huber_seq_00.mat', 'Estimate');