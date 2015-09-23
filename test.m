% Test

frame1 = single(imread('data/000000.png'));
frame2 = single(imread('data/000001.png'));

% Get interest points.

[point1, descr1] = vl_sift(frame1);
[point2, descr2] = vl_sift(frame2);

% Match points.

[matches, score] = vl_ubcmatch(descr1, descr2);

% Select 200 best matches.

[~, idx] = sort(score);
good_matches = matches(:,idx(1:200))';

m_1 = point1(1:2,good_matches(:,1))';
m_2 = point2(1:2,good_matches(:,2))';

% Get initial estimate using the best 200 scores and get the inliers.
[F, Inl] = estimateFundamentalMatrix(m_1, m_2);

% Compute essential matrix out of the fundamental.

P0 = [7.188560000000e+02 0.000000000000e+00 6.071928000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 7.188560000000e+02 1.852157000000e+02 0.000000000000e+00; ...
    0.000000000000e+00 0.000000000000e+00 1.000000000000e+00 0.000000000000e+00];

Calib = P0(1:3, 1:3);
Calib_inv = inv(Calib);

InGuess = Calib' * F * Calib;

% Choose Inliers and Calibrate them.
m_1_inl = (Calib_inv * [m_1(Inl,:), ones(sum(Inl),1)]')';
m_2_inl = (Calib_inv * [m_2(Inl,:), ones(sum(Inl),1)]')';

% Use algorithm.

E = helmke_alg(InGuess, 0.1, m_1_inl, m_2_inl);
