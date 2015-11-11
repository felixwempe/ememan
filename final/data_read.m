function [ M1, M2 ] = data_read(data_dir, Calib_mat, numb_frames)
% Reads all images (up to numb_frames) from a directory and computes there
% best (max. 200) matches.
%   Input parameters are:
%   data_dir:       directory containing the image data.
%   Calib_mat:      Calibration matrix of the image data (here for all
%                   images the same)
%   numb_frames:    maximum of images read.
%   Output parameters are:
%   M1:     Struct containing the matched image points of the first image.
%           At each iteration.
%   M2:     Struct containing the matched image points of the second image.
%           At each iteration.
M1(numb_frames) = struct;
M2(numb_frames) = struct;

filenames = dir(data_dir);

im1 = single(imread([data_dir, filenames(3).name]));
[point1, descr1] = vl_sift(im1);
for i=1:numb_frames
    im2 = single(imread([data_dir, filenames(3+i).name]));
    [point2, descr2] = vl_sift(im2);
    [matches, score] = vl_ubcmatch(descr1, descr2);
    % Select 200 best matches.
    [~, idx] = sort(score);
    good_matches = matches(:,idx(1:200))';
    
    m_1 = point1(1:2,good_matches(:,1))';
    m_2 = point2(1:2,good_matches(:,2))';
    
    % Use the calibration matrix to normalize the image points.
    M1(i).m = Calib_mat\[m_1, ones(length(m_1),1)]';
    M2(i).m = Calib_mat\[m_2, ones(length(m_2),1)]';
    
    M1(i).m = M1(i).m';
    M2(i).m = M2(i).m';
    
    % Iterate.
    point1 = point2;
    descr1 = descr2;
end


end

