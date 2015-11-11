function [ P ] = gt_read(gt_dir, sequence, numb_frames)
%The function reads (up to numb_frames) groundtruth rotation and
%translation informations into the struct P.
%   Input arguments:
%       gt_dir:     directory containing the groundtruth informations.
%       sequence:   The sequence that should be evaluated.
%       numb_frames:Number of frames that shall be evaluated.
%   Output arguments:
%       P:  struct containing at each entery the groundtruth rotation and
%           translation informations.
P(numb_frames) = struct;

tmp = textread([gt_dir, sequence, '.txt']);

for i=1:numb_frames
    P(i).P = reshape(tmp(i,:), [4,3])';
end

end

