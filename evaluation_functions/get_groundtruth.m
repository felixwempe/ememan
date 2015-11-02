% Read the groundtruth informations from the file.

directory = '/home/felix/Universität/Praktikum_Andreas_Neufeld/dataset/dataset/poses/';

tmp = textread([directory, '00.txt']);
P(length(tmp)) = struct;

for i=1:length(tmp)
    P(i).P = reshape(tmp(i,:), [4,3])';
    
end

clear tmp i directory