Evaluations of the sequences 00, 01 and 02.
For each sequence the three cost functions were evaluated. Once using RANSAC and once without.

RANSAC:
Input always 200 points of image 1 and image 2.
Takes random 10 points and computes the essential matrix, then computes m1'*E*m2 for all points to decides if a point is feasible to the estimation or not. Therefore accuracy chosen as 0.0001.
A set of feasible points is added to the consensus set if 100 points are in the set.

Huber parameter:
delta = 1.

