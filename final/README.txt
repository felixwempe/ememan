The algorithms:

calc_M.m: Is used for calculating gradient and Hessian of the cost function.

compare_to_groundtruth.m: Compares the result of the estimation to the groundtruh. Since no scale is estimated, the scale of the groundtruth vector is used to scale the estimated translation.

data_read.m: Reads data from a folder. Takes images computes the SIFT descriptors and matches them. Returns maximal 200 points.

evaluate_alg.m: Evaluates the algorithm on all frames given.

gt_read.m: Reads groundtruth data(in the form given in KITTI data set).

Helmke.m: Estimates the Essential Matrix for given points.

main.m: Script to read matched points and groundtruth that calls "evaluate_alg.m".

plot_rigid_motion.m: Plots the paths, estimated and groundtruth, the car drives.

ransacs.m: Random Sample Consensus algorithm for the given problem.
