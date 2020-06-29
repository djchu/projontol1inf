%% demo for the Euclidean projection onto the L_{1,Inf}-norm ball
% ===================================================================
% To run this demo, you first need to compile the mex file with Matlab:
%   >> mex -output myssnewton mex-ssnewton.c
% ===================================================================

clc,clear,close all;

dim = 1e4;
task = 1e4;

coef = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5];
width = length(coef);

A = randn(dim,task);
A_l1inf = sum(max(abs(A), [], 2));

for iter = 1:width
    C = coef(iter)*A_l1inf;
    tic
    AT = A';
    B_SNewton = myssnewton(reshape(AT(:),dim,task), C);
    BT = reshape(B_SNewton(:), task, dim);
    B_SNewton = BT';
    runningtime = toc;
    B_l1inf = sum(max(abs(B_SNewton), [], 2));
    
    fprintf('factor: %.2f, time: %f\n', coef(iter), runningtime);
end