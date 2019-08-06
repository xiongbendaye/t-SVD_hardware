clear;
clc;

load('A_50_50_31.mat');
A = A(:,:,1:16);
tic
[U,S,V] = tsvd(A, 0);
toc
