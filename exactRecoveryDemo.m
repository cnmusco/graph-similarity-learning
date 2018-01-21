% Example of how to exactly  recover a graph from its effective resistances
% or personalized PageRank scores

% Make sure that /utils is in the path
addpath(genpath('./'));

% Generate a random graph adjacency matrix
n = 50;
L = generateRandomLaplacian(n);
w = L2w(L);

% Compute exact resistances and personalized PageRank scores
effRes = getRes(w);
P = getPageRank(w,.01);

% Recover edge weights from these measures
wRes = exactRecover(effRes,0);
wPR = exactPageRankRecover(P,.01,0);

% Check that recovery is correct.
norm(wRes - w)
norm(wPR - w)