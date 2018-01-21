# Learning Networks from Random Walk-Based Node Similarities
This is a collection of algorithms for learning network structure from effective resistances and other random-walk-based similarities, as described in the paper [Learning Networks from Random Walk-Based Node Similarities](http://thePaper). 

Includes methods for exact graph recovery, heuristic methods, and optimization-based approaches (both convex and non-convex). See the paper for details on and comparision of these methods.

# Matlab Code

The code is in Matlab. A number of functions depend on files in the `/utils` folder. Ensure that this folder is added to your Matlab path.

## Full graph recovery from pairwise node similarities

**exactRecover.m**: Given a full set of `(n choose 2)` effective resistances, recovers the unique graph with these resistances. May also be used with regularization as a heuristic method to match a noisy or incomplete set of effective resistances. See Section 4.2 of [the paper](https://thePaper). 

**exactPageRankRecover.m**: Given an `n x n` matrix of all pairwise personalized PageRank scores, recovers the unique graph  matching these scores. As with `exactRecover.m`, this method may be used heuristically with a regularization parameter. 

**exactRecoveryDemo.m**: Demonstrates how to use `exactRecover.m` and `exactPageRankRecover.m` to recover a graph from a full set of pairwise node similarities.

## Heuristic recovery from incomplete pairwise effective  resistances
**recoverMissing.m**: Given an incomplete set of effective resistances, fills in the missing resistances via a shortest path heuristic described in Section 4.2 of [the paper](https://thePaper). The method then attempts to recover a set of edge weights using  `exactRecover.m`, possibly with regularization. Note that some of these edge weights may be negative. One option to clean them up is via `utils/noisyRecoveryCleanup.m`.

## Graph learning via least squares minimization

We provide two gradient/stochastic coordinate descent based methods which attempt to solve the non-convex problem of finding a graph whose effective resistances are as close as possible to the given resistances in `l2` norm. See Sections 3.3 and 4.2 of [the paper](http://thePaper) for details. 

Any effective resistance input as 0 is considered to be un-constrained. See comments in the code for details on tuning and optimization method options. Both methods allow for a regularization parameter `lambda`, and minimize the distance to the target resistances plus `lambda*tr(L)`.

**effResGDSmall.m**: Use this method for small graphs, where is is possible to fit an `(n choose 2) x n` edge-vertex incidence matrix in memory.

**effResGD.m**: Use this method for large graphs. It will avoid computing an `(n choose 2) x n` edge-vertex incidence matrix. For large graphs, it is generally advised to set `batchSize << (n choose 2)`.

Note that both the above methods make use of paralleism via `parfor` loops. You can set the number of paralell workers before running these methods, using a snippet like:
```
myCluster = parcluster('local');
myCluster.NumWorkers = 4;
parpool(4);
```

**leastSquaresDemo.m**: Demonstrates how to use `graphGDSmall.m` and `graphGD.m` on a small random k-nearest neighbors graph.


## Graph learning via convex relaxation

**sdpRecover.m**: Given a set of `(n choose 2)` effective resistance constraints, finds the graph with minimal total degree with all effective resistances below these constraints, by solving the SDP described in Section 4.3 of [the paper](http://thePaper). Any resistance input as 0 is considered to be un-constrained. Requires [CVX](http://cvxr.com/cvx/) convex programming system to be installed.

# Citation

> @article{muscostsourakakis2018similarities, 
> title={Learning Networks from Random Walk-Based Node Similarities},
> author={Hoskins, Jeremy and Musco, Cameron, and Musco, Christopher, and Tsourakakis, Charalampos},
> year={2018},
