# Learning Networks from Random Walk-Based Node Similarities
Algorithms for learning network structure from effective resistances and other random-walk-based similarities.

# Matlab Code

The code is in Matlab. A number of functions depend on files in the `/utils` folder. Ensure that this folder is added to your Matlab path.

## Full graph recovery from pairwise node similarities

**exactRecover.m**: Given a full set of (n choose 2) effective resistances, recovers the unique graph with these resistances. May also be used with regularization as a heuristic method to match a noisy or incomplete set of effective resistances. See Section 4.2 of [the paper](https://thePaper). 

**exactPageRankRecover.m**: Given an n x n matrix of all pairwise personalized PageRank scores, recovers the unique graph  matching these scores. As with `exactRecover.m`, may be used heuristically with a regularization parameter. 

**exactRecoveryDemo.m**: Demonstrates how to use `exactRecover.m` and `exactPageRankRecover.m` to recover a graph from a full set of pairwise node similarities.

## Heuristic recovery from incomplete pairwise effective  resistances
**recoverMissing.m**: Given an incomplete set of effective resistances, fills in the missing resistances via the shortest path heuristic described in Section 4.2 of [the paper](https://thePaper). Then attempts to recover a set of edge weights using  `exactRecover.m`, possibly with regularization. Note that some of these edge weights may be negative. One option to clean them up is via `utils/noisyRecoveryCleanup.m`.

## Graph recovery via convex relaxation

**sdpRecover.m**: Given a set of (n choose 2) effective resistance constraints, finds the graph with minimal total degree with all effective resistances below these constraints, by solving the SDP described in Section 4.3 of [the paper](http://thePaper). Any resistance constraint input as 0 is considered to be a non-constraint. Requires [CVX](http://cvxr.com/cvx/) convex programming system to be installed.

# Citation

> @article{muscostsourakakis2018similarities, 
> title={Learning Networks from Random Walk-Based Node Similarities},
> author={Hoskins, Jeremy and Musco, Cameron, and Musco, Christopher, and Tsourakakis, Charalampos},
> year={2018},
