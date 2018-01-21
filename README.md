# Learning Networks from Random Walk-Based Node Similarities
Algorithms for learning network structure from effective resistances and other random-walk-based similarities.

# Matlab Code

The code is in Matlab. A number of functions depend on files in the `/utils` folder. Ensure that this folder is added to your Matlab path.

## Full graph recovery from pairwise node similarities

**exactRecover.m**: Given a full set of (n choose 2) effective resistances, recovers the unique graph with these resistances. May also be used with regularization as a heuristic method to match a noisy or incomplete set of effective resistances. See Section 4.2 of [the paper](https://thePaper). 

**exactPageRankRecover.m**: Given an n x n matrix of all pairwise personalized PageRank scores, recovers the unique graph  matching these scores. As with `exactRecover.m`, may be used heuristically with a regularization parameter. 

**exactRecoveryDemo.m**: Demonstrates how to use `exactRecover.m` and `exactPageRankRecover.m` to recover a graph from a full set of pairwise node similarities.

# Citation

> @article{muscostsourakakis2018similarities, 
> title={Learning Networks from Random Walk-Based Node Similarities},
> author={Hoskins, Jeremy and Musco, Cameron, and Musco, Christopher, and Tsourakakis, Charalampos},
> year={2018},
