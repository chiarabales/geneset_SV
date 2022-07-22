# redundancy aware Shapley values

Demo code for the redundancy awareness integration into Shapley values for ranking sets within family of sets. 

## Introduction 
The well founded concept of Shapley values is frequently used to rank players in a fair manner with respect to their role into a cooperative game. However, the notion itself of Shapley values is missing awareness of redundancy among players. In particular, when two or more players are interchangeable they are assigned same Shapley values therefore using Shapley values to select players in the game can lead to a selection of redundant players.
The redundany-aware Shapley values based rankings proposed are a direct solution to this issue in the context of family of sets. A direct application of the proposed methods is gene set and pathways analysis. The rankings are obtained in two steps
* computation of the Shapley values for the sets 
* introduction of greedy punishments based on the Jaccard rate
We encourage you to read the [full paper](https://arxiv.org)

## Requirements
* numpy 1.18.5
* pandas 1.0.5

## Citation
If you found this work useful, please cite our paper:

```
@inproceedings{balestraUSV,
  author    = {Chiara Balestra  and
               Carlo Maj and
               Emmanuel M\"uller and
               Andreas Mayr},
  title     = {Redundancy-aware unsupervised ranking based on game theory - application to gene enrichment analysis},
  year      = {2022}
 }
```

## Use

The code consists of different files. Use the \_\_example.py file to test the code. 

The data can be downloaded from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp and should be saved in the folder '../data' using .csv extension with comma separation.
