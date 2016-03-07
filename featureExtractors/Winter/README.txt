Some Notes on the "Winter/Gene Rank Method" ~ Christine Staiger ~ Jan 2013
==============================================================

The method is based on the Page Rank algorithm.
Genes are ranked according to an initial node score, determined on the data,
and their connactivity and surrounding nodes in the network.

The different variationas only differ in the initial node score.
Winter: Correlation between gene's expression and patients' classlabels
WinterTime [1]: Correlation between gene's expression and patients' DMFS or RFS time
GeneRank [2]: absolute mean difference between gene's expression across the two patient classes.
GeneRankTscore: t-statistic of the gene's expression across the two patient classes.


REFERENCES
==========
[1] Winter C, Kristiansen G, Kersting S, Roy J, Aust D, et al. (2012) 
    Google Goes Cancer: Improving Outcome Prediction for Cancer Patients by 
    Network-Based Ranking of Marker Genes. 
    PLoS Comput Biol 8(5): e1002511. doi:10.1371/journal.pcbi.1002511

[2] GeneRank: Using search engine technology for the analysis of microarray experiments
    Julie L Morrison, Rainer Breitling, Desmond J Higham and David R Gilbert
    BMC Bioinformatics 2005 doi:10.1186/1471-2105-6-233
