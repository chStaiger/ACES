Some Notes on the "Dao Method" ~ Christine Staiger ~ Jan 2013
==============================================================

This is an exact method finding all heavy subnetworks that show deregulation in a subset of the poor outcome patients.
In this respect it differs from the other feature selection methods that assume genes or metagenes have to be deregulated
in all poor outcome patients. To score networks for the classification step, the average expression across all 
genes are calculated for each patient.

The input is an expression data set and a PPI network.

The method uses a lot of memory during execution, it is advisable to just let it run on a bigger cluster.

Prior to execute the feature extractor run
./make clean 
./make all


* The method is fully deterministic. Identical runs of the Dao feature extraction method
  will yield identical subnetworks.

* Genes that occurr in the network data but not in the expression data are removed before 
  the algorithm starts the search. This might disrupt the network.

REFERENCES
==========

[1] Inferring cancer subnetwork markers using density-constrained biclustering
    Phuong Dao*, Recep Colak*, Raheleh Salari, Flavia Moser, Elai Davicioni, Alexander Schonhuth**, Martin Ester** 
    9th European Conference on Computational Biology (ECCB 2010)
    Also: Bioinformatics 26(13): 1608-1615 (2010)

