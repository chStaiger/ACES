                            ACES, July 2013
===================================================================================
View also DOI: 10.3389/fgene.2013.00289

The software consists of 4 main parts that provide functionality to run benchmarking
experiments on feature extractors and classifiers:

    featureExtractors
    classifiers
    statistics
    datatypes

See the "README.txt" in these four subdirectories for more information.

To start the crossvalidation use the file:
    StartCrossValidations.py

The "experiments" directory contains scripts to analyse the performance and
stability of classifiers. See the README files contained in the 
'experiments' directory for more information.

The "PinnacleZ" directory contains a dump of the 'PinnacleZ' repository that we
used for implementing the 'Chuang' method. This directory also contains a few
patches that we made; these are described in a README file.

Grid computing:
The files

    pipeline.py
    SetUpGrid.py
    create_tokens.py

contain scripts and classes to run the experiments on the Life Sciences grid. They employ
a couchDB to handle tokens that specify the needed runs.
Usually one parameter combination (dataset, network, algorithm, fold, shuffle_of_database)
is one instance that will be executed on one node.
The code can be adopted to other HPC facilities.

Data:
The expression and secondary data lie in the folder experiments/data.

Dependencies:
===============

python dependencies:

* numpy version 1.6.1.
* scipy version 0.10.0
* sklearn version 0.11
* h5py version 2.0.0
* xlrd
* couchDB version 0.9
* rpy2 version 2.3.6 (we only need functions that are present in both R 2.1*.* and R 3.0.1)
* networkx 1.6
* picas: download from github and save folder in the CrossValidation folder
* libRmath.so.1: employed by the Dao feature selection method, provided in the Crossvalidation folder
                    execute:
                 export LD_LIBRARY_PATH=~/CrossValidationTool:$LD_LIBRARY_PATH
