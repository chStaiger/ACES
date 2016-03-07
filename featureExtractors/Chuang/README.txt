
Some Notes on the "Chuang Method" ~ Sidney Cadot ~ September 2011
=================================================================

The original code as used for the Chuang et al. paper [1] was irretrievably lost. We are
using a re-implementation of the method in Java called "PinnacleZ". This re-implementation
was developed in the Ideker group that authored the original article [1].

Our "PinnacleZ.py" file wraps a single PinnacleZ run as a Python function that calls the
Java-based PinnacleZ as a subprocess.

The "ChuangFeatureExtractor.py" file implements the Chuang feature-extraction method and
its score calculation. The "ChuangFeatureExtractor.py" file uses the runPinnacleZ()
function as defined in "PinnacleZ.py" for training.

ISSUES
======

There are several issues with the Chuang method that are worth mentioning.

* The PinnacleZ implementation is similar but not identical to the method described in [1].
  We know this, because if we run PinnacleZ on the same input (settings, expression data
  and protein-protein interaction network) as used in the original paper, we will get
  a larger number of significant modules than reported in [1].

  It is currently unknown what causes this.

* PinnacleZ is non-deterministic. It uses a random generator, the seed of which is not
  fixed; furthermore, it uses a single random generator from within multiple threads, which
  makes the code non-deterministic even if we initialize the seed to a known value.

  To make PinnacleZ deterministic, one would also have to force the number of threads to
  one. (We have in fact done this to confirm that this works, but when doing so the module
  search becomes prohibitively expensive.)

  Essentially, this means that PinnacleZ will return different outputs on identical runs;
  it is not a 'pure function'.

* As part of reading in the expression matrix, PinnacleZ performs a gene-wise normalization
  of expression values (to mu = 0, sigma = 1). However, this operation is undesirable in
  cases where we are looking at a subset of data that has previously been normalized.
  Precisely this happens during an n-fold cross-validation loop.

  It is unknown to us if, and how this issue was handled in the experiments described in [1].

  In order to skip the normalization step, we had to implement a small patch in the
  PinnacleZ source code. This patch adds a "-z" option that instructs PinnacleZ to *omit*
  its usual gene-wise z-normalization step.

* The paper [1] is silent on the issue of nodes within the interaction network for which
  no expression data is known. This should have been documented because it is an essential
  part of the feature extraction algorithm.

  As a minimum, it is clear that the [1] method does *something* with the proteins for
  which it does not have expression data: about 5% of the genes in the modules identified
  in [1] for the Wang and Vijver datasets contain nodes for which no expression data is
  known in the corresponding expression dataset. This is quite strange.

  This opens up the issue on how to handle features produced by PinnacleZ that contain
  genes for which we have no expression data. We handle this by filtering these 'unknown'
  genes from the modules; we therefore calculate the activation score as the sum of the
  activations of known genes, divided by the square root of known genes. This is consistent
  with the calculation performed in the third statistical test for significance ("ST3" in
  PinnacleZ).

* PinnacleZ includes the scores for the final modules in its output file. We can reproduce
  these numbers using the Python "scoreSubnetwork()" function as defined in "PinnacleZ.py";
  see that code for some comments.

  The "Mutual Information" scoring function used to score modules in PinnacleZ during network
  search appears to be buggy in cases where the (candidate) module contains nodes that have
  no corresponding expression data.

  Of course, this makes the modules generated suspect as well.

  We contacted the Ideker group about this issue but unfortunately we did not get an
  acknowledgement of the issue.

* The module activation is calculated as the sum of the normalized expressions of the
  genes, divided by the square root of the number of genes in the module. There appears
  to be no justification for dividing by sqrt(n). Also, it is noted that pathways where
  the expression of some constituent genes are anti-correlated will actually be scored
  as "less active".

REFERENCES
==========

[1] Network-based classification of breast cancer metastasis
    Han-Yu Chuang, Eunjung Lee, Yu-Tsueng Liu, Doheon Lee, Trey Ideker
    Molecular Systems Biology 3:140
