
Some Notes on the "Taylor Method" ~ Sidney Cadot ~ September 2011, ~ Christine Staiger ~, May 2013
===================================================================================================

The "Taylor" method hinges on the idea that some genes functions as regulatory
hubs. Hubs are postulated to regulate pathway activations, and the job of the
Taylor method is to identify those hubs.

First, potential hubs are identified by considering their network connectivity
in the protein-protein interaction network. The top 15% network nodes (genes) in
terms of number-of-neighbors are selected as "potential hubs".

For each potential hub, an "AverageHubDiff" value is calculated. This value's
calculation is hard to describe in anything other than code:

=================

def calculateAverageHubDiff(hub, interactors, cFalse, cTrue):

    """
    Given a hub index, a set of interactors, a "false" expression matrix, and a "true" expression matrix,
      calculate the AverageHubDiff measure.

    cFalse  corresponds roughly to "A" group ("Alive")
    cTrue   corresponds roughly to "D" group ("Deceased") in Taylor et al.
    """

    hubA = cFalse[:, hub]
    hubD = cTrue [:, hub]

    delta_r = []

    for interactor in interactors:

        interactorA = cFalse[:, interactor]
        interactorD = cTrue [:, interactor]

        pccA = PearsonCorrellationCoefficient(hubA, interactorA)
        pccD = PearsonCorrellationCoefficient(hubD, interactorD)

        delta_r.append(pccA - pccD)

    delta_r = numpy.array(delta_r)

    AverageHubDiff = numpy.mean(numpy.abs(delta_r))

    return AverageHubDiff

=================

The value thus obtained is a measure of how differently the hub correlates with its interactors
when considering the two outcome groups.

The AverageHubDiff value defined above is re-calculated 1000 times with the training data where
the outcome labels have been randomly permuted; only if the AverageHubDiff value for the
non-permuted outcome vector is in the top 5% of the randomized AverageHubDiff values, the potential
hub is accepted. Rather unsurprisingly, this step discards about 95% of the original hubs.

ISSUES
======

There is one issue with the Taylor method that is worth mentioning.

* The score as defined above only allows for ranking the hubs. Since for
  classification the hubs' edges are used as features, rather than the hubs themselves,
  there is no direct ranking of the features.

  In our implementation we rank the hubs and add all edges belonging to one hub. The number trained
  in CV is thus not the best number of features but the best number of fetaure sets.

* A second version of feature values is to determine the average difference between the hub and its interactors.
  In that case features are hubs and hubs are ranked in the previous step.

  average = True in the train function triggers the second version of features.


REFERENCES
==========

[1] Dynamic modularity in protein interaction networks predicts breast cancer outcome
    Ian W Taylor, Rune Linding, David Warde-Farley, Yongmei Liu, Catia Pesquita, Daniel Faria,
        Shelly Byll, Tony Pawson, Quaid Morris, Jeffrey L Wrana
    Nature Biotechnology ~ Volume 27 Number 2 ~ February 2009
