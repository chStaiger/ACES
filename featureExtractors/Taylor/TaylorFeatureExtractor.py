# @Author 
# Sidney Cadot, update Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import numpy, json
from statistics.Statistics import PearsonCorrellationCoefficient
import numpy.random

class TaylorFeatureExtractorFactory(object):

    productName = "TaylorFeatureExtractor"

    performSignificanceTest = True

    @staticmethod
    def selectHubs(network, hubFraction):

        assert 0.0 <= hubFraction <= 1.0

       # make a dictionary with all node-node connections

        interactions = {}

        for (v1, v2) in network.edges:
            if v1 not in interactions:
                interactions[v1] = set()
            interactions[v1].add(v2)

            if v2 not in interactions:
                interactions[v2] = set()
            interactions[v2].add(v1)

        # determine the 15% percentile

        number_of_interactors = map(len, interactions)

        MinimalNumberOfNeigborsToBeClassifiedAsHub = 0
        bestMinimalNumberOfNeigborsToBeClassifiedAsHub = None
        bestScore = None

        while True:

            hubs = [(hub, interactors) for (hub, interactors) in interactions.iteritems() if len(interactors) >= MinimalNumberOfNeigborsToBeClassifiedAsHub]
            nHubs = len(hubs)

            if nHubs == 0:
                print "NOTE: Threshold %d leads to zero hubs, ending search." % MinimalNumberOfNeigborsToBeClassifiedAsHub
                break

            fraction = float(nHubs) / len(number_of_interactors)
            score = -abs(numpy.log(fraction / hubFraction))

            if False: # make true for verbose output
                print "NOTE: MinimalNumberOfNeigborsToBeClassifiedAsHub == %d would lead to %d hubs out of %d (%.6f %%; score = %.6f)." % (
                    MinimalNumberOfNeigborsToBeClassifiedAsHub, nHubs, len(number_of_interactors), 100.0 * fraction, score)

            if bestScore is None or score > bestScore:
                bestMinimalNumberOfNeigborsToBeClassifiedAsHub = MinimalNumberOfNeigborsToBeClassifiedAsHub
                bestScore = score

            MinimalNumberOfNeigborsToBeClassifiedAsHub += 1

        MinimalNumberOfNeigborsToBeClassifiedAsHub = bestMinimalNumberOfNeigborsToBeClassifiedAsHub
        assert MinimalNumberOfNeigborsToBeClassifiedAsHub is not None

        print "NOTE: selected hub threshold value: %d" % MinimalNumberOfNeigborsToBeClassifiedAsHub

        # Pick highly connected nodes as "hubs". Cut-off is at 'MinimalNumberOfNeigborsToBeClassifiedAsHub'.

        hubs = [(hub, interactors) for (hub, interactors) in interactions.iteritems() if len(interactors) >= MinimalNumberOfNeigborsToBeClassifiedAsHub]

        print "NOTE: selected %d hubs (genes with at least %d neighbors) out of %d genes (%.2f %%; target is %.2f %%)" % \
            (len(hubs), MinimalNumberOfNeigborsToBeClassifiedAsHub, len(interactions), float(len(hubs)) / float(len(interactions)) * 100.0, 100.0 * hubFraction)

        return hubs

    @staticmethod
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

    def train(self, dataset, network, average = True):
        """
        average: if True employ AvgHubDiff as feature values, otherwise each edge will be a feature with the difference as value.
        """

        hubs = self.selectHubs(network, hubFraction = 0.15)

        geneLabelToIndex = dict(zip(dataset.geneLabels, xrange(len(dataset.geneLabels))))

        # CS: To save running time it would be wiser to first intersect the network with the expression data and then determine the hubs
        #     edges = [(v1, v2) for (v1, v2) in network.edges if v1 in geneLabelToIndex and v2 in geneLabelToIndex]
        #     network.edges = edges ?
        #
        # SC: Running time is not a big issue, but that could work as well. The main thing we need to decide is the following TODO,
        #     then we can implement it easily.
        #
        # TODO: decide and document PRECISELY how missing data is handled -- i.e., network nodes for which we have no expression data.
        #       This is TODO #4.

        # Discard hubs for which we do not have expression data
        hubs = [(hub, interactors) for (hub, interactors) in hubs if hub in geneLabelToIndex]
        print "NOTE: after pruning unknown hubs, %d hubs remain." % len(hubs)

        # Discard interactors for which we do not have expression data.
        hubs = [(hub, frozenset([interactor for interactor in interactors if interactor in geneLabelToIndex])) for (hub, interactors) in hubs]

        # Discard hubs for which the interactor set is empty.
        hubs = [(hub, interactors) for (hub, interactors) in hubs if len(interactors) > 0]
        print "NOTE: after pruning of hubs with zero interactors after considering available expression data, %d hubs remain." % len(hubs)

        # Translate hubs and interactors to be indices rather than names
        hubs = [(geneLabelToIndex[hub], frozenset([geneLabelToIndex[interactor] for interactor in interactors])) for (hub, interactors) in hubs]

        if self.performSignificanceTest:
            print "NOTE: determining statistical cut-off by randomization (this will take a mighty long time) ..."

            numberOfRandomizations = 1000
            scores = []
            for i in xrange(1, 1 + numberOfRandomizations):
                numpy.random.seed(i)
                if (i <= 10) or (i <= 100 and i % 5 == 0) or (i % 50 == 0):
                    print "NOTE: outcome randomization", i, "of", numberOfRandomizations, "..."

                patientClassLabelsShuffled = numpy.random.permutation(dataset.patientClassLabels)

                cTrue  = dataset.expressionData[                  patientClassLabelsShuffled ] # corresponds roughly to "D" group ("Deceased") in Taylor et al.
                cFalse = dataset.expressionData[numpy.logical_not(patientClassLabelsShuffled)] # corresponds roughly to "A" group ("Alive")

                score = []
                for (hub, interactors) in hubs:

                    AverageHubDiff = self.calculateAverageHubDiff(hub, interactors, cFalse, cTrue)
                    score.append(AverageHubDiff)

                scores.append(score)

            # scores is now a numberOfRandomizations x numberOfHubs array
            scores = numpy.array(scores)

            scores.sort(axis = 0) # sort in-place along the numberOfRandomizations axis

            cutOffIndex = int(round(0.95 * numberOfRandomizations))
            if cutOffIndex >= numberOfRandomizations:
                cutOffIndex = numberOfRandomizations - 1
            assert 0 <= cutOffIndex < numberOfRandomizations

            # The threshold; just 5% of the randomized AverageHubDifference scores are larger than this.
            cutOffPerHub = scores[cutOffIndex]

        else:

            # This cutoff vector will lead to acceptance of all hubs
            cutOffPerHub = [-1.0 for i in xrange(len(hubs))]


        print "NOTE: scoring all hubs without randomization ..."

        cTrue  = dataset.expressionData[                  dataset.patientClassLabels ] # corresponds roughly to "D" group ("Deceased") in Taylor et al.
        cFalse = dataset.expressionData[numpy.logical_not(dataset.patientClassLabels)] # corresponds roughly to "A" group ("Alive")

        scored_hubs = []
        for ((hub, interactors), cutOffValue) in zip(hubs, cutOffPerHub):

            AverageHubDiff = self.calculateAverageHubDiff(hub, interactors, cFalse, cTrue)

            HubIsSignificant = AverageHubDiff >= cutOffValue

            if HubIsSignificant:

                scored_hub = (AverageHubDiff, (hub, interactors))
                scored_hubs.append(scored_hub)

        # sort hubs from highest-scoring to lowest-scoring
        hubs = [(hub, interactors) for (AverageHubDiff, (hub, interactors)) in reversed(sorted(scored_hubs))]

        interactorCountPerHub = [len(interactors) for (hub, interactors) in hubs]

        print "NOTE: after pruning of insignificant hubs, %d hubs remain with a total of %d interactors (mean interactors/hub: %f))." % (len(hubs), sum(interactorCountPerHub), float(sum(interactorCountPerHub))/float(len(hubs)))
        print "NOTE: interactor count per hub:", interactorCountPerHub
        print "NOTE: done training Taylor feature-extractor."

        if average:
            print "NOTE: Employing average difference between hub and its interactors as feature values."
            return TaylorFeatureExtractorAvgHubDiff(dataset.geneLabels, hubs)

        print "NOTE: Employing each (hub, interactor)-edge as feature."
        return TaylorFeatureExtractor(dataset.geneLabels, hubs)


class TaylorFeatureExtractor(object):

    name       = "TaylorFeatureExtractor"

    def __init__(self, geneLabels, hubs):
        self.geneLabels = geneLabels
        self.hubs       = hubs

        # Features are added in groups. A group corresponds to all interactors of a hub.
        featureCounts = []
        nf = 0
        for (hub, interactors) in self.hubs:
            nf += len(interactors)
            featureCounts.append(nf) # number corresponds to sum_{h \in hubs} degree(h) ?

        self.validFeatureCounts = featureCounts

    def extract(self, dataset, k): # is k the number of edges or hubs?

        # Make sure we have the same idea of our source data
        assert all(dataset.geneLabels == self.geneLabels)

        # If k is the number of hubs --> assert k <= len(self.hubs) ?
        assert k in self.validFeatureCounts

        result = numpy.zeros(shape = (dataset.numPatients, k))

        col = 0

        for (hub, interactors) in self.hubs:

            for interactor in interactors:

                result[:, col] = dataset.expressionData[:, interactor] - dataset.expressionData[:, hub]

                col += 1

            assert col <= k # not if k is the number of hubs
            if col == k:
                return result

    def toJsonExpression(self):
        return json.dumps((self.__class__.__name__, [geneLabel for geneLabel in self.geneLabels], [(hub, sorted(interactors)) for (hub, interactors) in self.hubs]))


# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

class TaylorFeatureExtractorAvgHubDiff(object):

    name       = "TaylorFeatureExtractorAvgHubDiff"

    def __init__(self, geneLabels, hubs):
        self.geneLabels = geneLabels
        self.hubs       = hubs

        self.validFeatureCounts = range(1, len(hubs) + 1)

    def extract(self, dataset, k): # is k the number of edges or hubs?

        # Make sure we have the same idea of our source data
        assert all(dataset.geneLabels == self.geneLabels)

        # If k is the number of hubs --> assert k <= len(self.hubs) ?
        assert k in self.validFeatureCounts

        result = numpy.zeros(shape = (dataset.numPatients, k))

        col = 0

        for (hub, interactors) in self.hubs:
            #Uses the average difference between the hub and all of its interactors as feature values
            result[:, col] = (((dataset.expressionData[:, hub].T)*numpy.ones((len(interactors), dataset.expressionData.shape[0]))).T - dataset.expressionData[:, list(interactors)]).mean(axis = 1)
            col = col + 1 

            assert col <= k # not if k is the number of hubs
            if col == k:
                return result

    def toJsonExpression(self):
        return json.dumps((self.__class__.__name__, [geneLabel for geneLabel in self.geneLabels], [(hub, sorted(interactors)) for (hub, interactors) in self.hubs]))

