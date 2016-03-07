
# TODO #1 : Decide and document the handling of network nodes that have no corresponding expression data.

import numpy, json

class ChuangFeatureExtractorFactory(object):

    # CS: First intersect data and network to make sure we do not invoke the strange subnetwork score calculation in PinnacleZ ?
    #     (This is to be decided as per TODO #1)
    #
    # CS: If we first call PinnacleZ and then remove unknown genes, we will order the subnetworks according to a strange score.
    #     (This is to be decided as per TODO #1)

    productName = "ChuangFeatureExtractor"

    def __init__(self, pzJarFilename, pzNormalizeInputMatrix):

        self.pzJarFilename = pzJarFilename
        self.pzNormalizeInputMatrix = pzNormalizeInputMatrix

    def train(self, dataset, network):

        if self.pzNormalizeInputMatrix and not dataset.checkIfExpressionDataProperlyNormalized():
            print "WARNING: we are about to run PinnacleZ on non-normalized data. PinnacleZ will start by normalizing the expression data !!!"

        from PinnacleZ import runPinnacleZ

        geneModules = runPinnacleZ(
            self.pzJarFilename,
            self.pzNormalizeInputMatrix,
            dataset.patientClassLabels,
            dataset.geneLabels,
            dataset.patientLabels,
            dataset.expressionData,
            network.edges,
            deleteTemporaryDirectory = False
        )

        print "NOTE: number of networks as received from PinnacleZ:", len(geneModules)

        features = {}
        for (score, pValues, genes) in geneModules:
            genes = frozenset(genes)
            if genes in features:
                assert features[genes] == score
            else:
                features[genes] = score

        if len(features)  < len(geneModules):
            print "WARNING: removed duplicate modules (%d modules -> %d features)" % (len(geneModules), len(features))

        # Rank the gene-sets by score, and drop the score:

        features = [genes for (score, genes) in reversed(sorted([(score, genes) for (genes, score) in features.iteritems()]))]

        # Up until now, the genes are given as identifiers.
        # We now replace them by their indices.

        # NOTE: we DISCARD those genes within modules that are not represented in the expression data.
        #       This is consistent with how the ST3 test of PinnacleZ works, but it is up for
        #       debate if this is the proper way to handle this.

        geneLabelToIndex = dict(zip(dataset.geneLabels, xrange(len(dataset.geneLabels))))

        features = [frozenset([geneLabelToIndex[gene] for gene in genes if gene in geneLabelToIndex]) for genes in features]

        # The previous step may have yielded features that are identical. The next step ensures that duplicates
        # are discarded.

        filtered_features = []
        for feature in features:
            if feature in filtered_features:
                pass
            else:
                filtered_features.append(feature)

        if len(filtered_features) < len(features):
            print "WARNING: removed duplicate features [dupes after removing unknown genes] (%d features -> %d features)" % (len(features), len(filtered_features))

        return ChuangFeatureExtractor(dataset.geneLabels, filtered_features)


class ChuangFeatureExtractor(object):

    name               = "ChuangFeatureExtractor"
    def __init__(self, geneLabels, features):

        # The "features" as received here are an ordered list of frozensets.
        # Each frozenset represents a module that can be scored.
        self.geneLabels         = geneLabels
        self.features           = features
        self.validFeatureCounts = range(1, len(self.features) + 1)

    @staticmethod
    def score(expressionData, feature):
        # The division by sqrt(n) rather than simply by n follows the paper.
        # There does not seem to be a good reason for this.
        return numpy.sum(expressionData[:, list(feature)], axis = 1) / numpy.sqrt(len(feature))

    def extract(self, dataset, k):

        assert all(dataset.geneLabels == self.geneLabels)
        assert k in self.validFeatureCounts

        # Return the network scores for the k best subnetworks
        return numpy.transpose(numpy.array([self.score(dataset.expressionData, feature) for feature in self.features[:k]]))

    def toJsonExpression(self):
        return json.dumps((self.__class__.__name__, [geneLabel for geneLabel in self.geneLabels], [sorted(feature) for feature in self.features]))
