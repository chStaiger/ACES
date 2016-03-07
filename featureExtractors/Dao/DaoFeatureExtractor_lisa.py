# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import numpy, json
from wDCB import runWCDB

class DaoFeatureExtractorFactory(object):

    productName = "DaoFeatureExtractor"

    # The feature extractor needs either an executable or a path to an executable which does the feature extraction.
    # Please make sure that your method does NOT z-normalise the input expression matrix. This is decpreciated in the cross validation.
    # Z-normalising will take place when initialising a new instance of a dataset.
    def __init__(self, wDCBFilename):

        self.wDCBFilename = wDCBFilename

    def train(self, dataset, network):
               
        net = network.copy()
        # If the network does not provide any edge weights, initialise them here with 999
        if net.edgeweights == None:
            net.edgeweights = dict(zip(list(net.edges), [999]*len(net.edges)))
        
        #Intersect network with expression data.
        #genes = dataset.geneLabels 
        #nodes = network.getNodes()
        #intersect = set(list(genes)).intersection(nodes)
        #remove = nodes.difference(intersect)

        #network.removeNodes(remove)
        #intersect = set(list(genes)).intersection(network.getNodes())

        # NOTE: The dataset is not sorted according to the patients class labels.
        # Sort dataset, first good outcome (reference group, label=False), then poor outcome (cases, label = True)
        # attach to each classlabel the index of the array.
        LabelToIndex = zip(dataset.patientClassLabels.tolist(), range(0, len(dataset.patientClassLabels)))
        # sort the classlabels (first 'good', then 'poor') and get new order of patients
        indices = [index for (label, index) in sorted(LabelToIndex)] 
        #dataset with the right order of the patients. 
        sortedDataset = dataset.extractPatientsByIndices(dataset.name, indices, False)
        #sortedShrunkenDataset = sortedDataset.extractGenesByIndices(sortedDataset.name, 
        #    numpy.argwhere(numpy.in1d(sortedDataset.geneLabels, intersect))[:, 0], False)        
        # TODO: Call wrapper to your method here
        # geneModules = [(score, genes)]; a list of tupels. "score" gives the discriminatory power (e.g. t-test), 
        #               "genes" gives the genes for each subnetwork

        geneModules = runWCDB(
            self.wDCBFilename,
            sortedDataset, 
            net, # contains edges and edgeweights
            MinimumCases = dataset.numPatientsBadOutcome*5/100,
            deleteTemporaryDirectory = True
        )

        print "NOTE: number of networks as received from wDCB:", len(geneModules)

        features = {}
        for (score, genes) in geneModules:
            genes = frozenset(genes)
            #NOTE: If the algorithm finds the same subnetwork several times, we remove the duplicates and keep it only once as a feature.
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
	print "Features found:", len(filtered_features)
        return DaoFeatureExtractor(dataset.geneLabels, filtered_features)


class DaoFeatureExtractor(object):

    name               = "DaoFeatureExtractor"
    def __init__(self, geneLabels, features):

        # The "features" as received here are an ordered list of frozensets.
        # Each frozenset represents a module that can be scored.
        self.geneLabels         = geneLabels
        self.features           = features
        self.validFeatureCounts = range(1, len(self.features) + 1)

    # Score a subnetwork for one patient
    @staticmethod
    def score(expressionData, feature):
        # sum all expression values of genes in the subnetwork and divide by the number of genes in the subnetwork.
        return numpy.sum(expressionData[:, list(feature)], axis = 1) /len(feature)

    # Calculate the subnetwork scores of the k best subnetworks for each patient
    def extract(self, dataset, k):

        assert all(dataset.geneLabels == self.geneLabels)
        assert k in self.validFeatureCounts

        # Return the network scores for the k best subnetworks
        return numpy.transpose(numpy.array([self.score(dataset.expressionData, feature) for feature in self.features[:k]]))

    def toJsonExpression(self):
        return json.dumps((self.__class__.__name__, [geneLabel for geneLabel in self.geneLabels], [sorted(feature) for feature in self.features]))
