# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import numpy
import json
from datatypes.EdgeSet import removeNodes
from statistics.Statistics import tStatisticForUnequalSampleSizeAndUnequalVariance

# Gene rabking algorithm based on Google's page rank algorithm.
# Adaptation from the GeneRank algorithm. In the original paper
# only one value per experiment were available --> difference between 
# two experiments. Here we calculate the mean expression for each gene
# in the poor and in the good outcome group.
# The expression change of each gene is calculated as 
# difference betqween these two means.

# The algorithm works with an undirected network.

class GeneRankFeatureExtractorFactory(object):
    productName = "GeneRankFeatureExtractor"

    def train(self, dataset, networkAndDampingFactor):
        if len(networkAndDampingFactor) == 1: # only network is given
            print "Use damping factor 0.3"
            dampingFactor = 0.3
            network = networkAndDampingFactor[0]
        else:
            dampingFactor = networkAndDampingFactor[1]
            network = networkAndDampingFactor[0]

        name = self.productName+'_'+str(dampingFactor)

        # NOTE: Intersect network and dataset
        # In GeneRank by Morrison et al. singleton genes get degree 1 
        #   --> keep genes  in data that are not present in the network as singletons
        # Un measured genes should be removed since they cannot be initiated with
        # expression, p-values or correlation. One can not classify with them later on.

        intersect = set(dataset.geneLabels).intersection(network.getNodes())
        remove = set(network.getNodes()).difference(intersect)
        cropNetwork = removeNodes(network, remove)
        # If you want no additional singletons, get rid of genes not present in the network: 
        #datasetI = dataset.extractGenesByIndices(dataset.name, numpy.argwhere(numpy.in1d(dataset.geneLabels, 
	    #cropNetwork.getNodes())).flatten(), False, False)

        # the order of indices corresponds to datasetI.geneLabels
        #ADJ = numpy.array(networkx.to_numpy_matrix(G, nodelist = datasetI.geneLabels))
        GeneToIndex = dict(zip(dataset.geneLabels, range(len(dataset.geneLabels))))
        IndexToGene = dict(zip(range(len(dataset.geneLabels)), dataset.geneLabels))
        ADJ = numpy.zeros((len(dataset.geneLabels), len(dataset.geneLabels)))
        for edge in cropNetwork.edges:
            edge = list(edge)
            ADJ[GeneToIndex[edge[0]], GeneToIndex[edge[1]]] = 1
            ADJ[GeneToIndex[edge[1]], GeneToIndex[edge[0]]] = 1      
        
        assert sum(sum(ADJ)) == len(cropNetwork.edges)*2

        mTrue  = dataset.expressionData[                  dataset.patientClassLabels ].mean(axis = 0)
        mFalse = dataset.expressionData[numpy.logical_not(dataset.patientClassLabels)].mean(axis = 0)

        # Note: order is important!
        meanFoldChange  = mTrue - mFalse
 
        # Analytical solution: W- adjacency; c - vector of abs(PC), D - diag(degree) 
        # (I - dW.T*inv(D))r = (1-d)c
        # r = inv(I - dW.T*inv(D))*(1-d)c
        I = numpy.identity(len(ADJ))
        d = dampingFactor
        W = ADJ
        #If genes have no neighbour, give them rank one (Morrison et al.)
        #Problem: Singletons and genepairs get same influence 
        #Solution: add 1 to the whole vector 
        D = numpy.identity(len(ADJ))*(sum(ADJ)+1)

        assert all(D.diagonal() > 0)
        #Check if (I - numpy.dot(d*W.T, numpy.linalg.inv(D)) is under determined --> pinv
        #      otherwise --> inv
        #if numpy.linalg.det(I - numpy.dot(d*W.T, numpy.linalg.inv(D))) == 0:
        #    r = numpy.dot(numpy.linalg.pinv(I - numpy.dot(d*W.T, numpy.linalg.inv(D))), (1-d)*numpy.abs(meanFoldChange))
        #else:
        #    r = numpy.dot(numpy.linalg.inv(I - numpy.dot(d*W.T, numpy.linalg.inv(D))), (1-d)*numpy.abs(meanFoldChange))

        r = numpy.dot(numpy.linalg.pinv(I - numpy.dot(d*W.T, numpy.linalg.inv(D))), (1-d)*numpy.abs(meanFoldChange))
        #sort genes according to rank
        featureGeneIndices = numpy.argsort(-numpy.abs(r))
        featureGeneIndices = map(int, featureGeneIndices) # make sure that the elements are regular ints, not NumPy types.

        return GeneRankFeatureExtractor(dataset.geneLabels, featureGeneIndices, name)

class GeneRankFeatureExtractor(object):

    def __init__(self, geneLabels, featureGeneIndices, name):
        self.geneLabels         = geneLabels
        self.featureGeneIndices = featureGeneIndices
        self.validFeatureCounts = range(1, len(self.featureGeneIndices) + 1)
        self.name               = name

    # k - number of features
    def extract(self, dataset, k):
        # Make sure we have the same idea of our source data
        assert all(dataset.geneLabels == self.geneLabels)

        assert k in self.validFeatureCounts

        # return expression values for k-most-significant features
        return dataset.expressionData[:, self.featureGeneIndices[:k]]

    def toJsonExpression(self):
        return json.dumps((self.__class__.__name__, [geneLabel for geneLabel in self.geneLabels], [featureGeneIndex for featureGeneIndex in self.featureGeneIndices]))

 
