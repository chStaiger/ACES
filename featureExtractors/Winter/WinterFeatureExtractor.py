# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

from scipy import sparse
import numpy
import networkx
import json
from datatypes.EdgeSet import removeNodes
from statistics.Statistics import PearsonCorrellationCoefficient


# Gene rabking algorithm based on Google's page rank algorithm.
# Adaptation from Winter et al. Google goes cancer: Improving outcome prediction for cancer
# patients by network-based ranking of marker genes, to be published in Plos Comp Biol.

# The algorithm works with an undirected network.

class WinterFeatureExtractorFactory(object):
    productName = "WinterFeatureExtractor"

    def train(self, dataset, networkAndDampingFactor, survTime = None):
        if len(networkAndDampingFactor) == 1: # only network is given
            print "Use damping factor 0.3"
            dampingFactor = 0.3
            network = networkAndDampingFactor[0]
        else:
            dampingFactor = networkAndDampingFactor[1]
            network = networkAndDampingFactor[0]
        
        if survTime == None:
            print "Use classlabels to calculate correlations."
            survTime = dataset.patientClassLabels

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

        GeneCorrelation = []
        labels = numpy.cast[float](survTime)
        for gene in range(len(dataset.geneLabels)):
            expr = dataset.expressionData[:, gene]
            GeneCorrelation.append(PearsonCorrellationCoefficient(expr, labels))
 
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
        #NOTE: if you have time check if exact inverse is possible	
        #Check if (I - numpy.dot(d*W.T, numpy.linalg.inv(D)) is under determined --> pinv
        #      otherwise --> inv
        #if numpy.linalg.det(I - numpy.dot(d*W.T, numpy.linalg.inv(D))) == 0:
        #    r = numpy.dot(numpy.linalg.pinv(I - numpy.dot(d*W.T, numpy.linalg.inv(D))), (1-d)*numpy.abs(GeneCorrelation))
        #else:
        #r = numpy.dot(numpy.linalg.inv(I - numpy.dot(d*W.T, numpy.linalg.inv(D))), (1-d)*numpy.abs(GeneCorrelation))
        #NOTE: If not, use pseudo inverse
        r = numpy.dot(numpy.linalg.pinv(I - numpy.dot(d*W.T, numpy.linalg.inv(D))), (1-d)*numpy.abs(GeneCorrelation))

        #sort genes according to rank
        featureGeneIndices = numpy.argsort(-numpy.abs(r))
        featureGeneIndices = map(int, featureGeneIndices) # make sure that the elements are regular ints, not NumPy types.

        if survTime.dtype == 'bool':
            return WinterFeatureExtractor(dataset.geneLabels, featureGeneIndices, "WinterFeatureExtractor"+"_"+str(dampingFactor))
        else:
            return WinterFeatureExtractor(dataset.geneLabels, featureGeneIndices, "WinterFeatureExtractor_SurvTime"+"_"+str(dampingFactor))

class WinterFeatureExtractor(object):

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

 
