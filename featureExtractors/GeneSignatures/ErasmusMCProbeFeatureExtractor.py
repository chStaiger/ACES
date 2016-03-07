# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import numpy, json
from statistics.Statistics import tStatisticForUnequalSampleSizeAndUnequalVariance
import xlrd

# Difference between the two classes:
#
#     SingleGeneFeatureExtractorFactory: receives a training dataset and determines the possible features and their ranking
#     SingleGeneFeatureExtractor: receives an order of features and determines their feature values given a dataset

class ErasmusMCProbeFeatureExtractorFactory(object):

    productName = "ErasmusMCProbeSignature"
    
    def __init__(self, signatureFile):
        self.signatureFile = signatureFile


    def train(self, dataset):

        filename = self.signatureFile
        f = open(filename)
        lines = f.readlines()
        probes = []
        for line in lines:
            if line[0].isdigit():
                probes.append(line.split(',')[0])

        featureGeneIndices = numpy.argwhere(numpy.in1d(dataset.geneLabels, probes)).flatten()
        assert len(featureGeneIndices) > 0

        return ErasmusMCProbeFeatureExtractor(dataset.geneLabels, featureGeneIndices)


class ErasmusMCProbeFeatureExtractor(object):
    name               = "ErasmusMCProbeSignature"
    def __init__(self, geneLabels, featureGeneIndices):
        self.geneLabels         = geneLabels
        self.featureGeneIndices = featureGeneIndices
        self.validFeatureCounts = len(featureGeneIndices)

    # k - number of features
    def extract(self, dataset):
        # Make sure we have the same idea of our source data
        assert all(dataset.geneLabels == self.geneLabels)

        # return expression values for k-most-significant features
        return dataset.expressionData[:, self.featureGeneIndices]

    def toJsonExpression(self):
        return json.dumps((self.__class__.__name__, [geneLabel for geneLabel in self.geneLabels], [int(featureGeneIndex) for featureGeneIndex in self.featureGeneIndices]))
