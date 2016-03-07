# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import numpy, json
from statistics.Statistics import tStatisticForUnequalSampleSizeAndUnequalVariance

class VantVeerFeatureExtractorFactory(object):

    productName = "VantVeerGeneSignature"

    def train(self, dataset):

        sig = ['Entrez_2131', 'Entrez_55351', 'Entrez_10874', 'Entrez_51203', 'Entrez_10403', 'Entrez_56942', 'Entrez_57110', 'Entrez_51377',
            'Entrez_6515', 'Entrez_8840', 'Entrez_92140', 'Entrez_5019' 'Entrez_8833', 'Entrez_79888', 'Entrez_8817', 'Entrez_9833', 'Entrez_27113',
            'Entrez_5984', 'Entrez_51514', 'Entrez_7043', 'Entrez_23594', 'Entrez_4175', 'Entrez_51560', 'Entrez_4318', 'Entrez_1284', 'Entrez_81624',
            'Entrez_2321', 'Entrez_57593', 'Entrez_3488', 'Entrez_10531', 'Entrez_57211', 'Entrez_8476', 'Entrez_8659', 'Entrez_57758', 'Entrez_163', 
            'Entrez_2781', 'Entrez_1058', 'Entrez_445815', 'Entrez_58475', 'Entrez_2947', 'Entrez_10455', 'Entrez_8293', 'Entrez_1633', 'Entrez_9055', 
            'Entrez_11082', 'Entrez_55321', 'Entrez_85453']

        # get the index of the genes

        featureGeneIndices = numpy.argwhere(numpy.in1d(dataset.geneLabels, sig)).flatten()

        return VantVeerFeatureExtractor(dataset.geneLabels, featureGeneIndices)


class VantVeerFeatureExtractor(object):
    name               = "VantVeerGeneSignature"
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
