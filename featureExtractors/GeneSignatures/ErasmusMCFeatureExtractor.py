# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import numpy, json
from statistics.Statistics import tStatisticForUnequalSampleSizeAndUnequalVariance
import xlrd

class ErasmusMCFeatureExtractorFactory(object):

    productName = "ErasmusMCGeneSignature"

    def train(self, dataset):

        geneEntrez = [
            'Entrez_2055',            'Entrez_51280',            'Entrez_10897',
            'Entrez_3606',            'Entrez_11335',            'Entrez_1846',
            'Entrez_5501',            'Entrez_3983',            'Entrez_8741',
            'Entrez_58986',            'Entrez_718',            'Entrez_54892',
            'Entrez_745',            'Entrez_4134',            'Entrez_10579',
            'Entrez_10051',            'Entrez_26529',            'Entrez_23595',
            'Entrez_29028',            'Entrez_10051',            'Entrez_9840',
            'Entrez_1846',            'Entrez_960',            'Entrez_5347',
            'Entrez_10256',            'Entrez_8365',            'Entrez_2237',
            'Entrez_2286',            'Entrez_3838',            'Entrez_51093',
            'Entrez_81577',            'Entrez_161',            'Entrez_79682',
            'Entrez_9148',            'Entrez_11198',            'Entrez_678',
            'Entrez_397',            'Entrez_5701',            'Entrez_55596',
            'Entrez_960',            'Entrez_10559',            'Entrez_9134',
            'Entrez_51131',            'Entrez_10721',            'Entrez_100287076',
            'Entrez_824',            'Entrez_2116',            'Entrez_54963',
            'Entrez_7940',            'Entrez_32',            'Entrez_636',
            'Entrez_1917',            'Entrez_4747',            'Entrez_50810',
            'Entrez_50810',            'Entrez_51512',            'Entrez_64864',
            'Entrez_1280',            'Entrez_2620',            'Entrez_149076',
            'Entrez_2525',            'Entrez_4620',            'Entrez_8743',
            'Entrez_118433',            'Entrez_55879',            'Entrez_143',
            'Entrez_25906',            'Entrez_79370',            'Entrez_9702',
            'Entrez_26027']


        geneEntrez = list(set(geneEntrez))
        featureGeneIndices = numpy.argwhere(numpy.in1d(dataset.geneLabels, geneEntrez)).flatten()

        return ErasmusMCFeatureExtractor(dataset.geneLabels, featureGeneIndices)


class ErasmusMCFeatureExtractor(object):
    name               = "ErasmusMCGeneSignature"
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
