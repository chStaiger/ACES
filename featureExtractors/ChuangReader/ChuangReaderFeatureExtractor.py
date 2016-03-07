
import numpy

class ChuangReaderFeatureExtractorFactory(object):

    productName = "ChuangReaderFeatureExtractor"

    def __init__(self, directory):

        self.directory = directory

    def train(self, dataset, network):

        # We must do the ordering here

        return ChuangReaderFeatureExtractor(dataset.geneLabels, filtered_features)


class ChuangReaderFeatureExtractor(object):

    def __init__(self, geneLabels, features):

        self.geneLabels         = geneLabels
        self.features           = features
        self.validFeatureCounts = range(1, len(self.features) + 1)

    @staticmethod
    def score(expressionData, feature):
        return numpy.sum(expressionData[:, list(feature)], axis = 1) / numpy.sqrt(len(feature))

    def extract(self, dataset, k):

        assert all(dataset.geneLabels == self.geneLabels)
        assert k in self.validFeatureCounts

        return numpy.transpose(numpy.array([self.score(dataset.expressionData, feature) for feature in self.features[:k]]))
