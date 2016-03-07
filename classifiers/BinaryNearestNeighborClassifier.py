
import numpy

class BinaryNearestNeighborClassifierFactory(object):

    def __init__(self, k, alpha, beta):

        self.k           = k
        self.alpha       = alpha
        self.beta        = beta
        self.productName = "BinaryNearestNeighiborClassifier_%d" % k

    def train(self, featureValues, classLabels):

        # Check dimensionality of training data
        (ns, nf) = featureValues.shape
        assert classLabels.shape == (ns, )

        # Check that the class labels are all booleans
        assert frozenset(classLabels) == frozenset([False, True])

        # Find the mean values for the False and True classes

        return BinaryNearestNeighborClassifier(featureValues, classLabels, self.k, self.alpha, self.beta)


class BinaryNearestNeighborClassifier(object):

    def __init__(self, featureValues, classLabels, k, alpha, beta):

        self.k             = k
        self.alpha         = alpha
        self.beta          = beta
        self.featureValues = featureValues
        self.classLabels   = classLabels

    def scoreSample(self, sample):

        distances = numpy.apply_along_axis(lambda x : numpy.linalg.norm(x), 1, self.featureValues - sample)
        idx = numpy.argsort(distances)
        idx = idx[:self.k]

        signs = self.classLabels[idx]
        signs = 2 * signs.astype(numpy.int) - 1

        distances = distances[idx]
        #weights = numpy.exp(- self.alpha * distances ** self.beta)
        epsilon = 1e-6
        weights = 1.0 / (distances + epsilon)

        return numpy.dot(signs, weights)

    def score(self, samples):

        return numpy.apply_along_axis(lambda sample: self.scoreSample(sample), 1, samples)

class BinaryNearestNeighborClassifierAbsDiffFactory(object):

    def __init__(self, k):

        self.k           = k
        self.productName = "BinaryNearestNeighorClassifierAbsDiff_%d" % k

    def train(self, featureValues, classLabels):

        # Check dimensionality of training data
        (ns, nf) = featureValues.shape
        assert classLabels.shape == (ns, )

        # Check that the class labels are all booleans
        assert frozenset(classLabels) == frozenset([False, True])

        # Find the mean values for the False and True classes

        return BinaryNearestNeighborClassifierAbsDiff(featureValues, classLabels, self.k)

class BinaryNearestNeighborClassifierAbsDiff(object):

    def __init__(self, featureValues, classLabels, k):

        self.k             = k
        self.featureValues = featureValues
        self.classLabels   = classLabels

    def scoreAbsDiffSample(self, sample):

        distances = numpy.apply_along_axis(lambda x : sum(x), 1, numpy.abs(self.featureValues - sample))
        idx = numpy.argsort(distances)
        idx = idx[:self.k]

        signs = self.classLabels[idx]
        signs = 2 * signs.astype(numpy.int) - 1
        
        distances = distances[idx]
        #weights = numpy.exp(- self.alpha * distances ** self.beta)
        epsilon = 1e-6
        weights = 1.0 / (distances + epsilon)
        
        return numpy.dot(signs, weights)

    def score(self, samples):

        return numpy.apply_along_axis(lambda sample: self.scoreAbsDiffSample(sample), 1, samples)

