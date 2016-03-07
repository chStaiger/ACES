# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import numpy
import scipy.spatial.distance as distance

class FourClassNearestMeanClassifierFactory(object):

    def __init__(self, scoringFunction):

        self.scoringFunction = scoringFunction
        self.productName     = "FourClassNearestMeanClassifier_" + self.scoringFunction.__name__

    def train(self, featureValues, classLabels, ERlabels):
        """
        featureValues:  Numpy array samples x genes.
        classLabels:    Numpy array. Entries must have the same ordering as the patients in the data.
        ERlabels:       Numpy array. Entries must have the same ordering as the patients in the data.
        """

        # Check dimensionality of training data
        (ns, nf) = featureValues.shape
        assert classLabels.shape == (ns, )
        assert ERlabels.shape == (ns, )

        # Check that the class labels are all booleans
        assert frozenset(classLabels) == frozenset([False, True])
        assert frozenset(ERlabels) == frozenset([False, True])

        # Find the mean values for the False and True classes
        meanERposFalse = numpy.mean(featureValues[numpy.where(numpy.logical_not(classLabels) & ERlabels)[0]], axis = 0)
        meanERnegFalse = numpy.mean(featureValues[numpy.where(numpy.logical_not(classLabels) & numpy.logical_not(ERlabels))[0]], axis = 0)
        meanERposTrue  = numpy.mean(featureValues[numpy.where(classLabels & ERlabels)[0]], axis = 0)
        meanERnegTrue  = numpy.mean(featureValues[numpy.where(classLabels & numpy.logical_not(ERlabels))[0]], axis = 0)

        return FourClassNearestMeanClassifier(self.scoringFunction, meanERposFalse, meanERposTrue, meanERnegFalse, meanERnegTrue)

class FourClassNearestMeanClassifier(object):

    def __init__(self, scoringFunction, meanERposFalse, meanERposTrue, meanERnegFalse, meanERnegTrue):

        self.scoringFunction = scoringFunction
        self.meanERposFalse       = meanERposFalse
        self.meanERnegFalse       = meanERnegFalse
        self.meanERposTrue        = meanERposTrue
        self.meanERnegTrue        = meanERnegTrue

    def score(self, samples):
        #get closest mean
        scores = []
        for sample in range(samples.shape[0]):
            dist = [] #check with correlation which centroid is closest
            #close to ERpos
            dist.append(distance.pdist(numpy.row_stack([samples[sample, ], self.meanERposFalse]), 'correlation')[0])
            dist.append(distance.pdist(numpy.row_stack([samples[sample, ], self.meanERposTrue]), 'correlation')[0])
            dist.append(distance.pdist(numpy.row_stack([samples[sample, ], self.meanERnegFalse]), 'correlation')[0])
            dist.append(distance.pdist(numpy.row_stack([samples[sample, ], self.meanERnegTrue]), 'correlation')[0])
            #normal binary scoring of the classification, closest centroid determines the ER means
            argmin = numpy.argmin(dist)
            #score as in binary classifier but only with ER pos orER neg means, 
            #depending which ER mean is closest
            if argmin <=1: #ERpos
                scores.append(self.scoringFunction(self.meanERposFalse, self.meanERposTrue, samples[sample, ]))
            else: #ERneg
                scores.append(self.scoringFunction(self.meanERnegFalse, self.meanERnegTrue, samples[sample, ]))
    
        return numpy.array(scores)

# V1, V2a, V2b, and V3 implement different distance metrics for the NMC classifiers.
#
# Note that the calculations as given here are normalized in such a way that the scores are invariant under
# rotations, translations, and scalings of the vectors {sample, meanFalse, meanTrue}.

def V1(meanFalse, meanTrue, sample):
    """
    This NMC distance metric projects the sample onto the line from meanFalse -> meanTrue, and normalizes the value;
    Points that project to meanFalse are scored as 0 (zero), points that project to meanTrue are scored as 1 (one).
    """
    return numpy.inner(sample - meanFalse, meanTrue - meanFalse) / numpy.inner(meanTrue - meanFalse, meanTrue - meanFalse)

def V2a(meanFalse, meanTrue, sample):
    """
    This NMC distance metric scores a sample by considering the (meanFalse, sample, meanTrue) triangle;
    The score is calculated as the length (Euclidean distance) of the (meanFalse, sample) side, divided over the distance of meanFalse meanTrue via 'sample'.
    """
    return numpy.linalg.norm(sample - meanFalse) / (numpy.linalg.norm(sample - meanFalse) + numpy.linalg.norm(sample - meanTrue))


def V2b(meanFalse, meanTrue, sample):
    """
    This NMC distance metric scores a sample by subtracting its Euclidean distance to meanTrue from its distance to meanFalse.
    The score is normalized to the distance between {meanFalse, meanTrue}.
    The iso-score lines in this case are hyperbolas.
    """
    return (numpy.linalg.norm(sample - meanFalse) - numpy.linalg.norm(sample - meanTrue)) / numpy.linalg.norm(meanTrue - meanFalse)


def V3(meanFalse, meanTrue, sample):
    """
    This NMC distance metric scores samples by considering a point that is halfway {meanFalse, meanTrue}, then calculating the
    cosine of the angle {sample, halfway, meanTrue}.
    Points towards meanTrue get a score of close to +1, while points towards meanFalse get a score close to -1.
    """
    halfway = 0.5 * (meanFalse + meanTrue)
    return numpy.inner(sample - halfway, meanTrue - halfway) / (numpy.linalg.norm(sample - halfway) * numpy.linalg.norm(meanTrue - halfway))

