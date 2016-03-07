import numpy

# Make sure we can transparantly transfer NumPy types to and from R
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate() # In the new version of RPy2 you need to activate that feature.


class BinaryLogisticRegressionClassifierFactory(object):

    def __init__(self, trainingFunction):

        self.trainingFunction = trainingFunction
        self.productName      = "BinaryLogisticRegressionClassifier_" + self.trainingFunction.__name__

    def train(self, featureValues, classLabels):

        # Check dimensionality of training data
        (ns, nf) = featureValues.shape
        assert classLabels.shape == (ns, )

        # Check that the class labels are all booleans
        assert frozenset(classLabels) == frozenset([False, True])

        beta = self.trainingFunction(featureValues, classLabels)

        assert beta is None or isinstance(beta, numpy.ndarray)

        if isinstance(beta, numpy.ndarray) and not all(numpy.isfinite(beta)):
            print "WARNING: attempt to initialize a BinaryLogisticRegressionClassifier '%s' with a weird 'beta' vector; treating as beta = None (beta = %s)" % (name, beta)
            beta = None

        if beta is None:
            return None

        return BinaryLogisticRegressionClassifier(beta)


class BinaryLogisticRegressionClassifier(object):

    def __init__(self, beta):

        assert isinstance(beta, numpy.ndarray) and all(numpy.isfinite(beta))
        self.beta = beta

    def score(self, samples):

        (ns, nf) = samples.shape
        assert self.beta.shape == (nf + 1,)

        # What we return here is a number that is monotonously increasing with the predicted probability.
        #
        # This is all the score function needs to do, because it may only be used
        #   to rank samples according to probability.
        #
        # The advantage of not returning the probability is that we prevent trouble with
        #   large negative of z (overflow) and large positive values of z (underflow).
        #
        # If you do want the actual predicted probabilities, use the "probability()" method that is defined below.

        return numpy.inner(samples, self.beta[1:])

    def probability(self, samples):

        # If the training process that produced us did not work out, we cannot score the samples.
        if self.beta is None:
            return None

        (ns, nf) = samples.shape
        assert self.beta.shape == (nf + 1,)

        z = self.beta[0] + numpy.inner(samples, self.beta[1:])

        return 1.0 / (1.0 + numpy.exp(-z))


def R_GLM(featureValues, classLabels):

    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri

    (ns, nf) = featureValues.shape
    assert classLabels.shape == (ns, )

    # Suppress warnings
    robjects.r.options(warn = -1)

    robjects.globalenv["featureValues"] = featureValues
    robjects.globalenv["classLabels"] = classLabels

    fit = robjects.r.glm("classLabels ~ featureValues", family = "binomial")

    # Check convergence as reported by R's glm() function; abort if not converged.
    if not fit.rx("converged"):
        return None

    beta = numpy.array(fit.rx("coefficients"))

    # Reorder the coefficients
    beta = numpy.append(beta[-1], beta[:-1])

    if not all(numpy.isfinite(beta)):
        return None

    return beta


def R_NNET(featureValues, classLabels):

    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    from   rpy2.robjects.packages import importr

    # Suppress warnings
    robjects.r.options(warn = -1)

    robjects.packages.importr("nnet")

    robjects.globalenv["featureValues"] = featureValues
    robjects.globalenv["classLabels"] = classLabels

    fit = robjects.r.multinom("classLabels ~ featureValues", trace = False)

    # NOTA BENE (1): The "convergence" field, according to the documentation, is zero
    #   "if the maximum number of iterations is reached", i.e., when the search did NOT converge!
    #
    # NOTA BENE (2): The 'convenience' field is a vector of an int-vector, hence the [0][0]

    if fit.rx("convergence")[0][0] != 0:
        return None

    beta = numpy.array(fit.rx("wts"))

    assert beta[0, 0] == 0 # get rid of this dummy zero coefficient
    beta = beta[0, 1:]

    if not all(numpy.isfinite(beta)):
        return None

    return beta


def SciKits(featureValues, classLabels):

    # TODO: check how non-convergence is reported

    (ns, nf) = featureValues.shape
    assert classLabels.shape == (ns, )
    #NOTE: There is a version issue
    try:
        from scikits.learn.linear_model import LogisticRegression
    except:
        #scikits version >= 0.9
        from sklearn.linear_model       import LogisticRegression
    RegularisationParameter = 1.0e+30 # For us, higher is better (no regularisation!)
    Tolerance               = 1.0e-30 # Smaller is better

    # From the documentation page at
    #
    #     http://scikit-learn.sourceforge.net/modules/generated/scikits.learn.linear_model.LogisticRegression.html
    #
    # "The underlying C implementation uses a random number generator to select features when fitting the model.
    #  It is thus not uncommon, to have slightly different results for the same input data.
    #  If that happens, try with a smaller tol parameter."

    classifier = LogisticRegression(penalty = 'l1', C = RegularisationParameter, tol = Tolerance)
    classifier.fit(featureValues, classLabels)

    beta = -classifier.raw_coef_[0,:]
    beta = numpy.append(beta[-1], beta[:-1])

    if not all(numpy.isfinite(beta)):
        return None

    return beta


def BioPython(featureValues, classLabels):

    (ns, nf) = featureValues.shape
    assert classLabels.shape == (ns, )

    from Bio.LogisticRegression import train, calculate

    try:
        classifier = train(featureValues, classLabels)
    except:
        # BioPython will throw an exception if it detects non-convergence
        return None

    beta = numpy.array(classifier.beta)

    if not all(numpy.isfinite(beta)):
        return None

    return beta


def R_GLM_Verified(featureValues, classLabels):

    # Use multiple independent methods to solve for the beta vector.
    # If any of these fails, we fail, too.

    beta = R_GLM(featureValues, classLabels)
    if beta is None:
        return None

    # The value for TOLERANCE was determined by inspection of many cases.
    # If convergence of the methods occurs, it will usually occur to at least 10 decimal digits;
    #    if it doesn't occur, values can vary by 1 or more.
    # This value was chosen to clearly distinguish those cases.

    TOLERANCE = 1.0e-3

    for verificationTrainer in [R_NNET, SciKits]:

        betaCheck = verificationTrainer(featureValues, classLabels)

        if betaCheck is None:
            return None

        # If any of the beta components deviates more than TOLERANCE, report non-convergence.
        if any(numpy.abs(beta - betaCheck) > TOLERANCE):
            return None

    # All good. Return R_GLM value of beta.
    return beta
