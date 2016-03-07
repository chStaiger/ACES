# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013
from SetUpGrid import SetUpRun, RunInstance, splitData

def localProcess(token, db = None):
    """
    Processes the token. Number of folds is fixed to 5!
    db  :   specific for couchDB, leave empty or adopt the section after NOTE for a specific output.
    """

    dataset = token['input']['dataset']
    network = token['input']['network']
    method = token['input']['method']
    specific = token['input']['specific']
    repeat = token['input']['repeat']
    fold = token['input']['fold']
    shuffleNr =  token['input']['shuffleNr']
    
    print 'dataset:', dataset
    print 'network', network
    print 'method', method
    print 'specific', specific
    print 'repeat', repeat
    print 'fold', fold
    print 'shuffleNr', shuffleNr

    innerCV = False
    if 'innerFold' in token['input']:
        innerCV = True
        innerfold = token['input']['innerFold']
        innerrepeat = token['input']['innerRepeat']
    
    (data, net, featureSelector, classifiers, Dataset2Time) = SetUpRun(dataset, network, method)
    
    if not innerCV:
        if specific == True or specific == False:
            (dataName, featureExtractorproductName, netName, shuffle, featureExtractor, AucAndCi) =  RunInstance(data, net,
                featureSelector, specific, classifiers, repeat, 5, fold, shuffleNr, Dataset2Time, specific)
        else:
            (dataName, featureExtractorproductName, netName, shuffle, featureExtractor, AucAndCi) =  RunInstance(data, net,
                featureSelector, specific, classifiers, repeat, 5, fold, shuffleNr, Dataset2Time)
    else:
        dsOuterTraining, dsOuterTesting, _ = splitData(data, repeat, fold, 5)
        if specific == True or specific == False:
            (dataName, featureExtractorproductName, netName, shuffle, featureExtractor, AucAndCi) =  RunInstance(dsOuterTraining, net,
                featureSelector, specific, classifiers, innerrepeat, 5, innerfold, shuffleNr, Dataset2Time, specific)
        else:
            (dataName, featureExtractorproductName, netName, shuffle, featureExtractor, AucAndCi) =  RunInstance(dsOuterTraining, net,
                featureSelector, specific, classifiers, innerrepeat, 5, innerfold, shuffleNr, Dataset2Time)

    token['output'] = (dataName, featureExtractorproductName, netName, shuffleNr, shuffle, featureExtractor, AucAndCi)
    token['done'] = 1
    token['lock'] = 1
    #NOTE: specific for couchDB
    if db != None:
        db.update([token])
    
    return token

