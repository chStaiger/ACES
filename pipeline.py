# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013
from picas.clients import CouchClient
from couchdb import Server
from picas.modifiers import BasicTokenModifier
from picas.actors import RunActor
from picas.iterators import BasicViewIterator
from SetUpGrid import SetUpRun, RunInstance, splitData

class ExampleActor(RunActor):
    def __init__(self, iterator, modifier, nrFolds):
        self.iterator = iterator
        self.modifier = modifier
        self.client = iterator.client
        self.db = self.client.db
        self.nrFolds = nrFolds

    def prepare_env(self, *kargs, **kvargs):
        pass

    def prepare_run(self, *kargs, **kvargs):
        pass

    def process_token(self, ref, token):
       
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
                    featureSelector, specific, classifiers, repeat, self.nrFolds, fold, shuffleNr, Dataset2Time, specific)
            else: 
                (dataName, featureExtractorproductName, netName, shuffle, featureExtractor, AucAndCi) =  RunInstance(data, net, 
                    featureSelector, specific, classifiers, repeat, self.nrFolds, fold, shuffleNr, Dataset2Time) 
        else:
            dsOuterTraining, dsOuterTesting,_ = splitData(data, repeat, fold, self.nrFolds)
            print 'dsOuterTraining', dsOuterTraining
            print 'dsOuterTesting', dsOuterTesting
            if specific == True or specific == False:
                (dataName, featureExtractorproductName, netName, shuffle, featureExtractor, AucAndCi) =  RunInstance(dsOuterTraining, net,
                    featureSelector, specific, classifiers, innerrepeat, self.nrFolds, innerfold, shuffleNr, Dataset2Time, specific)
            else:
                (dataName, featureExtractorproductName, netName, shuffle, featureExtractor, AucAndCi) =  RunInstance(dsOuterTraining, net,
                    featureSelector, specific, classifiers, innerrepeat, self.nrFolds, innerfold, shuffleNr, Dataset2Time)
            
 
        token = self.modifier.close(token)
        token['output'] = (dataName, featureExtractorproductName, netName, shuffleNr, shuffle, featureExtractor, AucAndCi)
        self.db[token['_id']] = token

    def cleanup_run(self, *kargs, **kvargs):
        pass

    def cleanup_env(self, *kargs, **kvargs):
        pass


def main():
    nrFolds = 5
    client = CouchClient(url="http://<server>:<port>", db="<dbname>", username="<username>", password="***")
    modifier = BasicTokenModifier()
    iterator = BasicViewIterator(client, '<db-view>', modifier)
    actor = ExampleActor(iterator, modifier, nrFolds)
    actor.run()

if __name__ == '__main__':
    main()

