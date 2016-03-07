# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013
import sys, os, h5py, sqlite3, json, numpy, itertools

# Input of primary and secondary datasets

from datatypes.ExpressionDataset import HDF5GroupToExpressionDataset
from datatypes.ExpressionDataset import MakeRandomFoldMap
from datatypes.EdgeSet           import ReadSIF
from datatypes.GeneSetCollection import ReadGeneSetCollection
from experiments.AnnotationU133  import readAnnotation

# Feature-extraction algorithms

from featureExtractors.SingleGenes.SingleGeneFeatureExtractor import SingleGeneFeatureExtractorFactory
from featureExtractors.SingleGenes.RandomGeneFeatureExtractor import RandomGeneFeatureExtractorFactory
from featureExtractors.Chuang.ChuangFeatureExtractor          import ChuangFeatureExtractorFactory
from featureExtractors.Lee.LeeFeatureExtractor                import LeeFeatureExtractorFactory
from featureExtractors.Dao.DaoFeatureExtractor                import DaoFeatureExtractorFactory
from featureExtractors.Taylor.TaylorFeatureExtractor          import TaylorFeatureExtractorFactory
from featureExtractors.GeneSignatures.VantVeerFeatureExtractor      import VantVeerFeatureExtractorFactory
from featureExtractors.GeneSignatures.ErasmusMCFeatureExtractor     import ErasmusMCFeatureExtractorFactory

from featureExtractors.Winter.WinterFeatureExtractor          import WinterFeatureExtractorFactory
from featureExtractors.Winter.GeneRankFeatureExtractor        import GeneRankFeatureExtractorFactory
from featureExtractors.Winter.GeneRankTscoreFeatureExtractor  import GeneRankTscoreFeatureExtractorFactory

from classifiers.BinaryNearestMeanClassifier                  import BinaryNearestMeanClassifierFactory, V1, V2a, V2b, V3
from statistics.PerformanceCurve import CalculateFeatureCountDependentPerformanceCurve, CalculateFeatureCountDependentPerformance

import random


def CombineDataExperiment05(numShuffles):
    """
    Combines the methods, their parameters and the data for the "Shuffle Networks" experiment.
    """ 
    datasets = ['U133A_combat_DMFS', 'U133A_combat_RFS']
    networks = [
       'nwEdgesKEGG',  
       'nwEdgesHPRD9', 
       'nwEdgesI2D',   
       'nwEdgesIPP'
    ]
    pathways = ['nwGeneSetsKEGG', 'nwGeneSetsMsigDB'] 

    NetworkFeatureExtractorsSpecificParams = [
        ('Chuang', None),
        ('Taylor', True),
        ('Taylor', False),
        ('Dao', None),
    ]
    dampingF = [x * 0.1 for x in range(0, 11)]
    NetworkFeatureExtractorsSpecificParams.extend(list(itertools.product(['GeneRank', 'GeneRankTscore', 'Winter', 'WinterTime'], dampingF)))

    PathwayFeatureExtractorSpecificParams = [('Lee', None)]

    #Combine algos with sec data
    NetworkDataAndFeatureExtractors = list(itertools.product(NetworkFeatureExtractorsSpecificParams, networks))
    NetworkDataAndFeatureExtractors.extend(list(itertools.product(PathwayFeatureExtractorSpecificParams, pathways)))

    #Combine FEs with data and number of shuffles
    DataAndFeatureExtractors = list(itertools.product(datasets, NetworkDataAndFeatureExtractors, range(numShuffles)))
    
    return DataAndFeatureExtractors

def CombineDataExperiment01():
    """
    Combines the data, methods and parameters for the Performance and Overlap evaluation.
    """
    datasets = ['U133A_combat_DMFS', 'U133A_combat_RFS']
    networks = [
       'nwEdgesKEGG',
       'nwEdgesHPRD9',
       'nwEdgesI2D',
       'nwEdgesIPP'
    ]
    pathways = ['nwGeneSetsKEGG', 'nwGeneSetsMsigDB']

    NetworkFeatureExtractorsSpecificParams = [
        ('Chuang', None),
        ('Taylor', True),
        ('Taylor', False),
        ('Dao', None),
    ]
    dampingF = [x * 0.1 for x in range(0, 11)]
    NetworkFeatureExtractorsSpecificParams.extend(list(itertools.product(['GeneRank', 'GeneRankTscore', 'Winter', 'WinterTime'], dampingF)))

    PathwayFeatureExtractorSpecificParams = [('Lee', None)]
    SingleGene = [('SingleGenes', None), ('RandomGenes', None), ('NKI', None), ('Erasmus', None)]

    #Combine algos with sec data
    NetworkDataAndFeatureExtractors = list(itertools.product(NetworkFeatureExtractorsSpecificParams, networks))
    NetworkDataAndFeatureExtractors.extend(list(itertools.product(PathwayFeatureExtractorSpecificParams, pathways)))
    NetworkDataAndFeatureExtractors.extend(list(itertools.product(SingleGene, [None])))
    #Combine FEs with data and number of shuffles
    DataAndFeatureExtractors = list(itertools.product(datasets, NetworkDataAndFeatureExtractors, [None]))

    return DataAndFeatureExtractors

def CombineDataExperiment04():
    """
    Combines the methods, their parameters and the data for the "ER positive patients only" experiment.
    Data is not shown in the paper.
    """
    datasets = ['U133A_combat_DMFS_ERpos', 'U133A_combat_RFS_ERpos']
    networks = [
       'nwEdgesKEGG',
       'nwEdgesHPRD9',
       'nwEdgesI2D',
       'nwEdgesIPP'
    ]
    pathways = ['nwGeneSetsKEGG', 'nwGeneSetsMsigDB']

    NetworkFeatureExtractorsSpecificParams = [
        ('Chuang', None),
        ('Taylor', True),
        ('Taylor', False),
        ('Dao', None),
    ]
    dampingF = [x * 0.1 for x in range(0, 11)]
    NetworkFeatureExtractorsSpecificParams.extend(list(itertools.product(['GeneRank', 'GeneRankTscore', 'Winter', 'WinterTime'], dampingF)))

    PathwayFeatureExtractorSpecificParams = [('Lee', None)]
    SingleGene = [('SingleGenes', None), ('RandomGenes', None), ('NKI', None), ('Erasmus', None)]

    #Combine algos with sec data
    NetworkDataAndFeatureExtractors = list(itertools.product(NetworkFeatureExtractorsSpecificParams, networks))
    NetworkDataAndFeatureExtractors.extend(list(itertools.product(PathwayFeatureExtractorSpecificParams, pathways)))
    NetworkDataAndFeatureExtractors.extend(SingleGene)
    #Combine FEs with data and number of shuffles
    DataAndFeatureExtractors = list(itertools.product(datasets, NetworkDataAndFeatureExtractors, [None]))

    return DataAndFeatureExtractors
    
    
def SetUpRun(dataset, network, method):
    
    #get dataset
    if '_ERpos' in dataset:
        f = h5py.File("experiments/data/U133A_combat_ERpos.h5")
    else:
        f = h5py.File("experiments/data/U133A_combat.h5")
    
    data = [HDF5GroupToExpressionDataset(f[group]) for group in f.keys() if dataset in group][0]
    f.close()

    #get network
    if network == "nwGeneSetsKEGG":
        net   = ReadGeneSetCollection("KEGGpw" , "experiments/data/KEGG1210_PathwayGeneSets_Entrez.txt" , "Entrez_")
    elif network == "nwGeneSetsMsigDB":
        net = ReadGeneSetCollection("MsigDBpw" , "experiments/data/C2V3_PathwayGeneSets_Entrez.txt"     , "Entrez_")
    elif network == "nwEdgesKEGG":
        net = ReadSIF("KEGG"  , "experiments/data/KEGG_edges1210.sif" , "Entrez_")
    elif network == "nwEdgesHPRD9":
        net = ReadSIF("HPRD9" , "experiments/data/HPRD9.sif"          , "Entrez_")
    elif network == "nwEdgesI2D":
        net = ReadSIF("I2D"   , "experiments/data/I2D_edges_0411.sif" , "Entrez_")
    elif network == "nwEdgesIPP":
        net = ReadSIF("IPP"   , "experiments/data/ipp.sif"            , "Entrez_")
    elif network == None:
        net = None
        print "SG no network"
    else:
        raise Exception("Network not known. Add network in SetUpGrid line 167.")

    #get survival time for Wintertime
    Dataset2Time = None
    if method == 'WinterTime':
        print 'WinterTime'
        annotationfile = 'experiments/data/@expdesc_breast_3597_qc_110916.xls'        
        classlabelData = "experiments/data/U133A_combat_classlabels.pickle"
        Dataset2Time = readAnnotation(annotationfile, [data])

    #get featureselection method
    if method == "Lee":
        featureSelector = LeeFeatureExtractorFactory()
    elif method == "Chuang":
        featureSelector = ChuangFeatureExtractorFactory("featureExtractors/Chuang/pinnaclez-OPTIONAL_NORMALIZE.jar", False)
    elif method.startswith("Taylor"): 
        featureSelector = TaylorFeatureExtractorFactory()
    elif method == 'Dao':
        featureSelector = DaoFeatureExtractorFactory("featureExtractors/Dao/runwDCB")
    elif method == 'GeneRank':
        featureSelector = GeneRankFeatureExtractorFactory() 
    elif method == 'GeneRankTscore':
        featureSelector = GeneRankTscoreFeatureExtractorFactory()
    elif method.startswith('Winter'):
        featureSelector = WinterFeatureExtractorFactory()
    elif method.startswith('SingleGenes'):
        featureSelector = SingleGeneFeatureExtractorFactory()
    elif method.startswith('RandomGenes'):
        featureSelector = RandomGeneFeatureExtractorFactory()
    elif method == 'NKI':
        featureSelector = VantVeerFeatureExtractorFactory()
    elif method == 'Erasmus':
        featureSelector = ErasmusMCFeatureExtractorFactory()
    else:
        raise Exception("Method not known. Add method in SetUpGrid line 201.")

    #get classifier
    classifiers = [
        BinaryNearestMeanClassifierFactory(V1),
        BinaryNearestMeanClassifierFactory(V2a),
        BinaryNearestMeanClassifierFactory(V2b),
        BinaryNearestMeanClassifierFactory(V3),
    ]

    return (data, net, featureSelector, classifiers, Dataset2Time)

def RunInstance(data, net, featureSelector, special, classifiers, repeat, nrFolds, fold, shuffleNr, survTime = None, TaylorParam = True):

    #in case if the inner loop the dataset name has an addition tag.
    if "training" in data.name:
        dName = "_".join(data.name.split('_')[:len(data.name.split('_'))-2])
    else:
        dName = data.name

    #shuffle network
    shuffle = None
    if shuffleNr != None:
        seed = (shuffleNr, net.name)
        shuffle = net.makeShuffle(seed = seed)
        shuffledNetworkDataset = net.makeShuffledVersion("%s_shuffle_%d" % (net.name, shuffleNr), shuffle)
    else:
        shuffledNetworkDataset = net
    
    #split datasets
    dsTraining, dsTesting, foldMap = splitData(data, fold, repeat, nrFolds)    

    if featureSelector.productName in ["SingleGeneFeatureExtractor", "RandomGeneFeatureExtractor", "VantVeerGeneSignature", "ErasmusMCGeneSignature"]:
        featureExtractor = featureSelector.train(dsTraining)
    elif featureSelector.productName.startswith("Taylor"): #if True : average over edges between hub and interactors, if False take each edge as one feature.
        featureExtractor = featureSelector.train(dsTraining, shuffledNetworkDataset, average = TaylorParam)
    elif special == None: #Chuang, Dao, Lee
        featureExtractor = featureSelector.train(dsTraining, shuffledNetworkDataset)
    elif featureSelector.productName == "WinterFeatureExtractor" and survTime != None: #WinterTime
        survTimeTrain = survTime[dName][numpy.array(foldMap.foldAssignments)-1 != fold]
        print len(survTimeTrain), dsTraining.numPatients
        featureExtractor = featureSelector.train(dsTraining, (shuffledNetworkDataset,
            float(special)), survTime = survTimeTrain)
    elif featureSelector.productName == "WinterFeatureExtractor":
        featureExtractor = featureSelector.train(dsTraining, (shuffledNetworkDataset,
            float(special)), survTime = None)
    else: #GeneRank
        featureExtractor = featureSelector.train(dsTraining, (shuffledNetworkDataset, float(special)))

    #NOTE: not valid for gene signatures, SCP features and the 4-Class classifier
    maxFeatureCount = 400
    AucAndCi = {}
    for CF in classifiers:
        print "-->", CF.productName
        if featureExtractor.name in ["VantVeerGeneSignature", "ErasmusMCGeneSignature"]:
            #These methods have a fixed number of features that all go into the classifier, no ranking.
            nf_to_auc = CalculateFeatureCountDependentPerformance(
                featureExtractor,
                CF,
                (dsTraining, dsTesting)
            )
        else:
            featureCounts = [fc for fc in featureExtractor.validFeatureCounts if fc <= maxFeatureCount]
            nf_to_auc = CalculateFeatureCountDependentPerformanceCurve(
                featureExtractor,
                CF,
                (dsTraining, dsTesting),
                featureCounts
            )
        AucAndCi[CF.productName] = nf_to_auc
    
    if net == None:
        return (data.name, featureExtractor.name, None, shuffle, featureExtractor.toJsonExpression(), AucAndCi)    
    else:
        return (data.name, featureExtractor.name, net.name, shuffle, featureExtractor.toJsonExpression(), AucAndCi)

def splitData(data, fold, repeat, nrFolds):
    #RunInnerLoopInstance(data, net,
    #                featureSelector, specific, classifiers, repeat, self.nrFolds, fold, shuffleNr, innerfold, innerrepeat, Dataset2Time, specific)

    #split datasets
    foldMap = MakeRandomFoldMap(data, nrFolds, repeat)
    foldList = range(0, nrFolds)
    dsTraining = data.extractPatientsByIndices("%s_fold-%d-of-%d_training" % (data.name, fold,
        len(foldList)), numpy.array(foldMap.foldAssignments)-1 != fold, checkNormalization = False)
    dsTesting = data.extractPatientsByIndices("%s_fold-%d-of-%d_testing"  % (data.name, fold,
        len(foldList)), numpy.array(foldMap.foldAssignments)-1 == fold, checkNormalization = False)
    
    return dsTraining, dsTesting, foldMap

def NextItem(iterable):
    try:
        first = next(iterable)
    except StopIteration:
        return None
    it = itertools.chain([first], iterable)
    return first, it

def getDoneTokens(db, sqlPath, sqlName = ""):
    """
    Sort all done Tokens according to experiment.
    db  :   dictionary item['_id']: item
    """    

    writeDocsExp01 = []
    writeDocsExp04 = []
    writeDocsExp05 = []
    writeDocsExp01InnerLoop = []
    writeDocsExp04InnerLoop = []
    writeDocsExp05InnerLoop = []

    otherDocs = []
    pending = []
    for item in db:
        doc = db[item]
        if 'output' in doc.keys():
            if item.startswith('EXP01InnerLoop'):
                writeDocsExp01InnerLoop.append(doc)
            elif item.startswith('EXP01'):
                writeDocsExp01.append(doc)
            elif item.startswith('EXP04InnerLoop'):
                writeDocsExp04InnerLoop.append(doc)
            elif item.startswith('EXP04'):
                writeDocsExp04.append(doc)
            elif item.startswith('EXP05InnerLoop'):
                writeDocsExp05InnerLoop.append(doc)
            elif item.startswith('EXP05'):
                writeDocsExp05.append(doc)
            else:
                otherDocs.append(doc)
        elif "lock" in doc.keys():
            if doc["lock"] > 0 and doc["done"] == 0:
                pending.append(doc)

    if len(writeDocsExp01) > 0:
        TokenToSqliteExperiment01(writeDocsExp01, sqlPath+"EXP01"+sqlName+".sqlite3")
    else:
        print "No results for Exp01."
    if len(writeDocsExp01InnerLoop) > 0:
        TokenToSqliteExperiment01InnerLoop(writeDocsExp01InnerLoop, sqlPath+"EXP01InnerLoop"+sqlName+".sqlite3")
    else:
        print "No results for Exp01InnerLoop."

    if len(writeDocsExp04) > 0:
        TokenToSqliteExperiment01(writeDocsExp04, sqlPath+"EXP04"+sqlName+".sqlite3")
    else:
        print "No results for Exp04."
    if len(writeDocsExp04InnerLoop) > 0:
        TokenToSqliteExperiment01InnerLoop(writeDocsExp04InnerLoop, sqlPath+"EXP04InnerLoop"+sqlName+".sqlite3")
    else:
        print "No results for Exp04InnerLoop."

    if len(writeDocsExp05) > 0:
        TokenToSqliteExperiment05(writeDocsExp05, sqlPath+"EXP05"+sqlName+".sqlite3")
    else:
        print "No results for Exp05."
    if len(writeDocsExp05InnerLoop) > 0:
        TokenToSqliteExperiment05InnerLoop(writeDocsExp05InnerLoop, sqlPath+"EXP05InnerLoop"+sqlName+".sqlite3")
    else:
        print "No results for Exp05InnerLoop."


    print "Other Docs found:", len(otherDocs)

    return pending, otherDocs


def getDoneTokens_lowMemory(db, sqlPath, sqlName = ""):
    """
    Fetches the results from a couchdb and stores them in an sql database.
    sqlPath :   path, must end with '/'
    sqlName :   an additional filename, there are default names, but for storing results with eg. a time stamp, use this parameter
    db      :   dictionary item['_id']: item or a couch.db

    returns pending tokens (for resetting them)
    """

    #all done tokens
    writeDocsExp01InnerLoop = (db[item] for item in db if 'output' in db[item] and item.startswith('EXP01InnerLoop'))
    writeDocsExp01 = (db[item] for item in db 
        if 'output' in db[item] and item.startswith('EXP01') and not item.startswith('EXP01InnerLoop'))
    writeDocsExp04InnerLoop = (db[item] for item in db if 'output' in db[item] and item.startswith('EXP04InnerLoop'))
    writeDocsExp04 = (db[item] for item in db 
        if 'output' in db[item] and item.startswith('EXP04') and not item.startswith('EXP04InnerLoop')) 
    writeDocsExp05InnerLoop = (db[item] for item in db if 'output' in db[item] and item.startswith('EXP05InnerLoop'))
    writeDocsExp05 = (db[item] for item in db 
        if 'output' in db[item] and item.startswith('EXP05') and not item.startswith('EXP05InnerLoop')) 
    pending = (db[item] for item in db if db[item]['output']['lock'] > 0 and db[item]['output']['done'])
    
    first, writeDocsExp01InnerLoop = NextItem(writeDocsExp01InnerLoop)
    if first != None:
        TokenToSqliteExperiment01InnerLoop(writeDocsExp01InnerLoop, sqlPath+"EXP01InnerLoop"+sqlName+".sqlite3")
    else:
        print "No results for Exp01Inner."
    first, writeDocsExp01 = NextItem(writeDocsExp01)
    if first != None:
        TokenToSqliteExperiment01(writeDocsExp01, sqlPath+"EXP01"+sqlName+".sqlite3")
    else:
        print "No results for Exp01."
    first, writeDocsExp04InnerLoop = NextItem(writeDocsExp04InnerLoop)
    if first != None:
        TokenToSqliteExperiment01InnerLoop(writeDocsExp04InnerLoop, sqlPath+"EXP04InnerLoop"+sqlName+".sqlite3")
    else:
        print "No results for Exp04Inner."
    first, writeDocsExp04 = NextItem(writeDocsExp04)
    if first != None:
        TokenToSqliteExperiment01(writeDocsExp04, sqlPath+"EXP04"+sqlName+".sqlite3")
    else:
        print "No results for Exp04."
    first, writeDocsExp05InnerLoop = NextItem(writeDocsExp05InnerLoop)
    if first != None:
        TokenToSqliteExperiment05InnerLoop(writeDocsExp05InnerLoop, sqlPath+"EXP05InnerLoop"+sqlName+".sqlite3")
    else:
        print "No results for Exp05Inner."
    first, writeDocsExp05 = NextItem(writeDocsExp05)
    if first != None:
        TokenToSqliteExperiment05(writeDocsExp01InnerLoop, sqlPath+"EXP05InnerLoop"+sqlName+".sqlite3")
    else:
        print "No results for Exp05."

    return pending


def resetPendingTokens(ids, db):

    for ID in ids:
        token = db.get(ID)
        scrub = token['scrub_count']
        updateContent = {'lock': 0, 
            'scrub_count' : scrub+1}
        token.update(updateContent)
        db[token['_id']] = token

def TokenToSqliteExperiment01(tokens, sqlFilename):
    try:
        os.remove(sqlFilename)
        print "NOTE: Removed existing", sqlFilename
    except OSError:
        pass # This happens if the file did not yet exist, which is ok.

    dbConnection = sqlite3.connect(sqlFilename)
    dbCursor = dbConnection.cursor()

    dbCursor.execute("""
        CREATE TABLE FeatureExtractors (
            dataset            TEXT    NOT NULL , -- the name of the dataset used for training;
            repeatNr           INTEGER NOT NULL , -- the repeat number;
            foldNr             INTEGER NOT NULL , -- the number of features used for training and testing;
            feature_extractor  TEXT    NOT NULL , -- the name of the feature extractor used;
            network_dataset    TEXT        NULL , -- the name of the network dataset (if any). Null if no secondary dataset was used ("Single Genes");
            feature_definition TEXT    NOT NULL   -- the feature definition, as returned by the FeatureExtractor's "toJsonExpression" method.
        );
        """)

    dbCursor.execute("""
        CREATE TABLE Results (
            dataset            TEXT    NOT NULL , -- the name of the dataset used for training and testing;
            feature_extractor  TEXT    NOT NULL , -- the name of the feature extractor used;
            network_dataset    TEXT        NULL , -- the name of the network dataset (if any). NULL if no secondary dataset was used (e.g. "Single Genes");
            classifier         TEXT    NOT NULL , -- the name of the classifier used;
            repeatNr           INTEGER NOT NULL , -- the repeat number;
            foldNr             INTEGER NOT NULL , -- the number of features used for training and testing;
            number_of_features INTEGER NOT NULL , -- the number of features used for training and testing;
            auc_value          REAL        NULL , -- the resulting AUC value. NaN is represented as NULL;
            ci                 REAL    NOT NULL   -- the resulting concordance index. NaN is represented as NULL.
        );
        """)

    for token in tokens:
        (dataName, featureExtractorproductName, netName, shuffleNr, shuffle, featureExtractor, AucAndCi) = token['output']
        if shuffle != None:
            continue
        fold = token['input']['fold']
        repeat = token['input']['repeat']
        specific = token['input']['specific']
        assert token['input']['dataset'] == dataName
        
        entry = (dataName, repeat, fold, featureExtractorproductName, netName, featureExtractor)
        dbCursor.execute("INSERT INTO FeatureExtractors(dataset, repeatNr, foldNr, feature_extractor, network_dataset, \
            feature_definition) VALUES(?, ?, ?, ?, ?, ?);", entry)
        dbConnection.commit()

        for CL in AucAndCi:
            nf_to_auc = AucAndCi[CL]
            if specific != None:
                featureExtractorproductName+'_'+str(specific)
            entries = [(dataName, featureExtractorproductName, netName, CL, 
                repeat, fold, featureCount, AucCi[0], AucCi[1]) for (featureCount, AucCi) in nf_to_auc.iteritems()]
            dbCursor.executemany("INSERT INTO Results(dataset, feature_extractor, network_dataset, classifier, repeatNr, \
                foldNr, number_of_features, auc_value, ci) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);", entries)
            dbConnection.commit()
    dbConnection.close()

    print "saved in", sqlFilename

def TokenToSqliteExperiment01InnerLoop(tokens, sqlFilename):
    try:
        os.remove(sqlFilename)
        print "NOTE: Removed existing", sqlFilename
    except OSError:
        pass # This happens if the file did not yet exist, which is ok.

    dbConnection = sqlite3.connect(sqlFilename)
    dbCursor = dbConnection.cursor()

    dbCursor.execute("""
        CREATE TABLE FeatureExtractors (
            dataset            TEXT    NOT NULL , -- the name of the dataset used for training;
            repeatNr           INTEGER NOT NULL , -- the repeat number;
            repeatNrInnerloop  INTEGER NOT NULL , -- the repeat number inner loop;
            foldNr             INTEGER NOT NULL , -- the fold used for training;
            foldNrInnerloop    INTEGER NOT NULL , -- the fold used for training inner loop;
            feature_extractor  TEXT    NOT NULL , -- the name of the feature extractor used;
            network_dataset    TEXT        NULL , -- the name of the network dataset (if any). Null if no secondary dataset was used ("Single Genes");
            feature_definition TEXT    NOT NULL   -- the feature definition, as returned by the FeatureExtractor's "toJsonExpression" method.
        );
        """)

    dbCursor.execute("""
        CREATE TABLE Results (
            dataset            TEXT    NOT NULL , -- the name of the dataset used for training and testing;
            feature_extractor  TEXT    NOT NULL , -- the name of the feature extractor used;
            network_dataset    TEXT        NULL , -- the name of the network dataset (if any). NULL if no secondary dataset was used (e.g. "Single Genes");
            classifier         TEXT    NOT NULL , -- the name of the classifier used;
            repeatNr           INTEGER NOT NULL , -- the repeat number;
            repeatNrInnerloop  INTEGER NOT NULL , -- the repeat number;
            foldNr             INTEGER NOT NULL , -- the number of features used for training and testing;
            foldNrInnerloop    INTEGER NOT NULL , -- the number of features used for training and testing;
            number_of_features INTEGER NOT NULL , -- the number of features used for training and testing;
            auc_value          REAL        NULL , -- the resulting AUC value. NaN is represented as NULL;
            ci                 REAL        NULL   -- the resulting concordance index. NaN is represented as NULL.
        );
        """)
    for token in tokens:
        (dataName, featureExtractorproductName, netName, shuffleNr, shuffle, featureExtractor, AucAndCi) = token['output']
        if shuffle != None:
            continue
        fold = token['input']['fold']
        repeat = token['input']['repeat']
        specific = token['input']['specific']
        repeatInner = token['input']['innerRepeat']
        foldInner = token['input']['innerFold']

        entry = (dataName, repeat, repeatInner, fold, foldInner, featureExtractorproductName, netName, featureExtractor)
        dbCursor.execute("INSERT INTO FeatureExtractors(dataset, repeatNr, repeatNrInnerloop, foldNr, foldNrInnerloop, feature_extractor, network_dataset, \
            feature_definition) VALUES(?, ?, ?, ?, ?, ?, ?, ?);", entry)
        dbConnection.commit()

        for CL in AucAndCi:
            nf_to_auc = AucAndCi[CL]
            if specific != None:
                featureExtractorproductName+'_'+str(specific)
            entries = [(dataName, featureExtractorproductName, netName, CL,
                repeat, repeatInner, fold, foldInner, featureCount, AucCi[0], AucCi[1]) for (featureCount, AucCi) in nf_to_auc.iteritems()]
            dbCursor.executemany("INSERT INTO Results(dataset, feature_extractor, network_dataset, classifier, repeatNr, repeatNrInnerloop, \
                foldNr, foldNrInnerloop, number_of_features, auc_value, ci) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", entries)
            dbConnection.commit()
    dbConnection.close()

    print "saved in", sqlFilename



def TokenToSqliteExperiment05(tokens, sqlFilename):
    
    try:
        os.remove(sqlFilename)
        print "NOTE: Removed existing", sqlFilename
    except OSError:
        pass # This happens if the file did not yet exist, which is ok.


    dbConnection = sqlite3.connect(sqlFilename)
    dbCursor = dbConnection.cursor()

    dbCursor.execute("""
        CREATE TABLE Shuffles (

            network_dataset    TEXT    NOT NULL , -- the name of the network dataset;
            shuffleNr          INTEGER NOT NULL , -- the shuffle number;
            shuffle            TEXT    NOT NULL , -- the shuffle specification (Python dict as JSON expression).

            PRIMARY KEY (network_dataset, shuffleNr)
        );
        """)

    dbCursor.execute("""
        CREATE TABLE FeatureExtractors (
            training_dataset   TEXT    NOT NULL , -- the name of the dataset used for training;
            repeatNr           INTEGER NOT NULL , -- the repeat number;
            foldNr             INTEGER NOT NULL , -- the number of features used for training and testing;
            feature_extractor  TEXT    NOT NULL , -- the name of the feature extractor used;
            network_dataset    TEXT        NULL , -- the name of the network dataset (if any). Null if no secondary dataset was used ("Single Genes");
            shuffleNr          INTEGER NOT NULL , -- the shuffle number;
            feature_definition TEXT    NOT NULL , -- the feature definition, as returned by the FeatureExtractor's "toJsonExpression" method.


            FOREIGN KEY (network_dataset, shuffleNr) REFERENCES Shuffles (network_dataset, shuffleNr)
        );
        """)

    dbCursor.execute("""
        CREATE TABLE Results (
            dataset            TEXT    NOT NULL , -- the name of the dataset used for training;
            feature_extractor  TEXT    NOT NULL , -- the name of the feature extractor used;
            network_dataset    TEXT        NULL , -- the name of the network dataset (if any). Null if no secondary dataset was used ("Single Genes");
            shuffleNr          INTEGER NOT NULL , -- the shuffle number;
            classifier         TEXT    NOT NULL , -- the name of the classifier used;
            repeatNr           INTEGER NOT NULL , -- the repeat number;
            foldNr             INTEGER NOT NULL , -- the number of features used for training and testing;
            number_of_features INTEGER NOT NULL , -- the number of features used for training and testing;
            auc_value          REAL        NULL , -- the resulting AUC value. NaN is represented as NULL;
            ci                 REAL        NULL , -- the resulting AUC value. NaN is represented as NULL.


            FOREIGN KEY (network_dataset, shuffleNr) REFERENCES Shuffles (network_dataset, shuffleNr)
        );
        """)

    savedShuffles = set()
    for token in tokens:
        (dataName, featureExtractorproductName, netName, shuffleNr, shuffle, featureExtractor, AucAndCi) = token['output']
        if shuffle == None:
            continue
        if featureExtractor == None:
            print 'No features found for:', dataName, featureExtractorproductName, netName, shuffleNr
            continue
        fold = token['input']['fold']
        repeat = token['input']['repeat']
        specific = token['input']['specific']
        assert token['input']['dataset'] == dataName

        if (netName, shuffleNr, json.dumps(shuffle)) not in savedShuffles:
            dbCursor.execute("INSERT INTO Shuffles(network_dataset, shuffleNr, shuffle) VALUES(?, ?, ?)", [netName, shuffleNr, json.dumps(shuffle)])
            dbConnection.commit()
            savedShuffles.add((netName, shuffleNr, json.dumps(shuffle)))
        entry = (dataName, repeat, fold, featureExtractorproductName, netName, shuffleNr, featureExtractor)
        dbCursor.execute("INSERT INTO FeatureExtractors(training_dataset, repeatNr, foldNr, feature_extractor, network_dataset, shuffleNr,  \
            feature_definition) VALUES(?, ?, ?, ?, ?, ?, ?);", entry)
        dbConnection.commit()
        for CL in AucAndCi:
            nf_to_auc = AucAndCi[CL]
            if specific != None:
                featureExtractorproductName+'_'+str(specific)
            entries = [(dataName, featureExtractorproductName, netName, shuffleNr, CL, 
                repeat, fold, featureCount, AucCi[0], AucCi[1]) for ((featureCount, AucCi)) in nf_to_auc.iteritems()]
            dbCursor.executemany("INSERT INTO Results(dataset, feature_extractor, network_dataset, shuffleNr, classifier, repeatNr, \
                foldNr, number_of_features, auc_value, ci) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", entries)
            dbConnection.commit()
    dbConnection.close()
    print "saved in", sqlFilename 

def TokenToSqliteExperiment05Inner(tokens, sqlFilename):

    try:
        os.remove(sqlFilename)
        print "NOTE: Removed existing", sqlFilename
    except OSError:
        pass # This happens if the file did not yet exist, which is ok.


    dbConnection = sqlite3.connect(sqlFilename)
    dbCursor = dbConnection.cursor()

    dbCursor.execute("""
        CREATE TABLE Shuffles (

            network_dataset    TEXT    NOT NULL , -- the name of the network dataset;
            shuffleNr          INTEGER NOT NULL , -- the shuffle number;
            shuffle            TEXT    NOT NULL , -- the shuffle specification (Python dict as JSON expression).

            PRIMARY KEY (network_dataset, shuffleNr)
        );
        """)

    dbCursor.execute("""
        CREATE TABLE FeatureExtractors (
            training_dataset   TEXT    NOT NULL , -- the name of the dataset used for training;
            repeatNr           INTEGER NOT NULL , -- the repeat number;
            repeatNrInnerloop  INTEGER NOT NULL , -- the repeat number;
            foldNr             INTEGER NOT NULL , -- the number of features used for training and testing;
            foldNrInnerloop    INTEGER NOT NULL , -- the number of features used for training and testing;
            feature_extractor  TEXT    NOT NULL , -- the name of the feature extractor used;
            network_dataset    TEXT        NULL , -- the name of the network dataset (if any). Null if no secondary dataset was used ("Single Genes");
            shuffleNr          INTEGER NOT NULL , -- the shuffle number;
            feature_definition TEXT    NOT NULL , -- the feature definition, as returned by the FeatureExtractor's "toJsonExpression" method.


            FOREIGN KEY (network_dataset, shuffleNr) REFERENCES Shuffles (network_dataset, shuffleNr)
        );
        """)

    dbCursor.execute("""
        CREATE TABLE Results (
            dataset            TEXT    NOT NULL , -- the name of the dataset used for training;
            feature_extractor  TEXT    NOT NULL , -- the name of the feature extractor used;
            network_dataset    TEXT        NULL , -- the name of the network dataset (if any). Null if no secondary dataset was used ("Single Genes");
            shuffleNr          INTEGER NOT NULL , -- the shuffle number;
            classifier         TEXT    NOT NULL , -- the name of the classifier used;
            repeatNr           INTEGER NOT NULL , -- the repeat number;
            repeatNrInnerloop  INTEGER NOT NULL , -- the repeat number;
            foldNr             INTEGER NOT NULL , -- the number of features used for training and testing;
            foldNrInnerloop    INTEGER NOT NULL , -- the number of features used for training and testing;
            number_of_features INTEGER NOT NULL , -- the number of features used for training and testing;
            auc_value          REAL        NULL , -- the resulting AUC value. NaN is represented as NULL;
            ci                 REAL        NULL , -- the resulting AUC value. NaN is represented as NULL.


            FOREIGN KEY (network_dataset, shuffleNr) REFERENCES Shuffles (network_dataset, shuffleNr)
        );
        """)
    savedShuffles = set()
    for token in tokens:
        (dataName, featureExtractorproductName, netName, shuffleNr, shuffle, featureExtractor, AucAndCi) = token['output']
        if shuffle == None:
            continue
        if featureExtractor == None:
            print 'No features found for:', dataName, featureExtractorproductName, netName, shuffleNr
            continue
        fold = token['input']['fold']
        repeat = token['input']['repeat']
        specific = token['input']['specific']
        repeatInner = token['input']['innerRepeat']
        foldInner = token['input']['innerFold']

        if (netName, shuffleNr, json.dumps(shuffle)) not in savedShuffles:
            dbCursor.execute("INSERT INTO Shuffles(network_dataset, shuffleNr, shuffle) VALUES(?, ?, ?)", [netName, shuffleNr, json.dumps(shuffle)])
            dbConnection.commit()
            savedShuffles.add((netName, shuffleNr, json.dumps(shuffle)))
        entry = (dataName, repeat, repeatInner, fold, foldInner, featureExtractorproductName, netName, shuffleNr, featureExtractor)
        dbCursor.execute("INSERT INTO FeatureExtractors(training_dataset, repeatNr, repeatNrInnerloop, foldNr, foldNrInnerloop, \
            feature_extractor, network_dataset, shuffleNr,  \
            feature_definition) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?);", entry)
        dbConnection.commit()
        for CL in AucAndCi:
            nf_to_auc = AucAndCi[CL]
            if specific != None:
                featureExtractorproductName+'_'+str(specific)
            entries = [(dataName, featureExtractorproductName, netName, shuffleNr, CL,
                repeat, repeatInner, fold, foldInner, featureCount, AucCi[0], AucCi[1]) for ((featureCount, AucCi)) in nf_to_auc.iteritems()]
            dbCursor.executemany("INSERT INTO Results(dataset, feature_extractor, network_dataset, shuffleNr, classifier, repeatNr, \
                repeatNrInnerloop, foldNr, foldNrInnerloop, number_of_features, auc_value, ci) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);", entries)
            dbConnection.commit()
    dbConnection.close()
    print "saved in", sqlFilename

