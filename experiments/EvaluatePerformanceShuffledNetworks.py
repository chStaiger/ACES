# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013
import sqlite3
import itertools
from collections import Counter
import numpy as np

import matplotlib as mpl
mpl.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from statistics.Wilcoxon import Wilcoxon

def getResults_outerLoop(filename):

    dbConnection = sqlite3.connect(filename)
    dbCursor = dbConnection.cursor()
    dbCursor.execute("SELECT dataset, feature_extractor, network_dataset, shuffleNr, classifier, repeatNr, \
        foldNr, number_of_features, auc_value FROM Results WHERE classifier == 'BinaryNearestMeanClassifier_V1' AND repeatNr == 0;")
    dbResults = dbCursor.fetchall()
    dbConnection.close()

    return dbResults

def getResults_innerLoop(filename):
    
    dbConnection = sqlite3.connect(filename)
    dbCursor = dbConnection.cursor()

    # Average over repeats and fold numbers while extracting!
    dbCursor.execute("SELECT dataset, foldNr, network_dataset, shuffleNr, feature_extractor, classifier, number_of_features, \
        AVG(auc_value) AS avg_auc_value FROM Results \
        WHERE classifier == 'BinaryNearestMeanClassifier_V1' AND dataset LIKE '%RFS%' AND repeatNr == 0 GROUP BY dataset, feature_extractor, \
        network_dataset, shuffleNr, classifier, number_of_features, foldNr;")
    
    results = dbCursor.fetchall()

    dbCursor.execute("SELECT dataset, foldNr, network_dataset, shuffleNr, feature_extractor, classifier, number_of_features, \
        AVG(auc_value) AS avg_auc_value FROM Results \
        WHERE classifier == 'BinaryNearestMeanClassifier_V1' AND dataset LIKE '%DMFS%' AND repeatNr == 0 GROUP BY dataset, feature_extractor, \
        network_dataset, shuffleNr, classifier, number_of_features, foldNr;")

    results.extend(dbCursor.fetchall())

    return results

def formatMethodName(name):
    """
    name = method+'_'+network or just method (for Single genes etc.) or 
    """
    name = name.replace('ChuangFeatureExtractor', 'C').replace('LeeFeatureExtractor',
                    'L').replace('TaylorFeatureExtractor', 'T').replace('DaoFeatureExtractor', 'D').replace('IPP', 'NetC').replace('AvgHubDiff',
                    '').replace('pw', '').replace('SingleGeneFeatureExtractor', 'SG').replace('ErasmusMCGeneSignature',
                    'Erasmus').replace('VantVeerGeneSignature', 'NKI').replace('RandomGeneFeatureExtractor',
                    'Random').replace('GeneRankTscoreFeatureExtractor', 'GR-Tstat').replace('GeneRankFeatureExtractor',
                    'GR').replace('WinterFeatureExtractor_SurvTime', 'W-Time').replace('WinterFeatureExtractor',
                    'W')

    return name

def make_dot_boxplot(matrix, validXticks, title, xlabel, ylabel, pp_dot, sgUp = None, sgDown = None):
    """
    xticks: list of integer numbers ranging from 0 to 24, you can leave out numbers.
    """
    xticks = range(1, 26)
    validXticks = [tick+1 for tick in validXticks]
    #random.seed(0)
    #jitter = np.array([random.randint(-3, 3)/10. for r in range(len(matrix[0]))])
    #plot grey sqare
    plt.axis([0.5, len(xticks)+0.5, 0.4, 0.9])
    if sgUp != None and sgDown != None:
        print sgUp, sgDown
        rect = plt.Rectangle((0.5, sgDown), len(xticks)+0.5, sgUp-sgDown, facecolor="lightgrey", edgecolor = "lightgrey")
        plt.gca().add_patch(rect)
    #print "DotPlot"
    for i in range(len(validXticks)):
        plt.plot((validXticks[i])*np.ones(len(matrix[i]), dtype=np.int), matrix[i], 'o', color = 'blue')
        plt.plot([(validXticks[i])*np.ones(len(matrix[i]))-0.3, (validXticks[i])*np.ones(len(matrix[i]))+0.3], 
            [np.mean(matrix[i]), np.mean(matrix[i])], '-r')
    plt.ylabel(ylabel, fontsize=13)
    plt.xlabel(xlabel, fontsize=13)
    plt.title(title, fontsize=20)

    plt.grid()
    ax = plt.axes()
    ax.set_xticks(range(1,len(xticks)+1))
    #if matrix.shape[0] == len(xticks):
        #ax.set_xticklabels(xticks, rotation=90, ha = 'right', fontsize=6)
    #    ax.set_xticklabels(xticks, rotation=90, fontsize=10)
    fig = plt.gcf()
    ax1 = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

    fig.set_size_inches(max(5, len(xticks)),5)
    fig.tight_layout()
    pp_dot.savefig()


def evaluateInnerLoop(databases):
    """
    databases = [
        'EXP05Inner_Dao.sqlite3',
        'EXP05InnerLoopPart1Jul26.sqlite3',
        'EXP05InnerLoopPart1Jul28.sqlite3',
        'EXP05InnerLoopPart2Jul26.sqlite3',
        'EXP05InnerLoopPartJul19-Part1.sqlite3',
        'EXP05InnerLoopPartJul20-Part1.sqlite3',
        'EXP05InnerLoopPartJul24.sqlite3'
    ]
    """
    dbResults = []
    for f in databases:
        print f
        dbResults.extend(getResults_innerLoop(f))
    dbResults = list(set(dbResults))
    
    DatasetsAndFolds = set([(r[0], r[1]) for r in dbResults])
    featureExtractorsAndNetworks = set([(r[4], r[2]) for r in dbResults])
    classifiers = list(set([r[5] for r in dbResults]))
    classifier = classifiers[0]

    shuffles = set([shuffle for d, f, n, shuffle, fe, c, num, auc in dbResults])

    #All methods apart from Winter/GeneRank
    #1) get the data for one method, network
    #2) sort the items descending accoridng to the average auc
    #3) get number of features
    selectedF = [(f, n) for (f, n) in featureExtractorsAndNetworks if (not f.startswith("Winter") and not f.startswith("GeneRank"))]
    selectedFAndShuffles = list(itertools.product(selectedF, shuffles))
    MethodsToNumFeatures = {}
    for dataset, fold in DatasetsAndFolds:
        dataResults = [item for item in dbResults if dataset in item]
        for (method, network), shuffle in selectedFAndShuffles:
            if 'Lee' in method:
                continue
            data = [(auc_value, number_of_features, d, f) for (d, f, n, s, m, c,
                number_of_features, auc_value) in dataResults if d == dataset and method == m
                and c == classifier and network == n and s == shuffle and f == fold]
            if len(data) == 0:
                print 'No data for', dataset, fold , (method, network), shuffle
                continue
            sorted_by_first = sorted(data, key=lambda tup: tup[0], reverse=True)
            methodF = formatMethodName(method+'_'+network)
            dataName = dataset.split('_fold')[0]
            MethodsToNumFeatures[(dataName+'_'+str(fold), methodF, shuffle)] = sorted_by_first[0][1]
    #Winter and GeneRank:
    # 1) get the data for method, damping factor and network
    #2) sort the items descending accoridng to the average auc
    #3) get number of features and dampingfactor
    networks = set([n for m, n in featureExtractorsAndNetworks if m.startswith("Winter") or m.startswith("GeneRank")])
    prefixMethods = ['GeneRankFeatureExtractor', 'GeneRankTscore', 'WinterFeatureExtractor_', 'WinterFeatureExtractor_SurvTime']
    combi = list(itertools.product(prefixMethods, networks, shuffles))
    for dataset, fold in DatasetsAndFolds:
        dataResults = [item for item in dbResults if dataset in item]
        for method, network, shuffle in combi:
            if method == 'WinterFeatureExtractor_':
                data = [(auc_value, number_of_features, m, d, f) for (d, f, n, s, m, c,
                    number_of_features, auc_value) in dataResults if d == dataset and m.startswith(method)
                    and c == classifier and network == n and 'SurvTime' not in m and s == shuffle and fold == f]
            else:
                data = [(auc_value, number_of_features, m, d, f) for (d, f, n, s, m, c,
                    number_of_features, auc_value) in dataResults if d == dataset and m.startswith(method)
                    and c == classifier and network == n and s == shuffle and fold == f]
            if len(data) == 0:
                print 'No data for', dataset, fold , (method, network), shuffle
                continue
            sorted_by_first = sorted(data, key=lambda tup: tup[0], reverse=True)
            methodF = formatMethodName(sorted_by_first[0][2]+'_'+network)
            dataName = dataset.split('_fold')[0]
            MethodsToNumFeatures[(dataName+'_'+str(fold), methodF, shuffle)] = sorted_by_first[0][1]

    return MethodsToNumFeatures

def getAUCmatrixOuter(results, selectedF, datasets, MethodsToNumfeat):
    """
    Given the number of best performing or fixed features for each method
    return a matrix methods x  with the auc values, a list with the
    correctly sorted methods and the LatexTable.
    
    selectedF   :   (method, network) of the methods of interest
    """
    repeat = 0
    DatasetClassifierToAUCmatrix = {} #values: (matrix, list)
    DICT_methodToMethodF = {}
    for m in selectedF:
        DICT_methodToMethodF[formatMethodName(m[0]+'_'+m[1])] = m
    for dataset in datasets:
        dataResults = [item for item in results if dataset in item]
        classifier = 'BinaryNearestMeanClassifier_V1'
        matrix = []
        methods = []
        visited = []
        validMethods = [(tmpM, shuffle) for (d, tmpM, shuffle) in MethodsToNumfeat if d.startswith(dataset)]
        for (tmpM, shuffle) in validMethods:
            print tmpM, shuffle
            methodName = tmpM.split('_')[0]+'_'
            network = tmpM.split('_')[len(tmpM.split('_'))-1]
            if (dataset, methodName, network, shuffle) in visited:
                continue
            visited.append((dataset, methodName, network, shuffle))
            allFeaturekeys = [(d, m, s) for (d, m, s) in MethodsToNumfeat if d.startswith(dataset) 
                and m.startswith(methodName) and m.endswith(network) and s == shuffle]
            data = []
            if len(allFeaturekeys) < 0:
                print 'No data for:', dataset, methodName, network, shuffle        
                continue
            elif len(allFeaturekeys) < 5:
                print 'Entries for:', dataset, methodName, network, shuffle, ':', len(allFeaturekeys)
            for key in allFeaturekeys:
                print key, MethodsToNumfeat[key]
                # ds, fe, n, sh, c, repeat, fold, nf, auc
                maxNumFeat = max([number_of_features for (d, m, n, s, c, rn, fn,
                    number_of_features, auc_value) in dataResults if d == dataset and  m == DICT_methodToMethodF[key[1]][0]
                    and c == classifier and n == DICT_methodToMethodF[key[1]][1] and rn == repeat
                    and fn == int(key[0][len(key[0])-1]) and s == shuffle]+[0])
                print 'Maximum', maxNumFeat
                if maxNumFeat < MethodsToNumfeat[key]:
                    MethodsToNumfeat[key] = maxNumFeat
                if maxNumFeat == 0:
                    continue
                data.extend([(auc_value, fn, number_of_features) for (d, m, n, s, c, rn, fn,
                    number_of_features, auc_value) in dataResults if d == dataset and m == DICT_methodToMethodF[key[1]][0]
                    and c == classifier and rn == repeat and fn == int(key[0][len(key[0])-1]) and s == shuffle
                    and number_of_features == MethodsToNumfeat[key] and n == DICT_methodToMethodF[key[1]][1]])
            sorted_by_second = sorted(data, key=lambda tup: tup[1])
            data = [auc for auc, _, _ in sorted_by_second]
            matrix.append(data)
            methods.append((methodName+network, shuffle))
        DatasetClassifierToAUCmatrix[(dataset, classifier.replace('BinaryNearestMeanClassifier', 'NMC'))] = (matrix, methods)

    return DatasetClassifierToAUCmatrix


def evaluateOuterLoop(databases, num = None):
    """
    databases = [
        "EXP05_Dao_18-07.sqlite3",
        #"EXP05_backup_done_July1.sqlite3",
        "EXP05PartJul11.sqlite3",
        "EXP05PartJul15.sqlite3",
        "EXP05PartJul19-Part1.sqlite3",
        "EXP05PartJul20-Part1.sqlite3"
    ]    
    """

    dbResults = []
    for f in databases:
        print f
        dbResults.extend(getResults_outerLoop(f))

    dbResults  = list(set(dbResults))  

    # ds, fe, n, sh, c, repeat, fold, nf, auc
    datasets = set(r[0] for r in dbResults)
    featureExtractorsAndNetworks = set([(r[1], r[2]) for r in dbResults])
    shuffles = set([r[3] for r in dbResults]) 
    classifiers = set([r[4] for r in dbResults])
    repeats = set([r[5] for r in dbResults])

    if num != None:
        bestFeatures = False
        MethodsToNumfeat = makeDictMethodToNum(num, featureExtractorsAndNetworks, shuffles)
    else:
        bestFeatures = True
        num = 'best_features'
        databases = [
            'EXP05Inner_Dao.sqlite3',
            'EXP05InnerLoopPart1Jul26.sqlite3',
            'EXP05InnerLoopPart1Jul28.sqlite3',
            'EXP05InnerLoopPart2Jul26.sqlite3',
            'EXP05InnerLoopPartJul19-Part1.sqlite3',
            'EXP05InnerLoopPartJul20-Part1.sqlite3',
            'EXP05InnerLoopPartJul24.sqlite3'
        ]
        MethodsToNumfeat = evaluateInnerLoop(databases)

    pageRankF = [(f, n) for (f, n) in featureExtractorsAndNetworks if f.startswith("Winter") or f.startswith("GeneRank")]
    DatasetClassifierToAUCmatrix = getAUCmatrixOuter(dbResults, featureExtractorsAndNetworks, datasets, 
        MethodsToNumfeat)
    #We need DatasetClassifierToAUCmatrix also for the overlap, since it is costly to recalculate save it!
    import pickle
    pickle.dump((MethodsToNumfeat, DatasetClassifierToAUCmatrix), open('EXP05_bestFeatures', 'wb'))
 
    #get results from unshuffled networks to compare
    _, EXP01DatasetClassifierToAUCmatrix = pickle.load(open('EXP01_bestFeatures.pickle'))

    #Plot
    ttest = []
    for dataset, classifier in DatasetClassifierToAUCmatrix: 
        matrix, methodsAndShuffles = DatasetClassifierToAUCmatrix[(dataset, classifier)]
        matrix = np.array(matrix)
        methodsAndShuffles = np.array(methodsAndShuffles)
        matrixEXP01, methodsEXP01 = EXP01DatasetClassifierToAUCmatrix[(dataset, classifier)]
        matrixEXP01 = np.array(matrixEXP01)
        methodsEXP01 = np.array(methodsEXP01)
        ttest.extend([dataset+' '+str(num) + ' & & ', ' & Method & p-vals (no multiple testing correction)'])
        methods = list(set([m for m, s in methodsAndShuffles]))
        visitedPageRank = []
        for method in sorted(methods):
            #collect all entries for the method
            indices = []
            idxEXP01 = np.argwhere(methodsEXP01 == method).flatten()[0]
            for i in range(len(methodsAndShuffles)):
                m, s = methodsAndShuffles[i]
                if m == method:
                    indices.append(i)
            methodF = method
            up = max(matrixEXP01[idxEXP01, :])
            down = min(matrixEXP01[idxEXP01, :])
            #collect all entries for the shuffles
            validShuffles = [int(s) for m, s in methodsAndShuffles[indices]]           
            #sort shuffles
            indices = np.array(indices)[np.argsort(validShuffles)] 
            validShuffles = np.array(validShuffles)[np.argsort(validShuffles)]
            aucValues = matrix[indices, :]
            if [] in list(aucValues):
                indices = [i for i, x in enumerate(aucValues) if x == []]
                validShuffles = list((set(range(len(validShuffles))).difference(indices)))
                aucValues = aucValues[validShuffles]
            s = 'Results/Shuffles/Boxplot_Experiment_05_%s_%s_%s.pdf' % (dataset, methodF.replace('.', ''), str(num).replace(' ', ''))
            pp_box = PdfPages(s)
            make_dot_boxplot(aucValues, validShuffles, dataset+' '+methodF+' '+str(num), 'shuffles', 'AUC', pp_box, sgUp = up, sgDown = down)
            pp_box.close()
            plt.clf()
            #test difference with Wilcoxon
            ttest.append( ' & '+methodF+' & '+str("%0.4f" % (Wilcoxon(matrixEXP01[idxEXP01, :], np.array(list(itertools.chain(*aucValues))), False))))
    with open('Results/pVals_EXP05_'+str(num)+'.txt','w') as f:
        print >> f, "\\\\\n".join(ttest)


