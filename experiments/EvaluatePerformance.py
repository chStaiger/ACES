# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013
import sqlite3
import time
import numpy as np
import random
import pprint
import itertools
import pickle
from collections import Counter

import matplotlib as mpl
mpl.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from statistics.Wilcoxon import Wilcoxon

def getResults(sqliteFile):

    print "OUTER LOOP RESULTS: Reading database ..."
    
    db = sqlite3.connect(sqliteFile)
    dbCursor = db.cursor()
    
    resultsDict = dict()
    # Average over repeats and fold numbers while extracting!
    dbCursor.execute("SELECT dataset, network_dataset, feature_extractor, classifier, repeatNr, foldNr, number_of_features, auc_value FROM Results \
        WHERE classifier == 'BinaryNearestMeanClassifier_V1' AND repeatNr == 0 ;")
    results = dbCursor.fetchall()

    return results

def getResultsInnerLoop(sqliteFile):
    
    print "INNER LOOP RESULTS: Reading database ..."

    db = sqlite3.connect(sqliteFile)
    dbCursor = db.cursor()

    # Average over repeats and fold numbers while extracting!
    dbCursor.execute("SELECT dataset, foldNr, network_dataset, feature_extractor, classifier, number_of_features, \
        AVG(auc_value) AS avg_auc_value FROM Results \
        WHERE classifier == 'BinaryNearestMeanClassifier_V1' AND dataset LIKE '%RFS%' AND repeatNr == 0 GROUP BY dataset, foldNr, feature_extractor, \
        network_dataset, classifier, number_of_features;")

    results = dbCursor.fetchall()

    dbCursor.execute("SELECT dataset, foldNr, network_dataset, feature_extractor, classifier, number_of_features, \
        AVG(auc_value) AS avg_auc_value FROM Results \
        WHERE classifier == 'BinaryNearestMeanClassifier_V1' AND dataset LIKE '%DMFS%' AND repeatNr == 0 GROUP BY dataset, foldNr, feature_extractor, \
        network_dataset, classifier, number_of_features;")
    results.extend(dbCursor.fetchall())

    return results

def formatMethodName(name):
    """
    name = method+'_'+network or just method (for Single genes etc.) or 
    """
    name = name.replace('Method', 'FeatureExtractor').replace('ChuangFeatureExtractor', 'C').replace('LeeFeatureExtractor',
                    'L').replace('TaylorFeatureExtractor', 'T').replace('DaoFeatureExtractor', 'D').replace('IPP', 'NetC').replace('AvgHubDiff',
                    '').replace('pw', '').replace('SingleGeneFeatureExtractor', 'SG').replace('ErasmusMCGeneSignature',
                    'Erasmus').replace('VantVeerGeneSignature', 'NKI').replace('RandomGeneFeatureExtractor',
                    'Random').replace('GeneRankTscoreFeatureExtractor', 'GR-Tstat').replace('GeneRankFeatureExtractor',
                    'GR').replace('WinterFeatureExtractor_SurvTime', 'W-Time').replace('WinterFeatureExtractor',
                    'W').replace('BinaryNearestMeanClassifier', 'NMC')
    return name

def make_dot_boxplot(matrix, xticks, title, xlabel, ylabel, pp_dot, sgUp = None, sgDown = None, box = False):

    random.seed(0)
    jitter = np.array([random.randint(-3, 3)/10. for r in range(len(matrix[0]))])
    #plot grey sqare
    plt.axis([0.5, len(xticks)+0.5, 0.5, 1])
    if sgUp != None and sgDown != None:
        print sgUp, sgDown
        rect = plt.Rectangle((0.5, sgDown), len(xticks)+0.5, sgUp-sgDown, facecolor="lightgrey", edgecolor = "lightgrey")
        plt.gca().add_patch(rect)

    #print "DotPlot"
    for i in range(0, matrix.shape[0]):
        plt.plot((i+1)*np.ones(len(matrix[0]), dtype=np.int), matrix[i], 'o', color = 'blue')
        plt.plot([(i+1)*np.ones(len(matrix[i]))-0.3, (i+1)*np.ones(len(matrix[i]))+0.3], [np.mean(matrix[i]), np.mean(matrix[i])], '-r')
    if box:
        plt.boxplot(matrix.T)

    plt.ylabel(ylabel, fontsize=13)
    plt.xlabel(xlabel, fontsize=13)
    plt.title(title, fontsize=13)

    plt.grid()
    ax = plt.axes()
    ax.set_xticks(range(1,len(matrix)+1))
    if matrix.shape[0] == len(xticks):
        #ax.set_xticklabels(xticks, rotation=90, ha = 'right', fontsize=6)
        ax.set_xticklabels(xticks, rotation=90, fontsize=10)
    fig = plt.gcf()
    ax1 = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)

    fig.set_size_inches(max(5, len(xticks)),5)
    fig.tight_layout()
    pp_dot.savefig()

def makeDictMethodToNum(num, methods, folds, repeat = 0):
    """
    Initialises a dictionary with methods to number iof features
    with one specific number of features for all methods.
    
    """
    m = [formatMethodName(item[0]+'_'+item[1]) for item in featureExtractorsAndNetworks if item[1] != None]
    m.extend([formatMethodName(item[0]) for item in featureExtractorsAndNetworks if item[1] == None])    
    datasets = []
        for spec in list(itertools.product(['U133A_combat_RFS', 'U133A_combat_DMFS'], repeats, range(folds))):
        datasets.append(spec[0]+'_'+str(spec[1])+'_'+str(spec[2]))


    dsAndm = list(itertools.product(datasets, m))
    MethodsToNumfeat = dict(zip(dsAndm, len(dsAndm)*[num]))
    for d in datasets:
        MethodsToNumFeatures[(d, 'NKI')] = 41
        MethodsToNumFeatures[(d, 'Erasmus')] = 66
    
    return MethodsToNumfeat


def getAUCmatrixOuter(results, selectedF, datasets, MethodsToNumfeat, nrFolds, bestFeatures = False):
    """
    Given the number of best performing features for each method
    return a matrix methods x  with the auc values, a list with the
    correctly sorted methods. Note that each dataset split is treated
    as a single dataset. This means from the inner loop we get the
    best performing number of features for a dataset split. Thus we receive
    different numbers of features for each split.
    
    selectedF   :   (method, network) of the methods of interest
    nrFolds     :   number of folds used to generate the data
    datasets    :   names of the datasets for the outer loop, they need to be prefixes of
                    the dataset names from the inner loop
    MethodsToNumFeat    :   (dataset_innerloop, method) to numnber fo features
    
    """
    repeat = 0
    DatasetClassifierToAUCmatrix = {} #values: (matrix, list)
    DICT_methodToMethodF = {}
    for m in selectedF:
        if m[1] != None:
            DICT_methodToMethodF[formatMethodName(m[0]+'_'+m[1])] = m
        else:
            DICT_methodToMethodF[formatMethodName(m[0])] = m

    #for dataset, classifier in itertools.product(datasets, classifiers):
    for dataset in datasets:
        dataResults = [item for item in results if dataset in item]
        classifier = 'BinaryNearestMeanClassifier_V1'
        matrix = []
        methods = []
        visited = []
        validMethods = [tmpM for (d, tmpM) in MethodsToNumfeat if d.startswith(dataset)]
        for tmpM in validMethods:
            print tmpM
            #For winter and generank different damping factors apply for different folds
            if len(tmpM.split('_')) > 1:
                methodName = tmpM.split('_')[0]+'_'
                network = tmpM.split('_')[len(tmpM.split('_'))-1]
            else:
                methodName = tmpM
                network = ''
            if (dataset, methodName, network) in visited:
                continue
            visited.append((dataset, methodName, network))
            allFeaturekeys = [(d, m) for (d, m) in MethodsToNumfeat if d.startswith(dataset) and m.startswith(methodName) and m.endswith(network)]
            #NOTE: that features range from 1 to max, not like indices which range from 0 to max-1
            data = []
            for key in allFeaturekeys:
                print key, MethodsToNumfeat[key]
                maxNumFeat = max([number_of_features for (d, n, m, c, rn, fn,
                    number_of_features, auc_value) in dataResults if d == dataset and  m == DICT_methodToMethodF[key[1]][0]
                    and c == classifier and n == DICT_methodToMethodF[key[1]][1] and rn == repeat 
                    and fn == int(key[0][len(key[0])-1])])
                print 'Maximum', maxNumFeat
                if maxNumFeat < MethodsToNumfeat[key]:
                    MethodsToNumfeat[key] = maxNumFeat
                data.extend([(auc_value, fn, number_of_features) for (d, n, m, c, rn, fn,
                    number_of_features, auc_value) in dataResults if d == dataset and m == DICT_methodToMethodF[key[1]][0] 
                    and c == classifier and rn == repeat and fn == int(key[0][len(key[0])-1])
                    and number_of_features == MethodsToNumfeat[key] and n == DICT_methodToMethodF[key[1]][1]])
            if len(data) < 5:
                print "Not enough features for:", methodName
            #sort data according to foldNr --> important for paired t-test
            sorted_by_second = sorted(data, key=lambda tup: tup[1])
            data = [auc for auc, _, _ in sorted_by_second]
            matrix.append(data)
            if network != None:
                methods.append(methodName+network)
            else:
                methods.append(methodName)
                
        DatasetClassifierToAUCmatrix[(dataset, classifier.replace('BinaryNearestMeanClassifier', 'NMC'))] = (matrix, methods)

    return DatasetClassifierToAUCmatrix 

def AucVsFeaturesCurve(dbResults):
    selected_classifier = u'BinaryNearestMeanClassifier_V1'
    data = [(dataset, feature_extractor+'_'+network_dataset, repeatNr, foldNr, number_of_features, auc_value) for 
                (dataset, network_dataset, feature_extractor, classifier, repeatNr, foldNr, number_of_features, auc_value)
                in dbResults if classifier == selected_classifier and network_dataset != None and repeatNr == 0]
    data.extend([(dataset, feature_extractor, repeatNr, foldNr, number_of_features, auc_value) for
        (dataset, network_dataset, feature_extractor, classifier, repeatNr, foldNr, number_of_features, auc_value)
        in dbResults if classifier == selected_classifier and feature_extractor == 'SingleGeneFeatureExtractor' and repeatNr == 0])

    list_datasets = sorted(set([dataset for (dataset, method, repeatNr, foldNr, number_of_features, auc_value) in data]))
    list_methods = sorted(set([method for (dataset, method, repeatNr, foldNr, number_of_features, auc_value) in data]))
    list_repeats = sorted(set([repeatNr for (dataset, method, repeatNr, foldNr, number_of_features, auc_value) in data]))
    list_folds = sorted(set([foldNr for (dataset, method, repeatNr, foldNr, number_of_features, auc_value) in data]))
    
    print
    print "###########    DATA   ###############"
    print "DATASETS: ", list_datasets
    print "METHODS: ", list_methods
    print "TYPE:, ", list_repeats, "times", list_folds, "CV"

    #mold the aucs from folds and splits into mean(AUC) and std(AUC) for each number of features
    print "Calculating mean(auc)"
    data_to_mean = {}
    for method in list_methods:
        print method, selected_classifier
        for dataset in list_datasets:
            range_features = sorted(set([f for (DS, m, repeatNr, foldNr, f, auc_value) in data
                if m == method and dataset == DS]))
            for feature in range_features:
                aucs = [auc_value for (DS, m, repeatNr, foldNr, f, auc_value) in data
                    if DS == dataset and m == method and f  == feature]
                #if len(aucs) == repeats*folds:
                data_to_mean[dataset, method, feature] = (np.mean(aucs), np.std(aucs, ddof=1))
    pickle.dump(data_to_mean, open('EXP01outer_AUCvsFEAT_DataMean.pickle', 'wb')) 
    
    colors = ['r', 'b', 'c', 'g', 'm', 'y']
    dataset_color = dict(zip(list_datasets, colors[0:len(list_datasets)]))

    for selected_method in list_methods:
        s = 'Plots/CV_curve_%s_%s.pdf' %(formatMethodName(selected_classifier), formatMethodName(selected_method))
        pp = PdfPages(s)
        #data_to_plotdata = {}
        # save the dotplot part in plots to draw legend
        plots = []
        len_axis = 0
        for dataset in list_datasets:
            keys = [(d, m, feat) for (d, m, feat) in data_to_mean.keys() if d == dataset and m == selected_method]
            features = []
            means = []
            stds = []
            for key in keys:
                features.append(key[2])
                means.append(data_to_mean[key][0])
                stds.append(data_to_mean[key][1])
            if len(features) == 0:
                print "No data for", selected_method
                continue

            features = np.array(features)
            means = np.array(means)
            stds = np.array(stds)
            idx = np.argsort(features)
            features = features[idx]
            stds = stds[idx]
            means = means[idx]

            #plot cuve
            plt.errorbar(features[5:len(features):5],
                means[5:len(features):5],
                yerr = stds[5:len(features):5], fmt = '.'+dataset_color[dataset])
            p, = plt.plot(features, means, '.-'+dataset_color[dataset])
            plots.append(p)
            plt.legend(plots, list_datasets, loc=4)
        len_axis = max(len_axis, len(features))
        plt.title(formatMethodName(selected_classifier+" "+selected_method), fontsize=16)
        plt.axis([0, len_axis, 0, 1])
        plt.ylabel("AUC", fontsize=15)
        plt.xlabel("Number of features", fontsize=15)
        pp.savefig()
        pp.close()
        plt.clf()
    return data_to_mean

def AucVsDampingFactor(data_to_mean, numFeat = 50):
    """
    data_to_mean from AucVsFeaturesCurve
    Read in with data_to_mean = pickle.load(open('EXP01outer_AUCvsFEAT_DataMean.pickle') in ACES_results
    """
    methods = set(['_'.join(key[1].split('_')[:len(key[1].split('_'))-2]) for key in data_to_mean.keys() if key[1].startswith('GeneRank') or key[1].startswith('Winter')])
    networks = set([key[1].split('_')[len(key[1].split('_'))-1] for key in data_to_mean.keys() if key[1].startswith('GeneRank') or key[1].startswith('Winter')])
    datasets = set([key[0] for key in data_to_mean.keys()])
    colors = ['r', 'b', 'c', 'g', 'm', 'y']
    dataset_color = dict(zip(list(datasets), colors[0:len(datasets)]))

    for selected_method, selected_network in itertools.product(methods, networks):
        s = 'Plots/DampingVsAUC_curve_%s_%s_%s_%d.pdf' %(formatMethodName(selected_classifier), formatMethodName(selected_method), formatMethodName(selected_network), numFeat)
        pp = PdfPages(s)
        for dataset in datasets:
            keys = sorted([(d, m, feat) for (d, m, feat) in data_to_mean.keys() if d == dataset
                and (m.startswith(selected_method+'_0') or m.startswith(selected_method+'_1'))
                and m.endswith(selected_network) and d == dataset and feat == numFeat])
            dampings = []
            means = []
            stds = []
            for key in keys:
                dampings.append(key[1].split('_')[len(key[1].split('_'))-2])
                means.append(data_to_mean[key][0])
                stds.append(data_to_mean[key][1])
            plt.errorbar(range(0, len(dampings)), means, yerr = stds, fmt = '.'+dataset_color[dataset])
            p, = plt.plot(range(0, len(dampings)), means, '.-'+dataset_color[dataset])
            plots.append(p)
        plt.legend(plots, list_datasets, loc=4)
        plt.grid()
        ax = plt.axes()
        ax.set_xticks(range(0, len(dampings)))
        ax.set_xticklabels(dampings, fontsize=13)
        plt.title(formatMethodName(selected_classifier+" "+selected_method+" "+formatMethodName(selected_network))+' #features: '+str(numFeat), fontsize=16)
        plt.axis([-0.5, 10.5, 0, 1])
        plt.ylabel("AUC", fontsize=15)
        plt.xlabel("Damping factor", fontsize=15)
        pp.savefig()
        pp.close()
        plt.clf()


def evalInnerLoop(databases, folds, repeat = 0):

    """
    Determines for each method the best performing number of features across the
    innerloop evaluation.
    Parameters are set such that only the results for repeatNr = 0
    and classifier = NMC_V1 are selected.

    databases:
        'EXP01Inner_Dao.sqlite3',
        'EXP01InnerLoop_backup_done_July1.sqlite3',
        'EXP01InnerLoopPartJul11.sqlite3',
        'EXP01InnerLoopPartJul15.sqlite3',
        'EXP01InnerLoopPartJul19-Part1.sqlite3'
    """

    #Determine the best performing number of features. For Winter and GeneRank
    #get also the best performing damping factor
    MethodsToNumFeatures = {}
    datasets = []
    for spec in list(itertools.product(['U133A_combat_RFS', 'U133A_combat_DMFS'], repeats, range(folds))):
        datasets.append(spec[0]+'_'+str(spec[1])+'_'+str(spec[2]))
    for d in datasets:
        MethodsToNumFeatures[(d, 'NKI')] = 41
        MethodsToNumFeatures[(d, 'Erasmus')] = 66

    dbResults = []
    for f in databases:
        print f
        dbResults.extend(getResultsInnerLoop(f))
    dbResults = list(set(dbResults))
    
    datasets = set([r[0] for r in dbResults])
    folds = set([r[1] for r in dbResults])
    featureExtractorsAndNetworks = set([(r[3], r[2]) for r in dbResults])
    classifiers = list(set([r[4] for r in dbResults]))
    classifier = classifiers[0]

    #All methods apart from Winter/GeneRank
    #1) get the data for one method, network
    #2) sort the items descending accoridng to the average auc
    #3) get number of features
    selectedF = [(f, n) for (f, n) in featureExtractorsAndNetworks if (not f.startswith("Winter") and not f.startswith("GeneRank"))]
    for dataset, fold in itertools.product(datasets, folds):
        for method, network in selectedF:
            data = [(auc_value, number_of_features) for (d, f, n, m, c, 
                number_of_features, auc_value) in dbResults if d == dataset and method == m
                and c == classifier and network == n and f == fold]
            sorted_by_first = sorted(data, key=lambda tup: tup[0], reverse=True)
            if network != None:
                methodF = formatMethodName(method+'_'+network)
            else:
                methodF = formatMethodName(method)
            dataName = dataset.split('_fold')[0]+'_'+str(repeat)+'_'+str(fold)
            MethodsToNumFeatures[(dataName, methodF)] = sorted_by_first[0][1]
    #Winter and GeneRank:
    # 1) get the data for method, damping factor and network
    #2) sort the items descending accoridng to the average auc
    #3) get number of features and dampingfactor
    networks = set([n for m, n in featureExtractorsAndNetworks if m.startswith("Winter") or m.startswith("GeneRank")])
    prefixMethods = ['GeneRankFeatureExtractor', 'GeneRankTscore', 'WinterFeatureExtractor_', 'WinterFeatureExtractor_SurvTime']
    combi = list(itertools.product(prefixMethods, networks))
    for dataset, fold in itertools.product(datasets, folds):
        for method, network in combi:
            if method == 'WinterFeatureExtractor_':
                data = [(auc_value, number_of_features, m) for (d, f, n, m, c,
                    number_of_features, auc_value) in dbResults if d == dataset and m.startswith(method)
                    and c == classifier and network == n and 'SurvTime' not in m and f == fold]
            else:
                data = [(auc_value, number_of_features, m) for (d, f, n, m, c,
                    number_of_features, auc_value) in dbResults if d == dataset and m.startswith(method)
                    and c == classifier and network == n and f == fold]
            sorted_by_first = sorted(data, key=lambda tup: tup[0], reverse=True)
            methodF = formatMethodName(sorted_by_first[0][2]+'_'+network)
            dataName = dataset.split('_fold')[0]+'_'+str(repeat)+'_'+str(fold)
            MethodsToNumFeatures[(dataName, methodF)] = sorted_by_first[0][1]
    
    return MethodsToNumFeatures

def EvaluatePerformanceOuterLoop(sqlFiles, num = None):
    """
    MethodsToNumFeat: dictionary
    Parameters are set such that only the results for repeatNr = 0
    and classifier = NMC_V1 are selected.

    num: for fixed number of parameters, eg, 50, 100 or 150; if none the inner loop results will be
         evalauted to get the best number of features
    sqlFiles:
        'EXP01_FINAL.sqlite3'
        'EXP01_ChuangKEGG_DMFS.sqlite3'
        'Experiment_01_10xWithin_U133A_Signatures_Random.sqlite3'
    """
    #
    dbResults = []
    for f in sqlFiles:
        print f
        dbResults.extend(getResults(f))
    
    dbResults = list(set(dbResults))

    datasets = set([r[0] for r in dbResults])
    featureExtractorsAndNetworks = set([(r[2], r[1]) for r in dbResults])
    classifiers = set([r[3] for r in dbResults])
    repeats = set([r[4] for r in dbResults])
    
    folds = set([r[5] for r in dbResults])
    
    #either give fixed number of features or provide dictionary
    if num != None:
        bestFeatures = False
        MethodsToNumfeat = makeDictMethodToNum(num, featureExtractorsAndNetworks)
    else:
        bestFeatures = True
        num = 'best_Features'
        #NOTE: check path of the databases!!!!
        databases = [
            'EXP01Inner_Dao.sqlite3',
            'EXP01InnerLoop_backup_done_July1.sqlite3',
            'EXP01InnerLoopPartJul11.sqlite3',
            'EXP01InnerLoopPartJul15.sqlite3',
            'EXP01InnerLoopPartJul19-Part1.sqlite3',
            'EXP01InnerLoopPart1Jul28.sqlite3',
        ]
        MethodsToNumfeat = evalInnerLoop(databases, max(folds)+1)
            
    #All methods apart from Winter and GeneRank
    # 1) get the data for method, network
    # 2) get aucs corresponding to the number of features
    # 3) Plot
    # In the code there are some additional ifs, to handle incomplete runs of the whole Experiment 1
    selectedF = [(f, n) for (f, n) in featureExtractorsAndNetworks if (not f.startswith("Winter") and not f.startswith("GeneRank"))]
    #Winter and GeneRank:
    # 1) get the data for method, damping factor and network
    # 2) save in dictionary me:damping:numfeatures:aucs
    # 3) get aucs for number of features for each damping factor
    # 4) choose best performing damping factor
    # 4) Plot
    pageRankF = [(f, n) for (f, n) in featureExtractorsAndNetworks if f.startswith("Winter") or f.startswith("GeneRank")]
    DatasetClassifierToAUCmatrix = getAUCmatrixOuter(dbResults, selectedF+pageRankF, datasets,
        MethodsToNumfeat, max(folds)+1, bestFeatures = bestFeatures)

    #plot1 [SG, Random, NKI, Erasmus]
    methods = [formatMethodName(m+'_'+n) for m, n in featureExtractorsAndNetworks if n != None]
    methods.extend([formatMethodName(m) for m, n in featureExtractorsAndNetworks if n == None])
    plot1 = ['SG', 'Random', 'NKI', 'Erasmus']
    plot2 = sorted([m for m in methods if m.startswith('L') or m.startswith('C') or m.startswith('D') or m.startswith('T')])
    plot3 = sorted([m for m in methods if (m.startswith('GR') or m.startswith('W')) and "HPRD" in m])
    plot4 = sorted([m for m in methods if (m.startswith('GR') or m.startswith('W')) and "HPRD" not in m])
    p_values = []
    for dataset, classifier in DatasetClassifierToAUCmatrix:
        matrix, methods = DatasetClassifierToAUCmatrix[(dataset, classifier)]
        matrix = np.array(matrix)
        methods = np.array(methods)
        p_values.extend([dataset+' '+str(num) + ' & & ', ' & Method & p-vals (no multiple testing correction)'])
        idxSG = np.argwhere(methods == 'SG').flatten()[0]
        up = max(matrix[idxSG, :])
        down = min(matrix[idxSG, :])
        #plot1
        idxMethods = [list(methods).index(m) for m in plot1 if m in methods]
        aucValues = matrix[idxMethods, :]
        s = 'Results/Plots/Boxplot_Experiment_01_%s_%s_plot1_%s.pdf' % (dataset, classifier.replace('BinaryNearestMeanClassifier', 'NMC'), str(num))
        pp_box = PdfPages(s)
        make_dot_boxplot(aucValues, list(methods[idxMethods]), "", "", "AUC", pp_box, up, down)
        pp_box.close()
        plt.clf()
        #plot2
        idxMethods = [list(methods).index(m) for m in plot2 if m in methods]
        aucValues = matrix[idxMethods, :]
        formatMethods = [m.replace('_', '\n') for m in list(methods[idxMethods])]
        s = 'Results/Plots/Boxplot_Experiment_01_%s_%s_plot2_%s.pdf' % (dataset, classifier.replace('BinaryNearestMeanClassifier', 'NMC'), str(num))
        pp_box = PdfPages(s)
        make_dot_boxplot(aucValues, formatMethods, dataset+" "+classifier.replace('BinaryNearestMeanClassifier', 'NMC')+' #'+str(num), 
            "", "AUC", pp_box, up, down)
        pp_box.close()
        plt.clf()
        #plot3
        idxMethods = [list(methods).index(m) for m in plot3 if m in methods]
        aucValues = matrix[idxMethods, :]
        formatMethods = [m.replace('_', '\n') for m in list(methods[idxMethods])]
        s = 'Results/Plots/Boxplot_Experiment_01_%s_%s_plot3_%s.pdf' % (dataset, classifier.replace('BinaryNearestMeanClassifier', 'NMC'), str(num))
        pp_box = PdfPages(s)
        make_dot_boxplot(aucValues, formatMethods, "", "", "AUC", pp_box, up, down)
        pp_box.close()
        plt.clf()
        #plot4
        idxMethods = [list(methods).index(m) for m in plot4 if m in methods]
        aucValues = matrix[idxMethods, :]
        formatMethods = [m.replace('_', '\n') for m in list(methods[idxMethods])]
        s = 'Results/Plots/S_Boxplot_Experiment_01_%s_%s_plot4_%s.pdf' % (dataset, classifier.replace('BinaryNearestMeanClassifier', 'NMC'), str(num))
        pp_box = PdfPages(s)
        make_dot_boxplot(aucValues, formatMethods, "", "", "AUC", pp_box, up, down)
        pp_box.close()
        plt.clf()
        idx = np.argsort(methods)
        for i in idx:
            if methods[i] == 'SG':
                continue
            p_values.append( ' & '+methods[i]+' & '+str("%0.4f" % (Wilcoxon(matrix[idxSG], matrix[i, :], True))))

    with open('Results/pValsWilcoxon_EXP01_'+str(num)+'.txt','w') as f:
        print >> f, "\\\\\n".join(p_values)
  
if __name__ == "__main__":
    print "Plotting functions to analyse the pclassifiers' performance (no sec. data shuffling):"
    print """Functions: 
             EvaluatePerformanceOuterLoop: makes theplots for the outer cross validation loop
                sqlFile, calculates also the paired t-statistic and saves it as a txt.
                num - set to the required number of features if a fixed number of features should be tested
                      if num == None evalInner is called to determine the best number of features.
             evalInnerLoop: finds the best performing number of features for each method
                database - path to sql database
                returns dictionary

             For the Following two functions it is advisable to download the ACES_results, since compiling the
             the necessary information takes looooooooooooooooooooooooooong.
             AucVsFeaturesCurve: makes the AUC versus number of features plots in the supplement
                dbResults: list of tuples as received by getResults
                returns a dictionry
             AucVsDampingFactor: makes the AUC versus damping factor plots for the Winter and GeneRank method
                data_to_mean - dictionary as returned by AucVsFeaturesCurve
                numFeat - fixed number of features to be tested, default 50
             
    """        

