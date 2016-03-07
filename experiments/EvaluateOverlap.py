# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013
import numpy as np
import sqlite3
import json
import h5py
import itertools

from datatypes.ExpressionDataset import HDF5GroupToExpressionDataset
from datatypes.GeneSetCollection import ReadGeneSetCollection
from datatypes.EdgeSet import ReadSIF
from statistics.Fisher import FisherExactG # own implementation
from statistics.Fisher import LogarithmOfFraction

from rpy2 import *
import rpy2.robjects as robjects

import matplotlib
matplotlib.use("Agg") # make sure we can use matplotlib without X
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from matplotlib.patches import Polygon

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

def getFeatures(databases):
    print "Takes long when database is big."
    dbResults = []
    for database in databases:
        print database
        dbConnection = sqlite3.connect(database)
        dbCursor = dbConnection.cursor()
        #NOTE: for the shuffling experiment there is another variable to extract --> shuffle
        dbCursor.execute("SELECT dataset, foldNr, feature_extractor, network_dataset, feature_definition FROM FeatureExtractors \
            WHERE repeatNr == 0;")
        dbResults.extend(dbCursor.fetchall())
        dbConnection.close()
    dbResults = list(set(dbResults))
    #Sort the features
    #Make a dictionary Dataset:featureExtractorNetwork:rankedFeatures
    datasetSpec = list(set([(dataset, fold) for (dataset, fold, _, _, _) in dbResults]))
    methodsAndNetworks = list(set([(method, network) for (_, _, method, network, _) in dbResults]))
    MethodToDatasetToFeatures = {}
    for m, n in methodsAndNetworks:
        if n != None:
            methodName = formatMethodName(m+"_"+n)
        else:
            methodName = formatMethodName(m)
        if methodName not in ['Erasmus', 'NKI']:
            dataToFeatures = {}
            for dataset, fold in datasetSpec:
                featureDefGenerator = (featuredef for (d, foldNr, method, network, featuredef) in dbResults
                    if dataset == d  and fold == foldNr and method == m and network == n)
                try:
                    featuredef = next(featureDefGenerator)
                    dataName = dataset+'_'+str(0)+'_'+str(fold)
                    (fdef_name, fdef_lookup, fdef_indices) = json.loads(featuredef)
                    dataToFeatures[dataName] = fdef_indices
                except:
                    print "No features for:", dataset, fold, m, n
            MethodToDatasetToFeatures[methodName] = dataToFeatures
    return MethodToDatasetToFeatures

def getUniverses(hdf5File, secData, secDataNames):

    "SET UP UNIVERSE DictPWtoGenes"
    
    DictPWtoGenes = {} # PWname:Genelist

    #get array genes
    f = h5py.File(hdf5File)
    
    dataset = [HDF5GroupToExpressionDataset(f[group]) for group in f.keys()][0]
    f.close()
    DictPWtoGenes['SG'] = list(dataset.geneLabels)
   
    #get secData genes
    for i in range(len(secData)):
        try:
            nwEdges             = ReadSIF(secDataNames[i], secData[i] , "Entrez_")
            DictPWtoGenes[formatMethodName(secDataNames[i])]   = set(nwEdges.getNodes()).intersection(dataset.geneLabels)
        except:
            try:
                nwGeneSets        = ReadGeneSetCollection(secDataNames[i] , secData[i], "Entrez_")
                DictPWtoGenes[formatMethodName(secDataNames[i])] = set(nwGeneSets.getNodes()).intersection(dataset.geneLabels)
            except:
                print "File not found or invalid format:", secData[i]
                raise
    
    print 'Dict keys secData', DictPWtoGenes.keys()
    return DictPWtoGenes

def prefixDatasets(datasets):
    """
    This method is especially designed for the U133A cohorts.
    """
    return set(['_'.join(d.split('_')[:3]) for d in datasets])

def fisherOverlap(DatasetMethodNumber, MethodToDatasetToFeatures, DatasetToFeatures=None):
    if len(DatasetMethodNumber) == 5:
        if DatasetToFeatures == None:
            DatasetToFeatures = {}
        fisherG = []
        for dmn1, dmn2 in itertools.combinations(DatasetMethodNumber, 2):
            fdef_indices_kbest1 = MethodToDatasetToFeatures[dmn1[1]][dmn1[0]][:dmn1[2]]
            fdef_indices_kbest2 = MethodToDatasetToFeatures[dmn2[1]][dmn2[0]][:dmn2[2]]
            if dmn1[1] in ['Random'] or dmn1[1].startswith('GR') or dmn1[1].startswith('W'):
                geneset1 = frozenset(fdef_indices_kbest1)
                geneset2 = frozenset(fdef_indices_kbest2)
            elif dmn1[1].startswith('D') or dmn1[1].startswith('C') or dmn1[1].startswith('L'):
                geneset1 = frozenset.union(*[frozenset(f) for f in fdef_indices_kbest1])
                geneset2 = frozenset.union(*[frozenset(f) for f in fdef_indices_kbest2])
            elif dmn1[1].startswith('T'):
                geneset1 = []
                for gSet in fdef_indices_kbest1:
                    geneset1.extend([gSet[0]]+gSet[1])
                geneset2 = []
                for gSet in fdef_indices_kbest2:
                    geneset2.extend([gSet[0]]+gSet[1])
                geneset1 = frozenset(geneset1)
                geneset2 = frozenset(geneset2)
            elif dmn1[1] == 'SG':
                fdef_indices_kbest1 = MethodToDatasetToFeatures[dmn1[1]][dmn1[0]][:DatasetToFeatures[dmn1[0]]]
                fdef_indices_kbest2 = MethodToDatasetToFeatures[dmn2[1]][dmn2[0]][:DatasetToFeatures[dmn2[0]]]
                geneset1 = frozenset(fdef_indices_kbest1)
                geneset2 = frozenset(fdef_indices_kbest2)
            else:
                print "Method not defined", dmn1, dmn2
                raise Error
            #get Universe
            if dmn1[1] in ["SG", 'Random']:
                universe = DictUniverses['SG']
            elif dmn1[1].startswith('GR') or dmn1[1].startswith('W'):
                network = dmn1[1].split('_')[len(dmn1[1].split('_'))-1]
                universe = DictUniverses[network]
            else:
                #get network
                network = dmn1[1].split('_')[1]
                universe = DictUniverses[network]
            AandB = len(set(geneset1).intersection(geneset2))
            AandNotB = len(set(geneset1).difference(geneset2))
            BandNotA = len(set(geneset2).difference(geneset1))
            U = len(universe) - (AandB + AandNotB + BandNotA)
            fisherG.extend([-1*LogarithmOfFraction(FisherExactG(AandB, AandNotB, BandNotA, U))])
            DatasetToFeatures[dmn1[0]] = len(geneset1)
            DatasetToFeatures[dmn2[0]] = len(geneset2)
        return fisherG, DatasetToFeatures
    return None

def boxplotPairedDistribution(DATA, methods, numFeat, legend, pp, title):
    fig = plt.figure(figsize=(15, 6)) # numbers in inch
    #fig.set_size_inches(6.83, 5) #Plos format
    ax1 = fig.add_subplot(111)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15) # adjusting output figure
    bp = plt.boxplot(DATA.T)
    plt.setp(bp['boxes'], color='black')
    ax1.set_ylabel("Fisher exact test, -log10(p-value)", fontsize=15)
    numDists = DATA.shape[0]/2
    # Now fill the boxes with desired colors
    boxColors = ['orange','royalblue']
    numBoxes = numDists*2
    medians = range(numBoxes)

    for i in range(numBoxes):
        box = bp['boxes'][i]
        boxX = []
        boxY = []
        for j in range(0,5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = zip(boxX,boxY)
        k = i % 2
        boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
        ax1.add_patch(boxPolygon)
        # Now draw the median lines back over what we just filled in
        med = bp['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, 'k')
            medians[i] = medianY[0]
        #plt.plot([np.average(med.get_xdata())], markeredgecolor='k')
        plt.plot([np.average(med.get_xdata())], [np.average(DATA[i])],
            color='w', marker='*', markeredgecolor='k')
    plt.vlines([i+0.5 for i in range(2, numBoxes, 2)], -100, 5000, color='darkgrey', linestyles='dashed')
    ax = plt.axes()
    xtickNames = plt.setp(ax1, xticklabels = methods)
    plt.setp(xtickNames, rotation=90, fontsize=9)
    ax1.set_xticks(range(1, len(DATA), 2))
    if title == "Jaccard index":
        ax1.set_ylim(0, 0.3)
    else:
        ax1.set_ylim(DATA.min()-0.5, DATA.max()+10)
    ax1.set_xlim(0.5, numBoxes+0.5)
    ax1.yaxis.grid(True, linestyle='--', which='major', color='darkgrey',alpha=0.5, linewidth = 2)
    #ax1.xaxis.grid(True, linestyle='--', which='major', color='darkgrey',alpha=0.5)
    if numFeat != 'best_Features':
        s = 'Best '+str(numFeat)+'\n composite features'
    else:
        s = 'Best \n composite features'
    plt.figtext(0.7, 0.85,  s,
        backgroundcolor=boxColors[0], color='black')
    if legend.find("\n") == -1:
        plt.figtext(0.7, 0.79, legend , backgroundcolor=boxColors[1],
            color='white')
    else:
        plt.figtext(0.7, 0.715, legend , backgroundcolor=boxColors[1],
            color='white')
    fig = plt.gcf()
    fig.tight_layout()
    pp.savefig()


def makeDictMethodToNum(num, methods, datasets):
    """
    Initialises a dictionary with methods to number iof features
    with one specific number of features for all methods.
    
    """
    m = [formatMethodName(item[0]+'_'+item[1]) for item in featureExtractorsAndNetworks if item[1] != None]
    m.extend([formatMethodName(item[0]) for item in featureExtractorsAndNetworks if item[1] == None])
    MethodsToNumfeat = dict(zip(list(itertools.product(datasets, m)), len(m)*[num]))
    MethodsToNumfeat['NKI'] = 41
    MethodsToNumfeat['Erasmus'] = 66

    return MethodsToNumfeat



def FisherCntrlSize(MethodToNum, databases, hdf5File, secData, secDataNames, nrFolds = 5, repeats = [0]):
    """
    Calculates the overlap as Fisher Exact test from the Cross validation results.
    Since network markers are calculated more genes than features we 
    need to correct for this when examining the Single Genes and the Random gene markers.
    NOTE that Gene signatures as NKI (Van't veer) and Erasmus (Wang)
    are not altered when using a different training set, hence calculating an overlap
    is not useful.
    We assume that all datasets were measured on the same genes.

    
    Input
    MethodToNum:    The required number of features for each method, given as a dictionary.
                    MethodToNum,_ = pickle.load(open('EXP01_bestFeatures.pickle'))
    database:       filename for the sqlite database, e.g. 
                    database = "GridResults/EXP01_backup_done_June18.sqlite3"
                    databases = [
                    'EXP01_FINAL.sqlite3',
                    'EXP01_ChuangKEGG_DMFS.sqlite3',
                    'Experiment_01_10xWithin_U133A_Signatures_Random.sqlite3',
                    ]

    hdf5File:    Path to the dataset .h5 file (needed to extract the number of total genes in the Fisher test)
    secData:     list of paths to the secondary data sources employed in the features
    secDataNames: Names of the secData, same order as secData

    Output
    """
    
    MethodToDatasetToFeatures = getFeatures(databases)
    if secData == None:
        secData = ["experiments/data/KEGG_edges1210.sif", "experiments/data/HPRD9.sif", "experiments/data/I2D_edges_0411.sif", 
                        "experiments/data/ipp.sif", "experiments/data/C2V3_PathwayGeneSets_Entrez.txt"]
        secDataNames = ["KEGG", "HPRD9", "I2D", "IPP", "MsigDB"]

    DictUniverses = getUniverses(hdf5File, secData, secDataNames)
    datasets = MethodToDatasetToFeatures['SG'].keys()
    prefixDS = prefixDatasets(datasets)
    
    DictFisher = {}
    for prefix in prefixDS:
        validMethods = []
        overlap = []
        methods = sorted([m for (d, m) in MethodToNum if d.startswith(prefix)])
        visited = []
        for key in methods:
            print key
            if len(key.split('_')) > 1:
                methodName = key.split('_')[0]+'_'
                network = key.split('_')[len(key.split('_'))-1]
            else:
                methodName = key
                network = ''
            if (prefix, methodName, network) in visited or key in ['SG', 'NKI', 'Erasmus']:
                continue
            visited.append((prefix, methodName, network))
            DatasetMethodNumber = [(d, m, MethodToNum[d, m]) for (d, m) in MethodToNum if d.startswith(prefix) 
                and m.startswith(methodName) and m.endswith(network)]
            #get pairwise overlap between all features from the different datasets
            result, DatasetToFeatures = fisherOverlap(DatasetMethodNumber, MethodToDatasetToFeatures)
            if result == None:
                print "Not enough features for:", key, prefix
                continue
            if network != None:
                validMethods.append(methodName+network)
            else:
                validMethods.append(methodName)
            overlap.append(result)
            DatasetMethodNumber = [(d, m, MethodToNum[d, m]) for (d, m) in MethodToNum if d.startswith(prefix)
                and m.startswith('SG') and m.endswith('')]
            resultSG, _ = fisherOverlap(DatasetMethodNumber, MethodToDatasetToFeatures, DatasetToFeatures)
            validMethods.append('SG')
            overlap.append(resultSG)
        DictFisher[prefix] = (validMethods, overlap)

    return DictFisher

def makePlot(DictFisher, dataset, num):

    #dataset, classifier = DatasetToFeaturesToMethods.keys()[1]
    #methods = DatasetToFeaturesToMethods[(dataset, 'NMC_V1')][num]
    #methods = [m.replace('Tscore', 'Tstat') for m in methods]
    #DictFisher = DictFisher_fixedFeatures[0]
    #

    DATA = np.array(DictFisher[dataset][1])
    METHODS = np.array(DictFisher[dataset][0])
    methodsIDX = [list(METHODS).index(item) for item in METHODS if
        item != 'Erasmus' and item != 'NKI' and item != 'SG']
    sgIDX = list(np.array(methodsIDX)+1)
    #alternate the two index lists
    IDX = [None]*(len(methodsIDX)+len(sgIDX))
    IDX[::2] = methodsIDX
    IDX[1::2] = sgIDX

    #There are numerically zeros
    plotData = DATA[IDX]

    plotMethods = METHODS
    plotMethods = [m.replace('_', '\n') for m in plotMethods]

    s = 'Results/Plots/Fisher_'+dataset+'_'+str(num)+'.pdf'
    pp = pp_box = PdfPages(s)
    boxplotPairedDistribution(plotData, plotMethods, num, "Control for size \n single genes", pp, "Fisher Exact "+dataset)
    pp.close()
    plt.clf()


if __name__ == "__main__":
    print """ We recommend to use ipython for all the evaluation procedures.
              Calculate the overlap:
              start ipython
              from experiments.EvaluateOverlap import *
              Get all the necessary parameters from  the performance analysis ...
              MethodToNum - best performing number of features for each method
              MethodToNum = dict(zip(methods, num))
              or
              MethodToNum = dict(zip(methods, len(methods)*[50])) for the best 50 features
              methods - should only contain methods that you are interested in. 
                Eg for Winter and GeneRank we are just interested in the best performing 
                damping factor. That reduces the running time when sorting the ~1000000 items
              DictFisher = FisherCntrlSize(MethodToNum, database, methods, hdf5File, secData, secDataNames)
                conatins the methods for which enough features were present and a list of lists containing the pvalues
              DictFisher['dataset_prefix'] = (validMethods, list)
              convert the list to np.array and plot with xticklabels = validMethods
          """

