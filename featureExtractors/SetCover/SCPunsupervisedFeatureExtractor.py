import numpy as np
import scipy
from scipy import stats
import json
#from GeneSetCollection import GeneSetCollection
from SetCoverExperiments import WholeDatasetGeneCoverFromExpressionData
#from datatypes.GeneSetCollection import GeneSetCollection
import tempfile
import shutil

class SCPunsupervisedFeatureExtractorFactory(object):
    productName = "SCPunsupervisedFeatureExtractor"

    def train(self, dataset, K = 1, mul = 1):
        if K > 1:
            self.productName = str(K)+"_SCPunsupervisedFeatureExtractor"
        print self.productName
        #create bipartite graph, do not take classlabels into account
        #setup SCP instance
        #solve SCP
        #retrieve solution genes
        #get upper and lower threshold for each solution gene

        #creeate a tmp folder
        tempdir = tempfile.mkdtemp()
        experimentName = "tmp_SCP"
        prefixGenes = "Entrez_"
        mul = mul

        solutionGenes = WholeDatasetGeneCoverFromExpressionData([dataset], prefixGenes, tempdir, 
            experimentName, useCostFunction = "DeregulationAvg", mul = mul, K = K, timelimit = 1200)

        shutil.rmtree(tempdir)

        # --> features = dict(gene) = [upper, lower]
        solutionGenes = solutionGenes[dataset.name]
        sol, GRAPHS, COSTS, gap = solutionGenes
        print "SCP Solution:", dataset.name, len(sol), gap
        features = {}
    
        for gene in sol:
            idx = np.argwhere(dataset.geneLabels == gene)
            mean = dataset.expressionData[:, idx].mean()
            std = dataset.expressionData[:, idx].std()
            upper = mean + mul*std
            lower = mean - mul*std
            #upper, lower: cut-off values in the tail of the distribution.
            features[gene] = [upper, lower]
    
        return SCPunsupervisedFeatureExtractor(dataset.geneLabels, features)

class SCPunsupervisedFeatureExtractor(object):
    name               = "SCPunsupervisedFeatureExtractor"

    def __init__(self, geneLabels, features):
        self.geneLabels         = geneLabels
        self.features           = features
        self.validFeatureCounts = len(self.features)

    def extract(self, dataset):

        assert all(dataset.geneLabels == self.geneLabels)
    
        #Find which feature genes are deregulated for the samples
        #The derugulation is defined by the upper and lower thresholds
        idxGenes = np.argwhere(np.in1d(self.geneLabels, self.features.keys())).flatten()
        dereg = dataset.extractGenesByIndices(dataset.name+"dereg", idxGenes)
        #for gene in self.features:
        #    geneIdx = np.argwhere(dereg.geneLabels == gene).flatten()[0]
        #    upper, lower = self.features[gene]
        #    for pat in range(dereg.numPatients):
        #        if dereg.expressionData[pat, geneIdx] > upper or dereg.expressionData[pat, geneIdx] < lower:
        #            dereg.expressionData[pat, geneIdx] = 1
        #        else:
        #            dereg.expressionData[pat, geneIdx] = 0

        #return binary matrix patients x genes
        return dereg.expressionData
    
    def toJsonExpression(self):

        return json.dumps((self.__class__.__name__, [geneLabel for geneLabel in self.geneLabels], [sorted(feature) for feature in self.features.keys()], 
            [(feature, self.features[feature][0], self.features[feature][1]) for feature in self.features.keys()]))

