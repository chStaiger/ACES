import numpy as np
import scipy
from scipy import stats
import json
#from GeneSetCollection import GeneSetCollection
from SetCoverExperiments import GeneCoverFromExpressionData
#from datatypes.GeneSetCollection import GeneSetCollection
import tempfile
import shutil

class SCPFeatureExtractorFactory(object):
    productName = "SCPFeatureExtractor"

    def train(self, dataset, K = 1):
        if K > 1:
            self.productName = str(K)+"_SCPFeatureExtractor"
        print self.productName
        #create bipartite graph
        #setup SCP instance
        #solce SCP
        #retrieve solution genes
        #get upper and lower threshold for each solution gene

        #creeate a tmp folder
        tempdir = tempfile.mkdtemp()
        experimentName = "tmp_SCP"
        prefixGenes = "Entrez_"
        pVal = 0.01

        solutionGenes = GeneCoverFromExpressionData([dataset], prefixGenes, tempdir, experimentName, 
            useCostFunction = "DeregulationAvg", pVal = pVal, K = K, timelimit = 300)
        shutil.rmtree(tempdir)

        # --> features = dict(gene) = [upper, lower]
        solutionGenes = solutionGenes[dataset.name]
        sol, GRAPHS, COSTS, gap = solutionGenes
        print "SCP Solution:", dataset.name, len(sol), gap
        features = {}
    
        #divide into two patient groups
        GoodP = dataset.extractPatientsByIndices("Good", dataset.patientClassLabels==False, False, False)
        PoorP = dataset.extractPatientsByIndices("Poor", dataset.patientClassLabels==True, False, False)

        for gene in sol:
            geneIdx = np.argwhere(GoodP.geneLabels == gene).flatten()[0]
            dist = scipy.stats.norm.fit(GoodP.expressionData[:, geneIdx])
            #upper, lower: cut-off values in the tail of the distribution.
            upper = scipy.stats.norm.ppf(1-pVal/2, dist[0], dist[1])
            lower = scipy.stats.norm.ppf(pVal/2, dist[0], dist[1])
            features[gene] = [upper, lower]
    
        return SCPFeatureExtractor(dataset.geneLabels, features)

class SCPFeatureExtractor(object):
    name               = "SCPFeatureExtractor"

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

