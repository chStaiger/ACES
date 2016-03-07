#! /usr/bin/env python

import h5py, sys

# Input of primary and secondary datasets

import numpy as np
import networkx as nx
from networkx.algorithms import bipartite

import networkx as nx
import pylab
import pickle
from makeBiG import BipartiteGraphTwoSided, SampleGenePWDataGraph, SampleGeneInteractionGraph, SampleSecDataGraph, BipartiteGraphSTD
from SetCoverFunctions import SetCover, KCover_secData, KLCover_secData_connectedComponent

#defines costs for each gene in abipartite graph
#the smaller the costs the more important the gene
def geneCostFunction_Deregulation(B, prefixGenes):
    """
    Defines a cost function c:Genes \union Metagenes --> R+ based on the deregulation given in the graph B.

    B: bipartite graph with the two sets of nodes: Samples S, Genes G 
       Edges between samples and genes must carry an attribute 'deregulation' 
       (check edgeAnnotation = G[G.edges()[1][0]][G.edges()[1][1]])
    
    For each gene the 'deregulation' on the edges to the samples are summed 
        w'(g) = \sum_{s \in S} deregulation(s, g).
    The final weight for a gene g is then w(g) = 1/w'(g). If g has no neighbouring samples, its weight is set to 1.
    The costs of the metagenes is set to min_g(w(g)).
    """
    # get the two partitions
    Genes = [node for node in B.nodes() if node.startswith(prefixGenes)]
    Samples = [node for node in B.nodes() if not node.startswith(prefixGenes)]
    
    costs = {}
    for Gene in Genes:
        #get neighbours
        neighbours = nx.neighbors(B, Gene)
        neighbourMetagenes = []
        for n in neighbours:
            if n.startswith('META_'):
                neighbourMetagenes.append(n)
        neighbourSamples = set(neighbours).difference(neighbourMetagenes)

        c = 0
        for n in neighbourSamples:
            c = c + abs(B[n][Gene]['deregulation'])
        costs[Gene] = c
   
    #so far: the higher the costs, the further away gene's expression from the mean for the poor outcome patients.
    #NOTE: that some genes might score low because they are not deregulated in many poor patients.
    #Now we have to invert the score: the lower the costs, the further away gene's expression from the mean for the poor outcome patients.

    for Gene in costs:   
        if costs[Gene] > 0:
            costs[Gene] = 1./costs[Gene]
        else:
            costs[Gene] = 1

    #set the costs of the Netagenes to the costs of the best gene (only for Gene-Metagene graph) 
    if prefixGenes != "META_":
        minVal = min(costs.values())
        for node in B.nodes():
            if node.startswith("META_"):
                costs[node] = minVal

    assert max(costs.values()) <= 1

    return costs 

#defines costs for each gene in abipartite graph
#the smaller the costs the more important the gene
def geneCostFunction_DeregulationAvg(B, prefixGenes):
    """
    Defines a cost function c:Genes --> R+ based on the deregulation given in the graph B.
    
    B: bipartite graph with the two sets of nodes: Samples S, Genes G 
       Edges between samples and genes must carry an attribute 'deregulation' 
       (check edgeAnnotation = G[G.edges()[1][0]][G.edges()[1][1]])

    For each gene the 'deregulation' on the edges to the samples are summed and 
        w'(g) = 1/|S|*\sum_{s \in S} deregulation(s, g).
    The final weight for a gene g is then w(g) = 1/w'(g). If the gene is not connected to any patient the cost will be set to 1.
    The costs of the metagenes is set to min_g(w(g)).
    """
    # get the two partitions
    Genes = [node for node in B.nodes() if node.startswith(prefixGenes)]
    Samples = [node for node in B.nodes() if not node.startswith(prefixGenes)]

    costs = {}
    for Gene in Genes:
        #get neighbours
        neighbours = nx.neighbors(B, Gene)
        neighbourMetagenes = []
        for n in neighbours:
            if n.startswith('META_'):
                neighbourMetagenes.append(n)
        neighbourSamples = set(neighbours).difference(neighbourMetagenes)
        c = 0
        for n in neighbourSamples:
            c = c + abs(B[n][Gene]['deregulation']) #because of the construction of the graph
                                                    #all zscores between samples and genes are greater 1 (see pvalue in makeBig)
        if len(neighbourSamples) > 0: #there are genes in a triG without sample neighbours
            costs[Gene] = c/len(neighbourSamples)
        else:
            costs[Gene] = -1
    #so far: the higher the costs, the further away gene's expression from the mean for the poor outcome patients.
    #NOTE: that some genes might score low because they are not deregulated in many poor patients. That's why we took the average deregulation
    #Now we have to invert the score: the lower the costs, the further away gene's expression from the mean for the poor outcome patients.

    for Gene in costs:
        if costs[Gene] == -1:
            costs[Gene] = 1
        else:
            costs[Gene] = 1./costs[Gene]
        B.node[Gene]['score'] = costs[Gene] 
   
    #set the costs of the Netagenes to the costs of the best gene (only for Gene-Metagene graph) 
    if prefixGenes != "META_":
        minVal = min(costs.values())
        for node in B.nodes():
            if node.startswith("META_"):
                costs[node] = minVal
    
    assert max(costs.values()) <= 1

    return costs

def metageneCostFunction_sumGenes(B):
    """
    Defines the cost function c:Metagenes --> R+ 'normalised number of deregulated genes per sample' 
    based on the deregulation given a bipartite graph B = ((Samples, Metagenes), edges)
    
    B: bipartite graph with the two sets of nodes: Samples S, Metagenes M 
       Edges between samples and metagenes must carry an attribute 'cumDereg' 
       (check edgeAnnotation = B[B.edges()[1][0]][B.edges()[1][1]]) and the metagenes must have
       an attribute "size" (check e.g. B.node['META_Amoebiasis']['size'] for KEGG)

    For each metagene m the average number of deregulated genes per patient is calculated
       w'(m) = 1/(#neighbours(m)*#genes_in_m)*\sum_{n \in neighbours(m)} #deregulated_genes
    The final weight for a metagene m is then w(m) = 1/w'(m).
    """

    # get the two partitions
    Metagenes = [node for node in B.nodes() if node.startswith("META_")]
    Samples = [node for node in B.nodes() if not node.startswith("META_")]

    costs = {}
    for m in Metagenes:
        #get neighbours
        neighbours = B.neighbors(m)
        c = float(0)
        for n in neighbours:
            c = c + len(B[m][n]['cumDereg'])
        costs[m] = c/(len(neighbours)*B.node[m]['size'])

    for m in costs:
        costs[m] = 1/costs[m]
        B.node[m]['score'] = costs[m]

    return costs

def WholeDatasetGeneCoverFromExpressionData(datasets, prefixGenes, folder, experimentName, useCostFunction = None, mul = 1, K = 1, timelimit = None):
    """
    datasets:       list of expression datasets
    folder:         folder where results will be saved
    experimentName: name of the savefile(s)
    prefixGenes:    Since in the graph we need to distinguish between patient nodes and gene nodes
                    genes carry a prefix. For the datasets provided use "Entrez".
    K:              Coverage parameter.
    mul:            thresholdUp = mean(gene)+mul*std(gene) and thresholdDown = mean(gene)-mul*std(gene)
    timelimit:      timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.
    useCostfunction:one can attach costs to the genes. 
                    Possible functions are: "Deregulation" and "DeregulationAvg" for a Gene-Sample graph
 
    For each gene an expression profile (probability function on the expression values)
    is learned from all samples.
    Genes are defined "deregulated" for a patient if their expression values lie outside of the standard deviation.
    Which gene is deregulated for which poor outcome patient is
    stored as a bipartite graph with patients being one set of nodes and genes the second set of nodes.
    """

    print "CONSTRUCTING bipartite graphs"
    GRAPHS = {} #datasetName:BipartiteGraph
    COSTS = {}
    for dataset in datasets:
        print "STATUS", dataset.name
        GRAPHS[dataset.name] = BipartiteGraphSTD(dataset, mul = mul)
        #Determine costs
        if useCostFunction == "Deregulation":
            COSTS[dataset.name] = geneCostFunction_Deregulation(GRAPHS[dataset.name], prefixGenes)
        elif useCostFunction == "DeregulationAvg":
            COSTS[dataset.name] = geneCostFunction_DeregulationAvg(GRAPHS[dataset.name], prefixGenes)
        else:
            COSTS[dataset.name] = None

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print K, "- Set Cover for ", name
        sol, gap = SetCover(GRAPHS[name], prefixGenes, timelimit, COSTS[name], K)
        solutionGenes[name] = (sol, GRAPHS[name], COSTS[name], gap)

    #pickle dump solutions
    filename = folder+"/GeneCover_"+str(K)+"_"+experimentName+".pickle"
    print "Solution genes and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

def GeneCoverFromExpressionData(datasets, prefixGenes, folder, experimentName, useCostFunction = None, pVal = 0.01, K = 1, timelimit = None):
    """
    datasets:       list of expression datasets
    folder:         folder where results will be saved
    experimentName: name of the savefile(s)
    prefixGenes:    Since in the graph we need to distinguish between patient nodes and gene nodes
                    genes carry a prefix. For the datasets provided use "Entrez".
    timelimit:      timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.
    useCostfunction:one can attach costs to the genes. 
                    Possible functions are: "Deregulation" and "DeregulationAvg" for a Gene-Sample graph
    pVal:           determines the cut off for the upper and lower tail. 
 
    For each gene an expression profile (probability function on the expression values)
    is learned from the good outcome samples (classlabel = False).
    Genes are defined "deregulated" for a poor outcome patient if their expression values lie in the
    upper or lower tail of the profile. The cut-off value is determined by 1/2*pVal for each tail. At a pVal of 1%
    a poor outcome patient must show at least an as extreme low (high) expression than 0.5% percent of the good outcome patients. 
    Only then the gene is considered to be deregulated.
    Which gene is deregulated for which poor outcome patient is
    stored as a bipartite graph with patients being one set of nodes and genes the second set of nodes.
    """
    # determine bipartite graph from data, pVal determines the probability cut-off for the upper 
    # and lower tail.
    #pVal = 0.01
    print "CONSTRUCTING bipartite graphs with p-value = ", pVal
    GRAPHS = {} #datasetName:BipartiteGraph
    COSTS = {}
    for dataset in datasets:
        print "STATUS", dataset.name
        GRAPHS[dataset.name] = BipartiteGraphTwoSided(dataset, pVal)
        #Determine costs
        if useCostFunction == "Deregulation":
            COSTS[dataset.name] = geneCostFunction_Deregulation(GRAPHS[dataset.name], prefixGenes)
        elif useCostFunction == "DeregulationAvg":
            COSTS[dataset.name] = geneCostFunction_DeregulationAvg(GRAPHS[dataset.name], prefixGenes)
        else:
            COSTS[dataset.name] = None

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print K, "- Set Cover for ", name, "and p-value", pVal 
        sol, gap = SetCover(GRAPHS[name], prefixGenes, timelimit, COSTS[name], K)
        solutionGenes[name] = (sol, GRAPHS[name], COSTS[name], gap)

    print
    print "Solution:", solutionGenes
    print

    #pickle dump solutions
    filename = folder+"/GeneCover_"+str(K)+"_pValue_"+str(pVal).replace('.', '')+"_"+experimentName+".pickle"
    print "Solution genes and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

def MetageneCoverFromExpressionData(datasets, secData, prefixGenes, folder, experimentName, 
    numDeregulatedGenes = 1, useCostFunction = None, pVal = 0.01, timelimit = None, K = 1):
    """
    datasets:           list of expression datasets
    secData:            Some secondary data of type EdgeSet or GeneSetCollection. Either genes are directly linked to 
                        metagenes (GeneSetCollection) or pairwise relations on the genes are given (EdgeSet).
    folder:             folder where results will be saved
    experimentName:     name of the savefile(s)
    numDeregulatedGenes:number of deregulated genes in a pathway for one patient. 
    prefixGenes:        Since in the graph we need to distinguish between patient nodes and gene nodes
                        genes carry a prefix. For the datasets provided use "Entrez".
    timelimit:          timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.
    useCostfunction:    one can attach costs to the metagenes. 
                        Possible options: "sumGenes" for a Metagene-Sample graph.
    pVal:               determines the cut off for the upper and lower tail. 

    returns
    
    solutionGenes:      a dictionary [dataset.name] = (solution metagenes, bipartite graph, costs, gap) 

    Constructs the minimum metagene cover on the poor outcome samples.
    Firts a bipartite graph B = ((poor outcome patients, metagenes), edges) is constructed as follows:
    For each gene an expression profile (probability function on the expression values)
    is learned from the good outcome samples (classlabel = False).
    Genes are defined "deregulated" for a poor outcome patient if their expression values lie in the
    upper or lower tail of the profile. The cut-off value is determined by 1/2*pVal for each tail. At a pVal of 1%
    a poor outcome patient must show at least an as extreme low (high) expression than 0.5% percent of the good outcome patients. 
    Only then the gene is considered to be deregulated.
    Poor outcome patients are connected to metagenes if a certain number of genes in the metagene is deregulated.
    In case of pairwise relations on the genes a Metagene is introduced for each defined pair.
    The graph B is used to construct a cplex instance which is then solved. 
    """
    # determine bipartite graph from data, pVal determines the probability cut-off for the upper 
    # and lower tail.
    #pVal = 0.01
    print "CONSTRUCTING bipartite graphs with p-value = ", pVal
    GRAPHS = {} #datasetName:BipartiteGraph
    COSTS = {}
    for dataset in datasets:
        print "STATUS", dataset.name
        B = BipartiteGraphTwoSided(dataset, pVal)
        TriG = None
        if hasattr(secData, 'geneSets'):
            TriG = SampleGenePWDataGraph(B, secData, prefixGenes)
        elif hasattr(secData, 'edges'):
            TriG = SampleGeneInteractionGraph(B, secData, prefixGenes)
        else:
            print "ERROR secondary data not of type EdgeSet or GeneSetCollection."
            assert False
        GRAPHS[dataset.name] = SampleSecDataGraph(TriG, prefixGenes, numDeregulatedGenes)
        #Determine costs
        if useCostFunction == "sumGenes":
            COSTS[dataset.name] = metageneCostFunction_sumGenes(GRAPHS[dataset.name])
        else:
            COSTS[dataset.name] = None

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print "Set Cover for dataset", name, "and secondary data", secData.name
        sol, gap = SetCover(GRAPHS[name], "META_", timelimit, COSTS[name], K)
        solutionGenes[name] = (sol, GRAPHS[name], COSTS[name], gap)
    print
    print "Solution:", solutionGenes
    print
    #pickle dump solutions
    filename = folder+"/MetageneCover_"+secData.name+"_pValue_"+str(pVal).replace('.', '')+"_"+experimentName+".pickle"
    print "Solution metagenes and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

def KGeneMetageneCoverFromExpressionData(datasets, secData, K, prefixGenes, folder, experimentName,
    useCostFunction = None, pVal = 0.01, timelimit = None):
    """
    datasets:           list of expression datasets
    secData:            Some secondary data of type EdgeSet or GeneSetCollection. Either genes are directly linked to 
                        metagenes (GeneSetCollection) or pairwise relations on the genes are given (EdgeSet).
    K:                  minimum number of covering genes per pent
    folder:             folder where results will be saved
    experimentName:     name of the savefile(s)
    prefixGenes:        Since in the graph we need to distinguish between patient nodes and gene nodes
                        genes carry a prefix. For the datasets provided use "Entrez".
    timelimit:          timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.
    useCostfunction:    one can attach costs to the genes. 
                        Possible functions are: "Deregulation" and "DeregulationAvg" for a Gene-Sample graph
    pVal:               determines the cut off for the upper and lower tail. 

    returns
    
    solutionGenes:      a dictionary [dataset.name] = (solutiongenes, solution metagenes, solutionpatients, bipartite graph, costs, gap) 

    Constructs the minimum gene-metagene cover on the poor outcome samples.
    Firts a tripartite graph B = ((poor outcome patients, gene, metagenes), edges) is constructed as follows:
    For each gene an expression profile (probability function on the expression values)
    is learned from the good outcome samples (classlabel = False).
    Genes are defined "deregulated" for a poor outcome patient if their expression values lie in the
    upper or lower tail of the profile. The cut-off value is determined by 1/2*pVal for each tail. At a pVal of 1%
    a poor outcome patient must show at least an as extreme low (high) expression than 0.5% percent of the good outcome patients. 
    Only then the gene is considered to be deregulated.
    Genes are connected to metagenes accoring to the secondary data source.
    In case of pairwise relations on the genes a Metagene is introduced for each defined pair.
    """

    # determine tripartite graph from data, pVal determines the probability cut-off for the upper 
    # and lower tail.
    #pVal = 0.01
    print "CONSTRUCTING tripartite graphs with p-value = ", pVal
    GRAPHS = {} #datasetName:BipartiteGraph
    COSTS = {}
    for dataset in datasets:
        print "STATUS", dataset.name
        B = BipartiteGraphTwoSided(dataset, pVal)
        TriG = None
        if hasattr(secData, 'geneSets'):
            TriG = SampleGenePWDataGraph(B, secData, prefixGenes)
        elif hasattr(secData, 'edges'):
            TriG = SampleGeneInteractionGraph(B, secData, prefixGenes)
        else:
            print "ERROR secondary data not of type EdgeSet or GeneSetCollection."
            assert False
        GRAPHS[dataset.name] = TriG
        #Determine costs
        if useCostFunction == "Deregulation":
            COSTS[dataset.name] = geneCostFunction_Deregulation(GRAPHS[dataset.name], prefixGenes)
        elif useCostFunction == "DeregulationAvg":
            COSTS[dataset.name] = geneCostFunction_DeregulationAvg(GRAPHS[dataset.name], prefixGenes)
        else:
            COSTS[dataset.name] = None

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print "K Gene-Metagene Cover for dataset", name, "and secondary data", secData.name
        sol, solMetagenes, solPatients, gap = KCover_secData(GRAPHS[name], K, prefixGenes, "META_", timelimit = timelimit, costs = COSTS[name])
        solutionGenes[name] = (sol, solMetagenes, solPatients, GRAPHS[name], COSTS[name], gap)

    print
    print "Solution:", solutionGenes
    print

    #pickle dump solutions
    filename = folder+"/K_GeneMetageneCover_"+secData.name+"K_"+str(K)+"_pValue_"+str(pVal).replace('.', '')+"_"+experimentName+".pickle"
    print "Solution genes, metagenes, covered patients and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

def WholeDatasetKGeneMetageneCoverFromExpressionData(datasets, secData, K, prefixGenes, folder, experimentName,
    useCostFunction = None, mul = 3, timelimit = None):
    """
    datasets:           list of expression datasets
    secData:            Some secondary data of type EdgeSet or GeneSetCollection. Either genes are directly linked to 
                        metagenes (GeneSetCollection) or pairwise relations on the genes are given (EdgeSet).
    K:                  minimum number of covering genes per pent
    folder:             folder where results will be saved
    experimentName:     name of the savefile(s)
    prefixGenes:        Since in the graph we need to distinguish between patient nodes and gene nodes
                        genes carry a prefix. For the datasets provided use "Entrez".
    timelimit:          timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.
    useCostfunction:    one can attach costs to the genes. 
                        Possible functions are: "Deregulation" and "DeregulationAvg" for a Gene-Sample graph
    pVal:               determines the cut off for the upper and lower tail. 

    returns
    
    solutionGenes:      a dictionary [dataset.name] = (solutiongenes, solution metagenes, solutionpatients, bipartite graph, costs, gap) 

    Constructs the minimum gene-metagene cover on the poor outcome samples.
    Firts a tripartite graph B = ((poor outcome patients, gene, metagenes), edges) is constructed as follows:
    For each gene an expression profile (probability function on the expression values)
    is learned from the good outcome samples (classlabel = False).
    Genes are defined "deregulated" for a poor outcome patient if their expression values lie in the
    upper or lower tail of the profile. The cut-off value is determined by 1/2*pVal for each tail. At a pVal of 1%
    a poor outcome patient must show at least an as extreme low (high) expression than 0.5% percent of the good outcome patients. 
    Only then the gene is considered to be deregulated.
    Genes are connected to metagenes accoring to the secondary data source.
    In case of pairwise relations on the genes a Metagene is introduced for each defined pair.
    """

    # determine tripartite graph from data, pVal determines the probability cut-off for the upper 
    # and lower tail.
    #pVal = 0.01
    print "CONSTRUCTING tripartite graphs with thresholdUp/Down = mean(gene)+"+str(mul)+"*std(gene)"
    GRAPHS = {} #datasetName:BipartiteGraph
    COSTS = {}
    for dataset in datasets:
        print "STATUS", dataset.name
        B = BipartiteGraphSTD(dataset, mul = mul)
        TriG = None
        if hasattr(secData, 'geneSets'):
            TriG = SampleGenePWDataGraph(B, secData, prefixGenes)
        elif hasattr(secData, 'edges'):
            TriG = SampleGeneInteractionGraph(B, secData, prefixGenes)
        else:
            print "ERROR secondary data not of type EdgeSet or GeneSetCollection."
            assert False
        GRAPHS[dataset.name] = TriG
        #Determine costs
        if useCostFunction == "Deregulation":
            COSTS[dataset.name] = geneCostFunction_Deregulation(GRAPHS[dataset.name], prefixGenes)
        elif useCostFunction == "DeregulationAvg":
            COSTS[dataset.name] = geneCostFunction_DeregulationAvg(GRAPHS[dataset.name], prefixGenes)
        else:
            COSTS[dataset.name] = None

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print "K Gene-Metagene Cover for dataset", name, "and secondary data", secData.name
        sol, solMetagenes, solPatients, gap = KCover_secData(GRAPHS[name], K, prefixGenes, "META_", timelimit = timelimit, costs = COSTS[name])
        solutionGenes[name] = (sol, solMetagenes, solPatients, GRAPHS[name], COSTS[name], gap)

    #pickle dump solutions
    filename = folder+"/K_GeneMetageneCover_"+secData.name+"K_"+str(K)+"_mul_"+str(mul).replace('.', '')+"_"+experimentName+".pickle"
    print "Solution genes, metagenes, covered patients and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes


def KLConnectedGeneMetageneCoverFromExpressionData(datasets, secData, K, L, prefixGenes, folder, experimentName,
    useCostFunction = None, pVal = 0.01, timelimit = None):
    """
    datasets:           list of expression datasets
    secData:            Some secondary data of type EdgeSet or GeneSetCollection. Either genes are directly linked to 
                        metagenes (GeneSetCollection) or pairwise relations on the genes are given (EdgeSet).
    K:                  minimum number of covering genes per patient
    L:                  maximum number of uncovered patients
    folder:             folder where results will be saved
    experimentName:     name of the savefile(s)
    prefixGenes:        Since in the graph we need to distinguish between patient nodes and gene nodes
                        genes carry a prefix. For the datasets provided use "Entrez".
    timelimit:          timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.
    useCostfunction:    one can attach costs to the genes. 
                        Possible functions are: "Deregulation" and "DeregulationAvg" for a Gene-Sample graph
    pVal:               determines the cut off for the upper and lower tail. 

    returns
    
    solutionGenes:      a dictionary [dataset.name] = (solutiongenes, solution metagenes, solutionpatients, bipartite graph, costs, gap) 

    Constructs the minimum gene-metagene cover on the poor outcome samples such that the genes and metagenes induce a connected component.
    Firts a tripartite graph B = ((poor outcome patients, gene, metagenes), edges) is constructed as follows:
    For each gene an expression profile (probability function on the expression values)
    is learned from the good outcome samples (classlabel = False).
    Genes are defined "deregulated" for a poor outcome patient if their expression values lie in the
    upper or lower tail of the profile. The cut-off value is determined by 1/2*pVal for each tail. At a pVal of 1%
    a poor outcome patient must show at least an as extreme low (high) expression than 0.5% percent of the good outcome patients. 
    Only then the gene is considered to be deregulated.
    Genes are connected to metagenes accoring to the secondary data source.
    In case of pairwise relations on the genes a Metagene is introduced for each defined pair.
    """

    # determine tripartite graph from data, pVal determines the probability cut-off for the upper 
    # and lower tail.
    #pVal = 0.01
    print "CONSTRUCTING tripartite graphs with p-value = ", pVal
    GRAPHS = {} #datasetName:BipartiteGraph
    COSTS = {}
    for dataset in datasets:
        print "STATUS", dataset.name
        B = BipartiteGraphTwoSided(dataset, pVal)
        TriG = None
        if hasattr(secData, 'geneSets'):
            TriG = SampleGenePWDataGraph(B, secData, prefixGenes)
        elif hasattr(secData, 'edges'):
            TriG = SampleGeneInteractionGraph(B, secData, prefixGenes)
        else:
            print "ERROR secondary data not of type EdgeSet or GeneSetCollection."
            assert False
        GRAPHS[dataset.name] = TriG
        #Determine costs
        if useCostFunction == "Deregulation":
            COSTS[dataset.name] = geneCostFunction_Deregulation(GRAPHS[dataset.name], prefixGenes)
        elif useCostFunction == "DeregulationAvg":
            COSTS[dataset.name] = geneCostFunction_DeregulationAvg(GRAPHS[dataset.name], prefixGenes)
        else:
            COSTS[dataset.name] = None

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print "KL Gene-Metagene Cover for dataset", name, "and secondary data", secData.name
        sol, solMetagenes, solPatients, gap = KLCover_secData_connectedComponent(GRAPHS[name], K, L, prefixGenes, "META_", timelimit = timelimit, costs = COSTS[name])
        solutionGenes[name] = (sol, solMetagenes, solPatients, GRAPHS[name], COSTS[name], gap)

    print
    print "Solution:", solutionGenes
    print

    #pickle dump solutions
    filename = folder+"/KL_ConnectedGeneMetageneCover_"+secData.name+"K_"+str(K)+"_L_"+str(L)+"_pValue_"+str(pVal).replace('.', '')+"_"+experimentName+".pickle"
    print "Solution genes, metagenes, covered patients and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

 
def SetCoverFromMutagenisisData(datasets, prefixGenes, folder, experimentName, costs = None, timelimit = None, K = 1):
    """
    datasets: list of mutation datasets, entries in the matrix must be either 1 or 0
    folder: folder where results will be saved
    prefixGenes: Since in the graph we need to distinguish between patient nodes and gene nodes, genes carry a prefix. 
        For the datasets provided it is "Entrez".
    costs: a dictionary mapping fom genes to costs
  
    The bipartite graph is directly built from the mutation data for all poor outcome samples.
    """

    solutionGenes = {}
    for dataset in datasets:
        # The data itself already gives a bipartite graph
        B = nx.Graph()
        #select only the poor outcome samples
        for gene in dataset.geneLabels:
            if costs != None:
                B.add_node(gene, costs=costs[gene])
            else:
                B.add_node(gene)
        for sample in dataset.patientLabels[dataset.patientClassLabels == True]:
            mutDat = dataset.expressionData[dataset.patientLabels == sample, ].flatten()
            mutGenes = dataset.geneLabels[mutDat == 1]
            for gene in mutGenes:
                B.add_edge(sample, gene)
        assert sum(sum(dataset.expressionData)) == len(B.edges())

        sol = SetCover(B, prefixGenes, timelimit, costs, K)
        solutionGenes[dataset.name] = (sol, B)
    
    print "Solution genes and graphs are saved as pickle file:", folder+"/SimpleSetCover_"+str(K)+"_MutagenisisData_"+experimentName+".pickle"
    filename = folder+"/SimpleSetCover_"+str(K)+"_MutageneisisData_"+experimentName+".pickle"
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

def MetageneCoverFromMutagenisisData(datasets, secData, prefixGenes, folder, experimentName,
    numDeregulatedGenes = 1, timelimit = None, K = 1):
    """
    datasets:           list of expression datasets
    secData:            Some secondary data of type EdgeSet or GeneSetCollection. Either genes are directly linked to 
                        metagenes (GeneSetCollection) or pairwise relations on the genes are given (EdgeSet).
    folder:             folder where results will be saved
    experimentName:     name of the savefile(s)
    numDeregulatedGenes:number of deregulated genes in a pathway for one patient. 
    prefixGenes:        Since in the graph we need to distinguish between patient nodes and gene nodes
                        genes carry a prefix. For the datasets provided use "Entrez".
    timelimit:          timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.

    returns
    
    solutionGenes:      a dictionary [dataset.name] = (solution metagenes, bipartite graph, costs, gap) 

    Constructs the minimum metagene cover on the poor outcome samples.
    Firts a bipartite graph B = ((poor outcome patients, metagenes), edges) is constructed from the binary matrix.
    Poor outcome patients are connected to metagenes if a certain number of genes in the metagene is mutated.
    In case of pairwise relations on the genes a Metagene is introduced for each defined pair.
    The graph B is used to construct a cplex instance which is then solved. 
    """
    # determine bipartite graph from data, pVal determines the probability cut-off for the upper 
    # and lower tail.
    #pVal = 0.01
    print "CONSTRUCTING bipartite graphs"
    GRAPHS = {} #datasetName:BipartiteGraph
    for dataset in datasets:
        print "STATUS", dataset.name
        B = nx.Graph()
        #select only the poor outcome samples
        for sample in dataset.patientLabels[dataset.patientClassLabels == True]:
            mutDat = dataset.expressionData[dataset.patientLabels == sample, ].flatten()
            mutGenes = dataset.geneLabels[mutDat == 1]
            for gene in mutGenes:
                B.add_edge(sample, gene)
        TriG = None
        if hasattr(secData, 'geneSets'):
            TriG = SampleGenePWDataGraph(B, secData, prefixGenes)
        elif hasattr(secData, 'edges'):
            TriG = SampleGeneInteractionGraph(B, secData, prefixGenes)
        else:
            print "ERROR secondary data not of type EdgeSet or GeneSetCollection."
            assert False
        GRAPHS[dataset.name] = SampleSecDataGraph(TriG, prefixGenes, numDeregulatedGenes)

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print "Set Cover for dataset", name, "and secondary data", secData.name
        sol, gap = SetCover(GRAPHS[name], "META_", timelimit, None, K)
        solutionGenes[name] = (sol, GRAPHS[name], gap)
    print
    print "Solution:", solutionGenes
    print
    #pickle dump solutions
    filename = folder+"/MetageneCover_K_"+str(K)+"_"+experimentName+".pickle"
    print "Solution metagenes and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

def KGeneMetageneCoverFromMutagenisisData(datasets, secData, K, prefixGenes, folder, experimentName,
    costs = None, timelimit = None):
    """
    datasets:           list of expression datasets
    secData:            Some secondary data of type EdgeSet or GeneSetCollection. Either genes are directly linked to 
                        metagenes (GeneSetCollection) or pairwise relations on the genes are given (EdgeSet).
    K:                  minimum number of covering genes per patient
    folder:             folder where results will be saved
    experimentName:     name of the savefile(s)
    costs:              A dictionary mapping from Genes and Metagenes to costs
    prefixGenes:        Since in the graph we need to distinguish between patient nodes and gene nodes
                        genes carry a prefix. For the datasets provided use "Entrez".
    timelimit:          timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.

    returns
    
    solutionGenes:      a dictionary [dataset.name] = (solutiongenes, solution metagenes, solutionpatients, bipartite graph, costs, gap) 

    Constructs the minimum gene-metagene cover on the poor outcome samples.
    First a tripartite graph B = ((poor outcome patients, gene, metagenes), edges) is constructed as follows:
    For all samples with classlabel = True an edge is drawn to a gene, if this gene is mutated.
    Genes are connected to metagenes accoring to the secondary data source.
    In case of pairwise relations on the genes a Metagene is introduced for each defined pair.
    """
    
    print "CONSTRUCTING tripartite graphs"
    GRAPHS = {} #datasetName:BipartiteGraph
    for dataset in datasets:
        print "STATUS", dataset.name
        B = nx.Graph()
        if costs != None:
            for gene in costs:
                B.add_node(gene, costs=costs[gene])

        #select only the poor outcome samples
        for sample in dataset.patientLabels[dataset.patientClassLabels == True]:
            mutDat = dataset.expressionData[dataset.patientLabels == sample, ].flatten()
            mutGenes = dataset.geneLabels[mutDat == 1]
            B.add_nodes_from(dataset.geneLabels.tolist())
            for gene in mutGenes:
                B.add_edge(sample, gene)
        assert sum(sum(dataset.expressionData)) == len(B.edges())

        TriG = None
        if hasattr(secData, 'geneSets'):
            TriG = SampleGenePWDataGraph(B, secData, prefixGenes)
        elif hasattr(secData, 'edges'):
            TriG = SampleGeneInteractionGraph(B, secData, prefixGenes)
        else:
            print "ERROR secondary data not of type EdgeSet or GeneSetCollection."
            assert False
        GRAPHS[dataset.name] = TriG

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print "K Gene-Metagene Cover for dataset", name, "and secondary data", secData.name
        sol, solMetagenes, solPatients, gap = KCover_secData(GRAPHS[name], K, prefixGenes, "META_", None, costs)
        solutionGenes[name] = (sol, solMetagenes, solPatients, GRAPHS[name], gap)

        #pickle dump solutions
    filename = folder+"/K_GeneMetageneCover_K_"+str(K)+"_MutagenisisData_"+experimentName+".pickle"
    print "Solution genes, metagenes, covered patients and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

def KLConnectedGeneMetageneCoverFromMutatgenisisData(datasets, secData, K, L, prefixGenes, folder, experimentName,
    timelimit = None):
    """
    datasets:           list of expression datasets
    secData:            Some secondary data of type EdgeSet or GeneSetCollection. Either genes are directly linked to 
                        metagenes (GeneSetCollection) or pairwise relations on the genes are given (EdgeSet).
    K:                  minimum number of covering genes per patient
    L:                  maximum number of uncovered patients
    folder:             folder where results will be saved
    experimentName:     name of the savefile(s)
    prefixGenes:        Since in the graph we need to distinguish between patient nodes and gene nodes
                        genes carry a prefix. For the datasets provided use "Entrez".
    timelimit:          timelimit for each cplex instance (in seconds). Note that each dataset is one cplex instance.
    useCostfunction:    one can attach costs to the genes. 
                        Possible functions are: "Deregulation" and "DeregulationAvg" for a Gene-Sample graph
    pVal:               determines the cut off for the upper and lower tail. 

    returns
    
    solutionGenes:      a dictionary [dataset.name] = (solutiongenes, solution metagenes, solutionpatients, bipartite graph, costs, gap) 

    Constructs the minimum gene-metagene cover on the poor outcome samples such that the genes and metagenes induce a connected component.
    Firts a tripartite graph B = ((poor outcome patients, gene, metagenes), edges) is constructed as follows:
    For each gene an expression profile (probability function on the expression values)
    is learned from the good outcome samples (classlabel = False).
    Genes are defined "deregulated" for a poor outcome patient if their expression values lie in the
    upper or lower tail of the profile. The cut-off value is determined by 1/2*pVal for each tail. At a pVal of 1%
    a poor outcome patient must show at least an as extreme low (high) expression than 0.5% percent of the good outcome patients. 
    Only then the gene is considered to be deregulated.
    Genes are connected to metagenes accoring to the secondary data source.
    In case of pairwise relations on the genes a Metagene is introduced for each defined pair.
    """

    print "CONSTRUCTING tripartite graphs"
    GRAPHS = {} #datasetName:BipartiteGraph
    for dataset in datasets:
        print "STATUS", dataset.name
        B = nx.Graph()
        #select only the poor outcome samples
        for sample in dataset.patientLabels[dataset.patientClassLabels == True]:
            mutDat = dataset.expressionData[dataset.patientLabels == sample, ].flatten()
            mutGenes = dataset.geneLabels[mutDat == 1]
            for gene in mutGenes:
                B.add_edge(sample, gene)
        TriG = None
        if hasattr(secData, 'geneSets'):
            TriG = SampleGenePWDataGraph(B, secData, prefixGenes)
        elif hasattr(secData, 'edges'):
            TriG = SampleGeneInteractionGraph(B, secData, prefixGenes)
        else:
            print "ERROR secondary data not of type EdgeSet or GeneSetCollection."
            assert False
        GRAPHS[dataset.name] = TriG

    # set up cplex instance and solve it
    solutionGenes = {}
    for name in GRAPHS:
        print "KL connected Gene-Metagene Cover for dataset", name, "and secondary data", secData.name
        sol, solMetagenes, solPatients, gap = KLCover_secData_connectedComponent(GRAPHS[name], K, L, prefixGenes, "META_", timelimit, None)
        solutionGenes[name] = (sol, solMetagenes, solPatients, GRAPHS[name], gap)

        #pickle dump solutions
    filename = folder+"/KL_connectedGeneMetageneCover_K_"+str(K)+"_L_"+str(L)+"_MutagenisisData_"+experimentName+".pickle"
    print "Solution genes, metagenes, covered patients and graphs are saved as pickle file:", filename
    pickle.dump(solutionGenes, open(filename, "wb"))

    return solutionGenes

