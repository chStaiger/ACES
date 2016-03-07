import time
import scipy
from scipy import stats
import numpy as np
import networkx as nx

def CompleteWeightedBipartiteGraph(dataset, thresholdUp=None, thresholdDown=None, mul = 1, debug = False):
    """
    Constructs the complete bipartite graph.Edges are drawn between patients and genes. 
    One can give thresholds for the over expression and under expression. Edges will be scored as 'real' edges 
    if the expression in a patient lies above thresholdUp or below thresholdDown. 
    If these two parameters are not given, thresholdUp = mean(gene)+mul*std(gene) and thresholdDown = mean(gene)-mul*std(gene)

    The graph's edges are annotated with 'expression' (actual matrix entry), 'deregulation' (the z-score of the gene for this patient)
    and the 'score'.
    The 'score' is the distance to thresholdUp iff expression is higher than the mean expression or the distance to thresholdDown iff
    expression is lower than the mean expression of the gene. 
    The sign in front of 'score' indicates whether the expression value lies within [thresholdDown, thresholdUp] -- 
    indicated by a negative score; or outside of this interval -- positive score. 

    For binary data use thresholdUp = 0 and thresholdDown = -1. 
    If you work with binary data ignore the field 'deregulation' in the returned graph! 
    """

    edges = []
    #translate index of gene to Entrez GeneID
    dict_IdxToGeneNames = dict(zip(range(dataset.expressionData.shape[1]), dataset.geneLabels.tolist()))

    print "DETERMINE edges ..."

    tic = time.time()
    Tup = []
    Tdown = []
    for gene in range(0, dataset.expressionData.shape[1]):
        #upper, lower: cut-off values in the tail of the distribution.
        mean = dataset.expressionData[:, gene].mean()
        std = dataset.expressionData[:, gene].std()
        #means.append(mean)
        #stds.append(std)
        for pat in range(0, dataset.expressionData.shape[0]):
            try:
                patName = dataset.patientLabels[pat].split("/")[1] #if / in patient label only take the latter part
            except:
                patName = dataset.patientLabels[pat]
            #print pat, gene, dataset.expressionData[pat, gene]
            if thresholdUp != None and thresholdDown!=None:
                if dataset.expressionData[pat, gene]>thresholdUp:
                    edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene], 
                        (abs(dataset.expressionData[pat, gene]) - abs(thresholdUp))/std))
                elif dataset.expressionData[pat, gene]<thresholdDown:
                    edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene], 
                        (abs(dataset.expressionData[pat, gene]) - abs(thresholdDown))/std))
                elif dataset.expressionData[pat, gene] > mean:
                    #if the expr is close 0 then give a high score, the '-' indicates taht it is a negative edge
                    edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene], 
                        -((abs(thresholdUp) - abs(dataset.expressionData[pat, gene]))/std)))
                else:
                    edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene], 
                        -((abs(thresholdDown) - abs(dataset.expressionData[pat, gene]))/std)))
            else:
                Up = mean + mul*std
                Down = mean - mul*std
                Tup.append(Up)
                Tdown.append(Down)
                if dataset.expressionData[pat, gene]>Up:
                   edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene],
                        (abs(dataset.expressionData[pat, gene]) - abs(Up))/std))
                elif dataset.expressionData[pat, gene]<Down:
                    edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene],
                        (abs(dataset.expressionData[pat, gene]) - abs(Down))/std))
                elif dataset.expressionData[pat, gene] > mean:
                    #if the expr is close to the mean then give a high score, the '-' indicates taht it is a negative edge
                    edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene],
                        -((abs(Up) - abs(dataset.expressionData[pat, gene]))/std)))
                else:
                    edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene],
                        -((abs(Down) - abs(dataset.expressionData[pat, gene]))/std)))

    #bipartite graph
    print "CONSTRUCT bipartite graph"
    B = nx.Graph()
    for n1, n2, dereg, expr, score in edges:
        B.add_nodes_from([n1, n2])
        B.add_edge(n1, n2, deregulation = dereg, expression = expr, score = score)

    assert nx.is_bipartite(B)

    toc = time.time()
    print "Graph constructed in %.3f seconds." % (toc - tic)

    if debug:
        return (B, Tup, Tdown)
    else:
        return B

def BipartiteGraphTwoSided(dataset, pVal):
    """
    Construct bipartite graph B = (S, G)
    S are the poor outcome patients and G are the deregulated genes.
    For each gene a reference (gaussian) distribution is estimated from the good outcome patients. 
    Then for each poor outcome patient it is tested whether the gene's value (e.g. expression) lies in the top
    pVal/2 or bottom pVal/2 tails.

    dataset: An ExpressionDataset
    pVal:    cut off value for the tails
    """

    #divide into two patient groups
    GoodP = dataset.extractPatientsByIndices("Good", dataset.patientClassLabels==False, False, False)
    PoorP = dataset.extractPatientsByIndices("Poor", dataset.patientClassLabels==True, False, False)

    edges = []
    #translate index of gene to Entrez GeneID
    dict_IdxToGeneNames = dict(zip(range(dataset.expressionData.shape[1]), dataset.geneLabels.tolist()))

    print "DETERMINE edges ..."

    tic = time.time()
    for gene in range(0, GoodP.expressionData.shape[1]):
        dist = scipy.stats.norm.fit(GoodP.expressionData[:, gene])
        #upper, lower: cut-off values in the tail of the distribution.
        upper = scipy.stats.norm.ppf(1-pVal/2, dist[0], dist[1])
        lower = scipy.stats.norm.ppf(pVal/2, dist[0], dist[1])
        mean = GoodP.expressionData[:, gene].mean()
        std = GoodP.expressionData[:, gene].std()
        for pat in range(0, PoorP.expressionData.shape[0]):
            try:
                patName = PoorP.patientLabels[pat].split("/")[1] #if / in patient label only take the latter part
            except:
                patName = PoorP.patientLabels[pat]
            if PoorP.expressionData[pat, gene]>upper:
                edges.append((patName, dict_IdxToGeneNames[gene], (PoorP.expressionData[pat, gene]- mean)/std, PoorP.expressionData[pat, gene]))
            elif PoorP.expressionData[pat, gene]<lower:
                edges.append((patName, dict_IdxToGeneNames[gene], (PoorP.expressionData[pat, gene]- mean)/std, PoorP.expressionData[pat, gene]))

    #bipartite graph
    print "CONSTRUCT bipartite graph"
    B = nx.Graph()
    for n1, n2, dereg, expr in edges:
        B.add_edge(n1, n2, deregulation = dereg, expression = expr)

    assert nx.is_bipartite(B)

    toc = time.time()
    print "Graph constructed in %.3f seconds." % (toc - tic)

    return B

def BipartiteGraphSTD(dataset, mul = 1):
    """
    Construct bipartite graph B = (S, G)
    S are the patients and G are the deregulated genes.
    For each gene the expression profile is estimated from all patients.
    Then for each patient it is tested whether the gene's value (e.g. expression) lies outside of the standard deviation.
    If yes, an edge is drawn between the gene and the patient.

    dataset: An ExpressionDataset
    percentileUp, percentileLow: the cut-off percentiles
    usereferenceGood: if True the percentiles are determined on the profiles of good outcome patients only.
    """

    edges = []
    #translate index of gene to Entrez GeneID
    dict_IdxToGeneNames = dict(zip(range(dataset.expressionData.shape[1]), dataset.geneLabels.tolist()))

    print "DETERMINE edges ..."

    tic = time.time()
    for gene in range(0, dataset.expressionData.shape[1]):
        dist = scipy.stats.norm.fit(dataset.expressionData[:, gene])
        #upper, lower: cut-off values in the tail of the distribution.
        mean = dataset.expressionData[:, gene].mean()
        std = dataset.expressionData[:, gene].std()
        for pat in range(0, dataset.expressionData.shape[0]):
            try:
                patName = dataset.patientLabels[pat].split("/")[1] #if / in patient label only take the latter part
            except:
                patName = dataset.patientLabels[pat]
            if dataset.expressionData[pat, gene]>mean+mul*std:
                edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene]))
            elif dataset.expressionData[pat, gene]<mean-mul*std:
                edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene]))

    #bipartite graph
    print "CONSTRUCT bipartite graph"
    B = nx.Graph()
    for n1, n2, dereg, expr in edges:
        B.add_edge(n1, n2, deregulation = dereg, expression = expr)

    assert nx.is_bipartite(B)

    toc = time.time()
    print "Graph constructed in %.3f seconds." % (toc - tic)

    return B


def BipartiteGraphPercentiles(dataset, percentileUp, percentileLow, usereferenceGood = False):
    """
    Construct bipartite graph B = (S, G)
    S are the poor outcome patients and G are the deregulated genes.
    For each gene the corresponding values for the upper and lower percentiles are calculated from its expression profile. 
    Then for each patient it is tested whether the gene's value (e.g. expression) lies in the top or bottom percentile.

    dataset: An ExpressionDataset
    percentileUp, percentileLow: the cut-off percentiles
    usereferenceGood: if True the percentiles are determined on the profiles of good outcome patients only.
    """

    #divide into two patient groups
    if usereferenceGood:
        GoodP = dataset.extractPatientsByIndices("Good", dataset.patientClassLabels==False, False, False)
        PoorP = dataset.extractPatientsByIndices("Poor", dataset.patientClassLabels==True, False, False)
    else:
        GoodP = dataset
        PoorP = dataset

    edges = []
    #translate index of gene to Entrez GeneID
    dict_IdxToGeneNames = dict(zip(range(dataset.expressionData.shape[1]), dataset.geneLabels.tolist()))

    print "DETERMINE edges ..."

    tic = time.time()
    for gene in range(0, GoodP.expressionData.shape[1]):
        #upper, lower: cut-off values in the tail of the distribution.
        upper = scipy.stats.scoreatpercentile(GoodP.expressionData[:, gene], percentileUp)
        lower = scipy.stats.scoreatpercentile(GoodP.expressionData[:, gene], percentileLow)
        mean = GoodP.expressionData[:, gene].mean()
        std = GoodP.expressionData[:, gene].std()
        for pat in range(0, PoorP.expressionData.shape[0]):
            try:
                patName = PoorP.patientLabels[pat].split("/")[1] #if / in patient label only take the latter part
            except:
                patName = PoorP.patientLabels[pat]
            if PoorP.expressionData[pat, gene]>upper:
                edges.append((patName, dict_IdxToGeneNames[gene], (PoorP.expressionData[pat, gene]- mean)/std, PoorP.expressionData[pat, gene]))
            elif PoorP.expressionData[pat, gene]<lower:
                edges.append((patName, dict_IdxToGeneNames[gene], (PoorP.expressionData[pat, gene]- mean)/std, PoorP.expressionData[pat, gene]))

    #bipartite graph
    print "CONSTRUCT bipartite graph"
    B = nx.Graph()
    for n1, n2, dereg, expr in edges:
        B.add_edge(n1, n2, deregulation = dereg, expression = expr)

    assert nx.is_bipartite(B)

    toc = time.time()
    print "Graph constructed in %.3f seconds." % (toc - tic)

    return B

def BipartiteGraphExprThreshold(dataset, thresholdUp, thresholdDown):
    """
    Construct bipartite graph B = (S, G)
    S are the patients and G are the deregulated genes.
    A threshold on the expression values determines wether the gene is deregulated.

    dataset: An ExpressionDataset
    thresholdUp, thresholdDown: the cut-off values
    """

    edges = []
    #translate index of gene to Entrez GeneID
    dict_IdxToGeneNames = dict(zip(range(dataset.expressionData.shape[1]), dataset.geneLabels.tolist()))

    print "DETERMINE edges ..."

    tic = time.time()
    for gene in range(0, dataset.expressionData.shape[1]):
        #upper, lower: cut-off values in the tail of the distribution.
        mean = dataset.expressionData[:, gene].mean()
        std = dataset.expressionData[:, gene].std()
        for pat in range(0, dataset.expressionData.shape[0]):
            try:
                patName = dataset.patientLabels[pat].split("/")[1] #if / in patient label only take the latter part
            except:
                patName = dataset.patientLabels[pat]
            print pat, gene, dataset.expressionData[pat, gene]
            if dataset.expressionData[pat, gene]>thresholdUp:
                edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene]))
            elif dataset.expressionData[pat, gene]<thresholdDown:
                edges.append((patName, dict_IdxToGeneNames[gene], (dataset.expressionData[pat, gene]- mean)/std, dataset.expressionData[pat, gene]))

    #bipartite graph
    print "CONSTRUCT bipartite graph"
    B = nx.Graph()
    for n1, n2, dereg, expr in edges:
        B.add_edge(n1, n2, deregulation = dereg, expression = expr)

    assert nx.is_bipartite(B)

    toc = time.time()
    print "Graph constructed in %.3f seconds." % (toc - tic)

    return B

# As BipartiteGraphTeoSided. Here edges between the poor outcome patients
# and the genes are introduced when the value lies in the top of the distribution.
def BipartiteGraphTop(dataset, pVal):
    """
    As BipartiteGraphTwoSided. Here edges between the poor outcome patients
    and the genes are introduced when the value lies in the top of the distribution.
    """

    #divide into two patient groups
    GoodP = dataset.extractPatientsByIndices("Good", dataset.patientClassLabels==False, False, False)
    PoorP = dataset.extractPatientsByIndices("Poor", dataset.patientClassLabels==True, False, False)

    edges = []
    #translate index of gene to Entrez GeneID
    dict_IdxToGeneNames = dict(zip(range(dataset.expressionData.shape[1]), dataset.geneLabels.tolist()))

    print "DETERMINE edges ..."

    tic = time.time()
    for gene in range(0, GoodP.expressionData.shape[1]):
        dist = scipy.stats.norm.fit(GoodP.expressionData[:, gene])
        #upper, lower: cut-off values in the tail of the distribution.
        upper = scipy.stats.norm.ppf(1-pVal/2, dist[0], dist[1])
        mean = GoodP.expressionData[:, gene].mean()
        std = GoodP.expressionData[:, gene].std()
        for pat in range(0, PoorP.expressionData.shape[0]):
            patName = PoorP.patientLabels[pat].split("/")[1] #if / in patient label only take the latter part
            if PoorP.expressionData[pat, gene]>upper:
                edges.append((patName, dict_IdxToGeneNames[gene], (PoorP.expressionData[pat, gene]- mean)/std, PoorP.expressionData[pat, gene]))

    #bipartite graph
    print "CONSTRUCT bipartite graph"
    B = nx.Graph()
    for n1, n2, dereg, expr in edges:
        B.add_edge(n1, n2, deregulation = dereg, expression = expr)

    assert nx.is_bipartite(B)

    toc = time.time()
    print "Graph constructed in %.3f seconds." % (toc - tic)

    return B

# As BipartiteGraphTeoSided. Here edges between the poor outcome patients
# and the genes are introduced when the value lies on the bottom of the distribution.
def BipartiteGraphBottom(dataset, pVal):

    """
    As BipartiteGraphTeoSided. Here edges between the poor outcome patients
    and the genes are introduced when the value lies in the top of the distribution.
    """    
    #divide into two patient groups
    GoodP = dataset.extractPatientsByIndices("Good", dataset.patientClassLabels==False, False, False)
    PoorP = dataset.extractPatientsByIndices("Poor", dataset.patientClassLabels==True, False, False)

    edges = []
    #translate index of gene to GeneID
    dict_IdxToGeneNames = dict(zip(range(dataset.expressionData.shape[1]), dataset.geneLabels.tolist()))

    print "DETERMINE edges ..."

    tic = time.time()
    for gene in range(0, GoodP.expressionData.shape[1]):
        dist = scipy.stats.norm.fit(GoodP.expressionData[:, gene])
        #upper, lower: cut-off values in the tail of the distribution.
        lower = scipy.stats.norm.ppf(pVal/2, dist[0], dist[1])
        mean = GoodP.expressionData[:, gene].mean()
        std = GoodP.expressionData[:, gene].std()
        for pat in range(0, PoorP.expressionData.shape[0]):
            patName = PoorP.patientLabels[pat].split("/")[1] #if / in patient label only take the latter part
            if PoorP.expressionData[pat, gene]<lower:
                edges.append((patName, dict_IdxToGeneNames[gene], (PoorP.expressionData[pat, gene]- mean)/std, PoorP.expressionData[pat, gene]))

    #bipartite graph
    print "CONSTRUCT bipartite graph"
    B = nx.Graph()
    for n1, n2, dereg, expr in edges:
        B.add_edge(n1, n2, deregulation = dereg, expression = expr)

    assert nx.is_bipartite(B)

    toc = time.time()
    print "Graph constructed in %.3f seconds." % (toc - tic)

    return B

def SampleGenePWDataGraph(B, pwData, prefixGenes):
    """
    Constructs the tripartite graph between samples, genes and metagenes.
    Metagenes will receive the prefix 'META_' to distinguish them from the rest of the nodes.
    If the field geneSetsNames is empty in the secondary data, pathways will be called 'pw_number'.
    Note that in the tripartite graph genes that are not deregulated are also included.    
 
    B:                   a bipartite graph between samples and genes
    pwData:              gene sets as given by a pathway database e.g. Use the class GeneSetCollection
    prefixGenes:         A prefix that is common to all genes in B and distinguishes them from the patients.
    
    Output:
    Tripartite graph. Edges between samples and genes carry the annotation "deregulation" and "expression".
    Edges between genes and metagenes do not carry any annotation. The metagenes carry the annotation "size", 
    i.e. the number of genes in the pathway.
    """

    assert all(g.startswith(prefixGenes) for g in pwData.getNodes())
    
    Genes = [node for node in B.nodes() if node.startswith(prefixGenes)]
    Samples = [node for node in B.nodes() if not node.startswith(prefixGenes)]

    TriG = B.copy() #add all edges from samples to deregulated genes
    #add edges from genes to metagenes
    #Metagenes will receive the prefix 'META_'

    if len(pwData.geneSetsNames) > 0:
        pwNames = [pw.replace(' ', '_') for pw in pwData.geneSetsNames]
    else:
        pwNames = ['pw_'+str(num) for num in range(0, len(pwData))]

    for i in range(0, len(pwData.geneSets)):
        TriG.add_node("META_"+pwNames[i], size = len(pwData.geneSets[i]))
        for gene in pwData.geneSets[i]:
            TriG.add_edge("META_"+pwNames[i], gene)

    return TriG

def SampleGeneInteractionGraph(B, interactionData, prefixGenes):
    """
    Constructs the tripartite graph between samples, genes and their interactions.
    For each interaction a metagene is constructed connecting the two genes.
    The Metagenes will be called "META_number". Note that in the tripartite graph 
    genes that are not deregulated are also included.
    
    B:                   a bipartite graph between samples and genes
    interactionData:     pairwise interactions between genes. Use  the class EdgeSet.
    
    Output:
    Tripartite graph. Edges between samples and genes carry the annotation "deregulation" and "expression".
    Edges between genes and metagenes do naot carry any annotation. The metagenes carry the annotation "size", 
    i.e. the number of genes in the metagene (2).    
    """

    assert all(g.startswith(prefixGenes) for g in interactionData.getNodes())

    Genes = [node for node in B.nodes() if node.startswith(prefixGenes)]
    Samples = [node for node in B.nodes() if not node.startswith(prefixGenes)]

    TriG = B.copy() #add all edges from samples to deregulated genes

    #for each interaction (g1, g2) introduce a Metagene and add (g1, Metagene) and (g2, Metagene)
    #to the graph

    countMetagenes = 0
    for g1, g2 in interactionData.edges:
        TriG.add_node("META_"+str(countMetagenes), size = 2)
        TriG.add_edge("META_"+str(countMetagenes), g1)
        TriG.add_edge("META_"+str(countMetagenes), g2)
        countMetagenes = countMetagenes+1

    return TriG

def SampleSecDataGraph(TriG, prefixGenes, numDeregulatedGenes):
    """
    TriG:   a tripartite graph as constructed in SampleGenePWDataGraph or SampleGeneInteractionGraph.
            The edges must carry the annotation 'deregulation' and 'expression' or no annotation at all; and the metagene nodes
            must carry the annotation 'size'.
    numDeregulatedGenes: the value defines the minimum number of deregulated genes shared between a 
                         pathway and a sample. Only if there are at least numDeregulatedGenes genes shared between 
                         sample and pathway, the pathway will be counted as deregulated for the sample and the edge will 
                         be added to the graph.
 
    Removes the gene layer in the tripartite grapoh and connects the samples directly with the metagenes.
    
    The new edges will carry the annotaion 
        "cumDereg": a list with each gene's deregulation that also lies in te metagene
        "cumExpr":  a list with each gene's expression
    The metagenes will carry the annotaion 'size'
    """

    Genes = [node for node in TriG.nodes() if node.startswith(prefixGenes)]
    Metagenes = [node for node in TriG.nodes() if node.startswith("META_")]
    Samples = [node for node in TriG.nodes() if not node.startswith(prefixGenes) 
        and not node.startswith("META_")] 

    exampleEdge = None
    for edge in TriG.edges():
        if edge[0] in Samples or edge[1] in Samples:
            exampleEdge = edge
            break
    assert exampleEdge is not None

    exprData = False
    if len(TriG[exampleEdge[0]][exampleEdge[1]]) > 0:
        exprData = True
    print 'exprData', exprData
 
    edges = {} #here we will fill the new edges for the bipartite graph B = ((Samples, Metagenes), edges)
               #an edge will look like this edges[(sample, metagene)] = (cumDereg, cumExpr)
    for m in Metagenes:
        genes = TriG.neighbors(m)
        for gene in genes:
            samples = [s for s in TriG.neighbors(gene) if not s.startswith("META_")]
            for sample in samples:
                if exprData:
                    dereg = TriG[sample][gene]['deregulation']
                    expr  = TriG[sample][gene]['expression']
                if (sample, m) in edges:
                    if exprData:
                        cumDereg, cumExpr, deregGenes = edges[(sample, m)]
                        cumDereg.append(dereg)
                        cumExpr.append(expr)
                        deregGenes.append(gene)
                        edges[(sample, m)] = (cumDereg, cumExpr, deregGenes)
                    else:
                        deregGenes = edges[(sample, m)]
                        deregGenes.append(gene)
                        edges[(sample, m)] = deregGenes
                else:
                    if exprData:
                        edges[(sample, m)] = ([dereg], [expr], [gene])
                    else:
                        edges[(sample, m)] = ([gene])

    B = nx.Graph()
    #add edges with annotation
    for (sample, meta) in edges:
        if exprData:
            dereg, expr, deregGenes = edges[(sample, meta)]
            if len(dereg) >= numDeregulatedGenes:
                B.add_edge(sample, meta,  cumDereg = dereg, cumExpr = expr, deregGenes = deregGenes)
        else:
            deregGenes = edges[(sample, meta)]
            if len(deregGenes) >= numDeregulatedGenes:
                B.add_edge(sample, meta, deregGenes = deregGenes)

    #add Metagene annotation
    for m in Metagenes:
        if m in B:
            B.add_node(m, TriG.node[m])
    
    return B

def exportJena(G, outfile, NodeSet1 = None, NodeSet2 = None):
    """
    Exports a graph to a jena formatted file. If NodeSet1 and nodeSet2 are specified, 
    0-edges will be introduced between the elements of each set
    """

    numNodes = len(G.nodes())
    nodes = sorted(G.nodes())

    upperTriangle = []
    for i in range(len(nodes)-1):
        neighbours = []
        for j in range(i+1, len(nodes)):
            print nodes[i], nodes[j]
            if (nodes[i], nodes[j]) in B.edges():
                neighbours.append(np.sign(B[nodes[i]][nodes[j]]['score']))
            elif (nodes[j], nodes[i]) in B.edges():
                neighbours.append(np.sign(B[nodes[j]][nodes[i]]['score']))
            elif (nodes[i] in NodeSet1 and nodes[j] in NodeSet1) or (nodes[i] in NodeSet2 and nodes[j] in NodeSet2):
                neighbours.append(0)
            else:
                neighbours.append(-1)
        upperTriangle.append(neighbours)

    f = open(outfile, 'w')
    f.write(str(numNodes)+'\n')
    
    for n in nodes:
        f.write(n+'\n')
    
    for line in upperTriangle:
        l = str(line)
        l = l.replace('[', '').replace(']', '').replace(',', ' ')
        f.write(l+'\n')
    f.close() 



