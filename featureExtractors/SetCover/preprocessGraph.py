import networkx as nx

def preprocessTripartiteG(B, prefixGenes, prefixMetagenes, L):
    """
    Preprocesses a tripartite graph B = ((patients, genes, metagenes), E) where edges
    are only defined between patients and genes and genes and metagenes.
    Preprocessing:
    Find all connected components that cover all apart from L patients.
    Returns the connected components.

    B   : tripartite graph
    L   : number of outlier patients
    """

    print "Graph preprocessing: checking for connected components that do not cover all but L samples."
    Patients = [node for node in B.nodes() if not node.startswith(prefixGenes) and not node.startswith(prefixMetagenes)] 
    Genes = [node for node in B.nodes() if node.startswith(prefixGenes)]
    Metagenes = [node for node in B.nodes() if node.startswith(prefixMetagenes)]

    GMgraph = B.subgraph(Genes+Metagenes)
    
    CC = nx.connected_components(GMgraph)
    validCC = []
    otherCC = []
    for component in CC:
        #get number of patients
        genes = [n for n in component if n.startswith(prefixGenes)]
        coveredPatients = []
        for gene in genes:
            coveredPatients.extend(n for n in B.neighbors(gene) if not n.startswith(prefixMetagenes))

        if len(set(coveredPatients)) >= len(Patients)-L:
            validCC.append(component)
        else:
            otherCC.append(component)

    return validCC, otherCC
        
