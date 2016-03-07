#! /usr/bin/env python

#import time
#import scipy
#from scipy import stats
import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
import cplex
from cplex.callbacks import LazyConstraintCallback, UserCutCallback
import sys
import multiprocessing

from preprocessGraph import preprocessTripartiteG

class MyLazy(LazyConstraintCallback):
    """
    If the current LP solution (Genes and Metagenes) span a connected component in the tripartite
    graph B then the solution may also be valid for the ILP.
    If not we add a new constraint for all connected components in the current solution that do not contain 
    the root node. This constraint ensures that either the root node must be part of the connected component or 
    that the connected component must be expanded or the connected component is not part of the solution. 
    """

    def __call__(self):
        print "Adding lazy constraint"
        #get solution variables
        solVars = np.array(self.get_values()) #contains values for genes, metagenes and patients
        varNames = self.varNames
        varNames = np.array(varNames)
        B        = self.B
        prefixGenes = self.prefixGenes
        assert len(solVars) == len(varNames)       

        activeVars = varNames[solVars > 0]
        print "activeVars", activeVars 
        #get graph induced by the solution genes and metagenes
        # remove patients from G, we are only interested whether the solution genes and metagenes span a connected component
        #solution genes also contains the root variables and patient variables that are not nodes in B.
        G = B.subgraph(activeVars)
        print "Connected components:", nx.number_connected_components(G)
        ############################TESTING########################################
        #if 'g4' in activeVars:
        #    print 'bad gene g4 new try'
        #    idx = varNames == 'g4'
        #    self.add(cplex.SparsePair(['g4'], [1]), sense = "E", rhs = 0, use = 1)
        ###########################################################################
        if nx.number_connected_components(G) > 1:
            #add constraint
            root = [var.replace('Root', '') for var in activeVars if var.startswith('Root')][0]
            #########################################################
            # component without root cannot reach all patients -> sum var - sum rootvar - sum neighbours <= |C| - 2 
            # root component can also not reach  all patients -> sum var -sum rootvar - sum neighbours <= |C| - 2
            # --> each component may just be expanded if root and one other neighbour is present.
            components = nx.connected_components(G)
            for component in components:
                #print component
                allNeighbours = []
                for node in component:
                    neighbours = [neighbour for neighbour in B.neighbors(node) if neighbour not in component
                        and (neighbour.startswith(prefixGenes) or neighbour.startswith('META_'))]
                    allNeighbours.extend(neighbours)
                variables = component+['Root'+ n for n in component]+allNeighbours
                assert set(variables).issubset(varNames)
                values = list(np.ones(len(component)))+list(-1*np.ones(len(component)))+list(-1*np.ones(len(allNeighbours)))
                #print "New Constr", zip(variables, values), "<=", len(component)-1#sys.float_info.epsiloni
                #debug.append([cplex.SparsePair(variables, values), "L", len(component)-1])
                self.add(cplex.SparsePair(variables, values), sense = "L", rhs = len(component)-1)
        print


def SetCover(B, prefixGenes, timelimit, costs, K):
    """
    Given a bipartite graph with samples/patients as one set of nodes and genes as the second set of nodes 
    this function constructs the set cover formulation, which is then solved.

        OBJECTIVE   MIN \sum genes
        CONSTRAINTS for each patient: \sum genes >= K
    

    B:          bipartite graph between the samples and the genes or metagenes
    prefixGenes:genes or Metagenes must carry a prefix different from the prefixes for the samples.
    costs:      a dictionary mapping from Genes to costs. The objective function is minimised, thus genes with 
                low costs will be referred. costs: Genes -> [0, +Inf]. If undefined each gene will cost 1.
    timelimit:  huge instances will require a lot of memory and a lot of time. 
                The memory is by default 12,000 MB. To avoid long runningtimes you can also set the time limit. 
                NOTE: The solution returned then is not necessarily optimal.
    K: covereage parameter, default K = 1

    returns: 
    solutionGenes: The genes in the cover
    (gap, upper bound, solution value): If there is a solution the three values will be returned, otherwise the cplex Error code
    """
    assert bipartite.is_bipartite(B)    

    Patients = [node for node in B.nodes() if not node.startswith(prefixGenes)]
    Genes = [node for node in B.nodes() if node.startswith(prefixGenes)]
    Patients = sorted(list(Patients))
    Genes = sorted(list(Genes))
    print "Orig pat", len(Patients)
    
    # remove patients that have a lower degree than K
    Btmp = B.copy()
    removePatients = []
    for pat in Patients:
        if len(B.neighbors(pat)) < K:
            removePatients.append(pat)

    print "**********************Patients with low degree***********************************"
    print removePatients
    print "Number of outlier patients:", len(removePatients)
    print "*********************************************************************************"

    Btmp.remove_nodes_from(removePatients)
    Patients = [node for node in Btmp.nodes() if not node.startswith(prefixGenes)]
    Patients = sorted(list(Patients))

    # Initialise set cover
    # OBJECTIVE MIN \sum genes
    # CONSTRAINTS for each patient: \sum genes >= 1

    print "CPLEX ..."
    setCover = cplex.Cplex()
    setCover.parameters.workmem.set(12000)
    if timelimit is not None:
        setCover.parameters.timelimit.set(timelimit)
    setCover.parameters.threads.set(6)  

    setCover.parameters.mip.strategy.probe.set(3)
  
    print "Memory", setCover.parameters.workmem.get()
    print "Time", setCover.parameters.timelimit.get()
    print "Probing", setCover.parameters.mip.strategy.probe.get()
    print "Threads", setCover.parameters.threads.get()    

    setCover.objective.set_sense(setCover.objective.sense.minimize)
    
    #cost for each gene is 1
    if costs == None:
        my_obj = np.ones(len(Genes))
    else:
        c = np.ones(len(Genes))
        #fill costs vector
        for i in range(len(Genes)):
            c[i] = costs[Genes[i]]
        my_obj = np.array(c)
    # add a variable for each gene, variable name is "ENTREZ_number" 
    setCover.variables.add(obj = my_obj.tolist(), names = Genes, ub = np.ones(len(Genes)).tolist(), types = "B"*len(Genes))

    print "ADDING constraints"
    for pat in Patients:
        # get all neighbours of a patient (these are the deregulated genes), call by variable name, the weight for each gene is 1
        setCover.linear_constraints.add(lin_expr = [cplex.SparsePair(Btmp.neighbors(pat), np.ones(len(Btmp.neighbors(pat))).tolist())], \
            senses = ["G"], rhs = [K], names = [str(pat)])

    #solve
    print "SOLVING ... "
    setCover.solve()

    if setCover.solution.get_status() == 101 or setCover.solution.get_status() == 102: #Integer optimal
        #Retrieve solution genes
        print "RETRIEVING solution"
        solution = np.array(setCover.solution.get_values()) #binary solution vector of all variables, we just have genes as variables here
        print "Sum vars:", sum(solution)
        solutionGenes = np.array(setCover.variables.get_names())[solution > 0] # find the genes' (variables') names, that contribute to the solution
        print "Solution:", solutionGenes
        return (solutionGenes.tolist(), (None, setCover.solution.MIP.get_best_objective(), setCover.solution.MIP.get_best_objective()))
    elif setCover.solution.get_status() == 107: #timelimit exceeded
        print "Warning: Solution is suboptimal; check timelimit!" 
        print "Gap", setCover.solution.MIP.get_mip_relative_gap() 
        print "Lower bound", setCover.solution.MIP.get_best_objective()
        print "Actual cost", setCover.solution.get_objective_value()
        print "RETRIEVING solution"
        solution = np.array(setCover.solution.get_values()) #binary solution vector of all variables, we just have genes as variables here
        solutionGenes = np.array(setCover.variables.get_names())[solution == 1] # find the genes' (variables') names, that contribute to the solution
        print "Sum vars:", sum(solution)
        print "Solution:", solution        
        return (solutionGenes.tolist(), (setCover.solution.MIP.get_mip_relative_gap(), setCover.solution.MIP.get_best_objective(), 
            setCover.solution.get_objective_value()))
    else:
        print "ERROR", setCover.solution.get_status()
        return setCover.solution.get_status()

#Variant of the set cover.
# Combines the L-cover and the K-cover, employs secondary data
def KCover_secData(B, K, prefixGenes, prefixMetagenes, timelimit, costs):
    """
    B:               a tripartite graph B = ((Patients, Genes, Metagenes), (Patients x Genes, Genes x Metagenes))
    K:               number of genes required to cover a patient
    prefixGenes:     a prefix that is unique to all genes
    prefixMetagenes: a prefix that is unique to all metagenes
    costs:           a dictionary mapping from Genes \cup Metagenes -> R+

        OBJECTIVE   MIN \sum genes + \sum metagenes
        CONSTRAINTS for each patient: \sum genes >= K
                    for each gene   : \sum metagenes >= gene , if gene is chosen gene = 1 -> at lest one metagene must be chosen to cover gene

    Solves a set cover instance. Sets are the genes and the patients span the universe. The relation between patients, 
    genes and metagenes are given by a tripartite graph. We assume there are direct interctions between patients and genes
    (an edge translates to e.g. "is deregulated") and between genes and metagenes (an edge translates here to e.g. "belongs to")).
    You can use the parameters K to require that patients are covered by at least K genes.
    """
    
    Patients = [node for node in B.nodes() if not node.startswith(prefixGenes) and not node.startswith(prefixMetagenes)]
    Genes = [node for node in B.nodes() if node.startswith(prefixGenes)]
    Metagenes = [node for node in B.nodes() if node.startswith(prefixMetagenes)]
    GenesAndMetagenes = Genes+Metagenes

    # remove patients that have a lower degree than K
    Btmp = B.copy()
    removePatients = []
    for pat in Patients:
        if len(B.neighbors(pat)) < K:
            removePatients.append(pat)
    print "**********************Patients with low degree***********************************"
    print removePatients
    print "New number of outlier patients:", len(removePatients)
    print "*********************************************************************************"
    Btmp.remove_nodes_from(removePatients)
    Patients = [node for node in Btmp.nodes() if not node.startswith(prefixGenes) and not node.startswith(prefixMetagenes)]

    #initialise cplex instance
    KCover = cplex.Cplex()
    KCover.objective.set_sense(KCover.objective.sense.minimize)
    KCover.parameters.workmem.set(12000)
    if timelimit is not None:
        KCover.parameters.timelimit.set(timelimit)

    KCover.parameters.mip.strategy.probe.set(3)

    print "Memory", KCover.parameters.workmem.get()
    print "Time", KCover.parameters.timelimit.get()
    print "Probing", KCover.parameters.mip.strategy.probe.get()

    if costs == None:
        #genes cost 1, Metagenes cost 1
        c = np.hstack((np.ones(len(Genes)), np.ones(len(Metagenes))))
    else:
        print "costs", len(costs), "G&M", len(GenesAndMetagenes)
        assert len(costs) == len(GenesAndMetagenes)
        c = np.ones(len(GenesAndMetagenes))
        #fill costs vector
        for i in range(len(GenesAndMetagenes)):
            c[i] = costs[GenesAndMetagenes[i]]

    my_obj = np.concatenate([c, np.zeros(len(Patients))], 1)
    KCover.variables.add(obj = c.tolist()+np.zeros(len(Patients)).tolist(), names = Genes+Metagenes+Patients,
        ub = np.ones(len(Genes)+len(Patients)+len(Metagenes)).tolist(), types = "B"*(len(Genes)+len(Patients)+len(Metagenes)))

    #add constraints
    print "ADDING constraints"
    for pat in Patients:
        # get all neighbours of a patient (these are the deregulated genes), call by variable name
        KCover.linear_constraints.add(lin_expr = [cplex.SparsePair(Btmp.neighbors(pat), np.ones(len(Btmp.neighbors(pat))).tolist())], \
            senses = ["G"], rhs = [K], names = [str(pat)])

    for gene in Genes:
        # metagene1+metagene2+...+metageneN - gene >= 0
        metaGeneNeighbours = [m for m in Btmp.neighbors(gene) if m.startswith(prefixMetagenes)]
        if len(metaGeneNeighbours) > 0:
            KCover.linear_constraints.add(lin_expr = [cplex.SparsePair(metaGeneNeighbours+[gene], np.ones(len(metaGeneNeighbours)).tolist()+[-1.0])], 
                senses = ["G"], rhs = [0], names = ["constraint_"+str(gene)])

    print "SOLVING ... "
    KCover.solve()

    if KCover.solution.get_status() == 101 or KCover.solution.get_status() == 102: #Integer optimal
        #Retrieve solution genes
        print "RETRIEVING solution"
        solution = np.array(KCover.solution.get_values()) #binary solution vector of all variables
        solutionVars = np.array(KCover.variables.get_names())[solution > 0].tolist() #all variables that contribute to the solution
        print "Sum vars:", sum(solution)
        print "Solution:", solution
        solutionGenes = [item for item in solutionVars if item.startswith(prefixGenes)]
        solutionMetagenes = [item for item in solutionVars if item.startswith(prefixMetagenes)]
        solutionPatients = [item for item in solutionVars if not item.startswith(prefixGenes) and not item.startswith(prefixMetagenes)]

        return (solutionGenes, solutionMetagenes, solutionPatients, 
            (KCover.solution.MIP.get_mip_relative_gap(), KCover.solution.MIP.get_best_objective(), KCover.solution.get_objective_value()))
    elif KCover.solution.get_status() == 107: #timelimit exceeded
        print "Warning: Solution is suboptimal; check timelimit!"
        print "Gap", KCover.solution.MIP.get_mip_relative_gap()
        print "Lower bound", KCover.solution.MIP.get_best_objective()
        print "Actual cost", KCover.solution.get_objective_value()
        print "RETRIEVING solution"
        solution = np.array(KCover.solution.get_values()) #binary solution vector of all variables
        solutionVars = np.array(KCover.variables.get_names())[solution == 1].tolist() #all variables that contribute to the solution
        print "Sum vars:", sum(solution)
        print "Solution:", solution

        solutionGenes = [item for item in solutionVars if item.startswith(prefixGenes)]
        solutionMetagenes = [item for item in solutionVars if item.startswith(prefixMetagenes)]
        solutionPatients = [item for item in solutionVars if not item.startswith(prefixGenes) and not item.startswith(prefixMetagenes)]

        return (solutionGenes, solutionMetagenes, solutionPatients, 
            (KCover.solution.MIP.get_mip_relative_gap(), KCover.solution.MIP.get_best_objective(), KCover.solution.get_objective_value()))
    else:
        print "ERROR", KCover.solution.get_status()
        return KCover.solution.get_status()

#Variant of the set cover.
# Combines the L-cover and the K-cover, employs secondary data
def KLCover_secData_connectedComponent(B, K, L, prefixGenes, prefixMetagenes, timelimit, costs):
    """
    B:               a tripartite graph B = ((Patients, Genes, Metagenes), (Patients x Genes, Genes x Metagenes))
    K:               number of genes required to cover a patient
    L:               max number of patients not covered K-times
    prefixGenes:     a prefix that is unique to all genes
    prefixMetagenes: a prefix that is unique to all metagenes
    costs:           a dictionary mapping from Genes \cup Metagenes -> R+

        OBJECTIVE   MIN \sum genes 
        CONSTRAINTS \sum z_u >= len(Patients)-L, z_u is an indicator variable for the patients
                    for each patient: \sum genes >= K*z_u, if z_u is 1 for a patient then there must be k genes covering the patient
        Connectivity constraints:
                    \sum root <= 1, for each metagene and gene a root variable; choose one of them to be root
                    for each gene and metagene:     gene >= root; choose root node from covering genes or connecting metagenes
                    for each gene and metagene:     gene - root - \sum neighbours <= 0; gene is either root or it has neighbours in the cover
                    for each subsets of genes and metagenes:
                        \sum (gene - root) - \sum neighbour_of_subset < |subset| - 1; these are exponentially many constraints -> callback           

    Solves a set cover instance. Sets are the genes and the patients span the universe. The relation between patients, 
    genes and metagenes are given by a tripartite graph. We assume there are direct interctions between patients and genes
    (an edge translates to e.g. "is deregulated") and between genes and metagenes (an edge translates here to e.g. "belongs to")).
    You can use the parameters K, L and P to adjust the solution to certain requirements. E.g. if a pathway should only be part 
    of the solution only if a certain number of its genes are part of the cover.
    The solution is a connected component in the gene-metagene graph.
    """

    Btmp = B.copy()

    #Analyse input graph:
    # * find all connected components in the gene-metagene graph that cover all but L patients --> validCC
    # * find all connected components that cover less than n-L patients --> otherCC
    # * genes and metagenes from the otherCC can not be part of the solution, remove them from the graph
    validCC, otherCC = preprocessTripartiteG(B, prefixGenes, prefixMetagenes, L)
    
    if len(validCC) == 0:
        print "Graph does not contain a valid connected compoenent."
        return "ERROR", 103

    excludeGM = [item for sublist in otherCC for item in sublist]
    print "********************Exclude genes and metagenes******************************"
    print excludeGM
    print "*****************************************************************************"

    Btmp.remove_nodes_from(excludeGM)

    Genes = [node for node in Btmp.nodes() if node.startswith(prefixGenes)]
    Metagenes = [node for node in Btmp.nodes() if node.startswith(prefixMetagenes)]
    GenesAndMetagenes = Genes+Metagenes

    # remove patients that have a lower degree than K
    Patients = [node for node in Btmp.nodes() if not node.startswith(prefixGenes) and not node.startswith(prefixMetagenes)]
    removePatients = []
    for pat in Patients:
        if len(Btmp.neighbors(pat)) < K:
            removePatients.append(pat)

    print "**********************Patients with low degree***********************************"
    print removePatients
    print "New number of outlier patients:", L + len(removePatients)
    print "*********************************************************************************"

    Btmp.remove_nodes_from(removePatients)
    Patients = [node for node in Btmp.nodes() if not node.startswith(prefixGenes) and not node.startswith(prefixMetagenes)]

    Znames = ['Z'+str(pat) for pat in Patients]
    rootNames = ['Root'+element for element in GenesAndMetagenes]

    KLCover = cplex.Cplex()
    KLCover.objective.set_sense(KLCover.objective.sense.minimize)

    KLCover.parameters.workmem.set(12000)
    if timelimit is not None:
        KLCover.parameters.timelimit.set(timelimit)

    KLCover.parameters.mip.strategy.probe.set(3)
    KLCover.parameters.threads.set(min(multiprocessing.cpu_count(), 10))

    print "Memory", KLCover.parameters.workmem.get()
    print "Time", KLCover.parameters.timelimit.get()
    print "Probing", KLCover.parameters.mip.strategy.probe.get()
    print "Threads", KLCover.parameters.threads.get()

    if costs == None:
        #genes cost 1, Metagenes cost 0
        c = np.hstack((np.ones(len(Genes)), np.zeros(len(Metagenes))))
    else:
        assert len(costs) == len(GenesAndMetagenes)
        c = np.ones(len(GenesAndMetagenes))
        #fill costs vector
        for i in range(len(GenesAndMetagenes)):
            c[i] = costs[GenesAndMetagenes[i]]
    
    #introduce root variable for genes and metagenes
    my_obj = np.concatenate([c, np.zeros(len(Patients)), np.zeros(len(GenesAndMetagenes))], 1)

    KLCover.variables.add(obj = my_obj.tolist(), names = Genes+Metagenes+Znames+rootNames,
        ub = np.ones(len(Genes)+len(Patients)+len(Metagenes)+len(rootNames)).tolist(), types = "B"*(len(Genes)+len(Patients)+len(Metagenes)+len(rootNames)))

    #add constraints
    print "ADDING constraints"
    for pat in Patients:
        zvar = "Z"+str(pat)
        assert zvar in Znames
        # gene1+gene2+...+geneN - K*zvar >= 0
        # neighbours of patients are always genes
        KLCover.linear_constraints.add(lin_expr = [cplex.SparsePair(Btmp.neighbors(pat)+[zvar], np.ones(len(Btmp.neighbors(pat))).tolist()+[-1.0*K])],
            senses = ["G"], rhs = [0], names = ["constraint_"+str(pat)])

    #choose one root node
    KLCover.linear_constraints.add(lin_expr = [cplex.SparsePair(rootNames, np.ones(len(rootNames)).tolist())], senses = ["E"], rhs = [1],
        names = ["ROOT_CONSTR"])
   
    # x >= rootx 
    for rootVar in rootNames:
        #cut off 'Root' prefix
        nodeVar = rootVar.replace('Root', '')
        KLCover.linear_constraints.add([cplex.SparsePair([nodeVar, rootVar], [1, -1])],
            senses = ["G"], rhs = [0], names = ["VarGreaterRoot_"+str(nodeVar)])
    
    # sum z-u <= |U| - L
    KLCover.linear_constraints.add(lin_expr = [cplex.SparsePair(Znames, np.ones(len(Znames)).tolist())], senses = ["G"], rhs = [len(Patients)-L],
        names = ["Z_CONSTR"])

    #for each gene and metagene
    for node in rootNames:
        nodeVar = node.replace('Root', '')
        # x - rootx - sum (neighbours(x)) <= 0
        neighboursNode = [n for n in Btmp.neighbors(nodeVar) if not n in Patients]
        variables = [nodeVar, node]
        variables.extend(neighboursNode)
        values = [1, -1]
        values.extend(list((-1)*np.ones(len(neighboursNode))))
        KLCover.linear_constraints.add([cplex.SparsePair(variables, values)],
            senses = ["L"], rhs = [0], names = ["RootOrNeighbour_"+str(nodeVar)])

    lazyConstraint = KLCover.register_callback(MyLazy)
    lazyConstraint.varNames = Genes+Metagenes+Znames+rootNames
    lazyConstraint.B = Btmp
    lazyConstraint.prefixGenes = prefixGenes

    print "SOLVING ... "
    KLCover.solve()

    if KLCover.solution.get_status() == 101 or setCover.solution.get_status() == 102: #Integer optimal or very small gap
        #Retrieve solution genes
        print "RETRIEVING solution"
        solution = np.array(KLCover.solution.get_values()) #binary solution vector of all variables
        solutionVars = np.array(KLCover.variables.get_names())[solution > 0].tolist() #all variables that contribute to the solution
        print "Sum vars:", sum(solution)
        print "Solution:", solution
        solutionGenes = [item for item in solutionVars if item.startswith(prefixGenes)]
        solutionMetagenes = [item for item in solutionVars if item.startswith(prefixMetagenes)]
        solutionPatients = [item for item in solutionVars if item.startswith('Z')]

        return (solutionGenes, solutionMetagenes, solutionPatients, 
            (KLCover.solution.MIP.get_mip_relative_gap(), KLCover.solution.MIP.get_best_objective(), KLCover.solution.get_objective_value()))

    elif KLCover.solution.get_status() == 107: #timelimit exceeded
        print "Warning: Solution is suboptimal; check timelimit!"
        print "Gap", KLCover.solution.MIP.get_mip_relative_gap()
        print "Lower bound", KLCover.solution.MIP.get_best_objective()
        print "Actual cost", KLCover.solution.get_objective_value()
        print "RETRIEVING solution"
        solution = np.array(KLCover.solution.get_values()) #binary solution vector of all variables
        solutionVars = np.array(KLCover.variables.get_names())[solution == 1].tolist() #all variables that contribute to the solution
        print "Sum vars:", sum(solution)
        print "Solution:", solution

        solutionGenes = [item for item in solutionVars if item.startswith(prefixGenes)]
        solutionMetagenes = [item for item in solutionVars if item.startswith(prefixMetagenes)]
        solutionPatients = [item for item in solutionVars if item.startswith('Z')]

        return (solutionGenes, solutionMetagenes, solutionPatients,
            (KLCover.solution.MIP.get_mip_relative_gap(), KLCover.solution.MIP.get_best_objective(), KLCover.solution.get_objective_value()))

    else:
        print "ERROR", KLCover.solution.get_status()
        return KLCover.solution.get_status()

