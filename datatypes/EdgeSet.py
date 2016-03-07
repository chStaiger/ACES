# @Author 
# Sidney Cadot, update in 2012/13 Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import random

class EdgeSet:

    def __init__(self, name, edges, edgeweights):
        assert all([len (edge) == 2 for edge in edges])
        self.name = name
        self.edges = edges
        self.edgeweights = edgeweights

    def getNodes(self):

        return frozenset.union(*self.edges)

    def copy(self):
        return EdgeSet(self.name, self.edges, self.edgeweights)

    def makeShuffle(self, seed = 0):
        
        prng = random.Random(seed)

        nodes = sorted(frozenset.union(*self.edges))
        shuffled_nodes = nodes[:] # make a copy
        prng.shuffle(shuffled_nodes)

        return dict(zip(nodes, shuffled_nodes))

    def makeShuffledVersion(self, name, shuffle):

        edges = set()
        for (v1, v2) in self.edges:
            v1 = shuffle[v1]
            v2 = shuffle[v2]
            edge = frozenset([v1, v2])
            assert edge not in edges
            edges.add(edge)

        edges = frozenset(edges)
        #TODO: if edges are shuffled, what are the edgeweights?
        return EdgeSet(name, edges, self.edgeweights)
    
    def addEdgeWeights(self, weights):
    # edges is a list of sets (a, b) where a and b are nodes in the network.
        self.edgeweights = dict(zip(self.edges, weights))

    # This function is used for the Dao featutre extractor. The network needs to be written with its 
    # weights to a file 
    def writeEdgesPlusWeights(self, filename):
        """
        Writes edges and edge weights to a space-separated file.
        node1 node2 weight
        node3 node4 weight
        ...
        """
        if self.edgeweights == None:
            self.edgeweights = dict(zip(self.edges, [999]*len(self.edges)))

        f = open(filename, 'w')
        for key in self.edgeweights.keys():
            f.write(str(list(key)[0])+" "+str(list(key)[1])+" "+str(self.edgeweights[key])+"\n")
        f.close()

    def writeSIF(self, filename):
        """
        Writes the network to file in sif-format.
        """
        lines = []
        for edge in self.edges:
            lines.append(str(list(edge)[0])+' '+self.name+' '+str(list(edge)[1]))       
        f = open(filename, 'w')
        f.writelines(["%s\n" % item  for item in lines])
        f.close()


def removeNodes(network, nodes):
    print "Remove nodes from edgefile: copy edges"
    newEdges = set(network.edges.copy())
    edgeweights = None
    if network.edgeweights is not None:
        print "Remove nodes from edgefile: copy weights"
        edgeweights = network.edgeweights.copy()

    for e in network.edges:
        if len(e.intersection(nodes)) > 0:
            newEdges.remove(e)
            if network.edgeweights is not None:
                edgeweights.pop(e)
    newEdges = frozenset(newEdges)
    return EdgeSet(network.name, newEdges, edgeweights)
        
def ReadSIF(networkname, filename, prefix = None):

    # Read a SIF network file.

    edges = set()
    relationship = None
    numDuplicates = 0
    numEdgeToSelf = 0

    f = open(filename)
    for line in f:

        line = line.rstrip("\r\n")

        if len(line) == 0:
            continue # Ignore empty lines.

        line_parts = line.split()
        assert len(line_parts) == 3

        (x, rel, y) = line_parts

        if relationship is None:
            relationship = rel
        else:
            assert relationship == rel

        if prefix is not None and not x.startswith(prefix) and not y.startswith(prefix):
            x = prefix + x
            y = prefix + y

        edge = frozenset([x, y])

        if len(edge) == 1:
            numEdgeToSelf += 1
        elif edge in edges:
            numDuplicates += 1
        else:
            edges.add(edge)

    f.close()

    edges = frozenset(edges)

    numDistinctNodes = len(frozenset.union(*edges))

    print "NOTE: ReadSIF(\"%s\", \"%s\"): read %d good edges; ignored %d edges-to-self and %d duplicate edges; \
        %d distinct nodes." % (networkname, filename, len(edges), numEdgeToSelf, numDuplicates, numDistinctNodes)

    return EdgeSet(networkname, edges, None)

# The file should look like this:
# noe1 node2 weight
# delim specifies the delimiter.
# Prefix attaches a prefix to all nodes in the network. 
def ReadEdgesPlusWeights(networkname, filename, prefix, delim):

    edges = set()
    weights = {} 
    relationship = None
    numDuplicates = 0
    numEdgeToSelf = 0

    f = open(filename)
    for line in f:

        line = line.rstrip("\r\n")

        if len(line) == 0:
            continue # Ignore empty lines.

        line_parts = line.split(delim)
        assert len(line_parts) == 3

        (x, y, weight) = line_parts

        if prefix is not None:
            x = prefix + x
            y = prefix + y

        edge = frozenset([x, y])

        if len(edge) == 1:
            numEdgeToSelf += 1
        elif edge in edges:
            numDuplicates += 1
        else:
            edges.add(edge)
            weights[edge] = weight
    f.close()

    edges = frozenset(edges)

    numDistinctNodes = len(frozenset.union(*edges))

    print "NOTE: ReadSIF(\"%s\", \"%s\"): read %d good edges; ignored %d edges-to-self and %d duplicate edges; \
        %d distinct nodes." % (networkname, filename, len(edges), numEdgeToSelf, numDuplicates, numDistinctNodes)

    return EdgeSet(networkname, edges, weights)


