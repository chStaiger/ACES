# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

import subprocess
import tempfile
import os
import shutil
import numpy
import time

#The subprocess will be called from a temporary directory. Thus on one machine several calls can be done.
#wDCBFilename executable file

def runWCDB(wDCBFilename, dataset, network, ExpressedThreshold=0.1, 
    DensityThreshold=0.5, MinimumCases=4, deleteTemporaryDirectory = False):

    print 'MinimumCases', MinimumCases

    # Some properties of the expression data that you might need
    actualOutcome = dataset.patientClassLabels
    featureLabels = dataset.geneLabels
    sampleLabels = dataset.patientLabels
    numControls = dataset.numPatientsGoodOutcome
    numCases = dataset.numPatientsBadOutcome
    expressionData = dataset.expressionData
    
    # Show some stats about overlap between the genes in the expression dataset (primary data) and in the network dataset (secondary data).

    network_nodes = frozenset(network.getNodes())
    feature_nodes = frozenset(dataset.geneLabels)

    print "NOTE: Number of nodes in both EXPRESSION and NETWORK data:", len(network_nodes & feature_nodes)
    print "NOTE: Number of nodes in EXPRESSION data but not in NETWORK data:", len(feature_nodes - network_nodes)
    print "NOTE: Number of nodes in NETWORK data but not in EXPRESSION data:", len(network_nodes - feature_nodes)

    # Check dimensions of the input.

    (ns, nf) = expressionData.shape
    assert actualOutcome.shape == (ns, )
    assert featureLabels.shape == (nf, )
    assert sampleLabels.shape  == (ns, )

    #Create wDCB input-files in a temporary directory.

    print "Executable", wDCBFilename

    tempdir = tempfile.mkdtemp()
    
    print tempdir
    
    #Create inputfiles
    wDCBInputInitFilename     = os.path.join(tempdir, "init_file.txt")
    wDCBInputMappingFilename   = os.path.join(tempdir, "mapping_file.txt")
    wDCBInputMatrixFilename    = os.path.join(tempdir, "matrix_file.txt")
    wDCBInputNetworkFilename   = os.path.join(tempdir, "network_file.txt")
    wDCBOutputFilename         = os.path.join(tempdir, "output_file.txt")

    #write the network and the expression data to files
    dataset.writeToFile(wDCBInputMatrixFilename)
    network.writeEdgesPlusWeights(wDCBInputNetworkFilename)
    #write a mapping file nodes in network to gene in expression data
    f = open(wDCBInputMappingFilename, "wb")
    for i in network.getNodes():
        f.write(i+"\t"+i+"\n") 
    f.close()   

    print 'MinimumCases', MinimumCases
    f = open(wDCBInputInitFilename, "wb")
    f.write("WorkingFolder="+tempdir+"/data\n")
    f.write("NetworkFile="+wDCBInputNetworkFilename+"\n")
    f.write("IDFile="+wDCBInputMappingFilename+"\n")
    f.write("ExpressionFile="+wDCBInputMatrixFilename+"\n")
    f.write("#Controls="+str(dataset.numPatientsGoodOutcome)+"\n")
    f.write("#Cases="+str(dataset.numPatientsBadOutcome)+"\n")
    f.write("ExpressedThreshold=%.4f\n" % ExpressedThreshold)
    f.write("DensityThreshold=%.4f\n" % DensityThreshold)
    f.write("MinimumCases=%d\n" % MinimumCases)
    f.close()

    # Execute wDCB as a sub-process
    # subprocessWDCB calls your executable and writes the networks and their score back to a file. The name of the file is given as last parameter.
    subprocessWDCB(wDCBFilename, wDCBInputInitFilename)

    # Result lies in "tempdir/data/WDCB/alpha_"+str(0.5)+"_minGraph_"+str(MinimumCases)+"_minSize_4_distFunc_2.txt"
    # Mapping from node index to geneID lies in "tempdir/data/nodes.txt"

    resultfile = tempdir+"/data/WDCB/alpha_"+str(DensityThreshold)+"_minGraph_"+str(MinimumCases)+"_minSize_4_distFunc_2.txt"
    mappingfile = tempdir+"/data/nodes.txt" 

    # Parse the result file as produced by wDCB.
    sortedGeneModules = readWDCBOutputFile(resultfile, mappingfile)

    # We are done with the temporary directory. Delete it.

    if deleteTemporaryDirectory:
        shutil.rmtree(tempdir)

    return sortedGeneModules

def subprocessWDCB (
        wDCBFilename, 
        wDCBInputInitFilename
        # other parameters you might need
    ):

    #TODO: Here we create the cammand line call
    #execulatble lies in featureExtractors/Dao/runwDCB with respect to the experiemnts directory
    args = []
    #print wDCBFilename.split("../featureExtractors/Dao/")[1]
    print os.path.basename(wDCBFilename)
    #args.extend(["./"+wDCBFilename.split("../featureExtractors/Dao/")[1]])
    args.extend(["./"+os.path.basename(wDCBFilename)])
    args.extend([wDCBInputInitFilename])

    print "NOTE: Executing runwDCB; commandline: %s" % " ".join(args)
    print "Change to dir:", os.path.dirname(wDCBFilename)
    print "Init file", wDCBInputInitFilename   

    tic = time.time()
    proc = subprocess.Popen(args, cwd=os.path.dirname(wDCBFilename))
    exitcode = proc.wait()
    toc = time.time()
    
    print "NOTE: wDCB finished (exitcode = %d); runtime (wallclock time): %.3f seconds." % (exitcode, toc - tic)

    assert exitcode == 0

def readWDCBOutputFile(resultfile, mappingfile):

    geneModules = []

    f = open(resultfile)

    #Discard the first two lines
    lines = f.readlines()
    lines = lines[2:len(lines)]
    
    f.close()
    
    mapping = {}
    f = open(mappingfile)
    for line in f:
        key = line.split('\t')[0]
        value = line.split('\t')[1].strip('\n')
        mapping[key] = value
    f.close()
    

    for line in lines:
        if line.startswith("\n") or line.startswith("\t"):
            continue # Ignore
        if ' Weight = ' in line: 
            line = line.split(' Weight = ')
            module    = line[0].split(' ')
            for i in range(len(module)): #replace indices with geneID
                module[i] = mapping[module[i]]
            score     = float(line[1].split('\n')[0])
            geneModule = (score, module)
            geneModules.append(geneModule)

    # translate indices to geneIDs

    return geneModules
