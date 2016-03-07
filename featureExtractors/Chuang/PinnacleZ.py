# This is a Python-wrapper to the PinnacleZ Java program.

import subprocess
import tempfile
import os
import shutil
import numpy
import time


def runPinnacleZ(pzJarFilename, pzNormalizeInputMatrix, actualOutcome, featureLabels, sampleLabels, expressionData, network, deleteTemporaryDirectory = True):

    # Show some stats about overlap between the genes in the expression dataset (primary data) and in the network dataset (secondary data).

    network_nodes = frozenset.union(*network)
    feature_nodes = frozenset(featureLabels)

    print "NOTE: Number of nodes in both EXPRESSION and NETWORK data:", len(network_nodes & feature_nodes)
    print "NOTE: Number of nodes in EXPRESSION data but not in NETWORK data:", len(feature_nodes - network_nodes)
    print "NOTE: Number of nodes in NETWORK data but not in EXPRESSION data:", len(network_nodes - feature_nodes)

    # Check dimensions of the input.

    (ns, nf) = expressionData.shape
    assert actualOutcome.shape == (ns, )
    assert featureLabels.shape == (nf, )
    assert sampleLabels.shape  == (ns, )

    # Create PinnacleZ input-files in a temporary directory.

    tempdir = tempfile.mkdtemp()

    pzInputClassFilename   = os.path.join(tempdir, "class_file.txt")
    pzInputMatrixFilename  = os.path.join(tempdir, "matrix_file.txt")
    pzInputNetworkFilename = os.path.join(tempdir, "network_file.txt")
    pzOutputFilename       = os.path.join(tempdir, "output_file.txt")

    writePzClassFile(pzInputClassFilename, sampleLabels, actualOutcome)
    writePzMatrixFile(pzInputMatrixFilename, featureLabels, sampleLabels, expressionData)

    writePzNetworkFile(pzInputNetworkFilename, network)

    # Execute PinnacleZ as a sub-process

    subprocessPinnacleZ(pzJarFilename, pzInputClassFilename, pzInputMatrixFilename, pzInputNetworkFilename, pzOutputFilename, pzNormalizeInputMatrix)

    # Parse the result file as produced by PinnacleZ.

    geneModules = readPzOutputFile(pzOutputFilename)

    # We receive the scores with 6 digits of precision.
    # Recalculate the score to make sure we understand its algoritm properly.

    verifyModuleScores(geneModules, actualOutcome, featureLabels, expressionData)

    # We are done with the temporary directory. Delete it.

    if deleteTemporaryDirectory:
        shutil.rmtree(tempdir)

    return geneModules


def scoreSubnetwork(genes, actualOutcome, expressionData, geneLabelToIndex):

    # Determines the MI score for a subnetwork, following the implementation of PinnacleZ.
    #
    # This is rather confusing, because it EXACTLY replicates what PinnacleZ is doing,
    #   and the method used by PinnacleZ is strange -- probably buggy.

    # The scoring of a subnetwork is done in two stages:
    #
    #   (1) Calculate an activation value for each of the modules, for each of the patients
    #
    #   (2) Considering the activation scores and the classes, calculate the Mutual Information (MI) score.

    # Stage (1). Calculate an 'activation' value for each of the modules, for each of the patients

    n = len(actualOutcome)

    values = numpy.zeros(n)
    count = 0
    lastEmpty = True

    # NOTE: this is the loop over the genes; it counts nodes for which we have no expression data,
    # yet it (tries to) jump over them. We know that this code reproduces the module score as
    # calculated by PinnacleZ. We feel that the handling of missing expression data is buggy;
    # this makes it hard to trust the networks as produced by PinnacleZ.
    #
    # Unfortunately we didn't get a good response about this issue from the Ideker group.

    for node in genes:

        count += 1

        if node not in geneLabelToIndex:
            lastEmpty = True
        elif lastEmpty:
            values = expressionData[:,geneLabelToIndex[node]]
            lastEmpty = False
        else:
            values = (expressionData[:,geneLabelToIndex[node]] + values * numpy.sqrt(count - 1)) / numpy.sqrt(count)

    # The 'values' array now contains an activity score for each of the patients, used for determining the Mutual Information.

    # Stage (2). Calculate the Mutual Information (MI) score for the module.

    minv = min(values)
    maxv = max(values)

    # TODO: understand this widening of the bins. It is undocumented in the paper.
    dx = float(maxv - minv) / (n - 1) / 2.0

    minv -= dx
    maxv += dx

    def floor_log2(n):
        assert n > 0
        if n == 1:
            return 0
        else:
            return 1 + floor_log2(n // 2)

    nbins = floor_log2(n) + 1

    bincount = numpy.zeros(shape = (2, nbins), dtype = numpy.int)

    for i in range(n):
        c = int(actualOutcome[i])
        idx = numpy.round((values[i] - minv) / (maxv - minv) * nbins + 0.5)
        bincount[c, idx - 1] += 1

    # Check that all patients are accounted for
    assert numpy.sum(bincount) == n

    # Calculate the mutual information
    rowSum = numpy.sum(bincount, axis = 0)
    colSum = numpy.sum(bincount, axis = 1)

    mi = 0.0
    for i in range(nbins):
        for j in range(2):
            if bincount[j, i] > 0:
                mi += bincount[j, i] * numpy.log(float(bincount[j, i]) / rowSum[i] / colSum[j])

    # TODO: understand this bit.

    mi = mi / float(n) + numpy.log(n) - float(nbins - 1) / float(2 * n)

    return mi


def verifyModuleScores(geneModules, actualOutcome, featureLabels, expressionData):

    geneLabelToIndex = dict(zip(featureLabels, xrange(len(featureLabels))))

    ok = True
    for (score, stats, genes) in geneModules:

        mi = scoreSubnetwork(genes, actualOutcome, expressionData, geneLabelToIndex)

        if ("%.6f" % mi) != ("%.6f" % score):
            #print "PinnacleZ score difference: ", mi, score
            ok = False

    if ok:
        print "NOTE: PinnacleZ score verification PASSED."
    else:
        print "WARNING: PinnacleZ score verification FAILED (input data not normalized?)"


def writePzClassFile(filename, sampleLabels, actualOutcome):

    f = open(filename, "w")

    actualOutcome = map(int, actualOutcome)

    for i in xrange(len(actualOutcome)):
        print >> f, sampleLabels[i], actualOutcome[i]

    f.close()


def writePzMatrixFile(filename, featureLabels, sampleLabels, expressionData):

    f = open(filename, "w")

    # write header line: a single dummy string, followed by the sample (patient) labels

    print >> f, " ".join(["-"] + list(sampleLabels))

    # Write one line per feature (gene).
    # Each line consists of a feature label, followed by a number for each of the samples,
    #   denoting the expression.

    for fi in xrange(expressionData.shape[1]):
        print >> f, " ".join([str(featureLabels[fi])] + ["%f" % expression for expression in expressionData[:,fi]])

    f.close()


def writePzNetworkFile(filename, network):

    relationship = "is_connected_to"

    f = open(filename, "w")

    for (vFrom, vTo) in network:
        print >> f, "%s %s %s" % (vFrom, relationship, vTo)

    f.close()


def subprocessPinnacleZ (
        pzJarFilename,
        pzInputClassFilename,
        pzInputMatrixFilename,
        pzInputNetworkFilename,
        pzOutputFilename,
        pzNormalizeInputMatrix,
        st1cutoff      = 0.05,
        st2cutoff      = 0.05,
        st3cutoff      = 0.00005,
        st3trials      = 20000,
        maxNodeDegree  = 300,
        minImprovement = 0.05,
        maxModuleSize  = 20,
        maxRadius      = 2,
        scoreModel     = "MI",
        numberOfTrials = 100
    ):

    args = []
    args.extend(["java", "-jar", pzJarFilename])
    args.extend(["-v"])

    # NOTE: The "-z" option is introduced in our PATCHED version of PinnacleZ.
    # It suppresses the genewise z-score normalization usually performed by PinnacleZ on its expression matrix.
    if not pzNormalizeInputMatrix:
        args.extend(["-z"])

    args.extend(["-o", pzOutputFilename])
    args.extend(["--st1cutoff", str(st1cutoff)])
    args.extend(["--st2cutoff", str(st2cutoff)])
    args.extend(["--st3cutoff", str(st3cutoff)])
    args.extend(["--st3trials", str(st3trials)])
    args.extend(["--maxNodeDegree", str(maxNodeDegree)])
    args.extend(["--minImprovement", str(minImprovement)])
    args.extend(["--maxModuleSize", str(maxModuleSize)])
    args.extend(["--maxRadius", str(maxRadius)])
    args.extend(["--score", scoreModel])
    args.extend(["--trials", str(numberOfTrials)])
    args.extend([pzInputClassFilename, pzInputMatrixFilename, pzInputNetworkFilename])

    print "NOTE: Executing PinnacleZ; commandline: %s" % " ".join(args)

    tic = time.time()
    exitcode = subprocess.call(args)
    toc = time.time()
    print "NOTE: PinnacleZ finished (exitcode = %d); runtime (wallclock time): %.3f seconds." % (exitcode, toc - tic)
    assert exitcode == 0


def readPzOutputFile(filename):

    geneModules = []

    f = open(filename)

    for line in f:

        if line.startswith("#"):
            continue # Ignore comment lines

        line = line.split()

        seed    = line[0]
        score   = float(line[1])
        pValues = tuple(map(float, line[2:5]))
        genes   = line[5:]

        assert genes[0] == seed
        assert len(genes) == len(frozenset(genes))

        geneModule = (score, pValues, genes)

        geneModules.append(geneModule)

    f.close()

    return geneModules
