# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013
import xlrd
import numpy

def getReplicates(DictStudyToMeanToSamples):
    DictSampleToReplicates = {}
    allKeys = []
    for study in DictStudyToMeanToSamples:
        values = DictStudyToMeanToSamples[study].keys()
        allKeys.extend(values)
    allKeys = numpy.unique(allKeys)    

    for key in allKeys:
        samples = []
        for study in DictStudyToMeanToSamples:
            if key in DictStudyToMeanToSamples[study]:
                samples.append(DictStudyToMeanToSamples[study][key])
        for num in range(len(samples)):
            s = [sam.replace('gsm', 'GSM') for sam in samples]
            DictSampleToReplicates[samples[num].replace('gsm', 'GSM')] = s
 
    return DictSampleToReplicates

def readRecycledSamples():
    # reads the samples that were found by Gyorffy 
    # to be identical or rehybridised.
    # Samples with the same Avg expression value are identical.
    # It was checked for some examples that also the annotation data contained
    # the same information.
    # In Gyorffy et al. samples that were published first were 
    # taken into account and the later replicates were removed.

    wb = xlrd.open_workbook('experiments/data/Recycled breast arrays.xls')
    sh = wb.sheet_by_index(0)    

    #first column is empty
    DictStudyToMeanToSamples = {}
    #DictStudyToMeanToSamples[column[len(column)-1]] = []

    column = sh.col_values(1)

    Avg = sh.col_values(2)
    while '' in Avg:
        Avg.remove('')
    while u'Average' in Avg:
        Avg.remove(u'Average')

    pat = sh.col_values(3)
    while '' in pat:
        pat.remove('')
    while u'GEO sample' in pat:
        pat.remove(u'GEO sample')

    DictStudyToMeanToSamples[column[len(column)-1]] = dict(zip(Avg, pat))

    column = sh.col_values(5)
    while '' in column:
        column.remove('')

    Avg = sh.col_values(6)
    while '' in Avg:
        Avg.remove('')
    while u'Average' in Avg:
        Avg.remove(u'Average')

    pat = sh.col_values(7)
    while '' in pat:
        pat.remove('')
    while u'GEO sample' in pat:
        pat.remove(u'GEO sample')

    DictStudyToMeanToSamples[column[len(column)-1]] = dict(zip(Avg, pat))

    column = sh.col_values(9)
    while '' in column:
        column.remove('')

    Avg = sh.col_values(10)
    while '' in Avg:
        Avg.remove('')
    while u'Average' in Avg:
        Avg.remove(u'Average')

    pat = sh.col_values(11)
    while '' in pat:
        pat.remove('')
    while u'GEO sample' in pat:
        pat.remove(u'GEO sample')

    DictStudyToMeanToSamples[column[len(column)-1]] = dict(zip(Avg, pat))

    column = sh.col_values(13)
    while '' in column:
        column.remove('')

    Avg = sh.col_values(14)
    while '' in Avg:
        Avg.remove('')
    while u'Average' in Avg:
        Avg.remove(u'Average')

    pat = sh.col_values(15)
    while '' in pat:
        pat.remove('')
    while u'GEO sample' in pat:
        pat.remove(u'GEO sample')

    DictStudyToMeanToSamples[column[len(column)-1]] = dict(zip(Avg, pat))

    DictSampleToReplicates = getReplicates(DictStudyToMeanToSamples)

    return DictSampleToReplicates

def readAnnotation(annotationfile, datasets):
    """
    Input
    annotationfile: Excel file with annotation.

    Output
    Dataset2Time:   Dictionary, mapping from dataset.name to a numpy.array with the survival times
    """

    Dataset2Time = {}

    wb = xlrd.open_workbook(annotationfile)
    sh = wb.sheet_by_index(0)
    colNames = sh.row_values(0)
    Annotation = {} #Annotation[columnName] = {[gsmNumber] = value}
    for col in range(1, len(colNames)):
        name = colNames[col]
        D = dict(zip(sh.col_values(0), sh.col_values(col)))
        Annotation[name] = D

    #We are interested in RFS time and DMFS time
    IDXtimeKeys = [6, 19] #DMFS, RFS
    for d in datasets:
        #Samples in my data that are not present in Gyorffy's annotation
        if '/' in d.patientLabels[0]:
            samples = [l.split('/')[1] for l in d.patientLabels]
        d.patientLabels = numpy.array(samples)
        SamplesWithoutClasslabels = numpy.unique(numpy.setdiff1d(samples, Annotation[Annotation.keys()[0]].keys())).tolist()

        DictSampleToReplicates = readRecycledSamples()
        noRep = []
        for sample in SamplesWithoutClasslabels:
            if sample in DictSampleToReplicates:
                #check if corresponding sample is in Annotation
                for annoKey in Annotation.keys(): # for each keyword copy entries
                    #intersect pseudonyme with annotation keys
                    samplesInAnno = Annotation[annoKey].keys()
                    intersect = list(set(DictSampleToReplicates[sample]).intersection(samplesInAnno))
                    synonymKey = intersect[0]
                    if len(intersect) > 0:
                        Annotation[annoKey][sample] = Annotation[annoKey][synonymKey]
                    else:
                        noRep.append(sample)
        assert len(noRep) == 0
        survTime = []
        for patient in samples:
            if patient not in noRep:
                if 'DMFS' in d.name:
                    assert Annotation[Annotation.keys()[6]][patient] != ''
                    survTime.append(Annotation[Annotation.keys()[6]][patient])
                elif 'RFS' in d.name:
                    assert Annotation[Annotation.keys()[19]][patient] != ''
                    survTime.append(Annotation[Annotation.keys()[19]][patient])
                else:
                    break
        Dataset2Time[d.name] = numpy.array(survTime)
    return Dataset2Time

