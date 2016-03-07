# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com
# July 2013

from create_tokens import generate_tokens, generate_tokens_innerloop
from SetUpGrid import CombineDataExperiment01, CombineDataExperiment05 
from localProcessToken import localProcess

import datetime

def DoEXP01outer():
    """
    Executes the outer CV loop for all methods, networks and parameters and the two datasets.
    """
    #Do not repeat the CV
    nr_crossvals = 1
    #Split data in five parts
    nr_folds = 5

    #Combine data, networks methods and other parameters
    DataAndFeatureExtractors = CombineDataExperiment01()
    #Create tokens that also contain the info about the splits and number CVs we want to execute
    #For each possible combination of crossval, fold and parameters in DataAndFeatureExtractors
    #we receive one token. They all carry the Prefix EXP01.
    #If you work with a HPC you can store them as tokens in a database or other tool.
    tokens = generate_tokens(DataAndFeatureExtractors, nr_crossvals, nr_folds, 'EXP01')

    #Execute calculations and store results
    #1) Locally process tokens
    doneTokens = []
    for token in tokens:
        #db may be a couchDB db
        doneTokens.append(localProcess(token, db = None))
    #2) For a couchDB use the pipeline script

    #Write results to an sql database
    sqlName = "_"+datetime.datetime.now().strftime('%b-%d-%G')
    sqlPath = "Results/"
    sqlFilename = sqlPath+"EXP01"+sqlName+".sqlite3"
    TokenToSqliteExperiment01(tokens, sqlFilename)
    #If you work with a database like couchDB use getDoneTokens to sort the tokens according to the
    #experiments and safe them in the appropriate files.
    #getDoneTokens(db, sqlPath, sqlName = ""), db must be a couchDB  or a dictionary.

    #For further analysis and plots look into the folder experiments

def DoEXP01inner():
    """
    Executes the inner CV loop for all methods, networks and parameters and the two datasets.
    """
    
    #We need to know how the outer loop looks like
    nr_crossvals = 1
    nr_folds = 5
    #Parameters for the inner loop
    nr_crossvals_inner = 1
    nr_folds_inner = 5

    #Combine data, networks methods and other parameters
    DataAndFeatureExtractors = CombineDataExperiment01()
    #Create tokens.
    #For each possible combination of outer and inner crossval, outer and inner fold and parameters in DataAndFeatureExtractors
    #we receive one token. They all carry the Prefix EXP01.
    #If you work with a HPC you can store them as tokens in a database or other tool.
    tokens = generate_tokens_innerloop(DataAndFeatureExtractors, nr_crossvals, nr_folds, nr_crossvals_inner, nr_folds_inner, 'EXP01InnerLoop')
    
    #Execute calculations and store results
    #1) Locally process tokens
    doneTokens = []
    for token in tokens:
        #db may be a couchDB db
        doneTokens.append(localProcess(token, db = None))
    #2) For a couchDB use the pipeline script

    #Write results to an sql database
    sqlName = "_"+datetime.datetime.now().strftime('%b-%d-%G')
    sqlPath = "Results/"
    sqlFilename = sqlPath+"EXP01InnerLoop"+sqlName+".sqlite3"
    TokenToSqliteExperiment01InnerLoop(tokens, sqlFilename)
    #If you work with a database like couchDB use getDoneTokens to sort the tokens according to the
    #experiments and safe them in the appropriate files.
    #getDoneTokens(db, sqlPath, sqlName = ""), db must be a couchDB  or a dictionary.

    #For further analysis and plots look into the folder experiments

def DoEXP04outer():
    """
    Executes the outer CV loop for all methods, networks and parameters and the two datasets with the ER positive patients only.
    """
    #Comments as in DoEXP01outer
    nr_crossvals = 1
    nr_folds = 5

    DataAndFeatureExtractors = CombineDataExperiment01()
    tokens = generate_tokens(DataAndFeatureExtractors, nr_crossvals, nr_folds, 'EXP04')

    doneTokens = []
    for token in tokens:
        doneTokens.append(localProcess(token, db = None))
    
    sqlName = "_"+datetime.datetime.now().strftime('%b-%d-%G')
    sqlPath = "Results/"
    sqlFilename = sqlPath+"EXP04"+sqlName+".sqlite3"
    TokenToSqliteExperiment01(tokens, sqlFilename)

def DoEXP04inner():
    """
    Executes the inner CV loop for all methods, networks and parameters and the two datasets.
    """
    #Comments as in DoEXP01inner
    nr_crossvals = 1
    nr_folds = 5
    nr_crossvals_inner = 1
    nr_folds_inner = 5

    DataAndFeatureExtractors = CombineDataExperiment01()
    tokens = generate_tokens_innerloop(DataAndFeatureExtractors, nr_crossvals, nr_folds, nr_crossvals_inner, nr_folds_inner, 'EXP04InnerLoop')

    doneTokens = []
    for token in tokens:
        doneTokens.append(localProcess(token, db = None))

    sqlName = "_"+datetime.datetime.now().strftime('%b-%d-%G')
    sqlPath = "Results/"
    sqlFilename = sqlPath+"EXP04InnerLoop"+sqlName+".sqlite3"
    TokenToSqliteExperiment01InnerLoop(tokens, sqlFilename)

def DoEXP05outer():
    """
    As DoEXP01outer but also shuffling the network and pathway data.
    """
    #Comments as in DoEXP01outer
    nr_crossvals = 1
    nr_folds = 5

    #Shuffle networks and pathways 25 times
    num_shuffles = 25

    DataAndFeatureExtractors = CombineDataExperiment05()
    tokens = generate_tokens(DataAndFeatureExtractors, nr_crossvals, nr_folds, 'EXP05')

    doneTokens = []
    for token in tokens:
        doneTokens.append(localProcess(token, db = None))

    sqlName = "_"+datetime.datetime.now().strftime('%b-%d-%G')
    sqlPath = "Results/"
    sqlFilename = sqlPath+"EXP05"+sqlName+".sqlite3"
    TokenToSqliteExperiment05(tokens, sqlFilename)


def DoEXP05inner():
    """
    As DoEXP01inner, but does additionally shuffle the network data n times.
    Note, that methods, that do not employ networks are not considered here.
    """

    #We need to know how the outer loop looks like
    nr_crossvals = 1
    nr_folds = 5
    #Parameters for the inner loop
    nr_crossvals_inner = 1
    nr_folds_inner = 5

    #Shuffle the networks and pathways 25 times
    num_shuffles = 25

    #Combine data, networks methods and other parameters
    DataAndFeatureExtractors = CombineDataExperiment05(num_shuffles)
    #Create tokens.
    #For each possible combination of outer and inner crossval, outer and inner fold and parameters in DataAndFeatureExtractors
    #we receive one token. They all carry the Prefix EXP01.
    #If you work with a HPC you can store them as tokens in a database or other tool.
    tokens = generate_tokens_innerloop(DataAndFeatureExtractors, nr_crossvals, nr_folds, nr_crossvals_inner, nr_folds_inner, 'EXP05InnerLoop')

    #Execute calculations and store results
    #1) Locally process tokens
    doneTokens = []
    for token in tokens:
        #db may be a couchDB db
        doneTokens.append(localProcess(token, db = None))
    #2) For a couchDB use the pipeline script

    #Write results to an sql database
    sqlName = "_"+datetime.datetime.now().strftime('%b-%d-%G')
    sqlPath = "Results/"
    sqlFilename = sqlPath+"EXP05InnerLoop"+sqlName+".sqlite3"
    TokenToSqliteExperiment05InnerLoop(tokens, sqlFilename)
    #If you work with a database like couchDB use getDoneTokens to sort the tokens according to the
    #experiments and safe them in the appropriate files.
    #getDoneTokens(db, sqlPath, sqlName = ""), db must be a couchDB  or a dictionary.

    #For further analysis and plots look into the folder experiments

def EXPNewMethod():
    print """To add a new method you will have to
                1) add the method in SetUpGrid once in the CombineDataExperiment* functions
                2) give instructions on how to inititalise the method in SetUpGrid SetUpRun
                3) if the method needs some exotic parameters you might have to write a special case in
                    SetUpGrid RunInstance, eg. after 'if featureSelector.productName ... ' (line 233) you can check for the 
                    featureselctor.name and make your own case.
                4) execute only tokens with your method
                    --> selectedTokens = [token for token in tokens if "newMethod" in token['_id']]
          """


if __name__ == "__main__":
    import sys
    print "Execute the pipeline for the implemented methods "
    print """
            Call by ExecuteExperiments 
            DoEXP01outer
            DoEXP01inner
            DoEXP04outer
            DoEXP04inner
            DoEXP05outer
            DoEXP05inner
            EXPNewMethod
          """
    print "Results will be stores as sql data bases in the folder \"Results\". Make sure it exists."
    command= " ".join( sys.argv[1:] )
    if command == "DoEXP01outer":
        DoEXP01outer()
    elif command == "DoEXP01inner":
        DoEXP01inner()
    elif command == "DoEXP04outer":
        DoEXP04outer()
    elif command == "DoEXP04inner":
        DoEXP04inner()
    elif command == "DoEXP05outer":
        DoEXP05outer()
    elif command == "DoEXP05inner":
        DoEXP05inner()
    elif command == "EXPNewMethod":
        EXPNewMethod()
    else:
        print "Not a valid option"
