all: NCSC runwDCB createGraph createVector rankModules
 
NCSC : NCSC.cpp
		    g++ NCSC.cpp Loader.cpp Utility.cpp MultiHashIndex.cpp ConnectedComponent.cpp Node.cpp Pattern.cpp PatternComparator.cpp PatternID.cpp SimplePath.cpp -o NCSC

runwDCB : runwDCB.cpp
		    g++ runwDCB.cpp -o runwDCB

createGraph : createGraph.cpp
		    g++ createGraph.cpp -o createGraph

createVector : createVector.cpp
		    g++ createVector.cpp -o createVector

rankModules : rankModules.cpp
		    g++ -lRmath rankModules.cpp -o rankModules
		#for medoid g++ rankModules.cpp -o rankModules -lRmath
clean : 
		    rm NCSC runwDCB createVector rankModules *.o

