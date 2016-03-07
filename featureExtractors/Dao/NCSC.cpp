#include "NCSC.h"
#include "Node.h"
#include "Loader.h"
#include "Utility.h"
#include "ConnectedComponent.h"
#include "Pattern.h"
#include "PatternComparator.h"

map<short,string*> INDEX_DIMENSION_MAP;
set<short> UNIQE_GENES_IN_PATTERNS;
time_t START_TIME=0;
float MAX_ATTRIBUTE_VALUE=10000;
float MIN_ATTRIBUTE_VALUE=-10000;
float MISSING_VALUE=1111111;
int NUMBER_OF_EDGES=0;





//algorithm parameters
short NUMBER_OF_NODES=-1;
short NUMBER_OF_ATTRIBUTES=-1;
float MIN_DENSITY=-1;
short MIN_DIM=-1;
short MIN_SIZE=-1;
int DISTANCE_FUNCTION; //1- for over expression 2-for differential expression


// parallel proc parameters
int NUMBER_OF_PARTITIONS;
int NUMBER_OF_CPUs;

// custom size parameters
map<short,double> SIZE_MIN_DENSITY_MAP;
map<short,int> SIZE_MIN_DIMENSIONALITY_MAP;
short MAX_PATTERN_SIZE=-1;

Node** node_array;
vector<ConnectedComponent*> connectedComponents;


int main(int argc, char *argv[])
{
	Utility util;
	START_TIME=time (NULL);
    	string str;

	// PARSE PARAMETERS
	MIN_DENSITY=atof(argv[4]);
	MIN_DIM=atoi(argv[6]);
	MIN_SIZE=atoi(argv[8]);
	DISTANCE_FUNCTION=atoi(argv[10]);
	cout<<"Minimum density        : "<<MIN_DENSITY<<endl;
	cout<<"Minimum dim            : "<<MIN_DIM<<endl;
	cout<<"Minimum pattern size   : "<<MIN_SIZE<<endl;

	// REDIRECT IO
	std::ostringstream converter;
	converter<<argv[2]<<"experimental_results/alpha_"<<MIN_DENSITY<<"_minGraph_"<<MIN_DIM<<"_minSize_"<<MIN_SIZE<<"_distFunc_"<<DISTANCE_FUNCTION<<".txt";
	cout<<"Redirecting output to:"<<converter.str()<<endl;
	freopen (converter.str().c_str(),"w",stdout);

	// LOAD
	Loader loader;
	loader.loadDataset(argv[2]);

	// PREPROCESS  INPUT
	util.removeNonHomegenousEdges();
	//util.partitionGraph(NUMBER_OF_PARTITIONS);
	util.extractComponents();  // only in windows
	vector<ConnectedComponent*>::const_iterator iter=connectedComponents.begin();
	int numberOfComponents = (int)connectedComponents.size();
	for (int i=0; i<numberOfComponents;i++)
	{
		ConnectedComponent* ccPtr = connectedComponents[i];
		cout<<"\tComponent "<<i<< "  size:"<< ccPtr->getSize()<<endl;
	}
	cout<<"\n++++++++++++++++++++++++++++++++\n\n\n"<<endl;
	//cout<<"Number of processors:"<<omp_get_num_procs()<<endl;
	cout<<"Number of processors used:"<<NUMBER_OF_CPUs<<endl;
	//omp_set_num_threads(NUMBER_OF_CPUs);
	cout<<"\n\nStarting algorithm!"<<endl;

	//1- FIRST PHASE
	//#pragma omp parallel for
	for (int i=0; i<numberOfComponents;i++)
	{
		ConnectedComponent* ccPtr = connectedComponents[i];
		cout<<"\n\n\n\n======================\n\n\tStarting Connected Component: "<<ccPtr->getID()<<"\tSize: "<<ccPtr->getSize()<<endl;
		cout<<"Find seed patterns!"<<endl;
		ccPtr->findSeedPatterns();
		cout<<"Expand-by-one!"<<endl;
		ccPtr->expand_by_one(true);
		cout<<"\t\tFinished: "<<ccPtr->getID()<<"\tSize: "<<ccPtr->getSize()<<"  # of Found Patterns:"<<ccPtr->getTempMaximalPatterns()->size()<<endl;
	}
	cout<<"\n--------------------------------\n\n\n"<<endl;
	cout<<"\n\n\n\nRun time before merge (in minutes): "<<((float)time(NULL)-START_TIME)/60<<" "<<endl;
	cout<<"Run time before merge (in second): "<<(time(NULL)-START_TIME)<<" "<<endl;

	// all patterns will be merged to this guy
	ConnectedComponent* masterComponent = new ConnectedComponent(-1);


	// statistics
	cout<<"\n\n"<<endl;
	int patternCounter=0;
	for (int i=0; i<numberOfComponents;i++)
	{
		ConnectedComponent* ccPtr = connectedComponents[i];
		//cout<<"\tBefore Component "<<i<< " maximal pattern size:"<< ccPtr->getMaximalPatterns()->size()<<endl;
		//ccPtr->removeRedundancyInMaximalPatterns();
		patternCounter+=(int)ccPtr->getTempMaximalPatterns()->size();
		//cout<<"\tAfter Component "<<i<< " maximal pattern size:"<< ccPtr->getMaximalPatterns()->size()<<endl;
	}
	cout<<"Before redundancy elim, total patterns in components:"<<patternCounter<<endl;

	cout<<"\n\n\n\nRun time before redundancy (in minutes): "<<((float)time(NULL)-START_TIME)/60<<" "<<endl;
	cout<<"Run time before redundancy (in second): "<<(time(NULL)-START_TIME)<<" "<<endl;

	//3- REDUNDANCY ELIMINATION
	util.unifyPatternsOfAllComponents2(masterComponent);
	cout<<"After removing redundancy number of patterns:"<<masterComponent->getMaximalPatterns()->size()<<endl;


	// statistics
	cout<<"\n\n\n\n Final number of patterns:"<<masterComponent->getMaximalPatterns()->size()<<endl;

	masterComponent->printMaximalPatterns();


	//CLEAN MEMORY
	util.deallocateConnectedComponents();
	util.deallocateNodes();
	//cout<<"Pattern counter:"<<Pattern::counter<<endl;
	cout<<"Total run time(minutes): "<<((float)time(NULL)-START_TIME)/60<<" "<<endl;
	cout<<"Total run time(seconds): "<<(time(NULL)-START_TIME)<<" "<<endl;
	return 0;
}
