#include "ConnectedComponent.h"

typedef map<set<short>*,Pattern*,PatternComparator> PATTERN_MAP;
typedef map<set<short>*,SimplePath*,PatternComparator> PATH_MAP;

extern time_t START_TIME;
extern float MAX_ATTRIBUTE_VALUE;
extern float MIN_ATTRIBUTE_VALUE;
extern short NUMBER_OF_NODES;
extern short NUMBER_OF_ATTRIBUTES;
extern short MIN_DIM;
extern short MIN_SIZE;
extern float MIN_DENSITY;
extern Node** node_array;
extern map<short,double> SIZE_MIN_DENSITY_MAP;
extern map<short,int> SIZE_MIN_DIMENSIONALITY_MAP;

extern float ** weight_array;



void ConnectedComponent::removeRedundancyInMaximalPatterns()
{
	//cout<<"Removing redundancy in maximal patterns!"<<endl;
	PATTERN_MAP::iterator patternIter1;
	PATTERN_MAP::iterator patternIter2;
	PATTERN_MAP * non_redundant_maximal_patterns = new PATTERN_MAP();
	Pattern* ptr1;
	Pattern* ptr2;

	for (patternIter1=maximal_patterns->begin();patternIter1 != maximal_patterns->end();++patternIter1)
	{
		ptr1 = patternIter1->second;

		if(ptr1!=0) // if it is not deleted
		{

			patternIter2=patternIter1;
			patternIter2++;

			for (;patternIter2 != maximal_patterns->end();++patternIter2)
			{
				ptr2=patternIter2->second;

				if(ptr2!=0) // it may have been deleted before
				{

					if(ptr1->getSize()> ptr2->getSize())
					{
						if(ptr1->isSuperPatternOf(ptr2))
						{
							//cout<<"\nFound a non maximal:"<<endl;
							//ptr1->print();
							//cout<<"Is parent of:";
							//ptr2->print();
							delete ptr2;
							patternIter2->second=0;
						}
					}
					else
					{
						if(ptr1->getSize()< ptr2->getSize() && ptr2->isSuperPatternOf(ptr1))
						{
							//cout<<"\nFound a non maximal:"<<endl;
							//ptr2->print();
							//cout<<"Is parent of:";
							//ptr1->print();
							delete ptr1;
							ptr1=0;
							break;
						}
					}
				}
			}
		}

		if(ptr1!=0)
			(*non_redundant_maximal_patterns)[ptr1->getID()]=ptr1;
	}

	maximal_patterns->clear();
	maximal_patterns=non_redundant_maximal_patterns;
}

PATH_MAP * ConnectedComponent::findAllPaths(set<short>* startNodes,
									  set<short>* endNodes) const
{
	PATH_MAP * currentPaths = new PATH_MAP();
	PATH_MAP * extentedPaths = new PATH_MAP();
	PATH_MAP * maximal_paths = new PATH_MAP();
	PATH_MAP * tempPtr;
	PATH_MAP::iterator pathIter;
	set<short>::iterator nodeIter;
	set<short>* extNeighbors;
    SimplePath * curPathPtr;
	SimplePath * extentedPathPtr;
	SimplePath * cyclicPathPtr;

	// generate seed paths
	for (set<short>::iterator iter=startNodes->begin(); iter!=startNodes->end(); ++iter)
	{
		curPathPtr = new SimplePath(*iter);

		if(endNodes->find(*iter)!=endNodes->end())
		{
			//this node is both start and end node, so this is a maximal path
			(*maximal_paths)[curPathPtr->getElements()]=curPathPtr;

		}
		else
		{

			if(currentPaths->find(curPathPtr->getElements())==currentPaths->end())
			{
				(*currentPaths)[curPathPtr->getElements()]=curPathPtr;
			}
			else
			{
				delete curPathPtr;
			}
		}
	}

	//extend current paths
	while(currentPaths->size()>0)
	{
		//iterate over each path to extend
		for (pathIter=currentPaths->begin();pathIter != currentPaths->end();pathIter++)
		{
			curPathPtr=pathIter->second;
			extNeighbors = curPathPtr->getExtensibleNeighbors(startNodes,endNodes);

			// extend over each possible extension neighbor
			for (nodeIter=extNeighbors->begin();nodeIter != extNeighbors->end();nodeIter++)
			{
				extentedPathPtr = new SimplePath(curPathPtr,*nodeIter);

				if(endNodes->find(*nodeIter)!=endNodes->end())
				{
					(*maximal_paths)[extentedPathPtr->getElements()]=extentedPathPtr; // this is a maximal simple path
				}
				else
				{
					if(extentedPaths->find(extentedPathPtr->getElements())==extentedPaths->end())
					{
						(*extentedPaths)[extentedPathPtr->getElements()]=extentedPathPtr;
					}
					else
					{
						cyclicPathPtr = extentedPaths->find(extentedPathPtr->getElements())->second; // this path contains a cycle
						extentedPaths->erase(extentedPathPtr->getElements());
						delete extentedPathPtr;
						delete cyclicPathPtr;
					}
				}
			}

			delete curPathPtr;
			delete extNeighbors;
		}

		//swap the maps
		tempPtr = currentPaths;
		currentPaths->clear();
		currentPaths=extentedPaths;
		extentedPaths = tempPtr;
	}

	delete extentedPaths;
	delete currentPaths;

	return maximal_paths;
}

void ConnectedComponent::findPathStartAndEndNodes(const Pattern* islandPtr1,
									  const Pattern* islandPtr2,
									  const set<short>* overlapNodes,
									  set<short>* & startNodes,
									  set<short>* & endNodes)const
{
	short overlapNodeID;
	set<short>* neighborsPtr;
	startNodes = new set<short>();
	endNodes = new set<short>();

	for (set<short>::const_iterator iter=overlapNodes->begin(); iter!=overlapNodes->end(); ++iter)
	{
		overlapNodeID=*iter;
		if(islandPtr1->containsNode(overlapNodeID) || islandPtr2->containsNode(overlapNodeID)) // musn't be in the islands
			continue;

		neighborsPtr = node_array[overlapNodeID]->getNeigbors();
		for (set<short>::const_iterator neighIter=neighborsPtr->begin(); neighIter!=neighborsPtr->end(); ++neighIter)
		{
			if(islandPtr1->containsNode(*neighIter))
				startNodes->insert(overlapNodeID);

			if(islandPtr2->containsNode(*neighIter))
				endNodes->insert(overlapNodeID);
		}

	}
}


void ConnectedComponent::expand_by_one(bool firstPhase)
{
	PATTERN_MAP::iterator patternIter;
	Pattern* currentPatternPtr;
	Pattern* candidatePatternPtr;
	map<short, float>* extNeighbors;
	bool isCurrentPatterMaximal;
	short extensionNodeId;
	float weightOfExtensionEdges;
	float currentPatternWeight;
	short potential_critical_node;

	int numberOfProcessedPatterns=0;
	int numberOfGeneratedCandidates=0;
	int CURRENT_LEVEL=2;

	while(current_pattern_map->size()>0)
	{
		cout<<"---------------------\n\n\n";
		cout<<"Current level:"<<CURRENT_LEVEL<<endl;
		cout<<flush;
		CURRENT_LEVEL++;

		numberOfProcessedPatterns=0;
		numberOfGeneratedCandidates=0;

		//iterate over every single pattern


        //patternIter=current_pattern_map->begin();set<short>*,Pattern*,PatternComparator
        //#pragma omp parallel for
		//for(int i = 1; (i < 100000000)&&(patternIter != current_pattern_map->end()); i++)
		//if (patternIter != current_pattern_map->end())
		for (patternIter=current_pattern_map->begin();patternIter != current_pattern_map->end();patternIter++)
		{
			numberOfProcessedPatterns++;
			currentPatternPtr=patternIter->second;
			currentPatternWeight = currentPatternPtr->getWeight();
			//cout<<"\nExtending pattern: \n";
			currentPatternPtr->print();
			isCurrentPatterMaximal=true;

			//get extensible neighbors
			extNeighbors = currentPatternPtr->getExtensibleNeighbors();
			//generate candidates
			while( !extNeighbors->empty() )
			{
				//cout <<"\t Node:" <<node_array[(*extNeighbors->begin()).first]->getName() << "->" << (*extNeighbors->begin()).second <<"\t";

				//generateID
				numberOfGeneratedCandidates++;
				extensionNodeId = (*extNeighbors->begin()).first;
				weightOfExtensionEdges = (*extNeighbors->begin()).second;

				candidatePatternPtr = new Pattern(currentPatternPtr,extensionNodeId,currentPatternPtr->getRelevantDim(),weightOfExtensionEdges+currentPatternWeight);

				//candidatePatternPtr->print();

				if(candidatePatternPtr->isHomogenous())
				{
					//cout<<"\t\t\tHomogenous!"<<endl;

					if(candidatePatternPtr->getDensity()>=MIN_DENSITY)
					{
						isCurrentPatterMaximal=false;
						//cout<<"\t\t\t\t and dense!"<<endl;
						if(valid_candidate_map->find(candidatePatternPtr->getID())==valid_candidate_map->end())
						{
							(*valid_candidate_map)[candidatePatternPtr->getID()]=candidatePatternPtr;

						}
						else
						{
							delete candidatePatternPtr;
						}
						
					}
					else
					{
						//cout<<"\t\t\t\tbut not dense!"<<endl;
					}
				}
				else
				{
					//cout<<"\t\t\tNot homogenous!"<<endl;
					delete candidatePatternPtr;
				}
				extNeighbors->erase( extNeighbors->begin());
			}

			delete extNeighbors;

			if(isCurrentPatterMaximal && currentPatternPtr->getSize() >= MIN_SIZE
				  && isPatternModule(currentPatternPtr->getSize(),currentPatternPtr->getDensity(),currentPatternPtr->getNumberOfRelevantDimensions()))
			{
				//cout<<"\t\tmaximal!"<<endl;
				//(*maximal_patterns)[currentPatternPtr->getID()]=currentPatternPtr;
				if(firstPhase)
					temp_maximal_patterns->insert(pair<short,Pattern*>(currentPatternPtr->getSize(),currentPatternPtr));
				else
					(*maximal_patterns)[currentPatternPtr->getID()]=currentPatternPtr;

			}
			else // this guy has been extented, not maximal
			{
				delete currentPatternPtr;
			}
			//cout<<"\n*"<<endl;

			if(numberOfProcessedPatterns%100==0)
			{
				cout<<"\n\n\tNumber of processed patterns:"<<numberOfProcessedPatterns<<endl;
				cout<<"\tNumber of generated candidates patterns:"<<numberOfGeneratedCandidates<<endl;
				cout<<"\tNumber of valid candidates: "<<valid_candidate_map->size()<<endl;
			}
			//patternIter != current_pattern_map->end();
		}
		//else
            //break;
		cout<<"\tLevel: "<< CURRENT_LEVEL <<endl;
		cout<<"\tRuntime: "<<time(NULL)-START_TIME<<endl;
		cout<<"\tNumber of generated candidates patterns:"<<numberOfGeneratedCandidates<<endl;
		cout<<"\tNumber of valid candidates: "<<valid_candidate_map->size()<<endl;

		cleanMemoryAfterExpansion(CURRENT_LEVEL);
	}

	if(firstPhase)
		cout<<"\tNumber of potential critical patterns:"<<potential_critical_patterns->size()<<endl;

}

bool ConnectedComponent::isPatternModule(short size, double density, short dimensionality) const
{
	return density>=MIN_DENSITY;

}

void ConnectedComponent::cleanMemoryAfterExpansion(int currentLevel)
{
	PATTERN_MAP * temp_map;
	list<Pattern*>::iterator patternIter1;

	//cout<<"\n\nComponent:"<<this->getID()<<endl;
	//cout<<"\tCurrent level : "<<currentLevel<<endl;
	//cout<<"\tCurrent elapsed time: "<<((float)time(NULL)-START_TIME)/60<<endl;
	//cout<<"\t\tBefore expansion "<< current_pattern_map->size() << "  patterns"<<endl;
	//cout<<"\t\tDuring expansion "<<valid_candidate_map->size()<<"   valid patterns are generated"<<endl;
	//cout<<"\t\tUp to now "<<maximal_patterns->size()<<" maximal patterns are generated"<<endl;


	//clear current map
	current_pattern_map->clear();

	temp_map=current_pattern_map;
	current_pattern_map=valid_candidate_map;
	valid_candidate_map=temp_map;

	//cout<<"\t\tPattern counter: "<<Pattern::counter<<endl;
	//cout<<"\t\tMaximal patterns: "<<maximal_patterns->size()<<endl;
}

void ConnectedComponent::findSeedPatterns()
{
	float MIN_EDGE = 1 ;// (MIN_SIZE-1)*MIN_DENSITY;

	Node *currentNodePtr;
	map<short,Node*>::const_iterator iter= nodes.begin();
	short neighborID;
	set<short>* neighborsPtr;
	Pattern *seedPtr;
	int seedCounter=0;
	//cout<<"Seed generation for component "<<id<<" is started!"<<endl;

	for (int i=0;i<getSize();i++)
	{
        //cout <<iter->first<<" ";
		currentNodePtr=iter->second;
		neighborsPtr=currentNodePtr->getNeigbors();

		// add the neighbors to the queue
		for (set<short>::iterator it=neighborsPtr->begin(); it!=neighborsPtr->end(); ++it)
		{
			neighborID=*it;
			if(currentNodePtr->getID()<neighborID && weight_array[currentNodePtr->getID()][neighborID]>=MIN_DENSITY)
			{
				//cout<<"\n-------\t"<<currentNodePtr->getID()<<" - "<<neighborID<<endl;
				seedCounter++;
				seedPtr = new Pattern(currentNodePtr,node_array[neighborID]);
				(*current_pattern_map)[seedPtr->getID()]=seedPtr;
			}
		}
		iter++;
	}

	cout<<"Number of seeds generated: "<<seedCounter<<"\n------\n"<<endl;
}

void ConnectedComponent::printMaximalPatterns()
{
	multiset<Pattern> sortedPatterns;
	multiset<Pattern>::iterator it;
	cout<<"Maximal patterns: "<<endl;
	PATTERN_MAP::iterator patternIter1;


	for (patternIter1=maximal_patterns->begin();patternIter1 != maximal_patterns->end();++patternIter1)
	{
		sortedPatterns.insert(*(patternIter1->second));
	}

	for (it=sortedPatterns.begin(); it!=sortedPatterns.end(); it++)
	{
		(*it).print();
	}

}



int ConnectedComponent::getID()const
{
	return id;
}

Node * ConnectedComponent::getNode(short nodeID)
{
	return this->nodes[nodeID];
}



short ConnectedComponent::getSize()const
{
	return (short) nodes.size();
}

map<set<short>*,Pattern*,PatternComparator>* ConnectedComponent::getMaximalPatterns()const
{
	return maximal_patterns;
}

multimap<short,Pattern*>* ConnectedComponent::getTempMaximalPatterns()const
{
	return temp_maximal_patterns;
}

multimap<short,Pattern*>*  ConnectedComponent::getPotentialCriticalPatterns()const
{
	return potential_critical_patterns;
}

void ConnectedComponent::addNode(Node * newNodePtr)
{
	this->nodes[newNodePtr->getID()]=newNodePtr;
}

void ConnectedComponent::print()const
{
	map<short,Node*>::const_iterator iter;
	cout<<"\tI am connected component: "<<this->id<<"  (size:"<<this->nodes.size()<<" )";
	cout<<"My members are:\t";
	int counter=1;
	for (iter=nodes.begin(); iter != nodes.end(); ++iter)
	{
        cout <<iter->first<<" ";
		counter++;
		if(counter%20==0)
			cout<<endl;
    }
	cout<<endl;
}


ConnectedComponent::ConnectedComponent(int id)
{
	//cout<<"\tConstructor for connected component "<<id<<" is called"<<endl;
	nodes;
	this->id=id;

	current_pattern_map = new PATTERN_MAP();
    valid_candidate_map = new PATTERN_MAP();
	maximal_patterns = new PATTERN_MAP();
	potential_critical_patterns = new multimap<short,Pattern*>();
	temp_maximal_patterns = new multimap<short,Pattern*>();
}


ConnectedComponent::~ConnectedComponent(void)
{
	delete current_pattern_map;
	delete valid_candidate_map;
	delete maximal_patterns;
	delete temp_maximal_patterns;
	//cout<<"\tDestructor for connected component "<<id<<" is called"<<endl;
}
