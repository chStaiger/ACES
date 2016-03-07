#include "Utility.h"


const short NO_SUCH_NODE_FLAG = -1;

extern int DISTANCE_FUNCTION;
extern int NUMBER_OF_EDGES;
extern float MAX_ATTRIBUTE_VALUE;
extern float MIN_ATTRIBUTE_VALUE;
extern short NUMBER_OF_NODES;
extern short NUMBER_OF_ATTRIBUTES;
extern short MIN_DIM;
extern short MIN_SIZE;
extern short MAX_PATTERN_SIZE;
extern float MISSING_VALUE;
extern float MIN_DENSITY;
extern Node** node_array;
extern vector<ConnectedComponent *> connectedComponents;
extern set<short> UNIQE_GENES_IN_PATTERNS;

int * Utility::get_Node_PatternSetSize_Array() const
{
	int numberOfComponents = (int)connectedComponents.size();
	Pattern * patternPtr;
	map<set<short>*,Pattern*,PatternComparator>::iterator patternIter;
	multimap<short,Pattern*>* temp_maximal_patterns;
	set<short>::const_iterator elementsIter;
	ConnectedComponent* ccPtr;
	int* node_patSetSize = new int[NUMBER_OF_NODES];
	short sizeOfPattern=-1;


	set<short>* patternID;
	short nodeID=-1;

	// initalize everything to zero
	for(int j=0;j<NUMBER_OF_NODES;j++)
		node_patSetSize[j]=0;

	// component
	for (int i=0; i<numberOfComponents;i++)
	{
		ccPtr = connectedComponents[i];
		temp_maximal_patterns = ccPtr->getTempMaximalPatterns();
		//pattern
		for(multimap<short,Pattern*>::iterator iter = temp_maximal_patterns->begin(); iter != temp_maximal_patterns->end(); ++iter )
		{
			patternPtr=iter->second;
			sizeOfPattern=iter->first;
			patternID=patternPtr->getID();

			if(sizeOfPattern>MAX_PATTERN_SIZE)
				MAX_PATTERN_SIZE=sizeOfPattern;

			// node
			for (elementsIter=patternID->begin(); elementsIter!=patternID->end(); ++elementsIter)
			{
				node_patSetSize[*elementsIter]++;
			}
		}
	}
	cout<<"Maximum pattern size "<< MAX_PATTERN_SIZE<<endl;
	return node_patSetSize;
}

void Utility::findUniqueGenesInPatterns(ConnectedComponent* ccPtr) const
{
	map<set<short>*,Pattern*,PatternComparator>::iterator patternIter;
	map<set<short>*,Pattern*,PatternComparator>* maximal_patterns;
	Pattern* patternPtr;

	maximal_patterns = ccPtr->getMaximalPatterns();
	for (patternIter=maximal_patterns->begin();patternIter != maximal_patterns->end();++patternIter)
	{
		patternPtr=patternIter->second;

		for (set<short>::const_iterator elementsIter=patternPtr->getID()->begin(); elementsIter!=patternPtr->getID()->end(); ++elementsIter)
		{
			UNIQE_GENES_IN_PATTERNS.insert(*elementsIter);
		}
	}
}


void Utility::unifyCriticalPatternsOfAllComponents(ConnectedComponent* ccPtr) const
{
	multimap<short,Pattern*>::iterator it1,itlow,itup;
	multimap<short,Pattern*>*all_pot_crit_patterns = ccPtr->getPotentialCriticalPatterns();

	int numberOfComponents = (int)connectedComponents.size();

	Pattern* patternPtr;
	short critical_node_id=-1;;

	for (int i=0; i<numberOfComponents;i++)
	{
		ConnectedComponent* ccPtr = connectedComponents[i];
		multimap<short,Pattern*>* pot_crit_patterns = ccPtr->getPotentialCriticalPatterns();




		// for every pair
		for(short id=0;id<NUMBER_OF_NODES;id++)
		{

			itlow=pot_crit_patterns->lower_bound(id);  // itlow points to b
			itup=pot_crit_patterns->upper_bound(id);   // itup points to e (not d)
			for( it1=itlow ; it1 != itup; it1++ )
			{
				critical_node_id = (*it1).first;
				patternPtr=(*it1).second;

				//cout<<"Inserting pair:"<<endl;
				//cout<<node_array[critical_node_id]->getName()<<endl;
				//patternPtr->print();
				all_pot_crit_patterns->insert(pair<short,Pattern*>(critical_node_id,patternPtr));
			}
		}


		//for (patternIter=pot_crit_patterns->begin();patternIter != pot_crit_patterns->end();++patternIter)
		//{
		//	critical_node_id=patternIter->first;
		//	patternPtr=patternIter->second;

		//	cout<<"Inserting pair:"<<endl;
		//	cout<<node_array[critical_node_id]->getName()<<endl;
		//	patternPtr->print();

		//	// insert this critical pattern into global map of critical patterns
		//	all_pot_crit_patterns->insert(pair<short,Pattern*>(critical_node_id,patternPtr));
		//}
	}
	cout<<"Number of potential critical patterns for the merge phase:"<<all_pot_crit_patterns->size()<<endl;
}


void Utility::unifyPatternsOfAllComponents2(ConnectedComponent* masterComPtr) const
{
	multimap<short,Pattern*>::iterator it1,itlow,itup;
	multimap<short,Pattern*>* maximal_patterns2;
	map<set<short>*,Pattern*,PatternComparator>::iterator patternIter;
	map<set<short>*,Pattern*,PatternComparator>*all_patterns = masterComPtr->getMaximalPatterns();
	MultiHashIndex* index = new MultiHashIndex(NUMBER_OF_NODES);
	int numberOfComponents = (int)connectedComponents.size();
	int* node_patternSet_Size = get_Node_PatternSetSize_Array();

	Pattern* patternPtr;

	// STEP-2
	for (int i=0; i<numberOfComponents;i++)
	{
		//cout<<"Component "<<i<<endl;
		ConnectedComponent* ccPtr = connectedComponents[i];
		maximal_patterns2 = ccPtr->getTempMaximalPatterns();

		for(int size=MAX_PATTERN_SIZE;size>=MIN_SIZE;size--)
		{
			itlow=maximal_patterns2->lower_bound(size);  //
			itup=maximal_patterns2->upper_bound(size);   //
			for( it1=itlow ; it1 != itup; it1++ )
			{
				patternPtr=(*it1).second;
				if(index->insertPattern2(patternPtr,node_patternSet_Size)==false)
				{
				//if(index->insertPattern(patternPtr)==false)
				//{
				//		// delete this pattern
						delete patternPtr;
				}

			}
		}


		//map<set<short>*,Pattern*,PatternComparator>* maximal_patterns = ccPtr->getMaximalPatterns();
		//for (patternIter=maximal_patterns->begin();patternIter != maximal_patterns->end();++patternIter)
		//{
		//	patternPtr=patternIter->second;
		//	if(index->insertPattern2(patternPtr,node_patternSet_Size)==false)
		//	{
		//	//if(index->insertPattern(patternPtr)==false)
		//	//{
		//		// delete this pattern
		//		delete patternPtr;
		//	}
		//}
	}

	// fill in the all_patterns with unique patterns
	index->unifyPatterns(all_patterns);
}


void Utility::unifyPatternsOfAllComponents(ConnectedComponent* masterComponent) const
{
	map<set<short>*,Pattern*,PatternComparator>::iterator patternIter;
	map<set<short>*,Pattern*,PatternComparator>*all_patterns = masterComponent->getMaximalPatterns();

	int numberOfComponents = (int)connectedComponents.size();

	Pattern* patternPtr;

	for (int i=0; i<numberOfComponents;i++)
	{
		ConnectedComponent* ccPtr = connectedComponents[i];
		map<set<short>*,Pattern*,PatternComparator>* maximal_patterns = ccPtr->getMaximalPatterns();

		for (patternIter=maximal_patterns->begin();patternIter != maximal_patterns->end();++patternIter)
		{
			patternPtr=patternIter->second;

			if(all_patterns->find(patternPtr->getID())== all_patterns->end())
			{
				(*all_patterns)[patternPtr->getID()]=patternPtr;
			}
			else
			{
				delete patternPtr;
			}
		}
	}
}




void Utility::partitionGraph(int numberOfClusters)const
{
	char command[80],partitionFileName[200];
	std::ostringstream stm;
	stm << numberOfClusters;

	// write the graph
	writeGraphInGraclusFormat();

	// contruct the and run the command
	strcpy (command,"./graclus graclus_graph.txt ");
	strcat (command,stm.str().c_str());
	cout<<"Command :"<<command<<endl;
	system(command);

	// read the partitions
	strcpy (partitionFileName,"graclus_graph.txt");
	strcat (partitionFileName,".part.");
	strcat (partitionFileName,stm.str().c_str());
	cout<<"Reading partiotions from file:" << partitionFileName <<endl;
	readPartitions(partitionFileName, numberOfClusters);
}

void Utility::readPartitions(char* partitionFile, int numberOfPartitions) const
{
	ifstream partFileStream (partitionFile);
	short currentNodeID=0;
	short assignedPartition=-1;

	// create connected components(partitions)
	for(int i=0;i<numberOfPartitions;i++)
	{
		ConnectedComponent * ccPtr = new ConnectedComponent(i);
		connectedComponents.push_back(ccPtr);
	}


	// read partition assignment
	if (partFileStream.is_open())
	{
		while (! partFileStream.eof() && currentNodeID < NUMBER_OF_NODES)
		{
		  partFileStream>>assignedPartition;
		  //cout<<currentNodeID<<"->"<<assignedPartition<<endl;
		  connectedComponents[assignedPartition]->addNode(node_array[currentNodeID]);
		  currentNodeID++;
		}
		partFileStream.close();
		cout<<currentNodeID<<" nodes are loaded from partition file\n"<<endl;
	}
	else
	{
		cout << "Unable to open partition file";
		exit(1);
	}
}

void Utility::writeGraphInGraclusFormat() const
{
	ofstream graphFile;
	graphFile.open ("graclus_graph.txt");
	Node * node1Ptr;
	Node * node2Ptr;
	graphFile<<NUMBER_OF_NODES<<" "<<NUMBER_OF_EDGES<<endl;
	for(int id1=0;id1<NUMBER_OF_NODES;id1++)
	{
		node1Ptr=node_array[id1];

		for (set<short>::iterator it=node1Ptr->getNeigbors()->begin(); it!=node1Ptr->getNeigbors()->end();it++)
		{
			node2Ptr=node_array[*it];
			graphFile<<(node2Ptr->getID()+1)<<" ";
		}
		graphFile<<endl;
	}
	cout<<"Wrote graph in graclust format..!"<<endl;
}


void Utility::deallocateConnectedComponents(void) const
{
	for (int i=0;i<(int)connectedComponents.size();i++)
	{
		delete connectedComponents[i];
    }
}

void Utility::deallocateNodes(void) const
{
	for(int i=0;i<NUMBER_OF_NODES;i++)
	{
		delete node_array[i];
	}
}

void Utility::extractComponents(void) const
{
	cout<<"Extracting connected components."<<endl;
	Node * root;
	short numberOfConnectedComponents=0;
	bool* visitFlag = new bool[NUMBER_OF_NODES]; // deallocation is taken care
	short rootIndex=NO_SUCH_NODE_FLAG;

	for(int i=0;i<NUMBER_OF_NODES;i++)
		visitFlag[i]=false;

	rootIndex=getUnvisitedNode(visitFlag);
	while(rootIndex != NO_SUCH_NODE_FLAG)
	{
		numberOfConnectedComponents++;
		ConnectedComponent * ccPtr = new ConnectedComponent(numberOfConnectedComponents);
		root = node_array[rootIndex];
		extendComponent(ccPtr,root,visitFlag);
		//ccPtr->print();
		if(ccPtr->getSize()>=MIN_SIZE)
		{
			connectedComponents.push_back(ccPtr);
		}
		else
		{
			delete ccPtr;
		}
		rootIndex=getUnvisitedNode(visitFlag);
	}
	delete visitFlag;
	cout<<"Extracted "<<numberOfConnectedComponents<<"  connected components of which "<<connectedComponents.size()<<" are large enough \n\n"<<endl;
}


void Utility::extendComponent(ConnectedComponent* cc, Node * root,bool* visitFlag) const
{
	//cout<<"\tExtending root :"<<root->getID()<<endl;
	queue<int> q;
	set<short>* neighborsPtr;
	short neighborID;
	Node *neighborPtr;
	Node *currentNodePtr;

	visitFlag[root->getID()]=true;
	cc->addNode(root);
	q.push(root->getID());

	while(q.size()!=0)
	{
		int currentID = q.front();
		currentNodePtr = node_array[currentID];
		neighborsPtr = currentNodePtr->getNeigbors();
		q.pop();

		// add the neighbors to the queue
		for (set<short>::iterator it=neighborsPtr->begin(); it!=neighborsPtr->end(); ++it)
		{
			neighborID=*it;
			if(visitFlag[neighborID]!=true)
			{
				visitFlag[neighborID]=true;
				q.push(neighborID);

				neighborPtr=node_array[neighborID];
				cc->addNode(neighborPtr);
			}
		}
	}
}

short Utility::getUnvisitedNode(bool* visitFlag) const
{
	short index=NO_SUCH_NODE_FLAG;;
	for(short i=0;i<NUMBER_OF_NODES;i++)
	{
		if(visitFlag[i]==false)
		{
			index=i;
			break;
		}
	}
	//cout<<"Node :"<<index<<" is the root!"<<endl;
	return index;
}


void Utility::removeNonHomegenousEdges(void) const
{
	cout<<"Processing edges ..."<<endl;
	int deletedEdgeCounter=0;
	Node * node1Ptr;
	Node * node2Ptr;

	for(int id1=0;id1<NUMBER_OF_NODES;id1++)
	{
		node1Ptr=node_array[id1];

		for (set<short>::iterator it=node1Ptr->getNeigbors()->begin(); it!=node1Ptr->getNeigbors()->end();)
		{
			node2Ptr=node_array[*it];
			++it;
			if(node1Ptr->getID()<node2Ptr->getID())
			{
				cout<<"\tPreprocessing edge "<<node1Ptr->getID()<<"  "<<node2Ptr->getID()<<endl;
				if(isEdgeEdgeHomogeneous(node1Ptr,node2Ptr)==false)
				{
					cout<<"\t\tRemoved!!!"<<endl;
					removeEdge(node1Ptr,node2Ptr);
					deletedEdgeCounter++;

				}
			}
		}
	}
	NUMBER_OF_EDGES=NUMBER_OF_EDGES-deletedEdgeCounter;
	cout<<"\t\t"<<deletedEdgeCounter<<" edges are deleted!"<<endl;
	cout<<"\t\t"<<NUMBER_OF_EDGES<<" edges present!"<<endl;
	//cout<<"All edges processed ...\n-----\n\n"<<endl;
}

void Utility::removeEdge(Node *n1, Node *n2) const
{
	//cout<<"\t\tRemoving the edge: "<<n1->getID()<<"\t"<<n2->getID()<<endl;
	n1->removeNeighbor(n2->getID());
	n2->removeNeighbor(n1->getID());
}


bool Utility::isEdgeEdgeHomogeneous(const Node *n1, const Node *n2) const
{
	const float* attr1 = n1->getAttributes();
	const float* attr2 = n2->getAttributes();
	short numberOfValidDimensions=0;
	float difference=0;

	for(int i=0;i<NUMBER_OF_ATTRIBUTES;i++)
	{

		if(attr1[i]==MISSING_VALUE || attr2[i]==MISSING_VALUE)
		{
			continue;
		}
		else
		{
			if(DISTANCE_FUNCTION==1)// over expression
            {
                if(attr1[i]==1 && attr2[i]==1)
                    numberOfValidDimensions++;       
            }
            else // differential expression
            {
                if((attr1[i]==1 && attr2[i]==1) || (attr1[i]==-1 && attr2[i]==-1))
                    numberOfValidDimensions++;
            }

		}
	}

	return (numberOfValidDimensions>=MIN_DIM);
}


