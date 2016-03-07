#include "Pattern.h"

typedef set<short> PATTERN_ID_TYPE;
typedef map<PATTERN_ID_TYPE*,Pattern*,PatternComparator> PATTERN_MAP;

extern map<short,string*> INDEX_DIMENSION_MAP;
extern float MAX_ATTRIBUTE_VALUE;
extern float MIN_ATTRIBUTE_VALUE;
extern short NUMBER_OF_ATTRIBUTES;
extern float MIN_DENSITY;
extern short MIN_DIM;
extern float MISSING_VALUE;
extern Node** node_array;
extern float ** weight_array;
extern int DISTANCE_FUNCTION;

int Pattern::counter=0;


bool Pattern::operator < (const Pattern& refParam) const
{
	return (this->getSize() < refParam.getSize());
}


float Pattern::getScore() const
{
	return (((float)ID.size())*getDensity()*((float)getNumberOfRelevantDimensions()/NUMBER_OF_ATTRIBUTES));
}


void Pattern::calculateWE() 
{
	Node *nodePtr1;
	Node *nodePtr2;
	set<short>::iterator iterator1;
	set<short>::iterator iterator2;

	weight=0.0;
	numberOfEdges = 0;
	for ( iterator1=ID.begin(); iterator1!=ID.end(); ++iterator1)
	{
		nodePtr1=node_array[*iterator1];
		
		iterator2=iterator1;
		iterator2++;

		for(;iterator2 !=ID.end(); ++iterator2)
		{
			nodePtr2=node_array[*iterator2];

			if(nodePtr1->isNeighbor(nodePtr2->getID()))
			{
				//cout<<"\t"<<nodePtr1->getID()<<" - "<<nodePtr2->getID()<<endl;
				weight += weight_array[nodePtr1->getID()][nodePtr2->getID()];
				numberOfEdges++;
			}		
		}
	}
}


Pattern::Pattern(const Node *n1,const Node *n2)
{
	Pattern::counter++;

	assert(n1->getID()<n2->getID());

	//initalize ID
	ID;
	ID.insert(n1->getID());
	ID.insert(n2->getID());
	
	//set other fields
	weight = weight_array[n1->getID()][n2->getID()];
	numberOfEdges = 1;


	//initialize and calculate relevant dimensions
	relevantDim = new bool[NUMBER_OF_ATTRIBUTES];
	for(short i=0;i<NUMBER_OF_ATTRIBUTES;i++)
		relevantDim[i]=true;
	updateRelevantDimensions2();
	//print();
}


//member nodes would be dynamically created by pattern generator
//therefore it is already allocated. relevant dimensions however
//come from the old pattern, therefore we can use its value
Pattern::Pattern(Pattern * motherPatternPtr, short newNodeID, const bool* oldPatternRelDim, float newWeight)
{
	Pattern::counter++;

	//initalize ID
	ID = *(motherPatternPtr->getID());
	ID.insert(newNodeID);

	//initialize other fields
	//calculateWE();
	weight = newWeight;
	//numberOfEdges = newNumberOfEdges;*/


	//initialize and update relevant dimensions relevant dimensions
	relevantDim=new bool[NUMBER_OF_ATTRIBUTES];
	for(int i=0;i<NUMBER_OF_ATTRIBUTES;i++)
		relevantDim[i]=oldPatternRelDim[i];
	updateRelevantDimensions2();
	//print();
}



Pattern::~Pattern(void)
{
	Pattern::counter--;
	delete relevantDim;

}


short Pattern::getSize()const
{
	return (short)ID.size();
}

bool Pattern::isDense() const
{
	return getDensity() >= MIN_DENSITY;
}


float Pattern::getDensity() const
{
     return weight/(ID.size()*(ID.size()-1)/2);
}



bool Pattern::containsNode(short nodeID) const
{
	return ID.find(nodeID)!=ID.end();
}

short Pattern::getNumberOfRelevantDimensions() const
{
	return numberOfRelevantDims;
}

bool Pattern::isHomogenous() const
{
	return numberOfRelevantDims>=MIN_DIM;
}

//this methods returns a map with key as the id of the neighbor
//and the needed weight
map<short, float>* Pattern::getExtensibleNeighbors() const
{
	map<short, float>* extNeighbors = new map<short, float>(); // deallocation is done
	map<short, float>* validExtNeighbors = new map<short, float>(); //deallocation id done in Utility.cpp
	map<short, float>::iterator it;
	Node* currentNodePtr;
	set<short>* neighborsPtr;
	float minWeightRequired = static_cast<float> ((ID.size()/2)*(MIN_DENSITY*(ID.size()+1)- getDensity()* (ID.size()-1)));
	
	short neighborID;

	for (set<short>::const_iterator memberIter=ID.begin(); memberIter!=ID.end(); ++memberIter)
	{
		//cout<<"Node:\n\t"<<*memberIter<<endl;
		currentNodePtr=node_array[*memberIter];
		neighborsPtr = currentNodePtr->getNeigbors();
		
		// add the neighbors to the external neighbors map
		for (set<short>::iterator neighIter=neighborsPtr->begin(); neighIter!=neighborsPtr->end(); ++neighIter)
		{
			neighborID=*neighIter;
			//cout<<neighborID<<"-";
			if(!containsNode(neighborID))
			{
				it=extNeighbors->find(neighborID);
				if(it==extNeighbors->end())
				{
					(*extNeighbors)[neighborID] = weight_array[currentNodePtr->getID()][neighborID];
				}
				else
				{
					it->second+= weight_array[currentNodePtr->getID()][neighborID];
				}
			}
		}
		//cout<<endl;
	}

	//---- RECEP ---
	//remove nodes with too low weight to the current node
	//while( !extNeighbors->empty() ) 
	//{
	//	//cout << (*extNeighbors->begin()).first << "->" << (*extNeighbors->begin()).second <<"\t";
	//	if((*extNeighbors->begin()).second>=minWeightRequired)
	//	{
	//		(*validExtNeighbors)[(*extNeighbors->begin()).first]=(*extNeighbors->begin()).second;
	//		//cout<<"\t\tvalid extension!"<<endl;
	//	}
	//	extNeighbors->erase( extNeighbors->begin());
	//}

	//delete extNeighbors;
	return extNeighbors;
}


void Pattern::updateRelevantDimensions2()
{
	float currentValue;
	float prevValue;
	numberOfRelevantDims=0;
	int counter=0;
	bool isRelevant;

	for(int dim=0;dim<NUMBER_OF_ATTRIBUTES;dim++)
	{

	
		if(relevantDim[dim])//incremental update on relevant dimensions
		{

			isRelevant=true;
			
			
			currentValue=-999999;
			
			// go over each member
			for (set<short>::iterator iter=ID.begin(); iter!=ID.end(); ++iter)
			{
				prevValue=currentValue;

				currentValue = node_array[*iter]->getAttributes()[dim];
				

				if((DISTANCE_FUNCTION==1 && currentValue!=1)
					||(DISTANCE_FUNCTION==2 && currentValue==0)
					|| (DISTANCE_FUNCTION==2 && prevValue!=-999999 && currentValue!=prevValue))
				{
					isRelevant=false;
				}
				
			}

			// update relavant dimensions
			if(!isRelevant)
			{
				relevantDim[dim]=false;
			}
			else
			{
				numberOfRelevantDims++;
			}
		}

	}
}

//void Pattern::updateRelevantDimensions()
//{
//	float currentValue;
//	float max;
//	float min;
//	numberOfRelevantDims=0;
//	int nonMissingCounter=0;
//
//	for(int dim=0;dim<NUMBER_OF_ATTRIBUTES;dim++)
//	{
//
//
//		if(relevantDim[dim])//incremental update on relevant dimensions
//		{
//
//
//			// go over each member
//			for (set<short>::iterator iter=ID.begin(); iter!=ID.end(); ++iter)
//			{
//				currentValue = node_array[*iter]->getAttributes()[dim];
//
//				if(currentValue==MISSING_VALUE)
//				{
//					relevantDim[dim]=false;
//				}
//				else
//				{
//					nonMissingCounter++;
//
//					if(DISTANCE_FUNCTION==1)// over expression
//					{
//						if(currentValue!=1)
//							relevantDim[dim]=false;	
//					}
//					else // differential expression
//					{
//						if(currentValue!=1 && currentValue!=-1)
//							relevantDim[dim]=false;
//					}
//				}
//
//			}
//
//			if(relevantDim[dim]==true)
//			{
//				numberOfRelevantDims++;
//			}
//		}
//
//	}
//}


// finds the overlap nodes between two patterns
set<short> * Pattern::getOverlappingNodes(Pattern* otherPattern) const
{
	set<short> * commonNodes = new set<short>();
	for (set<short>::const_iterator elementsIter=ID.begin(); elementsIter!=ID.end(); ++elementsIter)
	{
		if(otherPattern->containsNode(*elementsIter))
		{
			commonNodes->insert(*elementsIter);
		}
	}	
	return commonNodes; //TODO delete not called yet
}

short Pattern::getNumberOfOverlappingNodes(Pattern* otherPattern) const
{
	//TODO this can be done in linear time
	short counter=0;
	for (set<short>::const_iterator elementsIter=ID.begin(); elementsIter!=ID.end(); ++elementsIter)
	{
		if(otherPattern->containsNode(*elementsIter))
		{
			counter++;
		}
	}	
	return counter;
}

short Pattern::getNumberOfCommonRelevantDimensions(Pattern* otherPattern) const
{
	int counter=0;
	bool * otherPatternRelDimensions=otherPattern->getRelevantDim();
	for(int dim=0;dim<NUMBER_OF_ATTRIBUTES;dim++)
	{
		if(relevantDim[dim] && otherPatternRelDimensions[dim])
		{
			counter++;
		}
	}
	return counter;
}


bool Pattern::isSuperPatternOf(Pattern* rhsPtr) const
{
	assert(getSize() >= rhsPtr->getSize());
	
	set<short>* rhsIDPtr = rhsPtr->getID();
	set<short>::const_iterator iter=rhsIDPtr->begin();

	while(iter!= rhsIDPtr->end())
	{
		
		if(containsNode(*iter))
		{
			iter++;
		}
		else
		{
			return false;
		}
	}

	return true;
}


void Pattern::print() const
{

	for (set<short>::const_iterator elementsIter=ID.begin(); elementsIter!=ID.end(); ++elementsIter)
	{
		cout<<node_array[*elementsIter]->getID()<<" ";
	}

	//cout<<"Number of edges = "<< numberOfEdges << " Weight = " << weight << endl;
	cout<<"Weight = " << weight << endl;

	cout<<endl;
	//Flaviacout<<"\tScore: "<<getScore()<<"\t Size: "<<ID.size()<<"\tDensity:"<<getDensity()<<" \tDimensionality:"<<getNumberOfRelevantDimensions()<<endl;	
	cout<<"\t Size: "<<ID.size()<<"\tWeighted Density:"<<getDensity() <<endl;	
	
	cout<<"\t";
	
	for(int dim=0;dim<NUMBER_OF_ATTRIBUTES;dim++)
	{
		if(relevantDim[dim])
		{
			cout<<dim<<"-";
		}
	}
	cout<<endl;
	cout<<endl;
}


//void Pattern::print() const
//{
//	cout<<"\t";
//	for (set<short>::const_iterator elementsIter=ID.begin(); elementsIter!=ID.end(); ++elementsIter)
//	{
//		//cout<<node_array[*elementsIter]->getID()<<" ";
//		cout<<node_array[*elementsIter]->getName()<<" ";
//	}
//
//	//cout<<"Number of edges = "<< numberOfEdges << " Weight = " << weight << endl;
//	//cout<<"Weight = " << weight << endl;
//
//	cout<<endl;
//	cout<<"\t\t\tSize: "<<ID.size()<<"\tWeight:"<<weight<<"\tWeighted Density:"<<getDensity()<<" \tDimensionality:"<<getNumberOfRelevantDimensions()<<endl;	
//	
//	
//	cout<<"\t\t\tDimensions:";
//	
//	for(int dim=0;dim<NUMBER_OF_ATTRIBUTES;dim++)
//	{
//		if(relevantDim[dim])
//		{
//			cout<<dim<<"-";
//		}
//	}
//	cout<<endl;
//	cout<<endl;
//}
