#include "MultiHashIndex.h"
extern Node** node_array;


MultiHashIndex::~MultiHashIndex()
{
	delete node_patterns_map;
	for (int i=0;i<numberOfBuckets;i++)
	{
		delete node_patterns_map[i];
	}
}

MultiHashIndex::MultiHashIndex(int numberOfNodes)
{
	numberOfBuckets=numberOfNodes;
	node_patterns_map = new map<set<short>*,Pattern*,PatternComparator>*[numberOfNodes];
	for (int i=0;i<numberOfBuckets;i++)
	{
		node_patterns_map[i]= new  map<set<short>*,Pattern*,PatternComparator>();
	}
}

// Martin's improvement
bool MultiHashIndex:: insertPattern2(Pattern * newPatptr,int* node_patternSize)
{
	
	set<short>* patternID = newPatptr->getID();
	map<set<short>*,Pattern*,PatternComparator>* bucketPtr;
	map<set<short>*,Pattern*,PatternComparator>* toBeDeleted;
	map<set<short>*,Pattern*,PatternComparator>::iterator bucketIterator;
	Pattern * oldPatptr;

	// find the element with the smalles number of patterns
	int desired_node=-1;
	int smallest=1000000;
	for (set<short>::const_iterator elementsIter=patternID->begin(); elementsIter!=patternID->end(); ++elementsIter)
	{
		if(node_patternSize[*elementsIter]<smallest)
		{
			smallest=node_patternSize[*elementsIter];
			desired_node=*elementsIter;
		}
	}

	bucketPtr = node_patterns_map[desired_node];

	//  for every pattern in this bucket, check for superset-subset relation
	for (bucketIterator=bucketPtr->begin();bucketIterator != bucketPtr->end();++bucketIterator)
	{
		oldPatptr=bucketIterator->second;

		// larger patterns were inserted earlier
		if(newPatptr->getSize()<= oldPatptr->getSize() && oldPatptr->isSuperPatternOf(newPatptr)) 
		{
			//nothing to do, a superset or equal pattern is already added, redundant pattern
			return false;
		}		
	}

	// new pattern, add this to all buckets of its elements
	for (set<short>::const_iterator elementsIter=patternID->begin(); elementsIter!=patternID->end(); ++elementsIter)
	{
		bucketPtr = node_patterns_map[*elementsIter];
		(*bucketPtr)[patternID]=newPatptr;
	}

	return true;
}

bool MultiHashIndex:: insertPattern(Pattern * newPatptr)
{
	
	set<short>* patternID = newPatptr->getID();
	map<set<short>*,Pattern*,PatternComparator>* bucketPtr;
	map<set<short>*,Pattern*,PatternComparator>* toBeDeleted;
	map<set<short>*,Pattern*,PatternComparator>::iterator bucketIterator;
	Pattern * oldPatptr;

	// for every member of the pattern
	for (set<short>::const_iterator elementsIter=patternID->begin(); elementsIter!=patternID->end(); ++elementsIter)
	{
		toBeDeleted = new map<set<short>*,Pattern*,PatternComparator>();
		bucketPtr = node_patterns_map[*elementsIter];

		//  for every pattern in this bucket, check for superset-subset relation
		for (bucketIterator=bucketPtr->begin();bucketIterator != bucketPtr->end();++bucketIterator)
		{
			oldPatptr=bucketIterator->second;

			if(newPatptr->getSize()> oldPatptr->getSize())
			{
				if(newPatptr->isSuperPatternOf(oldPatptr))
				{
					(*toBeDeleted)[oldPatptr->getID()]=oldPatptr;
				}			
			}
			else 
			{
				if(newPatptr->getSize()<= oldPatptr->getSize() && oldPatptr->isSuperPatternOf(newPatptr)) 
				{
					//nothing to do, a superset or equal pattern is already added, redundant pattern
					return false;
				}
			}
			
		}

		// remove subset patterns
		for (bucketIterator=toBeDeleted->begin();bucketIterator != toBeDeleted->end();++bucketIterator)
		{
			bucketPtr->erase(bucketIterator->first);
		}
		delete toBeDeleted;

		// add this pattern
		(*bucketPtr)[newPatptr->getID()]=newPatptr;
	}
	return true;
}


// this function generates a set from the patterns
// that exist in all buckets.
// note that, when this function is called, the
// there should not exist any non-maximal pattern in
// the buckets
void MultiHashIndex::unifyPatterns(map<set<short>*,Pattern*,PatternComparator>* maximal_patterns )
{
	map<set<short>*,Pattern*,PatternComparator>::iterator bucketIterator;
	map<set<short>*,Pattern*,PatternComparator>* bucketPtr;
	Pattern * oldPatptr;
	for (int i=0;i<numberOfBuckets;i++)
	{
		bucketPtr = node_patterns_map[i];
		bucketIterator=bucketPtr->begin();
		bucketIterator=bucketPtr->end();

		for (bucketIterator=bucketPtr->begin();bucketIterator != bucketPtr->end();++bucketIterator)
		{
			oldPatptr=bucketIterator->second;
			(*maximal_patterns)[oldPatptr->getID()]=oldPatptr;
			
		}
	}
}
