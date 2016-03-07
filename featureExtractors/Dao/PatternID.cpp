#include "PatternID.h"
#include "NCSC.h"

int PatternID::counter=0;

PatternID::PatternID(short id1,short id2)
{
	value=new short[2];
	value[0]=id1;
	value[1]=id2;
	size=2;
	PatternID::counter++;
	//cout<<"-------PatternID constructor for:"<<id1<<"-"<<id2<<" is called!"<<endl;
}

PatternID::PatternID(PatternID* existingPatternID,short newNodeID)
{
	PatternID::counter++;
	size = existingPatternID->getSize()+1;
	value = new short[size+1];
	int index=0;
	short * existingPatternIDValue=existingPatternID->getValue();

	//transfer smaller ones
	while(existingPatternIDValue[index]<newNodeID && index<size-1)
	{
		value[index]=existingPatternIDValue[index];
		index++;
	}
	//insert the new one
	value[index]=newNodeID;
	index++;
	//transfer the bigger ones
	while(index<size)
	{
		value[index]=existingPatternIDValue[index-1];
		index++;
	}

	//cout<<"-------PatternID constructor";
	//print();
}

PatternID::~PatternID()
{
	//cout<<"--------PatternID destructor is called for:";
	//print();
	PatternID::counter--;
	delete[] value;
}

short PatternID::getSize() const
{
	return size;
}

short * PatternID::getValue() const
{
	return value;
}

short PatternID::getIDOfNodeAt(short index) const
{
	return value[index];
}


bool PatternID::containsNode(short nodeID) const
{
	for(int i=0;i<size;i++)
	{
		if(value[i]==nodeID)
			return true;
	}
	return false;
}

short PatternID::getNumberOfOverlap(PatternID* otherPatternID) const
{
	int counter=0;
	for(short i=0;i<size;i++)
	{
		if(otherPatternID->containsNode(value[i]))
		{
			counter++;
		}
	}

	return counter;
}


void PatternID::print()const
{
	cout<<"\tID:";
	for(int i=0;i<size;i++)
		cout<<value[i]<<" ";
	cout<<endl;
}

