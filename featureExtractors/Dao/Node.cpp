#include "Node.h"

//
extern short NUMBER_OF_ATTRIBUTES;

// default constructor 
Node::Node(short id,char* name)
{
	//cout<<"\tDeafult constructor for "<< id<<" is called!"<<endl;
	this->id=id;
	strcpy(this->name,name);
	attributes = new float[NUMBER_OF_ATTRIBUTES];
	neighbors;
}


Node::~Node()
{
	delete[] attributes;
	//cout<<"\tDestructor for node:"<<id<<"is called!"<<endl;
}

void Node::addNeighbor(short id)
{
	neighbors.insert(id);
}

set<short>* Node::getNeigbors()
{
	return &neighbors;
}

char * Node::getName()
{
	return name;
}


unsigned int Node::getNumberOfNeighbors(void) const
{
	return (unsigned int)neighbors.size();
}

bool Node::isNeighbor(short id) const
{
	return neighbors.find(id)!=neighbors.end();
}

void Node::removeNeighbor(short id)
{
	neighbors.erase(id);
}


short Node::getID() const
{
	return id;
}

const float* Node::getAttributes() const
{
	return this->attributes;
}

void Node::setAttributes(const float * attr)
{
	for(int i=0;i<NUMBER_OF_ATTRIBUTES;i++)
	{
		attributes[i]=attr[i];
	}
}


void Node::print() const
{
	cout<<"\tI am node:"<<this->id<<" \n\t\tattributes are:\t";
	for(int i=0;i<NUMBER_OF_ATTRIBUTES;i++)
	{
		cout<<attributes[i]<<"  ";
	}
	cout<<"\n\t\tneighbors("<<this->neighbors.size()<<"):\t";


	for (set<short>::const_iterator neighIter=neighbors.begin(); neighIter!=neighbors.end(); ++neighIter)
	{
		cout << " " << *neighIter;
	}

	cout << endl;
}


