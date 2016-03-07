#include "SimplePath.h"
#include "Loader.h"
extern Node** node_array;


set<short>* SimplePath::getExtensibleNeighbors(set<short>* startNodes,set<short>* endNodes) const
{
	set<short>* extensibleNeighbors = new set<short>(); //deleted in merge function
	Node * lastNodePtr = node_array[IDOfLastNode];

	for (set<short>::iterator iter=startNodes->begin(); iter!=startNodes->end(); ++iter)
	{
		if(lastNodePtr->isNeighbor(*iter) && !willCreateCycle(*iter)) // must be a neighbor of last node but not neighbor of others
			extensibleNeighbors->insert(*iter);
	}

	for (set<short>::iterator iter=endNodes->begin(); iter!=endNodes->end(); ++iter)
	{
		if(lastNodePtr->isNeighbor(*iter) && !willCreateCycle(*iter)) // must be a neighbor of last node but not neighbor of others
			extensibleNeighbors->insert(*iter);
	}

	return extensibleNeighbors;
}

bool SimplePath::willCreateCycle(short newNodeID) const
{
	for (set<short>::iterator iter=elements->begin(); iter!=elements->end(); ++iter)
	{
		if(((*iter)!=IDOfLastNode) && node_array[*iter]->isNeighbor(newNodeID))
			return true;
	}
	return false;
}

set<short>* SimplePath::getElements() 
{
	return elements;
}

short SimplePath::getIDOfLastNode() const
{
	return IDOfLastNode;
}

bool SimplePath::containsNode(short nodeID) const
{
	return elements->find(nodeID)!=elements->end();
}

void SimplePath::addElement(short nodeID)
{
	elements->insert(nodeID);
	IDOfLastNode=nodeID;
}

SimplePath::SimplePath(short nodeID)
{
	elements = new set<short>();
	elements->insert(nodeID);
	IDOfLastNode=nodeID;
}

SimplePath::SimplePath(SimplePath * prefixPathPtr,short extensionNodeID)
{
	elements= new set<short>();

	for (set<short>::iterator iter=prefixPathPtr->getElements()->begin();
		iter!=prefixPathPtr->getElements()->end(); ++iter)
	{
		elements->insert(*iter);
	}
	//assert(elements->size()>0);
	elements->insert(extensionNodeID);
	IDOfLastNode=extensionNodeID;
}

SimplePath::~SimplePath(void)
{
	delete elements;
}

