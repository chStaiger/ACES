#pragma once
#include "NCSC.h"
#include "Node.h"

class SimplePath
{
public:
	SimplePath(short);
	SimplePath(SimplePath * prefixPathPtr,short extensionNodeID);
	~SimplePath(void);
	short getIDOfLastNode() const;
	bool containsNode(short nodeID) const;
	void addElement(short nodeID);
	set<short>* getExtensibleNeighbors(set<short>* startNodes,set<short>* endNodes) const;
	bool willCreateCycle(short newNodeID) const;
	set<short>* getElements();
private:
	set<short>* elements;
	short IDOfLastNode;
};
