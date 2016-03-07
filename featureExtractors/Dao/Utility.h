#pragma once
#include "NCSC.h"
#include "Node.h"
#include "ConnectedComponent.h"
#include "Pattern.h"
#include "SimplePath.h"
#include "MultiHashIndex.h"

class Utility
{
public:
	void findUniqueGenesInPatterns(ConnectedComponent*) const;
	bool isEdgeEdgeHomogeneous(const Node *, const Node *) const;
	void removeEdge(Node *,Node *) const;
	void removeNonHomegenousEdges(void) const;
	void extractComponents() const;
	short getUnvisitedNode(bool* visitFlag) const;
	void extendComponent(ConnectedComponent* cc, Node * root,bool* visitFlag) const;
	void deallocateNodes()const;
	void deallocateConnectedComponents(void) const;
	int getNumberOfOverlaps(const ConnectedComponent* ccPtr) const;
	void writeGraphInGraclusFormat() const;
	void partitionGraph(int)const;
	void readPartitions(char* partitionFile, int numberOfPartitions) const;
	void unifyPatternsOfAllComponents(ConnectedComponent* ) const;
	void unifyPatternsOfAllComponents2(ConnectedComponent* ) const; // new index
	void unifyCriticalPatternsOfAllComponents(ConnectedComponent* ccPtr) const; // new function based on new merge algo
	int * get_Node_PatternSetSize_Array() const;
};
