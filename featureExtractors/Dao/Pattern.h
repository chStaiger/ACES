#pragma once
#include "NCSC.h"
#include "Node.h"
#include "Pattern.h"
#include "PatternComparator.h"
#include "SimplePath.h"

extern short NUMBER_OF_ATTRIBUTES;
class Pattern
{
public:
	Pattern(const Node *n1,const Node *n2);
	Pattern(Pattern * motherPattern,short newNodeID, const bool* oldPatternRelDim, float newWeight);
	~Pattern(void);
	inline bool* getRelevantDim(){return relevantDim;};
	//void updateRelevantDimensions();
	void updateRelevantDimensions2();
	bool isHomogenous()const;
	void print() const;
	int getNumberOfEdgesTo(short * )const;
	bool containsNode(short) const;
	bool isCritical(short) const;
	bool isDense() const;
	int getNumberOfEdges(){return numberOfEdges;};
	float getWeight(){return weight;};
	map<short, float>* getExtensibleNeighbors() const;
	set<short> * getOverlappingNodes(Pattern* otherPattern) const;
	map<set<short>*,Pattern*,PatternComparator>* getIslands();
	short getNumberOfCommonRelevantDimensions(Pattern*) const;
	short getNumberOfRelevantDimensions() const;
	short getNumberOfOverlappingNodes(Pattern* otherPattern) const;
	float getDensity() const;
	set<short>* getID(){return &ID;};
	short getSize()const;
	void calculateWE();
	bool isSuperPatternOf(Pattern* rhsPtr) const;
	float getScore() const;
	bool operator < (const Pattern& refParam) const;
	static int counter;
private:
	set<short> ID;
	float weight;
	int numberOfEdges;
	bool* relevantDim;
	short numberOfRelevantDims;
};
