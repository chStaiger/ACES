#pragma once
#include "NCSC.h"
#include "Node.h"
#include "Pattern.h"
#include "PatternComparator.h"
#include "SimplePath.h"
#include <time.h>
#include <stack>

class ConnectedComponent
{
public:
	ConnectedComponent(int);
	~ConnectedComponent(void);
	Node * getNode(short);
	map<short,Node*>::const_iterator getNodeIter()const {return nodes.begin();};
	void addNode(Node * newNodePtr);
	void print(void)const;
	int getID()const;
	short getSize()const;
	void findSeedPatterns();
	void expand_by_one(bool);
	void cleanMemoryAfterExpansion(int);
	void printMaximalPatterns();
	void merge_patterns();
	void merge_patterns2();
	void findPathStartAndEndNodes(const Pattern* islandPtr1,const Pattern* islandPtr2, const set<short>* overlapNodes,set<short>* &startNodes,set<short>* &endNodes)const;
	map<set<short>*,SimplePath*,PatternComparator>* findAllPaths(set<short>* startNodes,set<short>* endNodes)const;
	void connectIslands(Pattern* islandPtr1,Pattern* islandPtr2,map<set<short>*,SimplePath*,PatternComparator>*paths);
	void insertNewlyFoundMergePattern(Pattern * ptr);
	void removeRedundancyInMaximalPatterns();
	map<set<short>*,Pattern*,PatternComparator>* getMaximalPatterns()const;
	multimap<short,Pattern*>* getPotentialCriticalPatterns()const;
	bool isPatternModule(short size, double density, short dimensionality) const;
	multimap<short,Pattern*>* getTempMaximalPatterns()const;
private:
	map<short, Node*> nodes;
	map<set<short>*,Pattern*,PatternComparator>* current_pattern_map;
	map<set<short>*,Pattern*,PatternComparator>* valid_candidate_map;
	map<set<short>*,Pattern*,PatternComparator>* maximal_patterns;
	multimap<short,Pattern*>* potential_critical_patterns;
	multimap<short,Pattern*>* temp_maximal_patterns;
	list<Pattern*> * critical_patterns;
	list<Pattern*> * final_patterns;
	int id;
};
