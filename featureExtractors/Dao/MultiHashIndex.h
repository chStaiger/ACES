#pragma once
#include "NCSC.h"
#include "Pattern.h"
#include "PatternComparator.h"

class MultiHashIndex
{
public:
	MultiHashIndex(int);
	~MultiHashIndex(void);
	bool insertPattern(Pattern *);
	void unifyPatterns(map<set<short>*,Pattern*,PatternComparator> *);
	bool insertPattern2(Pattern * ,int* );
	
private:
	map<set<short>*,Pattern*,PatternComparator>** node_patterns_map;
	int numberOfBuckets;
};
