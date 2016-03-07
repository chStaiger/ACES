#pragma once
#include "NCSC.h"

class PatternComparator
{
public:
	PatternComparator(void);
	~PatternComparator(void);
	bool operator() (const set<short>*,const set<short>*) const;
	void printKeys(const set<short>* ID1,const set<short>* ID2) const;
};
