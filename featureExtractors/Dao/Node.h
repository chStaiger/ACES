#pragma once
#include "NCSC.h"


class Node
{
public:
	Node(short,char*);
	Node(const Node &);
	~Node(void);
	//void initialize(short);
	void addNeighbor(short);
	void removeNeighbor(short);
	short getID() const;
	const float* getAttributes() const;
	char * getName();
	void setAttributes(const float *);
	void print() const;
	unsigned int getNumberOfNeighbors(void) const;
	set<short>* getNeigbors();
	bool isNeighbor(short) const;
private:
	char name[100];
	short id;
	float* attributes;
	set<short> neighbors;
};
