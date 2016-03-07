#pragma once

class Loader
{
public:
	Loader(void);
	~Loader(void);
	void loadDataset(const char* dir);
	void loadNodes(const char*);
	void loadEdges(const char*);
	void loadAttributes(const char * attributesFileName);
	void loadDimensions(const char * dimFileName);

private:
	void allocateNodeArray(const char * nodeFileName);
};
