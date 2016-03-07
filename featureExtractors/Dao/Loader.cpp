#include "Loader.h"
#include "Node.h"
#include "Pattern.h"
using namespace std;

extern map<short,double> SIZE_MIN_DENSITY_MAP;
extern map<short,int> SIZE_MIN_DIMENSIONALITY_MAP;

extern map<short,string*> INDEX_DIMENSION_MAP;
extern Node** node_array;
extern short NUMBER_OF_ATTRIBUTES;
extern short NUMBER_OF_NODES;
extern float MISSING_VALUE;
extern int NUMBER_OF_EDGES;

float ** weight_array;

void Loader::loadDataset(const char* dir)
{
	std::ostringstream converter1;
	converter1<<dir<<"nodes.txt";
	string nodeFile=converter1.str();
	
	std::ostringstream converter2;
	converter2<<dir<<"edges.txt";
	string edgeFile=converter2.str();
	
	std::ostringstream converter3;
	converter3.clear();
	converter3<<dir<<"attributes.txt";
	string attrFile=converter3.str(); 
	
	std::ostringstream converter4;
	converter4<<dir<<"dimensions.txt";
	string dimFile=converter4.str();

	std::ostringstream converter5;
	converter5<<dir<<"parameters.txt";
	string parFile=converter5.str();

	cout<<"==============================\nLoading the dataset"<<endl;
	loadDimensions(dimFile.c_str());
	loadNodes(nodeFile.c_str());
	loadEdges(edgeFile.c_str());
	loadAttributes(attrFile.c_str());
	cout<<"loaded all dataset\n========================\n\n\n\n\n"<<endl;
}

void Loader::loadDimensions(const char * dimFileName)
{
	cout<<"Loading dimensions from file: "<<dimFileName<<" to the index_dim map..."<<endl;
	ifstream dimFile (dimFileName);

	//map<int,string*> INDEX_DIM_MAP; // map of words and their frequencies

	short index = 0; 
	char dimName[200];
	string * dimNameStr;
	
	if (dimFile.is_open())
	{
		while (! dimFile.eof() )
		{		  
		  dimFile>>index;
		  dimFile>>dimName;

		  dimNameStr = new string(dimName);
		  cout<<index<<" --> "<<*dimNameStr<<endl;
		  INDEX_DIMENSION_MAP[index]=dimNameStr;

		}
		NUMBER_OF_ATTRIBUTES=index+1;
		cout<<"\tNumber of dimensions:"<<NUMBER_OF_ATTRIBUTES<<endl;
		dimFile.close();
	}
	else
	{
		cout << "Unable to open dimensions file"; 
		exit(1); 
	}
	cout<<"\tLoad dimensions finished!\n"<<endl;
}


void Loader::loadNodes(const char * nodeFileName)
{
	cout<<"Loading nodes from file: "<<nodeFileName<<" to the id_node map..."<<endl;
	allocateNodeArray(nodeFileName);
	ifstream nodeFile (nodeFileName);
	short id = 0; 
	char nodeName[100];
	if (nodeFile.is_open())
	{
		while (! nodeFile.eof() )
		{
		  nodeFile>>id;
		  nodeFile>>nodeName;
		  node_array[id] = new Node(id,nodeName);
		  node_array[id]->print();
		}
		nodeFile.close();
	}
	else
	{
		cout << "Unable to open file"; 
		exit(1); 
	}
	cout<<"\tLoad node finished!\n"<<endl;
}


void Loader::allocateNodeArray(const char * nodeFileName)
{
	cout<<"\t allocating node array ..."<<endl;
	ifstream nodeFile (nodeFileName);
	int nodeCounter = 0; 
	char name[100];
	if (nodeFile.is_open())
	{
		while (! nodeFile.eof() )
		{
		   nodeFile>>nodeCounter;
		   nodeFile>>name;
		   //cout<<nodeCounter++<<endl;
		}
		nodeFile.close();

		try 
		{ 
			NUMBER_OF_NODES=nodeCounter+1;
			node_array = new Node*[NUMBER_OF_NODES];
			cout<<"\t allocatated "<<NUMBER_OF_NODES<<" nodes!"<<endl;
		} 
		catch (bad_alloc xa) 
		{ 
			cout << "Node array allocation failed!"<<endl; 
			exit(1); 
		} 
		
	}
	else
	{
		cout << "Unable to open file"; 
		exit(1); 
	}
}



void Loader::loadEdges(const char * edgeFileName)
{
	cout<<"Loading edges from file: "<<edgeFileName<<endl;
	//short * edgeNumbers = new short[NUMBER_OF_NODES];
	ifstream edgeFile (edgeFileName);
	short id1 = 0; 
	short id2 = 0;
      float weight = 0.0;

	weight_array = (float **)malloc(NUMBER_OF_NODES*sizeof(float *));
	for (int counter=0; counter<NUMBER_OF_NODES; counter++)
	{
		weight_array[counter]=(float *) malloc(NUMBER_OF_NODES*sizeof(float));
	}

	if (edgeFile.is_open())
	{
		while (! edgeFile.eof() )
		{
		  edgeFile>>id1;
		  edgeFile>>id2;
		  edgeFile>>weight;
		  //cout<<"\tEdge:\t"<<id1<<" \t "<<id2<<endl;
		  node_array[id1]->addNeighbor(id2);
		  node_array[id2]->addNeighbor(id1);
              weight_array[id1][id2] = weight;
              weight_array[id2][id1] = weight;
		  NUMBER_OF_EDGES++;
		}
		edgeFile.close();
		cout<<"\t"<<NUMBER_OF_EDGES<<" edges are loaded \n\n\n"<<endl;
	}
	else
	{
		cout << "Unable to open edge file ";
		exit(1); 
	}
	cout<<"\tLoad edge finished!\n"<<endl;
	//delete [] edgeNumbers;
}


void Loader::loadAttributes(const char * attributesFileName)
{
	cout<<"Loading attributes from file: "<<attributesFileName<<endl;
	int missingValueCounter=0;
	ifstream attributesFile (attributesFileName);
	
	short id = 0; 
	float * attributes = new float[NUMBER_OF_ATTRIBUTES];
	char stringValue[100];

	if (attributesFile.is_open())
	{
		while (! attributesFile.eof() )
		{
		  attributesFile>>id;
		  //cout<<"\n\tLoading attributes of node:\t"<<id<<endl;
		  for(int i=0;i<NUMBER_OF_ATTRIBUTES;i++)
		  {			  
			  attributesFile>>stringValue;
			  //cout<<i<<"  "<<stringValue<<"\t"<<endl;
			  if(strcmp("NaN",stringValue)!=0)
			  {
				  attributes[i]=atof(stringValue);
			  }
			  else
			  {
				  //cout<<"missing value"<<endl;
				  attributes[i]=MISSING_VALUE;
				  missingValueCounter++;
			  }

		  }

		  node_array[id]->setAttributes(attributes);
		  //node_array[id]->print();
		}

		attributesFile.close();
		cout<<"\t"<<missingValueCounter<<" missing values!"<<endl;
		cout<<"\tLoaded all attributes...\n"<<endl;
	}
	else
	{
		cout << "Unable to open file attr"; 
		exit(1);
	}

	delete[] attributes;
}

Loader::Loader(void)
{
	cout<<"Loader constructor is callded..."<<endl;
}

Loader::~Loader(void)
{
}
