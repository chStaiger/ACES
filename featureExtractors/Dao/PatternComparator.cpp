#include "PatternComparator.h"


PatternComparator::PatternComparator(void)
{
}



PatternComparator::~PatternComparator(void)
{
}

bool PatternComparator::operator() (const set<short>* ID1,const set<short>* ID2) const 
{
	int min=-1;
	bool equal=true;
	set<short>::const_iterator iter1=ID1->begin();
	set<short>::const_iterator iter2=ID2->begin();

	//printKeys(ID1,ID2);

	if(ID1->size()<ID2->size()){min=(short)(ID1->size());}
	else{min=(short)(ID2->size());}

	for(int i=0;i<min; i++)
	{
		if(*iter1<*iter2){return true;} //p1 is for sure smaller
		if(*iter1>*iter2){return false;}//p1 is for sure greater	
		iter1++;
		iter2++;
	}
	
	if(ID1->size()< ID2->size())
	{
		return true;
	} //abc vs abcde
	
	return false; // abc<abc
}



void PatternComparator::printKeys(const set<short>* ID1,const set<short>* ID2) const
{
	cout<<"comparing keys:\n\t";
	for (set<short>::const_iterator iter=ID1->begin(); iter!=ID1->end(); ++iter)
	{
		cout <<*iter<<"-";
	}
	cout<<"\n\t";
	
	for (set<short>::const_iterator iter=ID2->begin(); iter!=ID2->end(); ++iter)
	{
		cout <<*iter<<"-";
	}
	cout<<endl;
}




void mainss()
{
	
	PatternComparator p;
	map<set<short>*,int,PatternComparator> map;

	set<short>* s1 = new set<short>;
	set<short>* s2 = new set<short>;
	set<short>* s3 = new set<short>;
	set<short>* s4 = new set<short>;

	s1->insert(1);
	
	s2->insert(1);
	s2->insert(2);
	
	s3->insert(1);
	s3->insert(2);

	s4->insert(4);

	map[s1]=1;
	cout<<"-----"<<endl;
	map[s2]=2;
	cout<<"-----"<<endl;
	map[s3]=3;

	cout<<"recep"<<endl;

}
