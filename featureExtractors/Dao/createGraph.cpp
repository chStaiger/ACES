#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <deque>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <time.h>

using namespace std;

const int LINELIMIT=100000;
double NULLMARKER=-999999.1;
const int MAXSAMPLES=2000;
const int MAXGENES=30000;
const char* DELIM=",\t\n";
const char* DELIM1="|";
const char* DELIM2=". \t\n";
const char* DELIM3=",";

int n;  // the number of genes
int m;  // the number of samples
int nn; // the number of control samples
int nu; // the number of case samples
int i,j,k,l;
map<string, int> g2id;
map<string, int> stringMap;
map<string, int> g2id2;
map<string, int> g2id1;
map<int, string> id2g;
char* line;
char cmdline[256];
char filename[256];
char line1[256];
bool* expressed[MAXSAMPLES];
double** eArr[MAXSAMPLES];
bool** deArr[MAXSAMPLES];
double meanArr[MAXGENES];
double stdArr[MAXGENES];
double meanArr1[MAXGENES];
double stdArr1[MAXGENES];
char* tempString1;
char* tempString2;
char* token;
char* token1;
int* present;
double **expressionArr;
double **ppi;
double **posGraph;
double **negGraph;
double **diffGraph;
double thres;
bool allPresent[MAXGENES];
char command[256];

FILE *f;
FILE *f1;
bool new1,new2;
int id1,id2;

void allocateGraph(double*** a,int n)
{
    *a=(double**)malloc(sizeof(double*)*n);
    int i,j;
    for (i=0;i<n;i++)
        (*a)[i]=(double*)malloc(sizeof(double)*n);
    for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    	(*a)[i][j]=0.0;
}

void deallocateGraph(double*** a,int n)
{
    int i;
    for (i=0;i<n;i++)
        free((*a)[i]);
    free(*a);
}

void resetExpressionArray()
{
	int i,j;
	for (i=0;i<n;i++)
		present[i]=0;
	for (i=0;i<n;i++)
	for (j=0;j<m;j++)
		expressionArr[i][j]=0.0;
}

bool diffExp(double value,double mean,double std,double thres)
{
	if ((value>=mean)&&((value-mean)/std>thres))
		return true;
	return false;
}

void processAllSamples()
{
	int i,j,k,l;
	for (i=0;i<n;i++)
	if (allPresent[i])
	{
	    k=0;

	    for (j=0;j<nn;j++)
	    if (eArr[j][i][0]>(NULLMARKER+1.0))
            k++;

        if (k>0)
        {

            for (j=0;j<nn;j++)
            if (eArr[j][i][0]>(NULLMARKER+1.0))
                meanArr[i]+=eArr[j][i][0];

            meanArr[i]=0.0;
            for (j=0;j<nn;j++)
            if (eArr[j][i][0]>(NULLMARKER+1.0))
                meanArr[i]+=eArr[j][i][0];
            meanArr[i]/=(k)*1.0;

            double var=0.0;
            for (j=0;j<nn;j++)
            if (eArr[j][i][0]>(NULLMARKER+1.0))
                var+=(eArr[j][i][0]-meanArr[i])*(eArr[j][i][0]-meanArr[i]);
            var/=(k)*1.0;
            stdArr[i]=sqrt(var);
        }
	}

	for (i=0;i<n;i++)
	if (allPresent[i])

	for (j=0;j<nn+nu;j++)
	if (eArr[j][i][0]>(NULLMARKER+1.0))
	{
        if (fabs(stdArr[i])<0.0001)
            eArr[j][i][0]=(eArr[j][i][0]-meanArr[i]);
        else
            eArr[j][i][0]=(eArr[j][i][0]-meanArr[i])/stdArr[i];
	}
}

int compare (const void * a, const void * b)
{
  if ((*(double*)a - *(double*)b)>=0.0)
    return 1;
  return -1;
}



void processSample(int index)
{
	int i,j;

    double a1[MAXGENES+1];

    int n11=0;

    for (i=0;i<n;i++)
    if (allPresent[i]>0)
    if (eArr[index][i][0]>NULLMARKER+1.0)
    {
        n11++;
        a1[n11-1]=eArr[index][i][0];
    }

    qsort(a1,n11,sizeof(double),compare);
    int pivot1=(int)round(n11*(1-thres));
    int pivot2=(int)round(n11*thres);

    double p1;
    int numDEs=0;
    	for (i=0;i<n;i++)
    	if (allPresent[i])
    	{
                if (eArr[index][i][0]>a1[pivot1])
                	expressed[index][i]=true;
                if ((eArr[index][i][0]<a1[pivot2])&&(fabs(eArr[index][i][0])>0.0001))
                    expressed[index][i]=true;

                if (expressed[index][i])
                    numDEs++;
    	}
    	cout<<"Number of differentially expressed genes in sample "<<index+1<<" : "<<numDEs<<endl;
}

double max(double a,double b)
{
	if (a>b)
		return a;
	return b;
}

int main(int argc,char* argv[])
{
    sprintf(command,"mkdir %s",argv[4]);
    system(command);
    sprintf(command,"mkdir %s/WDCB",argv[4]);
    system(command);
    sprintf(command,"mkdir %s/WDCB/modules",argv[4]);
    system(command);
    sprintf(command,"mkdir %s/experimental_results",argv[4]);
    system(command);


    sscanf(argv[3],"%lf",&thres);
    sscanf(argv[6],"%d",&nn);
    sscanf(argv[7],"%d",&nu);
    m=nn+nu;
    n=0;
    line=(char*)malloc(sizeof(char)*LINELIMIT);
    tempString1=(char*)malloc(sizeof(char)*LINELIMIT);
    tempString2=(char*)malloc(sizeof(char)*LINELIMIT);
    present=(int*)malloc(sizeof(int)*LINELIMIT);

    f=fopen(argv[5],"r");
    while (fgets(line, LINELIMIT, f)>0){
    	token = strtok(line, DELIM2);
        if (g2id[token]==0){
    		n++;
    		g2id[token]=n;
    		id2g[n]=token;
    	}
    	i=g2id[token];
    	token = strtok(NULL, DELIM2);
    	stringMap[token]=i;
    }
    fclose(f);

    allocateGraph(&ppi,n);
    for (i=0;i<MAXSAMPLES;i++)
    {
        eArr[i]=(double**)malloc(sizeof(double*)*n);

        for (j=0;j<n;j++)
        	eArr[i][j]=(double*)malloc(sizeof(double)*2);

        for (j=0;j<n;j++)
        for (k=0;k<2;k++)
        	eArr[i][j][k]=NULLMARKER;
    }

    for (i=0;i<MAXSAMPLES;i++)
    {
        expressed[i]=(bool*)malloc(sizeof(bool)*n);

        for (j=0;j<n;j++)
        	expressed[i][j]=false;
    }


    // Allocate the expression array
    expressionArr=(double**)malloc(sizeof(double*)*n);
    for (i=0;i<n;i++)
        expressionArr[i]=(double*)malloc(sizeof(double)*m);

    f=fopen(argv[1],"r");
    fgets(line, LINELIMIT, f);
    double weight;
    while (fgets(line, LINELIMIT, f)>0){
    	token = strtok(line, DELIM2);
    	//token = strtok(NULL, DELIM2);
    	id1=stringMap[token];
    	//token = strtok(NULL, DELIM2);
    	token = strtok(NULL, DELIM2);
    	id2=stringMap[token];
    	token = strtok(NULL, DELIM2);
    	sscanf(token,"%lf",&weight);
    	//cout<<id1<<" "<<id2<<" "<<weight<<endl;
    	if ((id1>0)&&(id2>0))
    	if ((ppi[id1-1][id2-1]<0.0001)||(ppi[id1-1][id2-1]>weight)){
            ppi[id1-1][id2-1]=weight/1000.0;
            ppi[id2-1][id1-1]=weight/1000.0;

    	}
    }
    fclose(f);


    f=fopen(argv[2],"r");

    int index;
    double v;
    fgets(line, LINELIMIT, f);
    //fgets(line, LINELIMIT, f);
    while (fgets(line, LINELIMIT, f)>0){
    		token = strtok(line, DELIM);
    		//cout<<token<<endl;
    		if ((g2id[token]>0)||(g2id2[token]>0)){
    		        if (g2id[token]>0)
    				id1=g2id[token]-1;
    			else
    				id1=g2id2[token]-1;
    			present[id1]++;
    			allPresent[id1]=true;

    			for (index=0;index<nn+nu;index++)
    			{
    			    token = strtok(NULL, DELIM);
                    sscanf(token,"%lf",&v);
                    if (eArr[index][id1][0]<(NULLMARKER+1.0))
                        eArr[index][id1][0]=v;
                    else
                        eArr[index][id1][0]+=v;
    			}
    		}
    }

    fclose(f);

    processAllSamples();
    for (i=0;i<nn+nu;i++){
    	processSample(i);
    }

    sprintf(filename,"%s/edges.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    if (present[i]>0)
    for (j=i+1;j<n;j++)
    if (present[j]>0)
    if (ppi[i][j]>0.00001)
    	fprintf(f,"%d\t%d\t%lf\n",i,j,ppi[i][j]);
    fclose(f);

    // REPLACE SPACES BY _s
    for (i=0;i<n;i++)
    if (id2g[i+1].find_first_of(" ")!=string::npos)
    {
    	id2g[i+1].erase(id2g[i+1].find_first_of(" "));
    	g2id[id2g[i+1]]=i+1;
    }

    int nms=0;
    int ms=0;

    sprintf(filename,"%s/nodes.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    	fprintf(f,"%d\t%s\n",i,id2g[i+1].c_str());
    fclose(f);

    sprintf(filename,"%s/population.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    	fprintf(f,"%s\n",id2g[i+1].c_str());
    fclose(f);

    sprintf(filename,"%s/dimensions.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<nu;i++)
    fprintf(f,"%d\tPATIENT%d\n",i,i);
    fclose(f);

    // Print out the differentially expressed genes in case samples
    sprintf(filename,"%s/attributes.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    {
    	fprintf(f,"%d",i);
    	    int count0=0;
    	    int count1=1;
    	    for (k=0;k<nn;k++)
    	    if (expressed[k][i])
                count1++;
            else
                count0++;

    	for (j=0;j<nu;j++)
    	{
            if (expressed[j+nn][i])
                fprintf(f,"\t%d",1);
            //else if ((!expressed[j+nn][i])&&(count1>=2*count0))
                //fprintf(f,"\t%d",1);
            else
                fprintf(f,"\t%d",0);

    	}
        fprintf(f,"\n");
    }
    fclose(f);

    sprintf(filename,"%s/expression.txt",argv[4]);
    f=fopen(filename,"w");
    for (i=0;i<n;i++)
    {
    	fprintf(f,"%d",i);
    	for (j=0;j<nn+nu;j++)
    	if (eArr[j][i][0]>(NULLMARKER+1.0))
    		fprintf(f,"\t%.4lf",eArr[j][i][0]);
        else
            fprintf(f,"\t0.0");
    	fprintf(f,"\n");
    }
    fclose(f);
    deallocateGraph(&ppi,n);

    for (i=0;i<n;i++)
        free(expressionArr[i]);
    free(expressionArr);

    free(line);
    free(tempString1);
    free(tempString2);
    free(present);
    return 0;
}
