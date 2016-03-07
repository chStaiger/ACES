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
#include <math.h>
#include <time.h>

using namespace std;

const int LINELIMIT=16384;
const int MAXMODULES=30000;
const int MAXPATIENTS=512;
const int MAXGENES=30000;
const char* DELIM="\t\n";

typedef set<int> module;

int n;  // the number of genes
int m;  // the number of modules
int nn,nu;
int i,j,k,l;
map<string, int> g2id;
map<int, string> id2g;
char cmdline[256];
char filename[256];
char filename1[256];
char line[LINELIMIT];
char line1[LINELIMIT];
vector< module> moduleArr;
vector< module> moduleArrDim;

char* token;
char* token1;
double **g;
double expG[MAXGENES][MAXPATIENTS];
double expMod[MAXMODULES][MAXPATIENTS];
double expTotalEdges=0.0;
double totalEdges=0.0;
double numEdges[MAXMODULES];
double expEdges[MAXMODULES];
double prob[MAXMODULES];
int moduleOrder[MAXMODULES];
bool choosen[MAXMODULES];
bool patientSpace[MAXPATIENTS];
int numPatient;


FILE *f;
FILE *f1;
FILE *f2;

void quickSort(int *a, int lo, int hi)
{
//  lo is the lower index, hi is the upper index
//  of the region of array a that is to be sorted
    int i=lo, j=hi,h;
    double x=prob[a[(lo+hi)/2]];

    //  partition
    //cout<<"fail where"<<endl;
    while (i<=j)
    {
        while (prob[a[i]]>x) i++;
        while (prob[a[j]]<x) j--;
        //cout<<"fail here"<<endl;
        if (i<=j)
        {
            h=a[i]; a[i]=a[j]; a[j]=h;
            i++; j--;
        }
    }

    //  recursion
    if (lo<j) quickSort(a, lo, j);
    if (i<hi) quickSort(a, i, hi);
}


double informationGain(double* v)
{
    int i,j,k,m;
    m=nn+nu;
    int numBins=ceil(log(m*1.0));
    //int numBins=4;
    double l=1000000.0;
    double r=-1000000.0;
    for (i=0;i<m;i++){
        if (v[i]<l)
            l=v[i];
        if (v[i]>r)
            r=v[i];
    }


    double binLength=(r-l)/(numBins*1.0);
    double prob0=nn/(m*1.0);
    double prob1=nu/(m*1.0);
    double entropy=-prob1*log(prob1)-prob0*log(prob0);

    //cout<<"ENTROPY "<<nn<<" "<<nu<<" "<<entropy;

    double maxR;
    double minThres=0.00001;
    int totalInBin;
    double pInBin;
    double prob1InBin;
    double prob0InBin;
    int num1InBin;
    int num0InBin;
    for (int i=1;i<=numBins;i++)
    {
        if (i==numBins)
            maxR=r+minThres;
        else
            maxR=l+i*binLength;

        totalInBin=0;
        for (j=0;j<m;j++)
        if ((v[j]>=l+(i-1)*binLength)&&(v[j]<maxR))
            totalInBin++;
        pInBin=totalInBin/(m*1.0);

        if (totalInBin>0)
        {
            num0InBin=0;
            for (j=0;j<nn;j++)
            if ((v[j]>=l+(i-1)*binLength)&&(v[j]<maxR))
                num0InBin++;
            prob0InBin=num0InBin/(totalInBin*1.0);

            num1InBin=0;
            for (j=nn;j<m;j++)
            if ((v[j]>=l+(i-1)*binLength)&&(v[j]<maxR))
                num1InBin++;
            prob1InBin=num1InBin/(totalInBin*1.0);

            if (prob1InBin>minThres)
                entropy=entropy+pInBin*prob1InBin*log(prob1InBin);
            if (prob0InBin>minThres)
                entropy=entropy+pInBin*prob0InBin*log(prob0InBin);

            //cout<<"ENTROPY "<<entropy;
        }
    }
    return entropy;
}


int main(int argc,char* argv[])
{

    double w;
    sprintf(filename,"%s/nodes.txt",argv[1]);
	n=0;
	f=fopen(filename,"r");
    while (fgets(line, LINELIMIT, f)>0){
        token = strtok(line, DELIM);
        sscanf(token,"%d",&i);
        if (i>n)
            n=i;
    	token = strtok(NULL, DELIM);
    	id2g[i]=token;
    	g2id[token]=i;
    }
	fclose(f);
	n++;

    nu=0;
    sprintf(filename,"%s/dimensions.txt",argv[1]);
    f=fopen(filename,"r");
    while (fgets(line, LINELIMIT, f)>0){
        token = strtok(line, DELIM);
    	nu++;    
    }
    fclose(f);

    sprintf(filename,"%s/expression.txt",argv[1]);
    f=fopen(filename,"r");
    while (fgets(line, LINELIMIT, f)>0){
        token = strtok(line, DELIM);
	sscanf(token,"%d",&i);
	j=0;
	token=strtok(NULL,DELIM);
    	while (token!=NULL){
		//token = strtok(NULL, DELIM);
    		sscanf(token,"%lf",&w);
		expG[i][j]=w;
		j++;	
		token = strtok(NULL, DELIM);	
	}
	nn=j-nu;   
    }
    fclose(f);


    g=(double**)malloc(sizeof(double*)*n);
    int i,j;
    for (i=0;i<n;i++)
        g[i]=(double*)malloc(sizeof(double)*n);
    for (i=0;i<n;i++)
    for (j=0;j<n;j++)
    	g[i][j]=0.0;
    sprintf(filename,"%s/edges.txt",argv[1]);
	f=fopen(filename,"r");
    while (fgets(line, LINELIMIT, f)>0){
        //cout<<"iii"<<endl;
        token = strtok(line, DELIM);
        sscanf(token,"%d",&i);
    	token = strtok(NULL, DELIM);
    	sscanf(token,"%d",&j);
    	token = strtok(NULL, DELIM);
    	sscanf(token,"%lf",&w);
    	totalEdges+=1.0;
    	expTotalEdges+=w;
    	g[i][j]=w;
    	g[j][i]=w;
    }
	fclose(f);

	//totalEdges=(n*(n-1))/2.0;

	sprintf(cmdline,"ls -1 %s/%s/modules/ > %s/%s_list.txt",argv[1],argv[2],argv[1],argv[2]);
	//cout<<cmdline<<endl;
	system(cmdline);

    vector< module >::iterator iter;
    set<int>::iterator iter1;
    set<int>::iterator iter2;
    module mod;
    module modDim;
    double expE;
    double numE;

    sprintf(filename,"%s/%s_list.txt",argv[1],argv[2]);
    f=fopen(filename,"r");
    m=0;
    while (fgets(line, LINELIMIT, f)>0){
        token = strtok(line, DELIM);
        sprintf(filename1,"%s/%s/modules/%s",argv[1],argv[2],token);
        //cout<<filename1<<endl;
        mod.clear();
        f1=fopen(filename1,"r");
        //cout<<"module";
        while (fgets(line1, LINELIMIT, f1)>0){
            token1 = strtok(line1, DELIM);
            mod.insert(g2id[token1]);
            //cout<<" "<<token1;
        }
        fclose(f1);
        //sprintf(filename1,"%s/%s/modulesDim/%s",argv[1],argv[2],token);
        //cout<<filename1<<endl;
        //modDim.clear();
        //f1=fopen(filename1,"r");
        //cout<<"module";
        //while (fgets(line1, LINELIMIT, f1)>0){
            //token1 = strtok(line1, DELIM);
            //if (token1!=NULL)
            //modDim.insert(atoi(token1));
            //cout<<" "<<token1;
        //}
        //fclose(f1);
        //cout<<endl;

            expE=0.0;
            numE=0.0;
            for (iter1=mod.begin();iter1 != mod.end();iter1++)
            for (iter2=mod.begin();iter2 != mod.end();iter2++)
            if (*iter1<*iter2){
                numE+=1.0;
                if (g[*iter1][*iter2])
                    expE+=g[*iter1][*iter2];
            }
            //if (expE/numE>0.95)

            moduleArr.push_back(mod);
            moduleArrDim.push_back(modDim);
            moduleOrder[m]=m;
            numEdges[m]=0.0;
            expEdges[m]=0.0;
	    for (i=0;i<nn+nu;i++){
		expMod[m][i]=0.0;
	    	for (iter1=mod.begin();iter1 != mod.end();iter1++)
			expMod[m][i]+=expG[*iter1][i];
		expMod[m][i]/=mod.size()*1.0;
	    }
            m++;


    }
    fclose(f);

    //cout<<m<<endl;

    //vector< module >::iterator iter;
    //set<int>::iterator iter1;
    //set<int>::iterator iter2;

    iter = moduleArr.begin();
    i=0;
    while( iter != moduleArr.end() ) {
        //cout<<"module";
        for (iter1=(*iter).begin();iter1 != (*iter).end();iter1++)
        for (iter2=(*iter).begin();iter2 != (*iter).end();iter2++)
        if (*iter1<*iter2){
            if (g[*iter1][*iter2]){
                numEdges[i]+=1.0;
                expEdges[i]+=g[*iter1][*iter2];
            }
        }
        //phyper(numBad1-1,colonCancerSize,n-colonCancerSize,modSize,0,0)
        //HGProbabilityCalc c((int)(round(totalEdges)),(int)(round(expTotalEdges)),(int)(round(numEdges[i])));
        //prob[i]=c.odds_against((int)(round(expEdges[i])));
        //prob[i]=phyper((int)(round(expEdges[i]))-1,(int)(round(expTotalEdges)),(int)(round(totalEdges))-(int)(round(expTotalEdges)),(int)(round(numEdges[i])),0,0);
	prob[i]=informationGain(expMod[i]);
	cout<<"Module "<<i<<" score: "<<prob[i]<<endl;
        //cout<<prob[i]<<" "<<(int)(round(expEdges[i]))<<" "<<(int)(round(numEdges[i]))<<" "<<(int)(round(expTotalEdges))<<" "<<(int)(round(totalEdges))<<endl;
        i++;
        iter++;
        //cout<<endl;
    }

    //cout<<(5+6)/2<<endl;
    //cout<<"recursive quick sort"<<endl;
    quickSort(moduleOrder,0,m-1);

    /*
    for (i=0;i<m;i++){
        cout<<"totalEdges= "<<totalEdges<<" expTotalEdges= "<<expTotalEdges<<" numEdges= "<<numEdges[moduleOrder[i]]<<" expEdges= "<<expEdges[moduleOrder[i]]<<" prob= "<<prob[moduleOrder[i]]<<endl;
    }
    */

    //cout<<"stupid down here"<<endl;
    sprintf(filename,"%s/%s/%s",argv[1],argv[2],argv[3]);
    f=fopen(filename,"w");
    fprintf(f," Final number of patterns:%d\n",m);
    fprintf(f,"Maximal patterns:\n");
    int top50;
    //if (m<=50)
        //top50=m;
    //else
        //top50=50;
    top50=0;
    set<int> allVertices;

    for (i=0;i<m;i++)
        choosen[i]=false;

    //choosen[0]=true;
    for (i=0;i<m;i++)
    if (top50<50)
    {
        int numNew=0;
	cout<<"Module "<<i<<" with order "<<moduleOrder[i]<<endl;
        for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++){
            if (allVertices.find(*iter1)==allVertices.end())
                numNew++;
        }
        //cout<<numNew<<"XXXXXXXXXXXXXXX"<<moduleArr[moduleOrder[i]].size()<<endl;
        if ((numNew*2>=moduleArr[moduleOrder[i]].size())||(i==0))
        //if (numNew>=1)
        //if (numNew==((int)moduleArr[moduleOrder[i]].size()))
        {
            choosen[i]=true;
            top50++;
            for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++)
                allVertices.insert(*iter1);
        }

        if (choosen[i]){
            for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++){
                fprintf(f,"%d ",*iter1);
            }
            fprintf(f,"Weight = %lf\n",expEdges[moduleOrder[i]]);
            fprintf(f,"\n");
            fprintf(f,"\t Size:%d\tWeighted Density:%lf\n",moduleArr[moduleOrder[i]].size(),(2.0*expEdges[moduleOrder[i]])/(1.0*moduleArr[moduleOrder[i]].size()*(moduleArr[moduleOrder[i]].size()-1)));
            fprintf(f,"\t");
            for (iter1=(moduleArrDim[moduleOrder[i]]).begin();iter1 != (moduleArrDim[moduleOrder[i]]).end();iter1++){
                fprintf(f,"%d-",*iter1);
            }
            fprintf(f,"\n");
            fprintf(f,"\n");
        }
    }
    fclose(f);

    //cout<<"stupid at the end"<<endl;

    sprintf(cmdline,"rm %s/%s/modules/*",argv[1],argv[2]);
    //cout<<cmdline<<endl;
    system(cmdline);

    //cout<<"stupid at the end"<<endl;

    j=0;
    for (i=0;i<m;i++)
    if (choosen[i]){
        j++;
        sprintf(filename,"%s/%s/modules/%d",argv[1],argv[2],j);
        //cout<<filename<<endl;
        f=fopen(filename,"w");
        for (iter1=(moduleArr[moduleOrder[i]]).begin();iter1 != (moduleArr[moduleOrder[i]]).end();iter1++)
            //fprintf(f,"%d\n",*iter1);
            fprintf(f,"%s\n",id2g[(*iter1)].c_str());
            fclose(f);
    }

    for (i=0;i<n;i++)
        free(g[i]);
    free(g);
    sprintf(cmdline,"%s/%s/modules/",argv[1],argv[2]);
    cout<<"\n\nTHE LIST OF MODULES ARE IN "<<cmdline<<endl;
}
