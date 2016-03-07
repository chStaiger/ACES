#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(int argc,char* argv[])
{
  int numControls;
  int numCases;
  int L;
  double alpha;
  double diffThres;
  char ExpressionFile[256];
  char PPIIDFile[256];
  char PPINetworkFile[256];
  char folder[256];
  char command[256];
  char line[256];
  //const char *DELIM="=\n \t";
  const char *DELIM="=\t \n";

  FILE *f;
  char *token;

  f=fopen(argv[1],"r");
  while (fgets(line, 256, f)>0){
    token = strtok(line, DELIM);
    //printf("%s\n", line);
    if (token!=NULL)
    {
        if (strcmp(token,"WorkingFolder")==0){
            token=strtok(NULL,DELIM);
            strcpy(folder,token);

        }
        else if (strcmp(token,"NetworkFile")==0){
            token=strtok(NULL,DELIM);

            strcpy(PPINetworkFile,token);
        }
        else if (strcmp(token,"IDFile")==0){
            token=strtok(NULL,DELIM);

            strcpy(PPIIDFile,token);
        }
        else if (strcmp(token,"ExpressionFile")==0){
            token=strtok(NULL,DELIM);

            strcpy(ExpressionFile,token);
        }
        else if (strcmp(token,"#Controls")==0){
            token=strtok(NULL,DELIM);

            sscanf(token,"%d",&numControls);
        }
        else if (strcmp(token,"#Cases")==0){
            token=strtok(NULL,DELIM);

            sscanf(token,"%d",&numCases);
        }
        else if (strcmp(token,"ExpressedThreshold")==0){
            token=strtok(NULL,DELIM);

            sscanf(token,"%lf",&diffThres);
        }
        else if (strcmp(token,"DensityThreshold")==0){
            token=strtok(NULL,DELIM);
            sscanf(token,"%lf",&alpha);
        }
        else if (strcmp(token,"MinimumCases")==0){
            token=strtok(NULL,DELIM);
            sscanf(token,"%d",&L);
        }
    }
  }
  fclose(f);

  printf("All read in \n");
  /*
  printf("Please enter the folder that contains the project:\n");
  scanf("%s",folder);

  printf("Please enter the file that contains the network:\n");
  scanf("%s",PPINetworkFile);
  printf("Please enter the file that contains mappings from IDs in %s to Gene Symbol IDs:\n",PPINetworkFile);
  scanf("%s",PPIIDFile);

  printf("Please enter the name of the file contain expression profiles:\n");
  scanf("%s", ExpressionFile);
  printf("Please enter the number of control samples:\n");
  scanf("%i", &numControls);
  printf("Please enter the number of case samples:\n");
  scanf("%i", &numCases);
  printf("Please enter the threshold to determine differentially expressed genes:\n");
  scanf("%lf", &diffThres);

  printf("Please enter the minimum density threshold:\n");
  scanf("%lf", &alpha);
  printf("Please enter the minimum number of case samples that a densely connected subnetwork should be present in:\n");
  scanf("%i", &L);
  */

  //./createGraph STRINGPPI.txt Normalized_10950.csv 0.1 TEST STRINGMAP1.txt 24 24
  printf("\n\nCREATING THE GRAPH STRUCTURES...\n");
  sprintf(command, "./createGraph %s %s %lf %s %s %d %d",PPINetworkFile,ExpressionFile,diffThres,folder,PPIIDFile,numControls,numCases);
  printf("%s\n",command);
  system(command);
  printf("\n\nEXTRACTING NETWORK MODULES...\n");
  //awk -f runWDCB.awk -v folder="TEST" -v density=0.5 -v minGraph=4 -v minSize=4
  sprintf(command, "awk -f runwDCB.awk -v folder=""%s"" -v density=%lf -v minGraph=%d -v minSize=%d",folder,alpha,L,4);
  printf("\n%s\n",command);
  system(command);

}
