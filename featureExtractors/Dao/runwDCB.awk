BEGIN{
		#print "Folder:",folder;
		#print "Density:",density;
		#print "Threshold:",threshold;
		#print "Min Patients:",minGraph;
		#print "Min Subnetwork Size:",minSize;

		while (match(density,"0$")){
            density=substr(density,1,length(density)-1);
            #print density"\n";
		}


		numGenes=0;
		mapFile=sprintf("%s/nodes.txt",folder);
		while ((getline line < mapFile)>0){
			split(line,c,"\t");
			if (length(c)>=2)
			{
				id2g[c[1]]=c[2];
				if (c[1]>numGenes)
				{
					numGenes=c[1];
					bad[c[1]]=0;
				}
			}
		}
		close(mapFile);

		totalNumEdges=0;

		edgeFile=sprintf("%s/edges.txt",folder);
		while ((getline line < edgeFile)>0){
			split(line,c,"\t");
			totalNumEdges++;
			g[c[1],c[2]]=c[3];
			g[c[2],c[1]]=c[3];
			numEdges+=c[3];
		}
		close(edgeFile);

		directory="experimental_results";
		#makingDirectory=sprintf("mkdir %s/%s",folder,directory);
		#system(makingDirectory);

		# run weighted DCB
		runningWeightedDCB=sprintf("./NCSC -dir %s/ -density %s -minGraph %d -minSize %d -distFunc 2",folder,density,minGraph,minSize);
		system(runningWeightedDCB);
		#print runningWeightedDCB;

		# create specific directory according to graphFile to store the results
		directory="WDCB";


		# moving the result file
		resultFile=sprintf("alpha_%s_minGraph_%d_minSize_%d_distFunc_2.txt",density,minGraph,minSize);
		movingResultFile=sprintf("mv %s/experimental_results/%s %s/%s/",folder,resultFile,folder,directory);
		system(movingResultFile);
		#print movingResultFile;

		# parse the result file to get the modules
		fullPathResultFile=sprintf("%s/%s/%s",folder,directory,resultFile);

		while ((getline line < fullPathResultFile)>0){
			if (line~" Final number of patterns:"){
				split(line,a,":");
				getline line < fullPathResultFile;
				#print a[2] > fullPathDensityFile;
				for (k=1;k<=a[2];k++)
				{
					getline line < fullPathResultFile;
					if (!(line ~ "Weight = "))
						break;
					moduleFile=sprintf("%s/%s/modules/%d",folder,directory,k);

					split(line,b," ");
					for (l=1;l<=length(b)-3;l++)
					{
						print id2g[b[l]] > moduleFile;
					}
					close(moduleFile);

					totalModuleEdges=0;
					numModuleEdges=0;
					for (i=1;i<length(b)-3;i++)
					for (j=i+1;j<=length(b)-3;j++)
					if (g[b[i],b[j]]>0)
					{
						totalModuleEdges++;
						numModuleEdges+=g[b[i],b[j]];
					}

					getline line < fullPathResultFile;
					getline line < fullPathResultFile;
					getline line < fullPathResultFile;
					split(line,b,"[-\t]+");

					getline line < fullPathResultFile;
				}

				break;
			}
		}
		close(fullPathResultFile);
		#close(fullPathDensityFile);

		# running ranking modules
		print "\n\nRANKING MODULES...\n";
		runRanking=sprintf("./rankModules %s WDCB %s",folder,resultFile);
		print runRanking;
		system(runRanking);
		#cleaningDirectory=sprintf("rm %s/%s/modulesDim/*",folder,directory);
		#system(cleaningDirectory);

		while ((getline line < fullPathResultFile)>0){
			#print line;
			if (line~" Final number of patterns:"){
				split(line,a,":");
				getline line < fullPathResultFile;
				for (k=1;k<=a[2];k++)
				{
					getline line < fullPathResultFile;
					if (!(line ~ "Weight = "))
						break;
					moduleFile=sprintf("%s/%s/modules/%d",folder,directory,k);
					split(line,b," ");
					for (l=1;l<=length(b)-3;l++)
					{
						#print b[l];
						print id2g[b[l]] > moduleFile;
						bad[b[l]]=1;
					}
					close(moduleFile);


					getline line < fullPathResultFile;
					getline line < fullPathResultFile;
					getline line < fullPathResultFile;
					split(line,b,"[-\t]+");
					getline line < fullPathResultFile;
				}

				break;
			}
		}
		close(fullPathResultFile);
}
