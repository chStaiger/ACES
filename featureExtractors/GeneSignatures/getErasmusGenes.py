# @Author 
# Christine Staiger
# staiger@cwi.nl; staigerchristine@gmail.com

filename = 'ErasmusMC_Wang_genes.csv'
f = open(filename)
lines = f.readlines()
probes = []
for line in lines:
    if line[0].isdigit():
        probes.append(line.split(',')[0])


