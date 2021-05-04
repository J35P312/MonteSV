import sys
import random
import os
import numpy
import copy

#Number of breakpoints
n_bp=6

#number of observed fusions
n_fusion=3

#observed average distance
d_obs=(160535951-16322795)/6

chromosome_min={}
chromosome_max={}

#chromosome of the rearrangement
autosomes=["6"]

p_dens=0
p_disp=0

for line in open("Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai"):
	content=line.strip().split()
	chromosome_min[content[0]]=0
	chromosome_max[content[0]]=int(content[1])

n=10000
p=0.0
n_ok=0
for m in range(0,n):
	pos_list=[]
	chromosomes=numpy.random.choice(autosomes,size=n_bp)
	k=0

	fragments=["A","B","C","D","E","F","G"]
	pos_to_gene={}
	genes_start={}
	genes_end={}

	pos_start={}
	pos_end={}

	i=0	
	for chromosome in chromosomes:
		pos=random.randint( chromosome_min[chromosome] , chromosome_max[chromosome] )
		os.system("tabix Homo_sapiens.GRCh37.87.clean.gtf.gz {}:{}-{} > tmpB.txt".format(chromosome,pos,pos+1 ))
		gene=False
		pos_list.append(pos)
		for line in open("tmpB.txt"):
			gene=line.split("gene_name")[-1].split(";")[0]

			break
		pos_to_gene[pos]=gene
		i+=1

	i=0
	for pos in sorted(pos_list):
		pos_end[fragments[i]]=pos
		pos_start[fragments[i+1]]=pos+1
		genes_end[fragments[i]]=pos_to_gene[pos] 
		genes_start[fragments[i+1]]=pos_to_gene[pos]

		i+=1
		

	genes_start["A"]=False
	genes_end["G"]=False

	pos_start["A"]=1
	pos_end["G"]=chromosome_max[chromosome]

	orientations=["+","-"]
	to_pick=set(["B","C","D","E","F"])
	rea=["A"]
	orientation=["+"]
	selected=set([])

	discard=False
	for i in range(0,len(to_pick)):
		tries=0
		while True:
			j=random.randint(0,5)
			k=random.randint(0,1)
			if fragments[j] == "A":
				continue

			#print("")
			#print(fragments[j])
			#print(rea[i])

			if (fragments[j-1] == rea[i]) and orientations[k] == "+" and orientation[i] == "+":
				continue

			if (fragments[j+1] == rea[i]) and orientations[k] == "-" and orientation[i] == "-":
				continue

			if not fragments[j] in selected:
				selected.add(fragments[j])
				rea.append(fragments[j])
				orientation.append(orientations[k])
				break
			tries+=1
			if tries > 30:
				discard=True
				break


	if discard:
		continue

	if rea[-1] == "F" and orientation[-1] == "+":
		continue

	rea.append("G")
	orientation.append("+")

	print("rea:" + str(m))

	fusions=[]
	for i in range(0,len(rea)-1):

		if orientation[i] == "+":
			A=genes_end[rea[i]]
		else:
			A=genes_start[rea[i]]

		if orientation[i+1] == "+":
			B=genes_start[rea[i+1]]
		else:
			B=genes_end[rea[i+1]]


		if A and B:
			fusions.append("{}-{}:{}-{}:".format(rea[i],rea[i+1],A,B) )

	print(" ".join(rea))
	print(" ".join(orientation))

	for i in range(len(rea)):
		if genes_start[rea[i]]:
			s=genes_start[rea[i]]
		else:
			s="none"

		if genes_end[rea[i]]:
			e=genes_end[rea[i]]
		else:
			e="none"
		#print(s)
		#print(e)
		#print(rea[i])
		print ("{}: {}_{} {}-{}".format(rea[i],s,e,pos_start[rea[i]],pos_end[rea[i]] ))

	print()		
	print("\n".join(fusions))
	print( (max(pos_list)-min(pos_list) )/len(pos_list) )
	density=(max(pos_list)-min(pos_list) )/len(pos_list)

	n_ok=n_ok+1
	if len(fusions) >= n_fusion:
		p+=1
	if density <= d_obs:
		p_dens+=1		
	else:
		p_disp+=1


print ("P fusion :" + str(float(p)/n_ok))
print ("P dens :" + str(float(p_dens)/n_ok))
print ("P disp :" + str(float(p_disp)/n_ok))
