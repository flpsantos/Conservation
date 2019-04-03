#!/usr/bin/env python

########################################################
# FILE: conbyproximity.py
# USAGE: conbyproximity.py
# DESCRIPTION: Cross orthologous ptcodinggenes to find ortholous lncRNAs
# OPTIONS: ---
# REQUIREMENTS: lncRNA-ptcoding genes pair (feelncRNA output) and ptcoding-ptcoding orthologous pairs by species
# BUGS: ---
# NOTES: ---
# AUTHOR: Felipe Rodolfo Camargo dos Santos
# COMPANY: Bioinformatics Lab - Hospital Sirio Libanes
# VERSION: 1.0
# CREATED: Tue Mar 12 14:30:15 -03 2019
# REVISION: ---
# LICENSE: GPL
########################################################

import argparse
import sys
import numpy as np  #Optimize speed
from multiprocessing import Process, Queue


#Getting  file names
filevars=["lncpairptAupstream", "lncpairptAdownstream", "lncpairptBupstream", "lncpairptBdownstream", "ptcpairptc"]  # variables that receive file names

parser = argparse.ArgumentParser(description='lncRNA conservation by proximity')

parser.add_argument('-aup', action = 'store', dest = filevars[0], help='lncRNA-ptcoding genes pair for organism A upstream')
parser.add_argument('-adown', action = 'store', dest = filevars[1], help='lncRNA-ptcoding genes pair for organism A downstream')
parser.add_argument('-bup', action = 'store', dest = filevars[2], help='lncRNA-ptcoding genes pair for organism B upstream')
parser.add_argument('-bdown', action = 'store', dest = filevars[3], help='lncRNA-ptcoding genes pair for organism B downstream')
parser.add_argument('-c', action = 'store', dest = filevars[4], help='ptcoding-ptcoding orthologous pairs by species')
parser.add_argument('--version')


print parser.parse_args()
entrada=parser.parse_args()

#################################Getting files###############################################


#One coding gene may have many closest lncRNAs associated to it. So dictionary key MUST be lncRNA
#Getting closest upstream/downstream lnc/protein coding proximity pair for organism A and B

corelist=[]
for i in range(0, len(filevars)-1):
	aux={}
	with open (eval("entrada."+filevars[i]), 'r') as lncpair:
	 	for line in lncpair:
			cols=line.split()
			if cols[3] not in aux:
				aux[cols[3]]=[]
			else:
				aux[cols[3]].append(cols[9]) 		#key:lncgene, ptcgene
	corelist.append(aux) 						#[dicupA, dicdownA, dicupB, dicdownB]

#Getting protein/protein orthologs for two organisms

corptc={}
with open (entrada.ptcpairptc, 'r') as ptcpair:
	for line in ptcpair:
		if line[0]!="G":
			cols=line.split()
			corptc[cols[0]]=cols[2] # key:lncgene, ptcgene



##################################Processing files##########################################

print ("Processing files...")


mergeA={}
for i in corelist[0]:
	for j in corelist[1]:
		if i==j:
			mergeA[i]=[corelist[0][i],corelist[1][j]]		#{lncA:[[upproteins],[downproteins]]


mergeB={}
for i in corelist[2]:
	for j in corelist[3]:
		if i==j:
			mergeB[i]=[corelist[2][i],corelist[3][j]]		#{lncB:[[upproteins],[downproteins]]


##########################Processing orthology protein as reference##################################

"""print ("Processing orthology...")


def procA(q,mergeA,corptc):
	ptcA={}
	for i in corptc:
		for j in mergeA:
			#upstream
			if i in mergeA[j][0] and i not in ptcA:
				ptcA[i]=()
				ptcA[i]=ptcA[i]+(j, mergeA[j][0].index(i))  #{Protein:(lncRNAA,distance)
			elif i in mergeA[j][0]:
				ptcA[i]=ptcA[i]+(j, mergeA[j][0].index(i))
			#downstream
			if i in mergeA[j][1] and i not in ptcA:
				ptcA[i]=()
				ptcA[i]=ptcA[i]+(j, -(mergeA[j][1].index(i)))
			elif i in mergeA[j][1]:
                        	ptcA[i]=ptcA[i]+(j, -(mergeA[j][1].index(i)))

	q.put(ptcA)


def procB(w,mergeB,corptc):
	ptcB={}
	for i in corptc:
		k=corptc[i]
		for j in mergeB:
                	#upstream
                	if k in mergeB[j][0] and k not in ptcB:
                        	ptcB[k]=()
                        	ptcB[k]=ptcB[k]+(j, mergeB[j][0].index(k))  #{Protein:(lncRNAA,distance)
                	elif i in mergeB[j][0]:
                        	ptcB[k]=ptcB[k]+(j, mergeB[j][0].index(k))
                	#downstream
                	if i in mergeB[j][1] and k not in ptcA:
                        	ptcB[k]=()
                        	ptcB[k]=ptcB[k]+(j, -(mergeB[j][1].index(k)))
                	elif i in mergeB[j][1]:
                        	ptcB[k]=ptcB[k]+(j, -(mergeB[j][1].index(k)))

	w.put(ptcB)


if __name__ == '__main__':
	q=Queue()
	p1=Process(target=procA,args=(q,mergeA,corptc))
	p1.start()
	w=Queue()
	p2=Process(target=procB,args=(w,mergeB,corptc))
	p2.start()
	p1.join()
	p2.join()

print w.get()
"""

#####################Processing orthology [lncRNA as reference]#################


for i in corptc:
	for j in mergeA:
		if i in mergeA[j][0]:
			for k in mergeB:
				if corptc[i] in mergeB[k][0]:
					score = mergeA[j][0].index(i) - mergeB[k][0].index(corptc[i])
					print j, k, score
				elif corptc[i] in mergeB[k][1]:
					score = mergeA[j][0].index(i) + mergeB[k][1].index(corptc[i])
					print j, k, score
		elif i in mergeA[j][1]:
			for k in mergeB:
				if corptc[i] in mergeB[k][0]:
					score = -(mergeA[j][1].index(i)) - mergeB[k][0].index(corptc[i])
                                        print j, k, score
                                elif corptc[i] in mergeB[k][1]:
					score = -(mergeA[j][1].index(i)) + mergeB[k][1].index(corptc[i])
                                        print j,k, score

#Getting file names


"""comporg=entrada.ptcpairptc.split('_')[0]
outname=comporg+".out"
print outname
#Compare
cons={}

dicup={}

for i in corptc: 										#iterate over {proteinA:proteinB}
	for j in corlncAup: 									#iterate over {lncRNAA:proteinA}
		if i == corlncAup[j]:
			for k in corlncBup:							#iterate over {lncRNAB:proteinB}
				if corptc[i] == corlncBup[k]:					#check if proteinA in {proteinA:proteinB} equals to proteinB in {lncRNAB:proteinB} 
					dicup[corlncAup[j]]=corlncBup[k]
#print dicup


dicdown={}

for i in corptc: 										#iterate over {proteinA:proteinB}
	for j in corlncAdown: 									#iterate over {lncRNAA:proteinA}
		if i == corlncAdown[j]:
			for k in corlncBdown:							#iterate over {lncRNAB:proteinB}
				if corptc[i] == corlncBdown[k]:					#check if proteinA in {proteinA:proteinB} equals to proteinB in {lncRNAB:proteinB} 
					dicdown[corlncAdown[j]]=corlncBdown[k]




#checking inversions

dicinv{}

for i in corptc:                                                                                #iterate over {proteinA:proteinB}
        for j in corlncAup:                                                                     #iterate over {lncRNAA:proteinA}
                if i == corlncAup[j]:
                        for k in corlncBdown:                                                   #iterate over {lncRNAB:proteinB}
                                if corptc[i] == corlncBdown[k]:                                 #check if proteinA in {proteinA:proteinB} equals to proteinB in {lncRNAB:proteinb}
                                        dicinv[corlncAup[j]]=corlncBup[k]


for i in dicinv:
	for j in corptc:
		if 

#with open (outname, "w") as out:
with open ("testedesaida", "w") as out:
	for i in dicup:
		out.write(i+"\t"+dicup[i]+"\n")
	for i in dicdown:
		out.write(i+"\t"+dicdown[i]+"\n")



def main():
	#Insert here the code#
	print ""

main();
"""

