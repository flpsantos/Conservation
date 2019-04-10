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
import multiprocessing



#Parameter of lncRNA distance checkout

######
step=4
######



#Getting  file names
filevars=["lncpairptAupstream", "lncpairptAdownstream", \
"lncpairptBupstream", "lncpairptBdownstream", \
"ptcpairptc","outname"]  # variables that receive file names

parser = argparse.ArgumentParser(description='lncRNA conservation by proximity')

parser.add_argument('-aup', action = 'store', dest = filevars[0], help='lncRNA-ptcoding genes pair for organism A upstream')
parser.add_argument('-adown', action = 'store', dest = filevars[1], help='lncRNA-ptcoding genes pair for organism A downstream')
parser.add_argument('-bup', action = 'store', dest = filevars[2], help='lncRNA-ptcoding genes pair for organism B upstream')
parser.add_argument('-bdown', action = 'store', dest = filevars[3], help='lncRNA-ptcoding genes pair for organism B downstream')
parser.add_argument('-c', action = 'store', dest = filevars[4], help='ptcoding-ptcoding orthologous pairs by species')
parser.add_argument('-out', action = 'store', dest = filevars[5], help='output file name')
parser.add_argument('--version')

entrada=parser.parse_args()

#################################Getting files###############################################

print "Processing files..."

#Getting protein/protein orthologs for two organisms

corptc={}
with open (entrada.ptcpairptc, 'r') as ptcpair:

        for line in ptcpair:

                if line[0]!="G":

                        cols=line.split()
                        corptc[cols[0]]=cols[2] # key:lncgene, ptcgene


#One coding gene may have many closest lncRNAs associated to it. So dictionary key MUST be lncRNA
#Getting closest upstream/downstream lnc/protein coding proximity pair for organism A and B

corelist=[]
aux={}
for i in range(0, len(filevars)-2):

	with open (eval("entrada."+filevars[i]), 'r') as lncpair:

	 	for line in lncpair:

			col=line.split()
			protein=col[9]
			lnc=col[3]
			if i<=1:

				if protein in corptc:

					if lnc not in aux:
						aux[lnc]=[[],[]]
					else:
						aux[lnc][i].append(protein) 		#dic={lncgene: ptcgeneup1, ptcgeneup2,...ptcgenedownN-1,ptcgenedowN]
			else:

				if protein in corptc.values():

                                        if lnc not in aux:
                                                aux[lnc]=[[],[]]
                                        else:
                                                aux[lnc][i-2].append(protein)            #dic={lncgene: ptcgeneup1, ptcgeneup2,...ptcgenedownN-1,ptcgenedowN]

	if i==1:
		dic1=aux
		aux={}

dic2=aux 								#[dicA, dicB] 


##################Optimizing orthology algorithm, test1#######################

"""print "Processing orthology..."


def lookupup(q, dic1,dic2):
	up=[]
	for i in dic1:
		for index1 in range(0,len(dic1[i][0])):
			protein1=dic1[i][0][index1]
			protein2=corptc[dic1[i][0][index1]]
			for j in dic2:
				if protein2 in dic2[j][0]:
							up.append([i,j,abs(index1-dic2[j][0].index(protein2))])
				elif protein2 in dic2[j][1]:
							up.append([i,j,abs(index1-dic2[j][1].index(protein2))])
#	return (up)
	q.put(up)

def lookupdown(w,dic1,dic2):
	down=[]
        for i in dic1:
                for index1 in range(0,len(dic1[i][1])):
			protein1=dic1[i][1][index1]
                        protein2=corptc[dic1[i][1][index1]]
                        for j in dic2:
                                if protein2 in dic2[j][0]:
                                                        down.append([i,j,abs(index1-dic2[j][0].index(protein2))])
                                elif protein2 in dic2[j][1]:
                                                        down.append([i,j,abs(index1-dic2[j][1].index(protein2))])
#	return (down)
	w.put(down)

dic12=dic1
dic22=dic2

if __name__ == '__main__':
        q=Queue()
        p1=Process(target=lookupup,args=(q,dic1,dic2))
        p1.start()
        w=Queue()
        p2=Process(target=lookupdown,args=(w,dic12,dic22))
        p2.start()
        p1.join()
        p2.join()


res1=q.get()
res2=w.get()

print res1
print res2
"""


##################Optimizing orthology algorithm, test2#######################

print "Processing orthology..."


def lookupup(dic1,dic2, out_dic):
        up=[]
        for i in dic1:

                for index1 in range(0,len(dic1[i][0])):

                        protein1=dic1[i][0][index1]
                        protein2=corptc[dic1[i][0][index1]]

                        for j in dic2:

                                if protein2 in dic2[j][0]:
                                                        up.extend([i,j,abs(index1-dic2[j][0].index(protein2))])
                                elif protein2 in dic2[j][1]:
                                                        up.extend([i,j,abs(index1-dic2[j][1].index(protein2))])
	out_dic["luup"]=up


def lookupdown(dic1,dic2, out_dic):
        down=[]
        for i in dic1:
                for index1 in range(0,len(dic1[i][1])):
                        protein1=dic1[i][1][index1]
                        protein2=corptc[dic1[i][1][index1]]
                        for j in dic2:
                                if protein2 in dic2[j][0]:
                                                        down.extend([i,j,abs(index1-dic2[j][0].index(protein2))])
                                elif protein2 in dic2[j][1]:
                                                        down.extend([i,j,abs(index1-dic2[j][1].index(protein2))])
	out_dic["ludown"]=down


if __name__ == '__main__':
	manager = multiprocessing.Manager()
	out_dic = manager.dict()
        p1=multiprocessing.Process(target=lookupup,args=(dic1,dic2, out_dic))
	p1.start()
        p2=multiprocessing.Process(target=lookupdown,args=(dic1,dic2, out_dic))
        p2.start()
        p1.join()
        p2.join()


##########################computing scores###################################

list1=out_dic["ludown"]
list2=out_dic["luup"]

score={}

for i in range (0,len(list1)-3,3):

	if list1[i] not in score:
		score[list1[i]]=[list1[i+1],list1[i+2]]
	elif list1[i+1] not in score[list1[i]]:
		score[list1[i]].extend([list1[i+1],list1[i+2]])
	else:
		index=score[list1[i]].index(list1[i+1])
		score[list1[i]][index+1]=score[list1[i]][index+1]+list1[i+2]

for i in range (0,len(list2)-3,3):

        if list2[i] not in score:
                score[list2[i]]=[list2[i+1],list2[i+2]]
        elif list2[i+1] not in score[list2[i]]:
                score[list2[i]].extend([list2[i+1],list1[i+2]])
        else:
                index=score[list2[i]].index(list2[i+1])
                score[list2[i]][index+1]=score[list2[i]][index+1]+list2[i+2]


##########################Getting bests####################################


filtered={}
for i in score:

	aux=[]
	for j in range(0,len(score[i])):

		if j%2!=0:
			aux.append(score[i][j])

	minimo=min(aux)
	key=i+"_"+str(minimo)

	filtered[key]=[]

	for t in range(0,len(score[i])):
		if score[i][t]==minimo:
			filtered[key].append(score[i][t-1])


#######################Returning results#################################


with open (entrada.outname, "w") as out:
	for i in filtered:
		for j in range(0, len(filtered[i])):
			out.write(i.split("_")[0]+"\t"+str(filtered[i][j])+"\t"+str(i.split("_")[1])+"\n")

