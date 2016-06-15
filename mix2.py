#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, regex

dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}

def inverComp (seq) :                              #return the inversed complementary sequence
	seq = seq.upper()                          #acgt -> ACGT
	seq_comp = ""
	for nuc in seq :                           #nuc stand for nucleotide
		seq_comp += dico_comp[nuc]         #use of the dictionnary to have the complementary nucleotide
	seq_comp = "".join(reversed(seq_comp))     #we reverse the sequence	
	return (seq_comp)

def positionsOfMatches (result,seq) : #get the matches positions in the fasta sequence
	pos = []  
	for res in result :
		pos.append([seq.find(res),res])
	return (pos)

def search_matches(nbmismatch,primer, seq) : #return all match(es) in the fasta sequence 
	return (regex.findall("("+primer+"){e<="+str(nbmismatch)+"}",seq,overlapped=True))

def mix (nb,primer,seq) : #function to find match(s) with a mismatch
	listFind =[]
	if nb==1 :
		invP=inverComp(primer)
		#use of findtiter first for one mismatch (faster)
		for i in range(len(primer)): #for each nucleotide of the primer
			reg1=primer[:i]
			reg2=primer[i+1:]
			reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2)
			searchseq = re.compile(reg)            #search request for finditer (Primer with a mismatch)
			tmpfinditer = searchseq.finditer(seq)  #tmpfinder : objects list with start() end() and group(0) functions
			tmp=[]
			for m in tmpfinditer :                 #for each result
				tmp.append([m.start(),"norm"]) #add the result in the tmp list
			if tmp != [] :                         #if results, we add the first result in listFind
				listFind.append(tmp[0])	
			else :                                 #same search with the complementary reversed primer
				reg1=invP[:i]
				reg2=invP[i+1:]
				reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2)
				searchseq = re.compile(reg)
				tmpfinditer = searchseq.finditer(seq)
				for m in tmpfinditer :
					tmp.append([m.start(),"inv"])
				if tmp != [] :
					listFind.append(tmp[0])	
			#if no match with one mismatch, use of regex.findall for two mismatch
			if listFind == [] :
				tmpfind_forward = search_matches(2,primer,seq) 		#get all matches (two mismatches)
				positions_f = positionsOfMatches(tmpfind_forward,seq)	#get positions of matches
				if tmpfind_forward == [] :
					tmpfind_reverse = search_matches(2,invP,seq)
					positions_r = positionsOfMatches(tmpfind_reverse,seq)
					if tmpfind_reverse != [] :
						for res in positions_r :
							listFind.append([res[0],"inv"])
				else :
					for res in positions_f :
						listFind.append([res[0],"norm"])
		
		
	elif nb==2 :					       #search only in the forward sense
		for i in range(len(primer)):
			reg1=primer[:i]
			reg2=primer[i+1:]
			reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2)
			searchseq = re.compile(reg)
			tmpfinditer = searchseq.finditer(seq) 
			if tmpfinditer !=[] :
				tmp=[]
				for m in tmpfinditer :
					tmp.append([m.start(),"norm"])
				if tmp != [] :
					listFind.append(tmp[0])	
			#if no match with one mismatch
			if listFind == [] :
				tmpfind_forward = search_matches(2,primer,seq) 	#get all matches
				positions_f = positionsOfMatches(tmpfind_forward,seq)		#get positions of matches
				if tmpfind_forward != [] :
					for res in positions_f :
						listFind.append([res[0],"norm"])
	return(listFind)

#Brucella_suis_NZ_CP009096.1_NZ_CP009097.1.fa
#Brucella_vulpis_NZ_LN997863.1_NZ_LN997864.1.fa
#Bruce11-211_63bp_257bp_2u	CTGTTGATCTGACCTTGCAACC	CCAGACAACAACCTACGTCCTG
#Bruce42-424_125bp_539bp_4u	CATCGCCTCAACTATACCGTCA	ACCGCAAAATTTACGCATCG
#Bruce22-322_8bp_158bp_6u	GATGAAGACGGCTATCGACTG	TAGGGGAGTATGTTTTGGTTGC

with open("/home/david/Documents/complete_genomes/brucellaceae/Brucella_suis_NZ_CP009096.1_NZ_CP009097.1.fa","r") as myfile :  #load of the fasta file
	myfile.readline()                                             #we don't keep the header
	tmp = myfile.read().split(">")
	chr1 = tmp[0].replace("\n","")
	chr2 = tmp[1].replace("\n","")
	fasta = chr1+chr2

print mix(1,"CCAGACAACAACCTACGTCCTG",fasta)
#print mix2(1,"CCAGACAACAACCTACGTCCTG",fasta)
#print search_matches(2,inverComp("ACCGCAAAATTTACGCATCG"),fasta)



