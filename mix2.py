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

def positionsOfMatches (result,seq,nbmismatch) : #get the matches positions in the fasta sequence
	pos = []  
	for res in result :
		pos.append([seq.find(res),res,nbmismatch])
	return (pos)

def search_matches(nbmismatch,primer, seq) : #return all match(es) in the fasta sequence 
	return (regex.findall("("+primer+"){e<="+str(nbmismatch)+"}",seq,overlapped=True))

def mix2 (nb,primer,seq) : #find matches with 0,1 or 2 mismatch for primer in forward and reverse sense
	listFind =[]
	if nb==1 :
		for nb_mismatch in xrange(3) :
			tmpfind_forward = search_matches(nb_mismatch,primer,seq) 		#get all matches
			positions_f = positionsOfMatches(tmpfind_forward,seq,nb_mismatch)	#get positions of matches
			if tmpfind_forward == [] :
				invP=inverComp(primer)
				tmpfind_reverse = search_matches(nb_mismatch,invP,seq)
				positions_r = positionsOfMatches(tmpfind_reverse,seq,nb_mismatch)
				if tmpfind_reverse != [] :
					for res in positions_r :
						listFind.append([res[0],"inv"])
					break
			else :
				for res in positions_f :
					listFind.append([res[0],"norm"])
				break
			
	elif nb==2 :
		for nb_mismatch in xrange(3) :
			tmpfind_forward = search_matches(nb_mismatch,primer,seq) 		#get all matches
			positions_f = positionsOfMatches(tmpfind_forward,seq,nb_mismatch)	#get positions of matches
			if tmpfind_forward != [] :
				for res in positions_f :
					listFind.append([res[0],"norm"])
				break
	return(listFind)

#Brucella_suis_NZ_CP009096.1_NZ_CP009097.1.fa
# CAGTATGATTTGTCGGAAGCCG -> 2 mismatch : CAGTATGATATGTCGGATGCCG
#Bruce42-424_125bp_539bp_4u	CATCGCCTCAACTATACCGTCA	ACCGCAAAATTTACGCATCG
with open("/home/david/Documents/complete_genomes/brucellaceae/Brucella_suis_NZ_CP009096.1_NZ_CP009097.1.fa","r") as myfile :  #load of the fasta file
	myfile.readline()                                             #we don't keep the header
	tmp = myfile.read().split(">")
	chr1 = tmp[0].replace("\n","")
	chr2 = tmp[1].replace("\n","")
	fasta = chr1+chr2

print mix2(1,inverComp("CAGTATGATAAGTCGGAAGCCG"),fasta)
#print search_matches(2,inverComp("ACCGCAAAATTTACGCATCG"),fasta)





