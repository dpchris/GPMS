#!/usr/bin/env python
# -*- coding: utf-8 -*-
import regex, re

# fonction trim = s.strip() #rstrip pour le coté droit #lstrip pour le coté gauche

#fasta test (from insilico site)

with open("/home/david/Documents/fasta_example.txt","r") as myfile :
	myfile.readline()                                             #on supprime la premiere ligne
	fasta = myfile.read().replace('\n', '')                       #on retire les \n 

#dictionnary to create complementary DNA sequence
dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}

def inverComp (seq) :                          #return the inversed complementary sequence
	seq = seq.upper() 
	seq_comp = ""
	for nuc in seq :                       #nuc stand for nucleotide
		seq_comp += dico_comp[nuc]
	seq_comp = "".join(reversed(seq_comp)) #we reverse the sequence	
	return (seq_comp)

def mix (nb,primer,seq) : #ne marche pas
	findlist = []
	if nb==1 :
		invP=inverComp(primer)
		for letter in "ACGT" :
			i=0
			while findlist == [] and i<len(primer) :
				reg1=primer[:i]
				reg2=primer[i+1:]
				reg=reg1+letter+reg2
				tmp = seq.find(reg)
				if tmp !=-1 : 
					findlist.append([tmp,'norm'])
				else :
					reg1=invP[0:i]
					reg2=invP[i:]
					reg=reg1+letter+reg2
					tmp = seq.find(reg)
					if tmp != [] :
						findlist.append([tmp,'inv'])
				i+=1
				print findlist
	elif nb==2 :
		i=0
		for letter in "ACGT" :
			i=0
			while findlist == [] and i<len(primer) :
				reg1=primer[:i]
				reg2=primer[i+1:]
				reg=reg1+letter+reg2
				tmp = seq.find(reg)
				if tmp !=-1 : 
					findlist.append([tmp,'norm'])
				i+=1
	return findlist
			
#searchSeqREStr = seqStr + '|' + '|'.join(seqStr[:i]+"[ACTGN]".replace(c,'') +seqStr[i+1:]


"""		tmpfind = regex.findall("("+primer+"){e<=1}",seq,overlapped=True) #autorise 1 mismatch en plus ou en moins
		tmpfind = re.finditer(primer,seq)
		for m in tmpfind:
   			print m.start(), m.end(), m.group(0)
		if tmpfind !=[] :
			for i,match in enumerate(tmpfind) :
				if len(match)>len(primer) : #on ne garde pas le(s) match avec un nucleotide en plus
					del tmpfind[i] 
			findlist.append([seq.index(tmpfind[0]),'norm'])
		else : #idem avec le brin complementaire
			tmpfind = regex.findall("("+invP+"){e<=1}",seq,overlapped=True) #autorise 1 mismatch en plus ou en moins
			for i,match in enumerate(tmpfind) :
				if len(match)>len(invP) : #on ne garde pas le(s) match avec un nucleotide en plus
					del tmpfind[i]
			findlist.append([seq.index(tmpfind[0]),'inv'])
		return findlist"""

def findFirst (primer,seq) :
	indLeft = []
	sense = []
	indlook=0
	resul=seq.find(primer,indlook)
	while (resul!=-1) :                             
		indLeft.append(resul)             
		indlook=indLeft[-1]+1  
		sense.append("norm")
		resul=seq.find(primer,indlook) 
	
	if len(indLeft)==0 :                                  #si on a pas trouvé de match avec le primer 
		indlook=0
		tmpP=inverComp(primer)                        #on recupere le brin complementaire du primer avec la fonction incerComp
		resul=seq.find(tmpP,indlook)
		while(resul!=-1):                             #on cherche des match avec le brin complementaire du primer
			indLeft.append(resul)                 
			indlook=indLeft[-1]+1   
			sense.append("inv")
			resul=seq.find(tmpP,indlook) 
		print len(indLeft)
		
		if (len(indLeft)==0):                         #si on toujours pas de match
			lastTest=mix(1,primer,seq)         #mix : recherche avec un mismatch #arg 1 : primer + primer complementaire
			print lastTest			
			for res in lastTest :
				indLeft.append(res[0]) #on recupere la position du match
				sense.append(res[1])   #on recupere le sens du primer (norm ou inv)
			
		
	
	return (indLeft,sense);

def findSec(primer,seq,indLeft,sense) :
	indRight = []
	if sense=="norm" :
		primer=inverComp(primer) #brin complementaire du primer
		indlook=indLeft
		resul=seq.find(primer,indlook)
		while(resul!=-1) :
			indRight.append(resul)
			indlook=indRight[-1]+1
			resul=seq.find(primer,indlook)
		
		if len(indRight) ==0 :
			lastTest=mix(2,primer,seq); #mix avec arg 2 : on cherche avec un mismatch seulement avec le primer donné
			for res in lastTest :
				indRight.append(res[0])
	else :
		indlook=0
		resul=seq.find(primer,indlook)
		while(resul!=-1) :
			indRight.append(resul)
			indlook=indRight[-1]+1
			resul=seq.find(primer,indlook)
		
		if len(indRight)==0 :
			lastTest=mix(2,primer,seq) 
			for res in lastTest :
				indRight.append(res[0])
	return indRight


print mix(1,"CACCTTCAACACCTT",fasta)

#print findFirst ("CACCTTCNACACCTT",fasta)

#print findSec ("AAGGTGTTGAAGGTG",fasta,0,"inv")




