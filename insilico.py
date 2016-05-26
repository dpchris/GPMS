#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, math

#fasta test (from insilico site)

with open("/home/david/Documents/complete_genomes/brucellaceae/Brucella_abortus_104M_NZ_CP009625.1_NZ_CP009626.1.fa","r") as myfile :  #load of the fasta file
	#myfile.readline()                                             #we don't keep the header
	fasta = myfile.read()

with open("/home/david/Documents/complete_genomes/brucellaceae/Brucella_abortus_104M_NZ_CP009625.1_NZ_CP009626.1.fa","r") as myfile :  #load of the fasta file
	myfile.readline()                                             #we don't keep the header
	fasta2 = myfile.read().split(">")[0].replace("\n","")

Primer ="Bruce06-1322_134bp_408bp_3u	ATGGGATGTGGTAGGGTAATCG	GCGTGACAATCGACTTTTTGTC"

#dictionnary to create complementary DNA sequences
dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}

def inverComp (seq) :                              #return the inversed complementary sequence
	seq = seq.upper()                          #acgt -> ACGT
	seq_comp = ""
	for nuc in seq :                           #nuc stand for nucleotide
		seq_comp += dico_comp[nuc]         #use of the dictionnary to have the complementary nucleotide
	seq_comp = "".join(reversed(seq_comp))     #we reverse the sequence	
	return (seq_comp)
"""
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
"""			

def mix (nb,primer,seq) : # return matches (and position) of the primer on the sequence (complementary matches only if none with not inversed primer) 
	if nb==1 :
		tmpfind =[]
		invP=inverComp(primer)                                       #get the complementary inversed primer
		#tmpfind = regex.findall("("+primer+"){e<=1}",seq,overlapped=True) #allow one mismatch (one more or less nucleotide)
		searchseq = primer + '|' + '|'.join(primer[:i]+"["+""\
			.join(dico_comp.keys()).replace("|","").replace(".","")\
			.replace(nuc,'')+"]"+primer[i+1:] for i,nuc in enumerate(primer)) #request for finditer which allow a mismatch
		searchseq = re.compile(searchseq)                            #you need to compil the request for finditer
		tmpfinditer = searchseq.finditer(seq)                        #get the matches and position of it

		if tmpfinditer !=[] :
			for m in tmpfinditer:
   				tmpfind.append([m.start(), m.end(), m.group(0),'norm'])	

		else : 	#same search with the complementary primer
			#tmpfind = regex.findall("("+invP+"){e<=1}",seq,overlapped=True) #allow one mismatch (one more or less nucleotide)
			searchseq = invP + '|' + '|'.join(invP[:i]+"["+"".join(dico_comp.keys()).replace("|","").replace(".","").replace(nuc,'')+"]"+invP[i+1:] for i,nuc in enumerate(invP))
			searchseq = re.compile(searchseq)
			tmpfinditer = searchseq.finditer(seq)
			if finditer != [] :
				for m in tmpfinditer:
   					tmpfind.append([m.start(), m.end(), m.group(0),'inv'])
		return tmpfind #return matches in a list (example : [0, 4, 'CACC', 'inv'] or [13, 16, 'TTA', 'norm'] )

	elif nb==2 :  # search only with the 
		tmpfind =[]                                    
		searchseq = primer + '|' + '|'.join(primer[:i]+"["+""\
			.join(dico_comp.keys()).replace("|","").replace(".","")\
			.replace(nuc,'')+"]"+primer[i+1:] for i,nuc in enumerate(primer)) #request for finditer which allow a mismatch
		searchseq = re.compile(searchseq)                            #you need to compil the request for finditer
		tmpfinditer = searchseq.finditer(seq)                        #get the matches and position of it

		if tmpfinditer !=[] :
			for m in tmpfinditer:
   				tmpfind.append([m.start(), m.end(), m.group(0),'norm'])	
		
		
		return  [tmpfind[0]] #we keep only the first match (as the javascript version)

def findFirst (primer,seq) : #first search of the primer on the sequence (use inverComp() and mix())
	indLeft = []
	sense = []           #to store the sense of search (normal or inversed) 
	indlook=0
	resul=seq.find(primer,indlook)         #get the first match (if it exists)
	while (resul!=-1) :                    #while the search has not been made on the entire sequence               
		indLeft.append(resul)             
		indlook=indLeft[-1]+1          #next search will start one nucleotide after the position of the last match
		sense.append("norm")
		resul=seq.find(primer,indlook) #next search, if no more match : resul = -1 -> end of the loop
	
	if len(indLeft)==0 :                   #if no match foud with the regular primer 
		indlook=0
		tmpP=inverComp(primer)         #get the inversed complementary primer with inverComp()
		resul=seq.find(tmpP,indlook)
		while(resul!=-1):              #same search with the converted primer
			indLeft.append(resul)                 
			indlook=indLeft[-1]+1   
			sense.append("inv")
			resul=seq.find(tmpP,indlook) 
		
		if (len(indLeft)==0):                  #if still no matches
			lastTest=mix(1,primer,seq)     #search with a mismatch with mix()		
			for res in lastTest :
				indLeft.append(res[0]) #get the position of the match
				sense.append(res[3])   #get the sens of the match
			
		
	
	return (indLeft,sense)

def findSec(primer,seq,indLeft,sense) : #search a match for the second primer
	indRight = []
	if sense=="norm" :
		primer=inverComp(primer) 
		indlook=indLeft                         #start the search at the match position of the first primer
		resul=seq.find(primer,indlook)
		while(resul!=-1) :
			indRight.append(resul)
			indlook=indRight[-1]+1
			resul=seq.find(primer,indlook)
		
		if len(indRight) ==0 :                  # if no match found : try with one mismatch
			lastTest=mix(2,primer,seq)     #mix with arg 2 : search only with the given primer, not the complementary
			if lastTest != [] :
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

def find(primers,text,round) : #main function : return the result of the matches 
	
	text = text.replace(" ","").replace("\t","")
	primers = primers.split("\n")
	
	if type(round) is str :
		round = round.replace(",",".") 
		round = float(round)    # parsefloat : transforme un string en float
	
	#primers=primers.filter(function(e){return e}) #??
	#primersSNP=primersSNP.filter(function(e){return e}) #??
	#var sizeSNP=parseInt(GetId('sizeSnp').value)
	
	files = text.split('>')  #split par > pour séparer les fasta
	del files[0]             # delete the '' created 
	tabRes = {}

	for f in xrange(len(files)) :                  # for each fasta sequence 
		text=">"+files[f]                      
		tmpSeq=text.split('\n')                
		titleSeq=tmpSeq[0]                     #get the name of the fasta file
		#titleSeq2=tmpSeq[0]
		del tmpSeq[0]		               #del the header of the fasta file to keep only the sequence
		text=''.join(tmpSeq).upper()           #merge the sequence list iinto a single string

		if primers :                 # si les primers on été rentré
			for i in xrange(len(primers)) :
				primer = primers[i].replace(" ",";").replace("\t",";").split(";")
				detPrimer = primer[0].split('_')
				#print primer[1]
				#print text[:100]
				tabIndLeft = findFirst(primer[1],text) #search match with findFirst() : primer[1]=primer, text=fasta, primer[0]= primer name
				tabIndRight = []
				#print tabIndLeft
				
				resul = []
				a = tabIndLeft
				for ind in xrange(len(tabIndLeft[0])) : #parse both lists of positions and sense of matches (primer1)
					tabIndRight.append(findSec(primer[2],text,tabIndLeft[0][ind],tabIndLeft[1][ind]))
					#print tabIndRight
					
					if tabIndRight != [] :
						for ind2 in xrange(len(tabIndRight[0])) : #parse both lists of positions and sense of matches (primer2)
							if tabIndLeft[1][ind]=="inv" :
								size = abs(tabIndLeft[0][ind]-tabIndRight[ind][ind2]+len(primer[1]))

							else :
								size = abs(int(tabIndRight[ind][ind2])-int(tabIndLeft[0][ind]+len(primer[2])))
							sizeU = abs(int(detPrimer[3].upper().replace("U",""))-((int(detPrimer[2].lower().replace("bp",""))-size)/int(detPrimer[1].lower().replace("bp",""))))
							resul.append([primer[0],tabIndLeft[0][ind],tabIndRight[ind][ind2],size,sizeU])
							#print resul
				#print detPrimer[0]

				if len(resul) == 0 and detPrimer[0] not in tabRes :
					tabRes[primer[0]]=[detPrimer[0],"","","",""]
				#print tabRes

				if (len(resul)>0) :
					ind=0
					for j in xrange(len(resul)) :
						if (resul[ind][4]>resul[j][4]) : ind=j

					if round !="" and round>0 :
						sizeU=resul[ind][4]
						#print sizeU
						if sizeU>=math.floor(sizeU) and sizeU<(math.floor(sizeU)+round) :
							sizeU = math.floor(sizeU)
						elif sizeU <= math.ceil(sizeU) and sizeU>(math.ceil(sizeU)-round) :
							sizeU=math.ceil(sizeU)
						else :
							sizeU=math.floor(sizeU)+0.5

						resul[ind][4]=sizeU
						#print resul
					tabRes[detPrimer[0]]=resul[ind]
					#print tabRes 
					#print primers[i]
					del primers[i]
	return tabRes


#print inverComp("TATAGAGACACA")

#print mix(1,"TGT","ACTGACGACCAYACAAGTTAC")

#print findFirst ("ATGGGATGTGGTAGGGTAATCG",fasta2) # -> ([2122861, 2122881, 2122892, 2122904], ['norm', 'norm', 'norm', 'norm'])

#print findSec ("GAC","ACTGACGACCAYACAAGTTAC",0,"norm") 
#print mix(2,"GAC","ACTGACGACCAACAAGTTAC")


print find(Primer,fasta,0.25)
#print findFirst ("ATGGGATGTGGTAGGGTAATCG",fasta2) # -> ([2122861, 2122881, 2122892, 2122904], ['norm', 'norm', 'norm', 'norm'])
#print findSec ("GCGTGACAATCGACTTTTTGTC",fasta2,684708,"norm") 



