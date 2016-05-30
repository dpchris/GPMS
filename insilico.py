#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, math

#fasta test (brucellaceae)
with open("/home/david/Documents/complete_genomes/brucellaceae/Brucella_abortus_104M_NZ_CP009625.1_NZ_CP009626.1.fa","r") as myfile :  #load of the fasta file
	#myfile.readline()                                             #we don't keep the header
	fasta = myfile.read()

with open("/home/david/Documents/complete_genomes/brucellaceae/Brucella_abortus_104M_NZ_CP009625.1_NZ_CP009626.1.fa","r") as myfile :  #load of the fasta file
	myfile.readline()                                             #we don't keep the header
	fasta2 = myfile.read().split(">")[0].replace("\n","")

Primers = open("/home/david/downloads/primers_brucella","r")
Primers = Primers.read()

fasta3="GTCTTGCGGCGTCCAGGGTTTGAGTTAGGAACACGCCATGACGGACAGCAGTGCAAATCCCCGTAATTAC\
GGCGAGCGTCCGGTGCACGAAGATGATCCGCTGATGGAACTTTCGCGGATTATGGACTTCGACACGCCTG\
CTGATGATAACGTCGCTCGGAATGAGCGCCGCCATGACAGTCAATTCGAGGATCAGGGCCGAGCAGAGCC\
GCGCTTTGATTCGGCGCAGGATGACCCCTCTTTTGATCCCGTTCTCGATCTTGAGCGCGAGCTTATGGGG\
CATTTCGACGATTATACGCAGCCTACAGCTCATAGCGAAACCGTTTCGACCGCAACCTTTGGCGTGGACG\
GGGAGCGGGCCTATGGCGAGCAGTCTCCGTTGGAGGAGGATGCGTTCGCTGCCGCGCTGGAAGAAGAATT\
CGATCTCGACCTTGGCTCGGCAGAAGCAGCACCGTCACCGGAAATTTTCGAGTTCGAGGACGTGGATCGC\
TCCGAGCCCGTTCCGCAATTCGACTATAACGATTATTCGCAGGCTCGTGATTCGAGCGAACACGCCCATA\
TGGCATCACCTGCGGCTCATGATGATCGCCAGCAGCCGATAGAGGCAGATGCGGGCGAATATGATCCGAC\
CGTCTATCAGCCTGCCATTGAGCCGGGTGAGCAAGCTTGGCGTGAAGATACCGGCGTGCAGGACAATTGG\
CAGGCACAGGATGAGTGGCCTGCACAGAATGACAGGCTGGAACAGAATGATTGGTCAGAACAGAATGACT\
GGCCAGCACAGAATAATTGGAATGCGCAGCAGCCTGTAGATGCTCCCGCTACGCATTTCGCGGCTGAGCC\
TGCACAGCAACCGCTTTCGCTTGAAGACGAGCTGGAAAACCTTCTGTTTGGCGATGAGCCGCAGCCTGTG"


#dictionnary to create complementary DNA sequences
dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}

def inverComp (seq) :                              #return the inversed complementary sequence
	seq = seq.upper()                          #acgt -> ACGT
	seq_comp = ""
	for nuc in seq :                           #nuc stand for nucleotide
		seq_comp += dico_comp[nuc]         #use of the dictionnary to have the complementary nucleotide
	seq_comp = "".join(reversed(seq_comp))     #we reverse the sequence	
	return (seq_comp)

def mix (nb,primer,seq) : # return best match (and position) of the primer on the sequence (complementary matches only if none with not inversed primer) 
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
		#print tmpfind
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
		if tmpfind == [] :
			return tmpfind
		else :
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
			#print "use of findfirst mix(1)",name,primer
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
		
		if len(indRight) ==0 :	                # if no match found : try with one mismatch
			#print "use of findSec mix(2)"
			lastTest=mix(2,primer,seq)      #mix with arg 2 : search only with the given primer, not the complementary
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
			#print "use of findSec mix(2)"
			lastTest=mix(2,primer,seq) 
			for res in lastTest :
				indRight.append(res[0])
	return indRight

def find(primers,text,round) : #main function : return the result of the matches 
	
	text = text.replace(" ","").replace("\t","") #delete spaces and tabulations
	primers = primers.split("\n") 
	del primers[-1]               
	#print primers
	if type(round) is str :
		round = round.replace(",",".") 
		round = float(round)    # float() : turn a string into a float 
	
	files = text.split('>')  #split par > pour séparer les fasta
	del files[0]             # delete the '' created
	tabRes = {}

	for f in xrange(len(files)) :                  # for each fasta sequence 
		text=">"+files[f]                      
		tmpSeq=text.split('\n')                
		titleSeq=tmpSeq[0]                     #get the name of the fasta file
		del tmpSeq[0]		               #del the header of the fasta file to keep only the sequence
		text=''.join(tmpSeq).upper()           #merge the sequence list into a single string

		if primers :                           # si les primers on été rentré
			for i in xrange(len(primers)) :
				print primers[i]
				primer = primers[i].replace(" ",";").replace("\t",";").split(";") #split the primers name, and primers into a list
				detPrimer = primer[0].split('_')        #split differents elements of the name into a list
				tabIndLeft = findFirst(primer[1],text,detPrimer[0])  #search match with findFirst() : primer[1]=primer, text=fasta, primer[0]= primer name
				tabIndRight = []
				resul = []
				print tabIndLeft[0]
				for ind in xrange(len(tabIndLeft[0])) : #parse result(s) of findfirst to do the second search with findsec()
					tabIndRight.append(findSec(primer[2],text,tabIndLeft[0][ind],tabIndLeft[1][ind])) #get result(s) of findsec with the second primer
					print len(tabIndRight)
					if tabIndRight != [] :
						for ind2 in xrange(len(tabIndRight[0])) : #parse both lists of positions and sense of matches (primer2)
							if tabIndLeft[1][ind]=="inv" :    #if the first match is made with the reversed complementary first primer
								size = abs(tabIndLeft[0][ind]-tabIndRight[ind][ind2]+len(primer[1])) #|position of first match minus position of second match|
								#print tabIndLeft[0]
												
							else :
								size = abs(int(tabIndRight[ind][ind2])-int(tabIndLeft[0][ind])+len(primer[2])) #|position of second match minus position of first match|
								#print tabIndLeft[0]
							sizeU = abs(int(detPrimer[3].upper().replace("U",""))-\
								((int(detPrimer[2].lower().replace("bp",""))-size)\
								/int(detPrimer[1].lower().replace("bp",""))))                            #computation of the sizeU value
							#print size,sizeU
							resul.append([primer[0],tabIndLeft[0][ind],tabIndRight[ind][ind2],size,sizeU])   

				if len(resul) == 0 and detPrimer[0] not in tabRes :    #if no result
					tabRes[detPrimer[0]]=[primer[0],"","","",""]   

				if (len(resul)>0) :
					ind=0
					for j in xrange(len(resul)) :                  #keep the result with the minimum sizeU value
						if (resul[ind][4]>resul[j][4]) : ind=j

					if round !="" and round>0 :                    #round of the sizeU value
						sizeU=resul[ind][4]
						if sizeU>=math.floor(sizeU) and sizeU<(math.floor(sizeU)+round) :
							sizeU = math.floor(sizeU)
							print "floor"
						elif sizeU <= math.ceil(sizeU) and sizeU>(math.ceil(sizeU)-round) :
							sizeU=math.ceil(sizeU)
							print "ceil"
						else :
							sizeU=math.floor(sizeU)+0.5
							print "0.5"
						resul[ind][4]=sizeU                    #set of the rounded sizeU value
					tabRes[detPrimer[0]]=resul[ind]                #set the best result as a new key : value in the dictionnary #replace the old dictionnary value if there is one
	return tabRes


#print inverComp("TATAGAGACACA")

#print mix(1,"TGT","ACTGACGACCAYACAAGTTAC")

#print findFirst ("GAGACGACGCTTGAGGTTTTT",fasta2) # -> ([2122861, 2122881, 2122892, 2122904], ['norm', 'norm', 'norm', 'norm'])

#print findSec ("GAC","ACTGACGACCAYACAAGTTAC",0,"norm") 
#print mix(1,"GCGGTGTTGTGTCTGTGGATA",fasta2)


#result = find(Primers,fasta,0.25)
print findFirst ("GCCGTCAGTATCCACGTCATAG",fasta2) # -> ([2122861, 2122881, 2122892, 2122904], ['norm', 'norm', 'norm', 'norm'])
#print findSec ("GCGTGACAATCGACTTTTTGTC",fasta2,684708,"norm") 
"""
Primers = Primers.split("\n")
Primers = [ pri.split("_")[0] for pri in Primers ]
del Primers[-1]

for Primer in Primers :
	print Primer,result[Primer]

res="\n"
for Primer in Primers :
	res+= "\t".join([Primer,str(result[Primer][4])])+"\n"

print res
"""
#different result from the site, surely due to the different mix function
#all different results are inferior -> better function mix ?

