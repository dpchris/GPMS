#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re

#fasta test (from insilico site)

with open("/home/david/Documents/fasta_example.txt","r") as myfile :  #load of the fasta file
	myfile.readline()                                             #we don't keep the header
	fasta = myfile.read().replace('\n', '')                       #erase \n to have a single string

primer="vrrB1_9bp_229bp_16U ATAGGTGGTTTTCCGCAAGTTATTC GATGAGTTTGATAAAGAATAGCCTGTG"

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
		print tmpfind
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

def find(text,primers,round) : #main function : return the result of the matches 
	errMult=""
	#var text=trim(GetId('textInput').value)                   # on recupere le fasta
	#var primers=GetId('primersInput').value.split('\n')       # on recupere les primers 
	#var primersSNP=GetId('primersInputSNP').value.split('\n') #on recupere les primers SNP
	
	text = text.strip()
	
	round=round.replace(',','.')
	round=float(round)    # parsefloat : transforme un string en float
	
	#primers=primers.filter(function(e){return e}) #??
	#primersSNP=primersSNP.filter(function(e){return e}) #??
	#var sizeSNP=parseInt(GetId('sizeSnp').value)
	
	files = text.split('>')              #split par > pour séparer les fasta
	tabRes = []
	for file in files :                  # pour chaque fasta 
		text=">"+file
		tmpSeq=text.split('\n') 
		titleSeq=tmpSeq[0]
		titleSeq2=tmpSeq[0]
		del tmpSeq[0]

		text=''.join(tmpSeq).upper()
		if primers :                 # si les primers on été rentré
			for primer in primers :
				primer = primer.replace(" ",";").replace("\t",";").split(";")
				detPrimer = primer[0].split('_')
				
				tabIndLeft = findFirst(primer[1],text,primer[0]) #recherche de match avec findFirst primer[1]=primer, text=fasta, primer[0]= nom primer
				tabIndRight = []
				
				resul = []
				for (var ind=0;ind<tabIndLeft[0].length;ind++) : 
					tabIndRight.push(findSec(primer[2],text,primer[0],tabIndLeft[0][ind],tabIndLeft[1][ind])) 
					if (tabIndRight!=""&&tabIndRight.length>0) :  # si resultat avec le deuxieme primer
						for (var ind2=0;ind2<tabIndRight[0].length;ind2++){
							// ## calcule de la sortie xxxbp - xxx recurence - xU
							if (tabIndLeft[1][ind]=="inv"){
								var size=Math.abs(tabIndLeft[0][ind]-tabIndRight[ind][ind2]+primer[1].length)
								var sizeU=Math.abs(detPrimer[3].toUpperCase().replace('U','')-((detPrimer[2].toLowerCase().replace('bp','')-size)/detPrimer[1].toLowerCase().replace('bp','')))
							}
							else{
								var size=Math.abs(tabIndRight[ind][ind2]-tabIndLeft[0][ind]+primer[2].length)
								var sizeU=Math.abs(detPrimer[3].toUpperCase().replace('U','')-((detPrimer[2].toLowerCase().replace('bp','')-size)/detPrimer[1].toLowerCase().replace('bp','')))
							}
							if (GetId('focusFull').checked)
								resul.push([detPrimer[0],tabIndLeft[0][ind],tabIndRight[ind][ind2],size,sizeU])
							else
								resul.push([primer[0],tabIndLeft[0][ind],tabIndRight[ind][ind2],size,sizeU])
						}
					}
				}
				
				if (resul.length==0&&!(detPrimer[0] in tabRes)) :
					if (GetId('focusFull').checked) :
						tabRes[detPrimer[0]]=[detPrimer[0],"","","",""];
					else
						tabRes[primer[0]]=[detPrimer[0],"","","",""];
				}
				#calcule de la valeur U (de son arrondi)
				if (resul.length>0) :
					var ind=0
					for (var j=1;j<resul.length;j++) :
						if (resul[ind][4]>resul[j][4])ind=j :
					}
					
					if (!isNaN(round)&&round>0){
						sizeU=resul[ind][4] :
						if (sizeU>=Math.floor(sizeU)&&sizeU<(Math.floor(sizeU)+round)) #math.floor : plus grand nombre entier inferieur
							sizeU=Math.floor(sizeU)
						else if(sizeU<=Math.ceil(sizeU)&&sizeU>(Math.ceil(sizeU)-round)) #math.ceil : arrondi a l'entier supperieur
							sizeU=Math.ceil(sizeU)
						else
							sizeU=Math.floor(sizeU)+0.5
						resul[ind][4]=sizeU
					}
					tabRes[detPrimer[0]]=resul[ind]
					delete primers[i]
				}
			}
		}



#print inverComp("TATAGAGACACA")

#print mix(1,"TGT","ACTGACGACCAYACAAGTTAC")

#print findFirst ("GAC","ACTGACGACCAYACAAGTTAC")

#print findSec ("GAC","ACTGACGACCAYACAAGTTAC",0,"norm")
#print mix(2,"GAC","ACTGACGACCAACAAGTTAC")






