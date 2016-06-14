#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, math, sys, os.path

#fasta = "/home/david/Documents/complete_genomes/brucellaceae"
#primers = "/home/david/downloads/primers_brucella"

if len(sys.argv)>1 :
	fasta_path = sys.argv[1]
else :
	fasta_path = input("Enter a fasta files directory : ")
files = os.listdir(fasta_path)

if len(sys.argv)>2 :
	primers_path = sys.argv[2]
else :
	primers_path = input("Enter a primers file : ")
Primers = open(primers_path,"r").read()

#dictionnary to create complementary DNA sequences
dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}

def inverComp (seq) :                              #return the inversed complementary sequence
	seq = seq.upper()                          #acgt -> ACGT
	seq_comp = ""
	for nuc in seq :                           #nuc stand for nucleotide
		seq_comp += dico_comp[nuc]         #use of the dictionnary to have the complementary nucleotide
	seq_comp = "".join(reversed(seq_comp))     #we reverse the sequence	
	return (seq_comp)

def mix (nb,primer,seq) : #function to find match(s) with a mismatch
	listFind =[]
	if nb==1 :
		invP=inverComp(primer)
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
	return(listFind)

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

def findFirst (primer,seq) :                         #first search of the primer on the sequence (use inverComp() and mix())
	indLeft = []
	sense = []                                   #to store the sense of search (normal or inversed) 
	indlook=0
	resul=seq.find(primer,indlook)               #get the first match (if it exists)
	while (resul!=-1) :                          #while the search has not been made on the entire sequence               
		indLeft.append(resul)             
		indlook=indLeft[-1]+1                #next search will start one nucleotide after the position of the last match
		sense.append("norm")
		resul=seq.find(primer,indlook)       #next search, if no more match : resul = -1 -> end of the loop
	
	if len(indLeft)==0 :                         #if no match foud with the regular primer 
		indlook=0
		tmpP=inverComp(primer)               #get the inversed complementary primer with inverComp()
		resul=seq.find(tmpP,indlook)
		while(resul!=-1):                    #same search with the converted primer
			indLeft.append(resul)                 
			indlook=indLeft[-1]+1   
			sense.append("inv")
			resul=seq.find(tmpP,indlook) 
		
		if (len(indLeft)==0):                    #if still no matches
			lastTest=mix(1,primer,seq)       #search with a mismatch with mix()		
			if lastTest != [] :
				for last in lastTest :
					indLeft.append(last[0])        #get the position of the match
					sense.append(last[1])          #get the sens of the match
			
		
	
	return (indLeft,sense)

def findSec(primer,seq,indLeft,sense) : #search a match for the second primer
	indRight = []
	if sense=="norm" :
		primer=inverComp(primer)                #the reversed complementary primer is used
		indlook=indLeft                         #start the search at the match position of the first primer
		resul=seq.find(primer,indlook)
		while(resul!=-1) :                      #while there's a result
			indRight.append(resul)          
			indlook=indRight[-1]+1		#indlook get the position of the following nucleotide for the next search
			resul=seq.find(primer,indlook)
		
		if len(indRight) ==0 :	                       #if no match found : try with one mismatch
			lastTest=mix(2,primer,seq)             #mix with arg 2 : search only with the given primer, not the complementary
			if lastTest != [] :
				indRight.append(lastTest[0][0])   #indRight get the position of the first match found by mix()
	else :						       #if sense=="inv" 
		indlook=0
		resul=seq.find(primer,indlook)
		while(resul!=-1) :
			indRight.append(resul)
			indlook=indRight[-1]+1
			resul=seq.find(primer,indlook)
		
		if len(indRight)==0 :
			lastTest=mix(2,primer,seq) 
			if lastTest != [] :
				for i in range(len(lastTest)) :
					indRight.append(lastTest[0][0])
	return indRight

def find(primers,text,round) : #main function : return the result of the matches 
	
	text = text.replace(" ","").replace("\t","")   #delete spaces and tabulations
	#primers = primers.split("\n") 
	#if primers[-1] =="" :
	#	del primers[-1]               

	if type(round) is str :
		round = round.replace(",",".") 
		round = float(round)                       
	
	files = text.split('>')                        #split the fasta files into a list of fasta file
	del files[0]                                   #delete the '' created
	dico_Res = {}

	for f in xrange(len(files)) :                  #for each fasta sequence 
		text=">"+files[f]                      
		tmpSeq=text.split('\n')                
		titleSeq=tmpSeq[0]                     #get the name of the fasta file
		del tmpSeq[0]		               #del the header of the fasta file to keep only the sequence
		text=''.join(tmpSeq).upper()           #merge the sequence list into a single string

		if primers :                           #if primers had been entered
			for i in range(len(primers)) :
				#print primers[i]
				primer = primers[i].replace(" ",";").replace("\t",";").split(";")        #split the primers name, and primers into a list
				detPrimer = primer[0].split('_')                                         #split differents elements of the name into a list
				tabIndLeft = findFirst(primer[1],text)                                   #search match with findFirst() : primer[1]=primer, text=fasta
				tabIndRight = []
				resul = []
				for ind in range(len(tabIndLeft[0])) :                                                      #parse result(s) of findfirst to do the second search with findsec()
					tabIndRight.append(findSec(primer[2],text,tabIndLeft[0][ind],tabIndLeft[1][ind]))   #get result(s) of findsec with the second primer
					if tabIndRight != [] :
						for ind2 in range(len(tabIndRight[0])) :                                  
							if tabIndLeft[1][ind]=="inv" :                                      #if the first match is made with the reversed complementary first primer
								size = abs(int(tabIndLeft[0][ind])-int(tabIndRight[ind][ind2])+len(primer[1])) 
							else :
								size = abs(int(tabIndRight[ind][ind2])-int(tabIndLeft[0][ind])+len(primer[2])) 
							sizeU = abs(float(detPrimer[3].upper().replace("U",""))-\
								((float(detPrimer[2].lower().replace("bp",""))-size)\
								/float(detPrimer[1].lower().replace("bp",""))))                            #computation of the sizeU value
							resul.append([primer[0],tabIndLeft[0][ind],tabIndRight[ind][ind2],size,sizeU])   

				if len(resul) == 0 and detPrimer[0] not in dico_Res :    #if no result
					dico_Res[detPrimer[0]]=[primer[0],"","","",""]   

				if (len(resul)>0) :                                    #if result(s)
					ind=0
					for j in range(len(resul)) :                   #keep the result with the minimum sizeU value
						if (resul[ind][4]>resul[j][4]) : ind=j

					if round !="" and round>0 :                    #round of the sizeU value
						sizeU=resul[ind][4]
						if sizeU>=math.floor(sizeU) and sizeU<(math.floor(sizeU)+round) :
							sizeU = math.floor(sizeU)
						elif sizeU <= math.ceil(sizeU) and sizeU>(math.ceil(sizeU)-round) :
							sizeU=math.ceil(sizeU)
						else :
							sizeU=math.floor(sizeU)+0.5
						resul[ind][4]=sizeU                    #set of the rounded sizeU value
					dico_Res[detPrimer[0]]=resul[ind]                #set the best result as a new key : value in the dictionnary #replace the old dictionnary value if there is one
	return dico_Res

Primers = Primers.split("\n")
if Primers[-1]=="" :
	del Primers[-1]
Primers_short = [ pri.split("_")[0] for pri in Primers ]

def main() :
	for i,file in enumerate(files) :
		print file
		pathfile = fasta_path+"/"+file
		fasta = open(pathfile,"r").read()
		fasta_names = []
		for line in fasta.split("\n") :
			if ">" in line :
				fasta_names.append(line.replace("\n","").split("|")) 

		result = find(Primers,fasta,0.25) 
					
		locus=[]
		mlva_score=[]
		for Primer in Primers_short :
			locus.append(Primer.split("_")[0])
			mlva_score.append(str(result[Primer][4]))
		if i==0 : 
			pathfile = "/home/david/Documents/mlva_results/MLVA_analysis_"+fasta_path.split("/")[-1]+".csv"	
			output = open(pathfile,"w") 
			output.write(";".join(["fasta_chr1","fasta_chr2","gi_chr1","gi_chr2","ref_chr1","ref_chr2"]+locus)+"\n")  #header
			output = open(pathfile,"a")
		output.write(";".join([fasta_names[0][4][1:],fasta_names[1][4][1:],fasta_names[0][1],fasta_names[1][1],fasta_names[0][3],fasta_names[1][3]]+mlva_score)+"\n")
	output.close()
	print "MLVA analysis finished for "+fasta_path.split("/")[-1]

main()

