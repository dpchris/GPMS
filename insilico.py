#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, regex, math, sys, os.path

#fasta = "/home/david/Documents/complete_genomes/brucellaceae/Brucella_abortus_104M_NZ_CP009625.1_NZ_CP009626.1.fa"
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

def positionsOfMatches (result,seq) : #get the matches positions in the fasta sequence
	pos = []  
	for res in result :
		pos.append([seq.find(res),res])
	return (pos)

def search_matches(nbmismatch,primer, seq) : #return all match(es) in the fasta sequence 
	return (regex.findall("("+primer+"){e<="+str(nbmismatch)+"}",seq,overlapped=False))

def mix (nb,primer,seq,nb_mismatch) : #function to find match(s) with a mismatch
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
				reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2) #search request for one mismatch
				searchseq = re.compile(reg)
				tmpfinditer = searchseq.finditer(seq)
				for m in tmpfinditer :
					tmp.append([m.start(),"inv"])
				if tmp != [] :
					listFind.append(tmp[0])	
			#if no match with one mismatch
		if listFind == [] and nb_mismatch == 2 :
			tmpfind_forward = search_matches(2,primer,seq) 		#get all matches (two mismatches)
			positions_f = positionsOfMatches(tmpfind_forward,seq)	#get positions of matches
			if tmpfind_forward == [] :
				tmpfind_reverse = search_matches(2,invP,seq)	#get the match(es) with two mismatch (use of regex.findall())
				positions_r = positionsOfMatches(tmpfind_reverse,seq)	#get the position(s) of match(es)
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
		#if no match with one mismatch, use of regex.findall for two mismatch
		if listFind == [] and nb_mismatch == 2 :
			tmpfind_forward = search_matches(2,primer,seq) 			#get all matches
			positions_f = positionsOfMatches(tmpfind_forward,seq)		#get positions of matches
			if tmpfind_forward != [] :
				for res in positions_f :
					listFind.append([res[0],"norm"])
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
			lastTest=mix(1,primer,seq,2)       #search with a mismatch with mix()		
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
			lastTest=mix(2,primer,seq,2)             #mix with arg 2 : search only with the given primer, not the complementary
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
			lastTest=mix(2,primer,seq,2) 
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
	tabRes = {}

	for f in xrange(len(files)) :                  #for each fasta sequence 
		text=">"+files[f]                      
		tmpSeq=text.split('\n')                
		titleSeq=tmpSeq[0]                     #get the name of the fasta file
		del tmpSeq[0]		               #del the header of the fasta file to keep only the sequence
		text=''.join(tmpSeq).upper()           #merge the sequence list into a single string

		if primers :                           #if primers had been entered
			for i in range(len(primers)) :
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
								size2 = int(tabIndLeft[0][ind])+len(primer[1])+(len(text)-tabIndRight[ind][ind2]) #if primers are separated by the splitted area in the sequence (circular chromosome)
								if size2 < size : size = size2 							  #if the pattern repetition is between the end and the beginning of the sequence 
								#print int(tabIndLeft[0][ind]),int(tabIndRight[ind][ind2]),size, size2, len(text)
							else :
								size = abs(int(tabIndRight[ind][ind2])-int(tabIndLeft[0][ind])+len(primer[2])) 
								size2 = int(tabIndRight[ind][ind2])+len(primer[2])+(len(text)-tabIndLeft[0][ind])
								#print int(tabIndLeft[0][ind]),int(tabIndRight[ind][ind2]),size, size2, len(text)								
								if size2 < size : size = size2
							sizeU = abs(float(detPrimer[3].upper().replace("U",""))-\
								((float(detPrimer[2].lower().replace("bp",""))-size)\
								/float(detPrimer[1].lower().replace("bp",""))))                            #computation of the sizeU value
							resul.append([primer[0],tabIndLeft[0][ind],tabIndRight[ind][ind2],size,sizeU])  

				if len(resul) == 0 and detPrimer[0] not in tabRes :    #if no result
					tabRes[detPrimer[0]]=[primer[0],"","","",""]   

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
					tabRes[detPrimer[0]]=resul[ind]                #set the best result as a new key : value in the dictionnary #replace the old dictionnary value if there is one
	return tabRes

Primers = Primers.split("\n")
if Primers[-1]=="" :
	del Primers[-1]
Primers_short = [ pri.split("_")[0] for pri in Primers ]

def main() : #run find() for each genome file in the directory with all primers in the primers file
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
			pathfile = "/home/david/Documents/MLVA/mlva_results/MLVA_analysis_"+fasta_path.split("/")[-1]+".csv"	
			if os.path.exists(pathfile) :
				pathfile = pathfile.split(".")[0]+"_bis.csv"
			output = open(pathfile,"w") 				#output is a csv file (delimiter=";")
			output.write(";".join(["fasta_chr1","fasta_chr2","gi_chr1","gi_chr2","ref_chr1","ref_chr2"]+locus)+"\n")  #header
			output = open(pathfile,"a")
		output.write(";".join([fasta_names[0][4][1:],fasta_names[1][4][1:],fasta_names[0][1],fasta_names[1][1],fasta_names[0][3],fasta_names[1][3]]+mlva_score)+"\n")
	output.close()
	print "MLVA analysis finished for "+fasta_path.split("/")[-1]

main()


#Ochrobactrum_anthropi_ATCC_49188_NC_009667.1_NC_009668.1.fa
#Bruce51-80_15bp_74bp_2u	TGACATGATGCAGAAAATGCA	GTCCCTTGCCGCCTTTCA
#Brucella_abortus_NZ_CP007700.1_NZ_CP007701.1.fa
#Bruce02-1923_339bp_787bp_3u	AACGCAGCATCACCAATGT	CCCAGATGTTCGGCTATAGTATG
"""
with open("/home/david/Documents/complete_genomes/brucellaceae/Brucella_abortus_NZ_CP007700.1_NZ_CP007701.1.fa","r") as myfile :  #load of the fasta file
	raw_fasta=myfile.read()
	myfile.seek(1)	
	tmp = myfile.read().split(">")
	chr1 = tmp[0].replace("\n","")
	chr2 = tmp[1].replace("\n","")
	fasta = chr1+chr2

#print mix(1,"TGACATGATGCAGAAAATGCA",chr1,2)
#print mix(1,"GTCCCTTGCCGCCTTTCA",chr1,2)
#print findFirst("TGACATGATGCAGAAAATGCA",chr1)
#print findSec("GTCCCTTGCCGCCTTTCA",chr1,"1031268","inv")
print find(["Bruce02-1923_339bp_787bp_3u	AACGCAGCATCACCAATGT	CCCAGATGTTCGGCTATAGTATG"],raw_fasta,0.25)
"""

