#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, regex, math, sys, os.path

#fasta = "/home/david/Documents/complete_genomes/brucellaceae"
#primers = "/home/david/downloads/primers_brucella"

if len(sys.argv)>1 :
	fasta_path = sys.argv[1]
else :
	#fasta_path = input("Enter a fasta files directory : ")
	fasta_path = "/home/david/Documents/complete_genomes/brucellaceae"
files = os.listdir(fasta_path)

if len(sys.argv)>2 :
	primers_path = sys.argv[2]
else :
	#primers_path = input("Enter a primers file : ")
	primers_path = "/home/david/downloads/primers_brucella"
Primers = open(primers_path,"r").read()

#dictionnary to create complementary DNA sequences
dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}

def inverComp (seq) :					#return the inversed complementary sequence
	seq = seq.upper()                  		#acgt -> ACGT
	seq_comp = ""
	for nuc in seq :                   		#nuc stand for nucleotide
		seq_comp += dico_comp[nuc]  		#use of the dictionnary to have the complementary nucleotide
	seq_comp = "".join(reversed(seq_comp))		#we reverse the sequence	
	return (seq_comp)

def positionsOfMatches (result,seq) : 			#get the matches positions in the fasta sequence
	pos = []  
	for res in result :
		pos.append([seq.find(res),res])
	return (pos)

def search_matches(nbmismatch,primer, seq) : 		#return all match(es) in the fasta sequence 
	return (regex.findall("("+primer+"){e<="+str(nbmismatch)+"}",seq,overlapped=True))

def mismatch (nb,primer,seq,nb_mismatch) : 		#function to find match(s) with a mismatch
	found = set()					#use of set to delete duplicate position
	if nb==1 :
		invP=inverComp(primer)
		for i in range(len(primer)): 		#for each nucleotide of the primer
			reg1=primer[:i]
			reg2=primer[i+1:]
			reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2)
			searchseq = re.compile(reg)           	 	#search request for finditer (Primer with a mismatch)
			tmpfinditer = searchseq.finditer(seq)  		#tmpfinder : objects list with start() end() and group(0) functions
			tmp=[]
			for m in tmpfinditer :                 		#for each result
				tmp.append([str(m.start()),"norm"]) 		#add the result in the tmp list
			if tmp != [] :                         		#if results, we add the first result in listFind
				found.add("_".join(tmp[0]))	
			else :                                 		#same search with the complementary reversed primer
				reg1=invP[:i]
				reg2=invP[i+1:]
				reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2) #search request for one mismatch
				searchseq = re.compile(reg)
				tmpfinditer = searchseq.finditer(seq)
				for m in tmpfinditer :
					tmp.append([str(m.start()),"inv"])
				if tmp != [] :
					found.add("_".join(tmp[0]))	
			#if no match with one mismatch
		if len(found) == 0 and nb_mismatch == 2 :
			tmpfind_forward = search_matches(2,primer,seq) 		#get all matches (two mismatches)
			positions_f = positionsOfMatches(tmpfind_forward,seq)	#get positions of matches
			if tmpfind_forward == [] :
				tmpfind_reverse = search_matches(2,invP,seq)	#get the match(es) with two mismatch (use of regex.findall())
				positions_r = positionsOfMatches(tmpfind_reverse,seq)	#get the position(s) of match(es)
				if tmpfind_reverse != [] :
					for res in positions_r :
						found.add("_".join([str(res[0]),"inv"]))
			else :
				for res in positions_f :
					found.add("_".join([str(res[0]),"norm"]))
		
		
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
					tmp.append([str(m.start()),"norm"])
				if tmp != [] :
					found.add("_".join(tmp[0]))	
		#if no match with one mismatch, use of regex.findall for two mismatch
		if len(found) == 0 and nb_mismatch == 2 :
			tmpfind_forward = search_matches(2,primer,seq) 			#get all matches
			positions_f = positionsOfMatches(tmpfind_forward,seq)		#get positions of matches
			if tmpfind_forward != [] :
				for res in positions_f :
					found.add("_".join([str(res[0]),"norm"]))
	foundlist = []	
	for e in found : 
		foundlist.append([int(e.split("_")[0]),e.split("_")[1]])
	del found
	return(foundlist)

def findFirst (primer,seq) :          			#first search of the primer on the sequence (use inverComp() and mismatch())

	match = []
	sense = []                    			#to store the sense of search (normal or inversed) 
	result=seq.find(primer)
	while (result!=-1) :          			#while the search has not been made on the entire sequence               
		match.append(result)             
		position=result+1     			#next search will start one nucleotide after the position of the last match
		sense.append("norm")
		result=seq.find(primer,position)	#next search, if no more match : resul = -1 -> end of the loop

	if match == [] :                         	#if no perfect match found with the regular primer 
		primer_inv=inverComp(primer)		#get the inversed complementary primer with inverComp()
		result=seq.find(primer_inv)
		while(result!=-1):			#same search with the converted primer
			match.append(result)                 
			position=result+1   
			sense.append("inv")
			result=seq.find(primer_inv,position) 

	if match == [] :                    		#if still no matches (perfect match)
		lastTest=mismatch(1,primer,seq,2)       #search with a mismatch first then with 2 mismatches		
		if lastTest != [] :
			for last in lastTest :
				match.append(last[0])	#get the position of the match
				sense.append(last[1])	#get the sens of the match
	return (match,sense)

def findSec(primer,seq,sense) : 			#search a match for the second primer

	if sense == "norm" : primer=inverComp(primer)
	match = []
	result=seq.find(primer)
	while(result!=-1) :                      	#while there's a result
			match.append(result)          
			position=result+1		#indlook get the position of the following nucleotide for the next search
			result=seq.find(primer,position)

	if match == [] :				#if no perfect match
		lastTest=mismatch(2,primer,seq,2)	
		if lastTest != [] :
			for last in lastTest :
				match.append(last[0])			
	return match

def find2(primers,fasta,round) : #main function : return the result of the matches 
	
	fasta = fasta.replace(" ","").replace("\t","")		#delete spaces and tabulations

	if type(round) is str :
		round = round.replace(",",".") 
		round = float(round)                       
	
	sequences = fasta.split('>')				#split the fasta files into a list of fasta file
	del sequences[0]					#delete the '' created
	dico_res = {}

	for seq in sequences :
		tmp_seq = seq.split("\n")
		title_seq = tmp_seq[0]
		del tmp_seq[0]
		seq = "".join(tmp_seq).upper()

		if primers :					#if primers had been entered
			for primer in primers :			#for each couple of primers 
				primer = primer.replace(" ",";").replace("\t",";").split(";")
				primer_info = primer[0].split('_')
				first_match = findFirst(primer[1],seq)		#search match(es) for the first primer 
				second_match = []
				result = []
				for i,pos_match in enumerate(first_match[0]) :
					second_match.extend(findSec(primer[2],seq,first_match[1][i]))	#search match(es) for the second primer
					if second_match != [] :						#if there is a match with the second primer on the complementary DNA sequence
						for pos_match2 in second_match :			#for each match found for the second primer
								if first_match[1][i] == "inv" :
									size = abs(int(pos_match)-int(pos_match2)+len(primer[1])) 	
									size2 = int(pos_match)+len(primer[1])+(len(seq)-pos_match2)	#if primers are separated by the splitted area in the sequence (circular chromosome)
									if size2 < size : size = size2
								else :
									size = abs(int(pos_match2)-int(pos_match)+len(primer[2])) 
									size2 =int(pos_match2)+len(primer[2])+(len(seq)-pos_match)
									if size2 < size : size = size2
								
								sizeU = abs(float(primer_info[3].upper().replace("U",""))-\
									((float(primer_info[2].lower().replace("bp",""))-size)\
									/float(primer_info[1].lower().replace("bp",""))))		#computation of sizeU
								result.append([primer[0],pos_match,pos_match2,size,sizeU])
								
				if len(result) == 0 and primer_info[0] not in dico_res.keys() :	#if no result
					dico_res[primer_info[0]]=[primer[0],"","","",""]	

				elif len(result) > 0 :						#if result(s)
					best_res = result[0]
					for res in result :					#keep the result with the minimum sizeU value
						if res[4]<best_res[4] : best_res=res	

					if round !="" and round>0 :                    		#round of the sizeU value
						sizeU=best_res[4]
						if sizeU>=math.floor(sizeU) and sizeU<(math.floor(sizeU)+round) :
							sizeU = math.floor(sizeU)
						elif sizeU <= math.ceil(sizeU) and sizeU>(math.ceil(sizeU)-round) :
							sizeU=math.ceil(sizeU)
						else :
							sizeU=math.floor(sizeU)+0.5
						best_res[4]=sizeU				#set of the rounded sizeU value
					dico_res[primer_info[0]]=best_res			#set the best result as a new key : value in the dictionnary #replace the old dictionnary value if there is one
	return dico_res

Primers = Primers.split("\n")
if Primers[-1]=="" :
	del Primers[-1]
Primers_short = [ pri.split("_")[0] for pri in Primers ]

def main() : #run find2() for each genome file in the directory with all primers in the primers file
	for i,file in enumerate(files) :
		print file
		pathfile = fasta_path+"/"+file
		fasta = open(pathfile,"r").read()
		fasta_names = []
		for line in fasta.split("\n") :
			if ">" in line :
				fasta_names.append(line.replace("\n","").split("|")) 

		result = find2(Primers,fasta,0.25) 
					
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

#main()



with open("/home/david/Documents/complete_genomes/brucellaceae/Brucella_ceti_TE10759-12_NC_022905.1_NC_022906.1.fa","r") as myfile :  #load of the fasta file
	raw_fasta=myfile.read()
	myfile.seek(1)	
	tmp = myfile.read().split(">")
	chr1 = tmp[0].replace("\n","")
	chr2 = tmp[1].replace("\n","")
	fasta = chr1+chr2

#Bruce11-211_63bp_257bp_2u	CTGTTGATCTGACCTTGCAACC	CCAGACAACAACCTACGTCCTG
#Bruce42-424_125bp_539bp_4u	CATCGCCTCAACTATACCGTCA	ACCGCAAAATTTACGCATCG
#Bruce22-322_8bp_158bp_6u	GATGAAGACGGCTATCGACTG	TAGGGGAGTATGTTTTGGTTGC
#Bruce02-1923_339bp_787bp_3u	AACGCAGCATCACCAATGT	CCCAGATGTTCGGCTATAGTATG
#Bruce21-329_8bp_148bp_6u	CTCATGCGCAACCAAAACA	GATCTCGTGGTCGATAATCTCATT

#print mismatch(1,"CTGTTGATCTGACCTTGCAACC",chr1,2)
#print findFirst("CTCATGCGCAACCAAAACA",chr2)
#print findSec ("GATCTCGTGGTCGATAATCTCATT",chr2,"inv")
dico = find2(Primers,raw_fasta,0.25)

for key in dico.keys() :
	print key,dico[key]




