#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, regex, math, sys, os.path, pickle

#fasta = "/home/david/Documents/spades_output/insilico_input"
#primers = "/home/david/downloads/primers_brucella"

if len(sys.argv)>1 :
	fasta_path = sys.argv[1]
else :
	#fasta_path = input("Enter a fasta files directory : ")
	fasta_path = "/home/david/Documents/complete_genomes/legionella_pneumophila"
files = os.listdir(fasta_path)

if len(sys.argv)>2 :
	primers_path = sys.argv[2]
else :
	#primers_path = input("Enter a primers file : ")
	primers_path = "/home/david/Documents/MLVA/primer_Legionella_pneumophila.csv"
Primers = open(primers_path,"r").read()

if len(sys.argv)>3 :
	nb_mismatch = int(sys.argv[3])
else :
	nb_mismatch = 2

log = open("/home/david/Documents/MLVA/mlva_results/"+fasta_path.split("/")[-1]+"_output.txt","w")
log = open("/home/david/Documents/MLVA/mlva_results/"+fasta_path.split("/")[-1]+"_output.txt","a")

#dictionnary to create complementary DNA sequences
dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}
dico_ref = pickle.load(open("/home/david/Documents/MLVA/dico_table_ref","r"))


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

def pretty_mismatch (primer, found) :

	if len(found) == len(primer) :
		tmp=[]
		for i,nuc in enumerate(primer) :
			if found[i] != nuc :
				tmp += found[i].lower()
			else : 
				tmp += found[i]
		found = "".join(tmp)
	else :
		diff = primer.find(found)		
		if diff == 0 :
			found = found+(len(primer)-len(found))*"."
		else :
			found = (len(primer)-len(found))*"."+found
	return(found)

def clean_mismatches (primer1,primer2,sense,found1,found2) :
	if sense == "norm" :
		primer2 = inverComp(primer2)
	else :
		primer1 = inverComp(primer1)

	res_found1 = []
	res_found2 = []
	for found in found1 :
		if len(found)<=len(primer1) :
			res_found1.append(pretty_mismatch(primer1,found))
	for found in found2 :
		if len(found)<=len(primer2) and found != primer2 :
			res_found2.append(pretty_mismatch(primer2,found))

	return (res_found1,res_found2)

def mismatch (nb,primer,seq) : 				#function to find match(s) with a mismatch
	match = []
	sense = []
	tmp_res = set()
	mismatches = []
	if nb==1 :
		invP=inverComp(primer)
		for i in range(len(primer)): 		#for each nucleotide of the primer
			reg1=primer[:i]
			reg2=primer[i+1:]
			reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2)
			searchseq = re.compile(reg)           	 	#search request for finditer (Primer with a mismatch)
			tmpfinditer = searchseq.finditer(seq)  		#tmpfinder : objects list with start() end() and group(0) functions
			tmp= set()
			for m in tmpfinditer :                 		#for each result
				#dico_mismatches[primer] = m.group()
				tmp.add("_".join([str(m.start()),"norm",m.group()]))	 		#add the result in the tmp list
			if tmp :                         		#if results, we add the first result in listFind
				for res in tmp :
					tmp_res.add(res)	
			else :                                 		#same search with the complementary reversed primer
				reg1=invP[:i]
				reg2=invP[i+1:]
				reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2) #search request for one mismatch
				searchseq = re.compile(reg)
				tmpfinditer = searchseq.finditer(seq)
				for m in tmpfinditer :
					#dico_mismatches[primer] = m.group()
					tmp.add("_".join([str(m.start()),"inv",m.group()]))
				if tmp :
					for res in tmp :
						tmp_res.add(res)
	elif nb==2 :					      		#search only in the forward sense
		for i in range(len(primer)):
			reg1=primer[:i]
			reg2=primer[i+1:]
			reg = str(reg1) + "["+"".join(dico_comp.keys()).replace("|","").replace(".","")+"]" + str(reg2)
			searchseq = re.compile(reg)
			tmpfinditer = searchseq.finditer(seq) 
			if tmpfinditer !=[] :
				tmp= set()
				for m in tmpfinditer :
					#dico_mismatches[primer] = m.group()
					tmp.add("_".join([str(m.start()),"norm",m.group()]))
				
				if tmp :
					for res in tmp :
						tmp_res.add(res)
	for e in tmp_res :
		match.append(int(e.split("_")[0]))
		sense.append(e.split("_")[1]) 	
		mismatches.append(e.split("_")[2])

	return(match,sense,mismatches)

def mismatches (nb,primer,seq,nbmismatch) :						#function to find match(es) with at least two mismatches 
	match = []
	sense = []
	invP=inverComp(primer)	
	dico_mismatches = {}
	if nb==1 :
		tmpfind = search_matches(nbmismatch,primer,seq) 		#get all matches (with mismatches)
		positions_f = positionsOfMatches(tmpfind,seq)			#get positions of matches
		if tmpfind == [] :
			tmpfind = search_matches(nbmismatch,invP,seq)		#get the match(es) with mismatches (use of regex.findall())
			positions_r = positionsOfMatches(tmpfind,seq)		#get the position(s) of match(es)
			if tmpfind != [] :
				for res in positions_r :
					match.append(res[0])
					sense.append("inv")
		else :
			for res in positions_f :
				match.append(res[0])
				sense.append("norm")

	elif nb==2 :
		tmpfind = search_matches(2,primer,seq) 			#get all matches
		positions_f = positionsOfMatches(tmpfind,seq)		#get positions of matches
		if tmpfind != [] :
			for res in positions_f :
				match.append(res[0])
				sense.append("norm")
	return (match,sense,tmpfind)							#return results (position of match + sense of primer)


def findFirst (primer,seq,nbmismatch) :          			#first search of the primer on the sequence (use inverComp() and mismatch())

	match = []
	sense = []                    					#to store the sense of search (normal or inversed) 
	mismatchs = []
	if nbmismatch == 0 :
		result=seq.find(primer)
		while (result!=-1) :          				#while the search has not been made on the entire sequence               
			match.append(result)             
			position=result+1     				#next search will start one nucleotide after the position of the last match
			sense.append("norm")
			result=seq.find(primer,position)		#next search, if no more match : resul = -1 -> end of the loop

		if match == [] :                         		#if no perfect match found with the regular primer 
			primer_inv=inverComp(primer)			#get the inversed complementary primer with inverComp()
			result=seq.find(primer_inv)
			while(result!=-1):				#same search with the converted primer
				match.append(result)                 
				position=result+1   
				sense.append("inv")
				result=seq.find(primer_inv,position) 

	if nbmismatch == 1 and match == [] :                    	#if still no matches (perfect match)
		match,sense,mismatchs = mismatch(1,primer,seq)       		#search with a mismatch 

	if nbmismatch >= 2 and match == [] :
		match,sense,mismatchs = mismatches(1,primer,seq,nbmismatch)				

	return (match,sense,mismatchs)

def findSec(primer,seq,sense,nbmismatch) : 				#search a match for the second primer
	match = []
	mismatchs = []
	if sense == "norm" : primer=inverComp(primer)
	if nbmismatch == 0 :
		result=seq.find(primer)
		while(result!=-1) :                      		#while there's a result
			match.append(result)          
			position=result+1				#indlook get the position of the following nucleotide for the next search
			result=seq.find(primer,position)

	if nbmismatch == 1 and match == [] :            		#if no perfect match
		match,trash,mismatchs = mismatch(2,primer,seq)	

	if nbmismatch >= 2 and match == [] :
		match,trash,mismatchs = mismatches(2,primer,seq,nbmismatch)
		
	return match,mismatchs

def table_ref (primer,size) : 						#return the sizeU value if the size is indexed in the table for brucella
	if primer in dico_ref.keys() :
		values = dico_ref[primer]
		sizes = []
		U = []
		for val in values :
			if "-" in val.split(" ")[0] :
				sizes.append(val.split(" ")[0].split("-")[0])
				sizes.append(val.split(" ")[0].split("-")[1])
				U.append(val.split(" ")[1].replace('(','').replace(")",""))
				U.append(val.split(" ")[1].replace('(','').replace(")",""))
			else :
				sizes.append(val.split(" ")[0]) 
		 		U.append(val.split(" ")[1].replace('(','').replace(")","")) 

	if str(size) in sizes :
		ind = sizes.index(str(size)) 
		sizeU = U[ind]
		return sizeU

def find2(primers,fasta,round,nbmismatch) : 			#return the result of the matches 
	
	fasta = fasta.replace(" ","").replace("\t","")		#delete spaces and tabulations
	mismatchs2 = []

	if type(round) is str :
		round = round.replace(",",".") 
		round = float(round)                       
	
	sequences = fasta.split('>')				#split the fasta files into a list of fasta file
	del sequences[0]					#delete the '' created
	dico_res = {}

	for s,seq in enumerate(sequences) :			#for each contig in the fasta file
		tmp_seq = seq.split("\n")
		title_seq = tmp_seq[0]
		del tmp_seq[0]
		seq = "".join(tmp_seq).upper()

		if primers :					#if primers had been entered
			for p,primer in enumerate(primers) :								#for each couple of primers 
				primer = primer.replace(" ",";").replace("\t",";").split(";")
				primer_info = primer[0].split('_')
				tmp1, tmp2, mismatchs = findFirst(primer[1],seq,nbmismatch)			#search match(es) for the first primer	
				first_match = tmp1,tmp2	
				second_match = []
				result = []
				for i,pos_match in enumerate(first_match[0]) :						#for each match of the first primer
					tmp,mismatchs2 = findSec(primer[2],seq,first_match[1][i],nbmismatch)	#search match(es) for the second primer
					second_match.extend(tmp)	
					if second_match != [] :								#if there is a match with the second primer on the complementary DNA sequence
						for pos_match2 in second_match :					#for each match found for the second primer
								if first_match[1][i] == "inv" :
									size = abs(int(pos_match)-int(pos_match2)+len(primer[1])) 	
								else :
									size = abs(int(pos_match2)-int(pos_match)+len(primer[2])) 

								sizeU = abs(float(primer_info[3].upper().replace("U",""))-\
									((float(primer_info[2].lower().replace("bp",""))-size)\
									/float(primer_info[1].lower().replace("bp",""))))		#computation of sizeU

								mismatchs, mismatchs2 = clean_mismatches(primer[1],primer[2],first_match[1][i],mismatchs,mismatchs2)
								result.append([primer[0],pos_match,pos_match2,size,sizeU,"contig"+str(s+1),nbmismatch,primer[1],mismatchs,primer[2],mismatchs2])
								
				if len(result) == 0 and primer_info[0] not in dico_res.keys() :	#if no result
					dico_res[primer_info[0]]=[["\t".join([primer[0],primer[1],primer[2]])],"","","","","contig"+str(s+1),nbmismatch,primer[1],mismatchs,primer[2],mismatchs2]	

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
					if primer_info[0] in dico_res.keys() and dico_res[primer_info[0]][4] != "" :
						best_res[5] = best_res[5]+", chr"+str(s+1) 	#if there's already a result with perfect matches
					dico_res[primer_info[0]]=best_res			#set the best result as a new key : value in the dictionnary #replace the old dictionnary value if there is one
	return dico_res

def get_empty_locus (dico_result) :
	tmpprimers = []
	for locus in dico_result.keys() :
		if dico_result[locus][4] == '' :
			tmpprimers.extend(dico_result[locus][0])
	return tmpprimers

def run (Primers,fasta,round,nbmismatch) :
	tmp = len(Primers)
	tmpPrimers = Primers
	result = {}
	for mismatch_allowed in range(int(nb_mismatch)+1) :
		tmp_dico = find2(tmpPrimers,fasta,round,mismatch_allowed) 			#no mismacth
		result = dict(result.items() + tmp_dico.items())
		tmpPrimers = get_empty_locus(result)
		nb_match = tmp -len(tmpPrimers) 
		print "results with",mismatch_allowed,"mismatch : ",nb_match,"/",len(Primers)
		log.write("".join(["results with ",str(mismatch_allowed)," mismatch : ",str(nb_match),"/",str(len(Primers))])+"\n")

		tmp = len(tmpPrimers)

		if len(tmpPrimers) == 0 :
			break 

	if len(tmpPrimers) != 0 : 
		print "no match : ", tmp
		log.write("".join(["no match : ",str(tmp)])+"\n")

	return result

def main() : #run find2() for each genome file in the directory with all primers in the primers file
	for i,file in enumerate(files) :
		print file, "\t strain ",i+1,"/",len(files)
		log.write("".join([file, "\t strain ",str(i+1),"/",str(len(files))])+"\n")
		pathfile = fasta_path+"/"+file
		fasta = open(pathfile,"r").read()
		fasta_names = []
		for line in fasta.split("\n") :
			if ">" in line :
				fasta_names.append(line.replace("\n","").split("|")) 

		dict_mismatches = {}
		result = run(Primers,fasta,0.25,nb_mismatch)			#use find for each number of mismatch	
				
		locus = []
		mlva_score = []
		mismatch = []
		header_mismatch = []
		for Primer in Primers_short :
			if result[Primer][1]=="" :
				result[Primer][6] = 99
			if result[Primer][6] not in [0,99] :
				print Primer, result[Primer][7:]
				log.write("\t".join([Primer]+[result[Primer][7]]+result[Primer][8]+["\n"]+[result[Primer][9]]+result[Primer][10])+"\n")
			locus.append(Primer.split("_")[0])
			mlva_score.append(str(result[Primer][4]))		#scores
			header_mismatch.append("nbmis_"+Primer.split("_")[0])
			mismatch.append(str(result[Primer][6]))			#number of mismatch allowed for each locus

		if i==0 :
			pathfile = "/home/david/Documents/MLVA/mlva_results/MLVA_analysis_test_"+fasta_path.split("/")[-1]+".csv"	
			output = open(pathfile,"w") 				#output is a csv file (delimiter=";")
			output.write(";".join(["Access_number"]+locus+header_mismatch)+"\n")  			#header
			output = open(pathfile,"a")
		output.write(file+";"+";".join(mlva_score+mismatch) + "\n")
	output.close()
	print "MLVA analysis finished for "+fasta_path.split("/")[-1]
	log.write("MLVA analysis finished for "+fasta_path.split("/")[-1])
	log.close

Primers = Primers.split("\n")
if Primers[-1]=="" :
	del Primers[-1]
Primers_short = [ pri.split("_")[0] for pri in Primers ]

main()

with open("/home/david/Documents/complete_genomes/legionella_pneumophila/Legionella_pneumophila_subsp._pneumophila_LPE509_NC_020521.1.fa","r") as myfile :  #load of the fasta file
	raw_fasta=myfile.read()
	myfile.seek(1)	
	tmp = myfile.read().split(">")
	chr1 = tmp[0].replace("\n","")
	#chr2 = tmp[1].replace("\n","")
	#fasta = chr1+chr2

#Lp13_24bp_428bp_11U;CAATAGCATCGGACTGAGCA;TGCCTGTGTATCTGGAAAAGC

#print mismatches (1,"CAATAGCATCGGACTGAGCA",chr1,1)
#print findFirst ("CAATAGCATCGGACTGAGCA",chr1,2)
#print findSec ("CAATAGCATCGGACTGAGCA",chr1,"inv",2)
#print Primers[2]
#print find2([Primers[2]],raw_fasta,0.25,2)

#CTGAAACAGTTGAGGATGTGA', ['TCACATCCTCAACTGTTTCAG'], 'TTATCAACCTCATCATCCCTG', ['TTATCAACCTCATCATCCCTG']
#print pretty_mismatch ("TTATCAACCTCATCATCCCTG", 'ATCAACCTCATCATCCCA')
#print clean_mismatches ("CTGAAACAGTTGAGGATGTGA",'TTATCAACCTCATCATCCCTG',"inv",'ACATCCTCAACTGTTTCAG','TTATCAAACTCATCATCCCTG')

