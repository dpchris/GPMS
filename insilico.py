#!/usr/bin/env python
# -*- coding: utf-8 -*-
import re, regex, math, sys, os.path, pickle, getopt, csv

#dictionnary to create complementary DNA sequences
dico_comp = {'A':'T','C':'G',"G":"C","T":"A","M":"K","R":"Y","W":"W","S":"S","Y":"R","K":"M","V":"B","H":"D","D":"H","B":"V","X":"X","N":"X",".":".","|":"|"}
#dico_ref = pickle.load(open("/home/david/Documents/MLVA/dico_table_ref","r"))

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

def pretty_mismatch (primer, found) : 			#lower the mismatches nucleotides 
	if len(found) == len(primer) :
		tmp=""
		for i,nuc in enumerate(primer) :
			if found[i] != nuc :
				tmp += found[i].lower()
			else : 
				tmp += found[i]
		found = tmp
	else :
		if found in primer :  
			diff = primer.find(found)		
			if diff == 0 :
				found = found+(len(primer)-len(found))*"."
			else :
				found = (len(primer)-len(found))*"."+found
			found = pretty_mismatch(primer,found)
	return(found)

def clean_mismatches (nbprimer,primer,sense_list,found_list) : #return the mismatches with differences in lower character from findfirst and findsec result

	res_found = []
	for sense,found in zip(sense_list,found_list) :
 
		if sense == "norm" :
			if nbprimer==2 : found = inverComp(found)
		else :
			if nbprimer==1 : found = inverComp(found)

		if len(found) <= len(primer) : 
			if found == primer : break
			res_found.append(pretty_mismatch(primer,found))

		if len(primer) in [len(res) for res in res_found] :				#if there's a mismatch with the same length as the primer
			res_found = [res for res in res_found if len(res) == len(primer)] 	#only keep mismatchs with the same length

	return (res_found)



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
					tmp.add("_".join([str(m.start()),"norm",m.group()]))
				
				if tmp :
					for res in tmp :
						tmp_res.add(res)
	for e in tmp_res :
		match.append(int(e.split("_")[0]))
		sense.append(e.split("_")[1]) 	
		mismatches.append(e.split("_")[2])

	return(match,sense,mismatches)

def mismatches (nb,primer,seq,nbmismatch) :					#function to find match(es) with at least two mismatches 
	match = []
	sense = []
	invP=inverComp(primer)	
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
	return (match,sense,tmpfind)					#return results (position of match + sense of primer)


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

	elif nbmismatch == 1 and match == [] :                    	#if still no matches (perfect match)
		match,sense,mismatchs = mismatch(1,primer,seq)       	#search with a mismatch 

	elif nbmismatch >= 2 and match == [] :
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
			position=result+1				#get the position of the following nucleotide for the next search
			result=seq.find(primer,position)

	elif nbmismatch == 1 and match == [] :            		#if no perfect match
		match,trash,mismatchs = mismatch(2,primer,seq)	

	elif nbmismatch >= 2 and match == [] :
		match,trash,mismatchs = mismatches(2,primer,seq,nbmismatch)

	return match,mismatchs

def find(primers,fasta,round,nbmismatch) : 			#return the result of the matches 
	
	fasta = fasta.replace(" ","").replace("\t","")		#delete spaces and tabulations
	mismatchs2 = []

	if type(round) is str :
		round = round.replace(",",".") 
		round = float(round)                       
	
	sequences = fasta.split('>')				#split the fasta files into a list of fasta file
	del sequences[0]					#delete the '' created
	dico_res = {}

	for s,seq in enumerate(sequences) :			#for each chromosome in the fasta file
		tmp_seq = seq.split("\n")
		title_seq = tmp_seq[0]
		del tmp_seq[0]
		seq = "".join(tmp_seq).upper()

		if primers :					#if primers had been entered
			for p,primer in enumerate(primers) :								#for each couple of primers 
				primer = primer.replace(" ",";").replace("\t",";").split(";")
				primer_info = primer[0].split('_')
				tmp1, tmp2, mismatchs = findFirst(primer[1],seq,nbmismatch)
				if mismatchs != [] : mismatchs = clean_mismatches(1,primer[1],tmp2,mismatchs)
				first_match = tmp1, tmp2								#search match(es) for the first primer 
				second_match = []
				result = []
				insert=""
				mismatchs2 = []
				for i,pos_match in enumerate(first_match[0]) :						#for each match of the first primer
					tmp, mismatchs2 = findSec(primer[2],seq,first_match[1][i],nbmismatch)		#search match(es) for the second primer
					if mismatchs2 != [] : mismatchs2 = clean_mismatches(2,primer[2],tmp2,mismatchs2) #lower the mismatched nucleotides 
					second_match.extend(tmp)							
					if second_match != [] :								#if there is a match with the second primer on the complementary DNA sequence
						for pos_match2 in second_match :					#for each match found for the second primer
								if first_match[1][i] == "inv" :
									size = abs(int(pos_match)-int(pos_match2)+len(primer[1])) 	
									size2 = int(pos_match)+len(primer[1])+(len(seq)-pos_match2)	#if primers are separated by the splitted area in the sequence (circular chromosome)
									if size2 < size and contig is False : size = size2
								else :
									size = abs(int(pos_match2)-int(pos_match)+len(primer[2])) 
									size2 =int(pos_match2)+len(primer[2])+(len(seq)-pos_match)
									if size2 < size and contig is False: size = size2					

								if size2==size and contig is False :					#if insert is separated by the splitted area in the sequence (circular chromosome)
									if pos_match < pos_match2 :insert = seq[int(pos_match2):]+seq[:int(pos_match)+len(primer[1])]	
									else : insert = seq[int(pos_match):]+seq[:int(pos_match2)+len(primer[2])]
								else :		
									if pos_match < pos_match2 : insert = seq[int(pos_match):int(pos_match2)+len(primer[2])]
									else : insert = seq[int(pos_match2):int(pos_match)+len(primer[1])]			

								sizeU = abs(float(primer_info[3].upper().replace("U",""))-\
									((float(primer_info[2].lower().replace("bp",""))-size)\
									/float(primer_info[1].lower().replace("bp",""))))		#computation of sizeU
								
								result.append([primer[0],pos_match,pos_match2,size,sizeU,sequence+str(s+1),nbmismatch,primer[1],mismatchs,primer[2],mismatchs2,insert])
								
				if len(result) == 0 and primer_info[0] not in dico_res.keys() :	#if no result
					dico_res[primer_info[0]]=["\t".join([primer[0],primer[1],primer[2]]),"","","","",sequence+str(s+1),nbmismatch,primer[1],"",primer[2],"",insert]	

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
						best_res[5] = best_res[5]+", "+sequence+str(s+1) 	#if there's already a result with perfect matches
					dico_res[primer_info[0]]=best_res			#set the best result as a new key : value in the dictionnary #replace the old dictionnary value if there is one
	return dico_res

def get_empty_locus (dico_result) :
	tmpprimers = []
	for locus in dico_result.keys() :
		if dico_result[locus][4] == '' :
			tmpprimers.append(dico_result[locus][0])
	return tmpprimers

def run (Primers,fasta,round,nbmismatch) :
	tmp = len(Primers)
	tmpPrimers = Primers
	result = {}
	for mismatch_allowed in range(int(nbmismatch)+1) :
		tmp_dico = find(tmpPrimers,fasta,round,mismatch_allowed) 			#no mismacth
		result = dict(result.items() + tmp_dico.items())
		tmpPrimers = get_empty_locus(result)
		nb_match = tmp -len(tmpPrimers) 
		print "results with",mismatch_allowed,"mismatch : ",nb_match,"/",len(Primers)

		tmp = len(tmpPrimers)

		if len(tmpPrimers) == 0 :
			break 

	if len(tmpPrimers) != 0 : 
		print "no match : ", tmp

	return result

def usage() :
	print "./insilico2.py -i <input_directory> -o <output_directory> -p <primers_file> [option -c for contigs]"  

def main() : #run find() for each genome file in the directory with all primers in the primers file
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hm:i:o:p:c", ["help", "mismatch=", "input=", "output=", "primer=", "contig"])
	except getopt.GetoptError as err:
		usage()
		sys.exit(2)
	nb_mismatch = 2
	global sequence
	sequence = "chr"
	global contig 
	contig = False
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-i", "--input"):
			fasta_path = arg
			if fasta_path[-1]=="/" : fasta_path=fasta_path[:-1]
			files = os.listdir(fasta_path)
			files=sorted(files)
		elif opt in ("-p", "--primer"): 
			Primers = open(arg,"r").read()
		elif opt in ("-o", "--output"):
			output_path = arg
			if output_path[-1] != "/" : output_path=output_path+"/"
		elif opt in ("-m","--mismatch"):
			nb_mismatch = int(arg)
		elif opt in ("-c", "--contig"):
			contig = True
			sequence = "contig"
		else:
			assert False, "unhandled option"

	header = ["strain","primer","position1","position2","size","score",sequence
				,"nb_mismatch","primer1","mismatch1","primer2","mismatch2","insert"]
	cr=[header]

	Primers = Primers.split("\n")
	if Primers[-1]=="" : del Primers[-1] 
	Primers_short = [ pri.split("_")[0] for pri in Primers ]

	for i,file in enumerate(files) :
		print file, "\t strain ",i+1,"/",len(files)
		pathfile = fasta_path+"/"+file
		fasta = open(pathfile,"r").read()
		fasta_names = []
		for line in fasta.split("\n") :
			if ">" in line :
				fasta_names.append(line.replace("\n","").split("|")) 
		
		result = run(Primers,fasta,0.25,nb_mismatch)			#use find for each number of mismatch
		
		locus = []
		mlva_score = []
		ch = []
		header_ch = []
		mismatch = []
		header_mismatch = []
		for Primer in Primers_short :
			tmp = ";".join([str(res) for res in result[Primer][1:]]).replace("[","").replace("]","").replace("'","")
			cr.append([file]+[result[Primer][0].split("\t")[0]]+tmp.split(";"))
			if result[Primer][1]=="" :
				result[Primer][6] = 99
			if result[Primer][6] not in [0,99] :
				print " ".join([Primer]+[" P1 :"]+[result[Primer][7]]+["\nmmatch(s):"]+result[Primer][8]+["\n\n"]+[Primer]+["P2 :"]+[result[Primer][9]]+["\nmmatch(s):"]+result[Primer][10]+["\n"])
			locus.append(Primer.split("_")[0])
			mlva_score.append(str(result[Primer][4]))		#scores 
			header_mismatch.append("nbmis_"+Primer.split("_")[0])
			mismatch.append(str(result[Primer][6]))			#number of mismatch allowed for each locus
			if contig is False :
				header_ch.append("ch_"+Primer.split("_")[0])	
				ch.append(result[Primer][5])			#chromosome of the result 

		if len(fasta_names)>1 and contig is False : 			#make firsts column of file depending on the number of chromosomes
			strain_chr=[]
			schr = []
			ref_chr = []
			rchr = []
			for r in range(len(fasta_names)) :
				strain_chr.append("strain_chr"+str(r+1))
				schr.append(fasta_names[r][4][1:].split(",")[0])
				ref_chr.append("ref_chr"+str(r+1))
				rchr.append(fasta_names[r][3].split(".")[0])
			header = ["key"]+strain_chr+ref_chr
			infos = [str(i+1).zfill(3)]+schr+rchr

		if contig is False :											#if it is running on complete genome
			if i==0 : 
				pathfile = output_path+"MLVA_analysis_"+fasta_path.split("/")[-1]+".csv"	
				output = open(pathfile,"w") 								#output is a csv file (delimiter=";")
				if len(fasta_names)>1 :
					output.write(";".join(header+locus+header_ch+header_mismatch)+"\n")  		#header
				else :
					output.write(";".join(["key","strain","ref"]+locus+header_mismatch)+"\n") 	#header
				output = open(pathfile,"a")

			if len(fasta_names) >1 :
				output.write(";".join(infos+mlva_score+ch+mismatch)+"\n")
			else :
				output.write(";".join([str(i+1).zfill(3),fasta_names[0][4][1:].split(",")[0],fasta_names[0][3].split(".")[0]]+mlva_score+mismatch)+"\n")
		
		else :													#if it is running for contigs
			if i==0 :
				pathfile = output_path+"MLVA_analysis_"+fasta_path.split("/")[-1]+".csv"
				output = open(pathfile,"w") 								#output is a csv file (delimiter=";")
				output.write(";".join(["key","Access_number"]+locus+header_mismatch)+"\n")  			#header
				output = open(pathfile,"a")
			output.write(str(i+1).zfill(3)+";"+file+";".join(mlva_score+mismatch) + "\n")
	out = csv.writer(open(output_path+fasta_path.split("/")[-1]+"_output.csv","w"), delimiter=';',quoting=csv.QUOTE_ALL)
	for row in cr :
		out.writerow(row)
	output.close()
	print "MLVA analysis finished for "+fasta_path.split("/")[-1]

	##### creation of mismatch summary txt file #####

	dico_mismatch = {}
	Primers = [p.replace("\t",";") for p in Primers]
	for primer in Primers :
		dico_mismatch[primer.split(";")[0]+"_FOR"] = set([])
		dico_mismatch[primer.split(";")[0]+"_REV"] = set([])

	for line in cr[1:] :
		if line[9] != "" : dico_mismatch[line[1]+"_FOR"].add(line[9])
		if line[11] != "" : dico_mismatch[line[1]+"_REV"].add(line[11])

	dico_mismatch = {key: value for key, value in dico_mismatch.items() if value != set() } 	#delete keys without value(s)
	
	output_mismatch = open(output_path+fasta_path.split("/")[-1]+"_mismatchs.txt","w")
	tmp_file =""
	for primer in Primers :
		if primer.split(";")[0]+"_FOR" in dico_mismatch.keys() :
			tmp_file += primer.split(";")[0]+"_FOR\n"+primer.split(";")[1]+"\n"+"\n".join(list(dico_mismatch[primer.split(";")[0]+"_FOR"]))+"\n\n"
		if primer.split(";")[0]+"_REV" in dico_mismatch.keys() :
			tmp_file += primer.split(";")[0]+"_REV\n"+primer.split(";")[2]+"\n"+"\n".join(list(dico_mismatch[primer.split(";")[0]+"_REV"]))+"\n\n"
	output_mismatch.write(tmp_file)
	output_mismatch.close()

if __name__ == "__main__" :
	main()

