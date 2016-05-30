import csv,sys,os.path,pickle

########## traitement du fichier suplementary data de Roberta Creti ##################
## mise en forme du csv pyogenes a partir du csv extrait du suplementary data de Roberta Creti
## silico## pour les souches silico qui n'ont pas de numero

cr = csv.reader(open("/home/david/Documents/csv_mlva/Supplementary_table_pyogenes.csv","r"),delimiter=";")

tmp = cr.next()
header = ["key","Strain ID","source"]
header.extend(tmp[2:])

output = open("/home/david/Documents/csv_mlva/Supplementary_table_pyogenes_corrected.csv","w")
output.write(";".join(header)+"\n")
output = open("/home/david/Documents/csv_mlva/Supplementary_table_pyogenes_corrected.csv","a")

i=1
for row in cr :
	if row[0] == '' :
		key = "silico"+str(i).zfill(3)
		i+=1
	else :
		key = "Imperi"+row[0].zfill(3)
	 
	if "(" in row[1] :
		source = row[1].split(" ")[1].split("(")[1].split(")")[0]
	else : 
		source =''	
	output.write(";".join([key,row[1].split(" ")[0],source,";".join(row[2:])] )+"\n")
output.close()


######## ajout de la colonne key en premiere position #############

if len(sys.argv) != 1 :
	path = sys.argv[1]
else :
	path = "/home/david/Documents/csv_mlva"

files = os.listdir(path)
if os.path.exists(path+"/corrected") is False :
	os.mkdir(path+"/corrected")

for file in files :
	exceptions = ["Supplementary_table_pyogenes.csv","Propionibacterium_access_number.csv","Propionibacterium_acnes_complement.csv"]
	if os.path.isfile(path+"/"+file) and file not in exceptions :
		print file
		input = open(path+"/"+file,"r")
		output = open(path+"/corrected/"+file,"w")
		header = input.readline()
		output.write("key;study;"+";".join(header.split(";")[1:]))
		output = open(path+"/corrected/"+file,"a")

		for i,line in enumerate(input) :
			if len(line)!=0 :
				output.write(str(i+1).zfill(3)+";"+line)
		input.close()
		output.close()

######## suppression de la colonne baseView #######################

path = "/home/david/Documents/csv_mlva/corrected"
files = os.listdir(path)


for file in files :
	input = open(path+"/"+file,"r")
	tmp = []
	cr = csv.reader(input,delimiter=";")
	header = cr.next()
	header_low = [e.lower() for e in header]
	position = -1
	if "baseview" in header_low :
		position = header_low.index("baseview")
		del header[position]
	tmp = ";".join(header) + "\n"
	for row in cr :
		if position != -1 :
			del row[position]
		tmp += ";".join(row) + "\n"

	with open(path+"/"+file,"w") as f :
		f.write(tmp)




################# ajout hypertext pour les contacts et publis ###############

dico_contact = pickle.load(open("/home/david/Documents/scripts/dico_contact","r"))
dico_publi = pickle.load(open("/home/david/Documents/scripts/dico_publi","r"))

path = "/home/david/Documents/csv_mlva/corrected"
files = os.listdir(path)

#on ajoute les mails pour les contact
for file in files :
	cr = csv.reader(open(path+"/"+file,"r"),delimiter=";")
	header = cr.next()
	tmp = ";".join(header)+"\n"
	if "contact" in header :
		position = header.index("contact")
		for line in cr :
			line[position]
			if line[position] in dico_contact.keys() :
				line [position] = dico_contact[line[position]] 
			tmp += ";".join(line)+"\n" 
		output = open(path+"/"+file,"w")
		output.write(tmp)
		output.close()

#on ajoute les liens pubmed 

for file in files :
	cr = csv.reader(open(path+"/"+file,"r"),delimiter=";")
	header = cr.next()
	tmp = ";".join(header)+"\n"
	if "publication" in header :
		position = header.index("publication")
		for line in cr :
			line[position].replace(" ","")
			if line[position] in dico_publi.keys() :
				line [position] = dico_publi[line[position]] 
			tmp += ";".join(line)+"\n" 
		output = open(path+"/"+file,"w")
		output.write(tmp)
		output.close()


##########traitement du fichier Acinetobacter.csv ##############

## compacte le nom de genus avec le nom de species en X. xxxxxx 

path = "/home/david/Documents/csv_mlva/corrected/Acinetobacter.csv"

cr = csv.reader(open(path,"r"),delimiter=";")
tmp = cr.next()
del tmp[4]
tmp2 = ";".join(tmp)+"\n"

for row in cr : 
	row [4] = row[4][0]+". "+row[5]
	del row[5]
	tmp2 += ";".join(row)+"\n"

with open(path,"w") as f :
	f.write(tmp2)
	f.close()

########## traitement du fichier Streptococcus pneumoniae.csv ###########

cr = csv.reader(open("/home/david/Documents/csv_mlva/corrected/Streptococcus_pneumoniae.csv","r"),delimiter=";")
tmp = ";".join(cr.next())+"\n"

for row in cr :
	row[7] = row[7].split(" ")[0]
	tmp += ";".join(row)+"\n"

with open("/home/david/Documents/csv_mlva/corrected/Streptococcus_pneumoniae.csv","w") as f :
	f.write(tmp)
	f.close()

############# Mycobacterium tuberculosis #################

cr = csv.reader(open("/home/david/Documents/csv_mlva/corrected/Mycobacterium_tuberculosis.csv","r"),delimiter=";")
tmp = cr.next()
del tmp[2]
del tmp[12:18]
tmp = ";".join(tmp)+"\n"

for row in cr :
	if row[4] and row[4][0] == "," :
		row[4] = row[4][2:]

	row[12] = str(row[12]).zfill(15) #on ajoute les 0 manquants

	del row[2]
	del row[12:19]

	tmp+= ";".join(row)+"\n"

output = open("/home/david/Documents/csv_mlva/corrected/Mycobacterium_tuberculosis.csv","w")
output.write(tmp)
output.close()

############# Pseudomonas_syringae_pv_actinidiae #################

cr = csv.reader(open("/home/david/Documents/csv_mlva/corrected/Pseudomonas_syringae_pv_actinidiae.csv","r"),delimiter=";")
header = cr.next()
tmp = ";".join(header[0:2])+";Strain_ID;"+";".join(header[2:4])+";publication;contact;isolated_in;"+";".join(header[7:])+"\n"


i=1
for row in cr :
	#on remplace les colonnes Country, Region et Locality par isolated in
	row[4] = ", ".join([row[6],row[5],row[4]])
	if row[4][0] == ',' :
		row[4] = row[4].replace(", , ","")
		row[4] = row[4].replace(", ","",1)
	del row[5:7]
	#ajout des colonnes publication et contact
	row = row[0:4]+["Ciarroni2015(http://www.ncbi.nlm.nih.gov/pubmed/26262683)","Angelo Mazzaglia(angmazza@unitus.it)"]+row[4:]

	#on modifie l'accession number obselete #ajout de la colonne Strain ID
	tmp2 = row[1]
	if tmp2 == "NZ_CM002752.1" :
		row[1] = "silico001"
		tmp2 = "ICMP18884"
		row[4] = "NZ_CP011972(http://www.ncbi.nlm.nih.gov/nuccore/NZ_CP011972)"
	else :
		row[1] = "2015Ciarroni"+str(i).zfill(3)
		i+=1
	row = row[0:2]+[tmp2]+row[2:]
	tmp += ";".join(row)+"\n"

with open("/home/david/Documents/csv_mlva/corrected/Pseudomonas_syringae_pv_actinidiae.csv","w") as f :
	f.write(tmp)
	f.close()
	


################ Propionibacterium ###############

cr = csv.reader(open("/home/david/Documents/csv_mlva/corrected/Propionibacterium_acnes.csv","r"),delimiter=";")
header = cr.next()
del header[9]
tmp = ";".join(header[0:3])+";accession_number;publication;contact;"+header[3]+";site;hospital;year;"+";".join(header[4:6])+";Belfast;"+";".join(header[7:])+"\n"

cr2 = csv.reader(open("/home/david/Documents/csv_mlva/Propionibacterium_access_number.csv","r"),delimiter=";")
cr2.next()
access = []
for row in cr2 :
	access.append(row[2])

nbline=0
for i,row in enumerate(cr) :
	del row[9] # on supprime la colonne MLVA7
	
	#on met silico dans la colonne collection
	row[1] = "silico"+str(i+1).zfill(3)
	
	tmp+= ";".join(row[0:3]+[access[i].split(".")[0]+"(http://www.ncbi.nlm.nih.gov/nuccore/"+access[i]+")"]+["Hauck2015(http://www.ncbi.nlm.nih.gov/pubmed/25965840)","Christine Pourcel(christine.pourcel@i2bc.paris-saclay.fr)"]+[row[3]]+["","",""]+row[4:])+"\n"
	nbline+=1

with open("/home/david/Documents/csv_mlva/corrected/Propionibacterium_acnes.csv","w") as f :
	f.write(tmp)
	f.close()

cr3 = csv.reader(open("/home/david/Documents/csv_mlva/Propionibacterium_acnes_complement.csv","r"),delimiter=";")
cr3.next()

f = open("/home/david/Documents/csv_mlva/corrected/Propionibacterium_acnes.csv","a")

for line in cr3 :
	nbline+=1
	line = [str(nbline).zfill(3)]+["Hauck2015_"+str(line[2][len(line[2])-2:]).zfill(3)]+[line[2]]+[" "]+["Hauck2015(http://www.ncbi.nlm.nih.gov/pubmed/25965840)","Christine Pourcel(christine.pourcel@i2bc.paris-saclay.fr)"]+line[3:7]+[";;;;"]+line[7:]
	f.write(";".join(line)+"\n")
f.close()

################# Chlamydia_psittaci.csv ###################

file = open("/home/david/Documents/csv_mlva/corrected/Chlamydophila_psittacii.csv","r")
cr = csv.reader(file,delimiter=";")
header = cr.next()
position = header.index("MLVA8") #get the position of the MLVA8 column 

mlva = [int(row[position]) for row in cr if 'temp' not in row[position] and row[position] != '' ] #get all the genome number of the MLVA8 column
next_gn = max(mlva) +1             
del mlva

file.seek(0) #return to the beginning of file (after the header)
cr.next()
 
tmp = ";".join(header) + "\n"
for row in cr :
	if 'temp' in row[position] :
		row[position] = str(next_gn)
		next_gn += 1
	if row[8] == "" :
		row[7] = row[1].split(".")[0]+"(http://www.ncbi.nlm.nih.gov/nuccore/"+row[1].split(".")[0]+")"
		row[1] = "in silico"
		row[6] = "David Christiany(David.christiany@i2bc.paris-saclay.fr)"
	
	tmp += ";".join(row) +"\n"

output = open("/home/david/Documents/csv_mlva/corrected/Chlamydophila_psittacii.csv","w")
output.write(tmp)
output.close()


#################### Coxiella burnetii ###########################

file = open("/home/david/Documents/csv_mlva/corrected/Coxiella_burnetii.csv","r")
cr = csv.reader(file,delimiter=";")
tmp = ";".join(cr.next()) +"\n"

for row in cr :
	if row[1][0:2] in ["NC","HG"] :
		row[7] = row[1].split(".")[0]+"(http://www.ncbi.nlm.nih.gov/nuccore/"+row[1].split(".")[0]+")"
		row[1] = "in silico"
		row[8] = "David Christiany(David.christiany@i2bc.paris-saclay.fr)"
	tmp += ";".join(row) + "\n"

with open("/home/david/Documents/csv_mlva/corrected/Coxiella_burnetii.csv","w") as output :
	output.write(tmp)








