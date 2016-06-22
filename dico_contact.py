import pickle, csv
# -*-coding:Latin-1 -*

dico_publi = {}
dico_publi["Hauck2012"] ="Hauck2012(http://www.ncbi.nlm.nih.gov/pubmed/22984530)"
dico_publi["Pourcel2004"] ="Pourcel2004(http://www.ncbi.nlm.nih.gov/pubmed/15186506)"
dico_publi["Pourcel2011"] ="Pourcel2011(http://www.ncbi.nlm.nih.gov/pubmed/21147956)"

dico_contact = {}

dico_contact["Christine Pourcel"] ="Christine Pourcel(christine.pourcel@u-psud.fr)"
dico_contact["Lenie Dijkshoorn"] ="Lenie Dijkshoorn(l.dijkshoorn@lumc.nl)"
dico_contact["Kirnpal Kaur Banga Singh"] ="Kirnpal Kaur Banga Singh(kiren@kck.usm.my)"

cr = csv.reader(open("/home/david/Documents/MLVA/contact.csv","r"),delimiter=";")

for row in cr :
	key_publi = row[0]
	key_contact= row[1]
	publi = row[2]
	contact = row[3]

	dico_contact[key_contact] = key_contact+"("+contact.replace("mailto:","")+")"
	dico_publi[key_publi] = key_publi+"("+publi+")"
	
del dico_contact[""]
del dico_publi[""]

filehandler = open("/home/david/Documents/MLVA/dico_publi","wb")
pickle.dump(dico_publi,filehandler)
filehandler.close

filehandler = open("/home/david/Documents/MLVA/dico_contact","wb")
pickle.dump(dico_contact,filehandler)
filehandler.close


