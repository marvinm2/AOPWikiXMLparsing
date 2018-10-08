import xml.etree.ElementTree as ET
tree = ET.parse('aop-wiki-xml-2018-04-01') #This is the function to parse the XML file.
root = tree.getroot()

writtenoutput = open('writtenoutput.txt', 'w') #The file to save all results to.

list = ['Molecular','Cellular','Tissue','Organ']

bpfdict = {} #A dictionary to link Biological Process IDs to their sources
for bpf in root.findall('{http://www.aopkb.org/aop-xml}biological-process'):
	bpfdict[bpf.get('id')] = bpf.find('{http://www.aopkb.org/aop-xml}source').text

bodict = {} #A dictionary to link Biological Object IDs to their sources
for bpf in root.findall('{http://www.aopkb.org/aop-xml}biological-object'):
	bodict[bpf.get('id')] = bpf.find('{http://www.aopkb.org/aop-xml}source').text


#Create a dictionary with all gene names to map to KE and KER descriptions
HGNC= open('HGNCgenelist.txt', 'r') #The HGNC gene name file from the HGNC database
genedict = {} #A dictionary to contain all variants of all genes provided by the HGNC file
for line in HGNC:
	if not 'HGNC ID	Approved Symbol	Approved Name	Previous Symbols	Synonyms	Ensembl ID(supplied by Ensembl)'in line: #To avoid the first line to be added to the dictionary
		a = line.split('\t')
		genedict[a[0]]=[]
		if not a[5][:-1] == '':
			genedict[a[0]].append(a[5][:-1])
		genedict[a[0]].extend((' '+a[1]+' ','('+a[1]+',',' '+a[1]+',','('+a[1]+')',' '+a[1]+')','['+a[1]+']','['+a[1]+',',' '+a[1]+']',' '+a[1]+'.','('+a[1]+' ','['+a[1]+' ',' '+a[2]+' ','('+a[2]+',',' '+a[2]+',','('+a[2]+')',' '+a[2]+')','['+a[2]+']','['+a[2]+',',' '+a[2]+']',' '+a[2]+'.','('+a[2]+' ','['+a[2]+' ')) #From all variants of gene names and synonyms, a range of written varieties are stored in the dictionary
		for item in a[3].split(', '):
			if not item == '':
				genedict[a[0]].extend((' '+item+' ','('+item+',',' '+item+',','('+item+')',' '+item+')','['+item+']','['+item+',',' '+item+']',' '+item+'.','('+item+' ','['+item+' '))
		for item in a[4].split(', '):
			if not item == '':
				genedict[a[0]].extend((' '+item+' ','('+item+',',' '+item+',','('+item+')',' '+item+')','['+item+']','['+item+',',' '+item+']',' '+item+'.','('+item+' ','['+item+' '))
HGNC.close()

KEIDlist=[] #A list to save all KE IDs on the molecular, cellular, tissue or organ level
for l in list:
	q = {} #A dictionary for all KE elements
	q['nodescription'] = 0
	q['hasdescription'] = 0
	q['nomeasurement'] = 0
	q['hasmeasurement'] = 0
	q['nocelterm'] = 0
	q['hascelterm'] = 0
	q['noorganterm'] = 0
	q['hasorganterm'] = 0
	q['noapplicabilitytax'] = 0
	q['hasapplicabilitytax'] = 0
	q['nobepr'] = 0 #be = biological event, pr = biological process
	q['hasbepr'] = 0
	q['nobeob'] = 0 #ob = biological object
	q['hasbeob'] = 0
	q['nobeac'] = 0 #ac = biological action
	q['hasbeac'] = 0
	amountofKEs = 0
	ensembllist = []
	overlapdict = {} #Dictionary to capture all genes in Key Event descriptions
	ode = {} #ode = OverlapDictEnsemble
	for KE in root.findall('{http://www.aopkb.org/aop-xml}key-event'):
		level = KE.find('{http://www.aopkb.org/aop-xml}biological-organization-level').text
		if level == l:
			if not KE.get('id') in KEIDlist:
				KEIDlist.append(KE.get('id'))
			amountofKEs += 1
			if KE.find('{http://www.aopkb.org/aop-xml}description').text == None:
				q['nodescription'] += 1
			else:
				q['hasdescription'] += 1
				overlapdict[KE.get('id')]=[]
				ode[KE.get('id')]=[]
				for key in genedict:
					for item in genedict[key]:
						if item in KE.find('{http://www.aopkb.org/aop-xml}description').text:
							overlapdict[KE.get('id')].append(item)
							ode[KE.get('id')].append(genedict[key][0])
							if not genedict[key][0] in ensembllist:
								ensembllist.append(genedict[key][0])
			if KE.find('{http://www.aopkb.org/aop-xml}cell-term') == None:
				q['nocelterm'] += 1
			else:
				q['hascelterm'] += 1
			if KE.find('{http://www.aopkb.org/aop-xml}organ-term') == None:
				q['noorganterm'] += 1
			else:
				q['hasorganterm'] += 1

			if not KE.find('{http://www.aopkb.org/aop-xml}biological-events') == None:
				if KE.find('{http://www.aopkb.org/aop-xml}biological-events').find('{http://www.aopkb.org/aop-xml}biological-event').get('process-id') == None:
					q['nobepr'] += 1
				else:
					q['hasbepr'] += 1
				if KE.find('{http://www.aopkb.org/aop-xml}biological-events').find('{http://www.aopkb.org/aop-xml}biological-event').get('object-id') == None:
					q['nobeob'] += 1
				else:
					q['hasbeob'] += 1
				if KE.find('{http://www.aopkb.org/aop-xml}biological-events').find('{http://www.aopkb.org/aop-xml}biological-event').get('action-id') == None:
					q['nobeac'] += 1
				else:
					q['hasbeac'] += 1
			else:
				q['nobepr'] += 1
				q['nobeob'] += 1
				q['nobeac'] += 1


	writtenoutput.write('Total amount of KEs at '+ l+ 'level: '+str(amountofKEs)+'\n'+str(q['nodescription'])+' KEs do not have description\n'+str(q['hasdescription'])+' KEs have description, which is '+str(float(100*q['hasdescription'])/float(amountofKEs))+' percent'+'\n'+str(q['nocelterm'])+' KEs do not have celterm\n'+str(q['hascelterm'])+' KEs have celterm\n'+str(q['noorganterm'])+' KEs do not have organterm\n'+str(q['hasorganterm'])+' KEs have organterm\n'+str(q['nobepr']) +' KEs do not have biological process\n'+str(q['hasbepr'])+' KEs have biological process\n'+str(q['nobeob'])+' KEs do not have biological object\n'+str(q['hasbeob'])+' KEs have biological object\n'+str(q['nobeac'])+' KEs do not have biological action\n'+str(q['hasbeac'])+' KEs have biological action'+'\n'+'Total amount of unique genes found in KE descriptions at '+l+' level: '+str(len(ensembllist))+'\n')

	kebp = {} #A dictionary for Key Event Biological Processes (kebp)
	kebp['GO'] =0
	kebp['HP'] = 0
	kebp['MESH']=0
	kebp['MI'] = 0
	kebp['MP'] = 0
	kebp['NBO'] = 0
	kebp['PCO'] = 0
	kebp['VT'] = 0
	kebo = {} #A dictionary for Key Event Biological Objects (kebo)
	kebo['CHEBI'] = 0
	kebo['CL'] =0
	kebo['FMA'] = 0
	kebo['GO'] =0
	kebo['MESH']=0
	kebo['MP'] = 0
	kebo['PCO'] = 0
	kebo['PR'] = 0
	kebo['TAIR'] = 0
	kebo['UBERON'] = 0
	kect={} #A dictionary for Key Event Cell-Terms (kect)
	kect['CL']=0
	kect['WIKI']=0
	keot={} #A dictionary for Key Event Organ-Terms (keot)
	keot['UBERON']=0
	keot['WIKI']=0

	for KE in root.findall('{http://www.aopkb.org/aop-xml}key-event'):
		level = KE.find('{http://www.aopkb.org/aop-xml}biological-organization-level').text
		if level == l:
			if not KE.find('{http://www.aopkb.org/aop-xml}biological-events') == None:
				KEBP=[]
				KEBO=[]
				if not KE.find('{http://www.aopkb.org/aop-xml}biological-events').find('{http://www.aopkb.org/aop-xml}biological-event').get('process-id') == None:
					for bptag in KE.find('{http://www.aopkb.org/aop-xml}biological-events').findall('{http://www.aopkb.org/aop-xml}biological-event'):
						if not bptag.get('process-id') == None:
							if not bpfdict[bptag.get('process-id')] in KEBP:
								KEBP.append(bpfdict[bptag.get('process-id')])
								kebp[bpfdict[bptag.get('process-id')]]+=1
				if not KE.find('{http://www.aopkb.org/aop-xml}biological-events').find('{http://www.aopkb.org/aop-xml}biological-event').get('object-id') == None:
					for botag in KE.find('{http://www.aopkb.org/aop-xml}biological-events').findall('{http://www.aopkb.org/aop-xml}biological-event'):
						if not botag.get('object-id') == None:
							if not bodict[botag.get('object-id')] in KEBO:
								KEBO.append(bodict[botag.get('object-id')])
								kebo[bodict[botag.get('object-id')]]+=1
			if not KE.find('{http://www.aopkb.org/aop-xml}cell-term') == None:
				KECT=[]
				for cterm in KE.findall('{http://www.aopkb.org/aop-xml}cell-term'):
					if not cterm.find('{http://www.aopkb.org/aop-xml}source').text == None:
						if not cterm.find('{http://www.aopkb.org/aop-xml}source').text in KECT:
							KECT.append(cterm.find('{http://www.aopkb.org/aop-xml}source').text)
							kect[cterm.find('{http://www.aopkb.org/aop-xml}source').text]+=1
			if not KE.find('{http://www.aopkb.org/aop-xml}organ-term') == None:
				KEOT=[]
				for oterm in KE.findall('{http://www.aopkb.org/aop-xml}organ-term'):
					if not oterm.find('{http://www.aopkb.org/aop-xml}source').text == None:
						if not oterm.find('{http://www.aopkb.org/aop-xml}source').text in KEOT:
							KEOT.append(oterm.find('{http://www.aopkb.org/aop-xml}source').text)
							keot[oterm.find('{http://www.aopkb.org/aop-xml}source').text]+=1
	writtenoutput.write('\nOntologies for Biological Process: '+str(kebp)+'\nOntologies for Biological Object: '+str(kebo)+'\nOntologies for cell-term: '+str(kect)+'\nOntologies for Organ-term: '+str(keot)+'\n')
	
	writtenoutput.write('\n'+'Key Event Biological Process percentages:\n')
	for key in kebp:
		if kebp[key] != 0:
			percentage = float(100*kebp[key])/float(amountofKEs)
			writtenoutput.write(key+' present for '+str(percentage)+' percent'+'\n')
	writtenoutput.write('\n'+'Key Event Biological Object percentages:\n')
	for key in kebo:
		if kebo[key] != 0:
			percentage = float(100*kebo[key])/float(amountofKEs)
			writtenoutput.write(key+' present for '+str(percentage)+' percent'+'\n')
	writtenoutput.write('\n'+'Key Event Cell Term percentages:\n')
	for key in kect:
		if kect[key] != 0:
			percentage = float(100*kect[key])/float(amountofKEs)
			writtenoutput.write(key+' present for '+str(percentage)+' percent'+'\n')
	writtenoutput.write('\n'+'Key Event Organ Term percentages:\n')
	for key in keot:
		if keot[key] != 0:
			percentage = float(100*keot[key])/float(amountofKEs)
			writtenoutput.write(key+' present for '+str(percentage)+' percent'+'\n')
			
	writtenoutput.write('\n')


	#Create a WikiPathways SPARQL query for Pathways by genes found in KEs, for each level of biological organization
	g = open('SPARQLqueryENSIDsafterHGNCmap'+l+'.txt', 'w')
	g.write('PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>'+'\n'+'PREFIX rdfs:    <http://www.w3.org/2000/01/rdf-schema#>'+'\n'+'PREFIX dcterms: <http://purl.org/dc/terms/>'+'\n'+'SELECT DISTINCT ?pathway str(?ensId) as ?geneProduct'+'\n'+'WHERE {'+'\n'+'    ?geneProduct a wp:GeneProduct . '+'\n'+'    ?geneProduct wp:bdbEnsembl ?ensId .'+'\n'+'    ?geneProduct dcterms:isPartOf ?pathway .'+'\n'+'    ?pathway a wp:Pathway .'+'\n'+'    FILTER (')
	n = 0
	for item in ensembllist:
		if not n==0:
			g.write('") || regex(str(?ensId), "'+item)
		else:
			g.write('regex(str(?ensId), "'+item)
		n+=1
	g.write('")).' + '\n'+'}')
	g.close()

writtenoutput.write('Total amount of KEs from molecular, cellular, tissue and organ level of organization: '+ str(len(KEIDlist))+'\n'+'\n')

#Key Event Relationships (KERs)
KRd = {} #Discionary to capture all info about KERs
KRd['nodescription'] = 0
KRd['hasdescription'] =0
KRd['nowoeb'] = 0 #woeb = Weight of Evidence Biological Plausibility
KRd['haswoeb'] =0
KRd['nowoee'] = 0 #woee = Weight of Evidence Empirical Support
KRd['haswoee'] =0

ode = {}
overlapdict = {}
ensembllist = []
totalamountofKERs = 0
amountofassessedKERs = 0
n=0

ndes=0
nbiopl=0
nemps=0
for KER in root.findall('{http://www.aopkb.org/aop-xml}key-event-relationship'):
	totalamountofKERs += 1
	if KER.find('{http://www.aopkb.org/aop-xml}title').find('{http://www.aopkb.org/aop-xml}downstream-id').text in KEIDlist:
		amountofassessedKERs += 1
		ode[KER.get('id')]=[]
		overlapdict[KER.get('id')]=[]
		if not KER.find('{http://www.aopkb.org/aop-xml}description').text == None:
			KRd['hasdescription'] += 1
			n+=1
			ndes+=1
			for key in genedict:
				for item in genedict[key]:
					if item in KER.find('{http://www.aopkb.org/aop-xml}description').text:
						overlapdict[KER.get('id')].append(item)
						ode[KER.get('id')].append(genedict[key][0])
						if not genedict[key][0] in ensembllist:
							ensembllist.append(genedict[key][0])
		else:
			KRd['nodescription'] += 1
		if not KER.find('{http://www.aopkb.org/aop-xml}weight-of-evidence').find('{http://www.aopkb.org/aop-xml}biological-plausibility').text == None:
			KRd['haswoeb'] += 1
			n+=1
			nbiopl+=1
			for key in genedict:
				for item in genedict[key]:
					if item in KER.find('{http://www.aopkb.org/aop-xml}weight-of-evidence').find('{http://www.aopkb.org/aop-xml}biological-plausibility').text:
						overlapdict[KER.get('id')].append(item)
						ode[KER.get('id')].append(genedict[key][0])
						if not genedict[key][0] in ensembllist:
							ensembllist.append(genedict[key][0])
		else:
			KRd['nowoeb'] += 1
		if not KER.find('{http://www.aopkb.org/aop-xml}weight-of-evidence').find('{http://www.aopkb.org/aop-xml}emperical-support-linkage').text == None:
			KRd['haswoee'] += 1
			n+=1
			nemps+=1
			for key in genedict:
				for item in genedict[key]:
					if item in KER.find('{http://www.aopkb.org/aop-xml}weight-of-evidence').find('{http://www.aopkb.org/aop-xml}emperical-support-linkage').text:
						overlapdict[KER.get('id')].append(item)
						ode[KER.get('id')].append(genedict[key][0])
						if not genedict[key][0] in ensembllist:
							ensembllist.append(genedict[key][0])
		else:
			KRd['nowoee'] += 1

#Write KER data to file
writtenoutput.write('\nKEY EVENT RELATIONSHIPS'+'\n'+'The total amount of KERs: '+str(totalamountofKERs)+'\n'+'The total amount of assessed KERs: '+str(amountofassessedKERs)+'\n'+str(KRd['nodescription'])+ ' KERs do not have description\n'+str(KRd['hasdescription'])+ ' KERs have description, which is '+str(float(100*KRd['hasdescription'])/float(amountofassessedKERs))+' percent'+'\n'+str(KRd['nowoeb']) + ' KERs do not have WoE biological plausibility\n'+str(KRd['haswoeb'])+ ' KERs have WoE biological plausibility\n'+str(KRd['nowoee']) + ' KERs do not have WoE emperical-support-linkage\n'+str(KRd['haswoee'])+ ' KERs have WoE emperical-support-linkage\n'+'\n'+'Total amount of texts for text mapping: '+str(n)+'\n'+'Total amount of descriptions for text mapping: '+str(ndes)+'\n'+'Total amount of biological plausibility for text mapping: '+str(nbiopl)+'\n'+'Total amount of emperical support for text mapping: '+str(nemps)+'\n')



g = open('SPARQLqueryENSIDsafterHGNCmapKER.txt', 'w')
g.write('PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>'+'\n'+'PREFIX rdfs:    <http://www.w3.org/2000/01/rdf-schema#>'+'\n'+'PREFIX dcterms: <http://purl.org/dc/terms/>'+'\n'+'SELECT DISTINCT ?pathway str(?ensId) as ?geneProduct'+'\n'+'WHERE {'+'\n'+'    ?geneProduct a wp:GeneProduct . '+'\n'+'    ?geneProduct wp:bdbEnsembl ?ensId .'+'\n'+'    ?geneProduct dcterms:isPartOf ?pathway .'+'\n'+'    ?pathway a wp:Pathway .'+'\n'+'    FILTER (')
n = 0
for item in ensembllist:
	if not n==0:
		g.write('") || regex(str(?ensId), "'+item)
	else:
		g.write('regex(str(?ensId), "'+item)
	n+=1
g.write('")).' + '\n'+'}')
g.close()

strd = {}
strd['nochem'] = 0
strd['haschem'] = 0
chemIDlist = []
amountofStressors = 0
totalchemicals = 0
for st in root.findall('{http://www.aopkb.org/aop-xml}stressor'):
	amountofStressors += 1
	if not st.find('{http://www.aopkb.org/aop-xml}chemicals')==None:
		strd['haschem'] += 1
		for chem in st.find('{http://www.aopkb.org/aop-xml}chemicals').findall('{http://www.aopkb.org/aop-xml}chemical-initiator'):
			totalchemicals+=1
			if not chem.get('chemical-id') == None:
				if not chem.get('chemical-id') in chemIDlist:
					chemIDlist.append(chem.get('chemical-id'))
	else:
		strd['nochem'] += 1

#Part to store all CAS RNs related to stressors in a file
amountCAS=0
g = open('ListofCasrns.txt', 'w')
for chemical in root.findall('{http://www.aopkb.org/aop-xml}chemical'):
	if chemical.get('id') in chemIDlist:
		if not 'NOCAS' in chemical.find('{http://www.aopkb.org/aop-xml}casrn').text:
			g.write(chemical.find('{http://www.aopkb.org/aop-xml}casrn').text + '\n')
			amountCAS+=1
g.close()

#Write stressor data in file
writtenoutput.write('\n'+'STRESSORS'+'\n'+'Total amount of stressors: '+ str(amountofStressors)+'\n'+ 'Stressors without chemical: '+str(strd['nochem'])+'\n'+ 'Stressors with chemical: '+str(strd['haschem'])+'\n'+'Total amount of chemicals linked to stressors: '+str(totalchemicals)+'\n'+'Amount of unique chemicals linked to stressors: '+str(len(chemIDlist))+'\n')
writtenoutput.write('\n'+ 'Amount of CAS RNs related to stressors: '+str(amountCAS)+'\n')
writtenoutput.close()




