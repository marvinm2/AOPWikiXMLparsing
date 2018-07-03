import xml.etree.ElementTree as ET
tree = ET.parse('aop-wiki-xml-2018-04-01')
root = tree.getroot()

writtenoutput = open('writtenoutput.txt', 'w')

list = ['Molecular','Cellular','Tissue','Organ']
godict = {}
bpfdict = {}
for bpf in root.findall('{http://www.aopkb.org/aop-xml}biological-process'):
	bpfdict[bpf.get('id')] = bpf.find('{http://www.aopkb.org/aop-xml}source').text
	if bpf.find('{http://www.aopkb.org/aop-xml}source').text == 'GO':
		godict[bpf.get('id')]=bpf.find('{http://www.aopkb.org/aop-xml}source-id').text
GOquery = open('GOsforKEsquery.txt', 'w')
GOquery.write('SELECT * FROM term \n WHERE ')
gonumber=0
for key in godict:
	if not gonumber==0:
		GOquery.write("or acc='"+godict[key]+"';")
	else:
		GOquery.write("acc='"+godict[key]+"';")
		gonumber+=1
GOquery.close()

bodict = {}
for bpf in root.findall('{http://www.aopkb.org/aop-xml}biological-object'):
	bodict[bpf.get('id')] = bpf.find('{http://www.aopkb.org/aop-xml}source').text
justGOs=open('GOlist.txt', 'w')
fileforGOs = open('GOsforKEs.txt', 'w')

HGNC= open('HGNCgenelist.txt', 'r')
genedict = {}

for line in HGNC:
	if not 'HGNC ID	Approved Symbol	Approved Name	Previous Symbols	Synonyms	Ensembl ID(supplied by Ensembl)'in line:
		a = line.split('\t')
		genedict[a[0]]=[]
		if not a[5][:-1] == '':
			genedict[a[0]].append(a[5][:-1])
		genedict[a[0]].append(' '+a[1]+' ')
		genedict[a[0]].append('('+a[1]+',')
		genedict[a[0]].append(' '+a[1]+',')
		genedict[a[0]].append('('+a[1]+')')
		genedict[a[0]].append(' '+a[1]+')')
		genedict[a[0]].append('['+a[1]+']')
		genedict[a[0]].append('['+a[1]+',')
		genedict[a[0]].append(' '+a[1]+']')
		genedict[a[0]].append(' '+a[1]+'.')
		genedict[a[0]].append('('+a[1]+' ')
		genedict[a[0]].append('['+a[1]+' ')
		genedict[a[0]].append(' '+a[2]+' ')
		genedict[a[0]].append('('+a[2]+',')
		genedict[a[0]].append(' '+a[2]+',')
		genedict[a[0]].append('('+a[2]+')')
		genedict[a[0]].append(' '+a[2]+')')
		genedict[a[0]].append('['+a[2]+']')
		genedict[a[0]].append('['+a[2]+',')
		genedict[a[0]].append(' '+a[2]+']')
		genedict[a[0]].append(' '+a[2]+'.')
		genedict[a[0]].append('('+a[2]+' ')
		genedict[a[0]].append('['+a[2]+' ')
		for item in a[3].split(', '):
			if not item == '':
				genedict[a[0]].append(' '+item+' ')
				genedict[a[0]].append('('+item+',')
				genedict[a[0]].append(' '+item+',')
				genedict[a[0]].append('('+item+')')
				genedict[a[0]].append(' '+item+')')
				genedict[a[0]].append('['+item+']')
				genedict[a[0]].append('['+item+',')
				genedict[a[0]].append(' '+item+']')
				genedict[a[0]].append(' '+item+'.')
				genedict[a[0]].append('('+item+' ')
				genedict[a[0]].append('['+item+' ')
		for item in a[4].split(', '):
			if not item == '':
				genedict[a[0]].append(' '+item+' ')
				genedict[a[0]].append('('+item+',')
				genedict[a[0]].append(' '+item+',')
				genedict[a[0]].append('('+item+')')
				genedict[a[0]].append(' '+item+')')
				genedict[a[0]].append('['+item+']')
				genedict[a[0]].append('['+item+',')
				genedict[a[0]].append(' '+item+']')
				genedict[a[0]].append(' '+item+'.')
				genedict[a[0]].append('('+item+' ')
				genedict[a[0]].append('['+item+' ')
HGNC.close()

KEIDlist=[]




for l in list:
	q = 'KE'+l+'dict'
	q = {}
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
	q['nobepr'] = 0
	#be = biological event, pr = biological process
	q['hasbepr'] = 0
	q['nobeob'] = 0
	#ob = biological object
	q['hasbeob'] = 0
	q['nobeac'] = 0
	#ac = biological action
	q['hasbeac'] = 0

	amountofKEs = 0
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
			if KE.find('{http://www.aopkb.org/aop-xml}measurement-methodology').text == None:
				q['nomeasurement'] += 1
			else:
				q['hasmeasurement'] += 1
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


	writtenoutput.write('Total amount of KEs for KEs at '+ l+ 'level:')
	writtenoutput.write(str(amountofKEs))
	writtenoutput.write('\n')
	writtenoutput.write(str(q['nodescription'])+'no description\n')
	writtenoutput.write(str(q['hasdescription'])+'has description\n')
	writtenoutput.write(str(q['nomeasurement'])+'no measurement\n')
	writtenoutput.write(str(q['hasmeasurement'])+'has measurement\n')
	writtenoutput.write(str(q['nocelterm'])+'no celterm\n')
	writtenoutput.write(str(q['hascelterm'])+'has celterm\n')
	writtenoutput.write(str(q['noorganterm'])+'no organterm\n')
	writtenoutput.write(str(q['hasorganterm'])+'has organterm\n')
	writtenoutput.write(str(q['nobepr']) +'no biological process\n')
	writtenoutput.write(str(q['hasbepr'])+'has biological process\n')
	writtenoutput.write(str(q['nobeob'])+'no biological object\n')
	writtenoutput.write(str(q['hasbeob'])+'has biological object\n')
	writtenoutput.write(str(q['nobeac'])+'no biological action\n')
	writtenoutput.write(str(q['hasbeac'])+'has biological action\n')
	writtenoutput.write(str(q['noapplicabilitytax']) +'no taxonomy\n')
	writtenoutput.write(str(q['hasapplicabilitytax']) +'has taxonomy\n')

	writtenoutput.write('\n')

	total=0
	kebp = {}
	kebp['GO'] =0
	kebp['MI'] = 0
	kebp['MESH']=0
	kebp['MP'] = 0
	kebp['MI'] = 0
	kebp['VT'] = 0
	kebp['HP'] = 0
	kebp['NBO'] = 0
	kebp['PCO'] = 0
	kebo = {}
	kebo['CL'] =0
	kebo['GO'] =0
	kebo['FMA'] = 0
	kebo['MESH']=0
	kebo['MP'] = 0
	kebo['PR'] = 0
	kebo['CHEBI'] = 0
	kebo['UBERON'] = 0
	kebo['PCO'] = 0

	kect={}
	kect['CL']=0
	kect['WIKI']=0
	keot={}
	keot['UBERON']=0
	keot['WIKI']=0

	for KE in root.findall('{http://www.aopkb.org/aop-xml}key-event'):
		level = KE.find('{http://www.aopkb.org/aop-xml}biological-organization-level').text
		if level == l:
			if not KE.find('{http://www.aopkb.org/aop-xml}biological-events') == None:
				if not KE.find('{http://www.aopkb.org/aop-xml}biological-events').find('{http://www.aopkb.org/aop-xml}biological-event').get('process-id') == None:
					for gotag in KE.find('{http://www.aopkb.org/aop-xml}biological-events').findall('{http://www.aopkb.org/aop-xml}biological-event'):
						if not gotag.get('process-id') == None:
							kebp[bpfdict[gotag.get('process-id')]]+=1
							total += 1
							if gotag.get('process-id') in godict:
								fileforGOs.write(str(KE.find('{http://www.aopkb.org/aop-xml}title').text))
								fileforGOs.write( '\t'+l+'\t')
								fileforGOs.write(str(godict[gotag.get('process-id')])+'\n')
								justGOs.write(str(godict[gotag.get('process-id')])+'\n')
				if not KE.find('{http://www.aopkb.org/aop-xml}biological-events').find('{http://www.aopkb.org/aop-xml}biological-event').get('object-id') == None: 
					kebo[bodict[KE.find('{http://www.aopkb.org/aop-xml}biological-events').find('{http://www.aopkb.org/aop-xml}biological-event').get('object-id')]]+=1
			if not KE.find('{http://www.aopkb.org/aop-xml}cell-term') == None:
				kect[KE.find('{http://www.aopkb.org/aop-xml}cell-term').find('{http://www.aopkb.org/aop-xml}source').text]+=1
			if not KE.find('{http://www.aopkb.org/aop-xml}organ-term') == None:
				keot[KE.find('{http://www.aopkb.org/aop-xml}organ-term').find('{http://www.aopkb.org/aop-xml}source').text]+=1
	writtenoutput.write('\nOntologies for Biological Process')
	writtenoutput.write(kebp)
	writtenoutput.write('\nOntologies for Biological Object')
	writtenoutput.write(kebo)
	writtenoutput.write('\nOntologies for cell-term')
	writtenoutput.write(kect)
	writtenoutput.write('\nOntologies for Organ-term')
	writtenoutput.write(keot)
	writtenoutput.write('\n')
	writtenoutput.write('\ntotal biological process:'+str(total))
	for key in kebp:
		if kebp[key] != 0:
			percentage = float(100*kebp[key])/float(total)
			writtenoutput.write(key+' present for '+str(percentage))
			
			
			
	ensembllist = []
	n=0
	o=0
	overlapdict = {}
	ode = {}
	for KE in root.findall('{http://www.aopkb.org/aop-xml}key-event'):
		level = KE.find('{http://www.aopkb.org/aop-xml}biological-organization-level').text
		if level == l:
			n+=1
			if not KE.find('{http://www.aopkb.org/aop-xml}description').text == None:
				overlapdict[KE.get('id')]=[]
				ode[KE.get('id')]=[]
				o+=1
				for key in genedict:
					for item in genedict[key]:
						if item in KE.find('{http://www.aopkb.org/aop-xml}description').text:
							overlapdict[KE.get('id')].append(item)
							ode[KE.get('id')].append(genedict[key][0])
							if not genedict[key][0] in ensembllist:
								ensembllist.append(genedict[key][0])

	writtenoutput.write(l)
	writtenoutput.write('\nTotal KEs')
	writtenoutput.write(n)
	writtenoutput.write('\nTotal with description')
	writtenoutput.write(o)
	writtenoutput.write('\nTotal amount of unique genes:')
	writtenoutput.write(len(ensembllist))
	
	g = open('HGNCmappedgenesonKEs'+l+'.txt', 'w')
	for item in overlapdict:
		g.write(item )
		for gene in overlapdict[item]:
			g.write( '\t'+gene)
		g.write('\n')
	g.close()

	g = open('SPARQLqueryENSIDsafterHGNCmap'+l+'.txt', 'w')
	g.write('PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>'+'\n'+'PREFIX rdfs:    <http://www.w3.org/2000/01/rdf-schema#>'+'\n'+'PREFIX dcterms: <http://purl.org/dc/terms/>'+'\n'+'SELECT DISTINCT ?pathway str(?ensId) as ?geneProduct'+'\n'+'WHERE {'+'\n'+'    ?geneProduct a wp:GeneProduct . '+'\n'+'    #?geneProduct rdfs:label ?label .'+'\n'+'    ?geneProduct wp:bdbEnsembl ?ensId .'+'\n'+'    ?geneProduct dcterms:isPartOf ?pathway .'+'\n'+'    ?pathway a wp:Pathway .'+'\n'+'    FILTER (')
	n = 0
	for item in ensembllist:
		if not n==0:
			g.write('") || regex(str(?ensId), "'+item)
		else:
			g.write('regex(str(?ensId), "'+item)
		n+=1
	g.write('")).' + '\n'+'}')
	g.close()
	writtenoutput.write('\n')
			
			
			
			
			
			
			
			
justGOs.close()
fileforGOs.close()
writtenoutput.write('total amount of KEs from molecular, cellular, tissue and organ:'+ len(KEIDlist))
KRd = {}
KRd['nodescription'] = 0
KRd['hasdescription'] =0
KRd['nowoeb'] = 0
KRd['haswoeb'] =0
KRd['nowoee'] = 0
KRd['haswoee'] =0
KRd['nowoeu'] = 0
KRd['haswoeu'] =0
KRd['notaxapp'] = 0
KRd['hastaxapp'] =0
KRd['notaxappev'] = 0
KRd['hastaxappev'] =0


ode = {}
overlapdict = {}
ensembllist = []
KERwithDes = []
amountofKERs = 0
n=0

ndes=0
nbiopl=0
nemps=0
for KER in root.findall('{http://www.aopkb.org/aop-xml}key-event-relationship'):
	amountofKERs += 1
	ode[KER.get('id')]=[]
	overlapdict[KER.get('id')]=[]
	if not KER.find('{http://www.aopkb.org/aop-xml}description').text == None:
		KRd['hasdescription'] += 1
		if KER.find('{http://www.aopkb.org/aop-xml}title').find('{http://www.aopkb.org/aop-xml}downstream-id').text in KEIDlist:
			n+=1
			ndes+=1
			o+=1
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
		if KER.find('{http://www.aopkb.org/aop-xml}title').find('{http://www.aopkb.org/aop-xml}downstream-id').text in KEIDlist:
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
		if KER.find('{http://www.aopkb.org/aop-xml}title').find('{http://www.aopkb.org/aop-xml}downstream-id').text in KEIDlist:
			n+=1
			nemps+=1
			o+=1
			for key in genedict:
				for item in genedict[key]:
					if item in KER.find('{http://www.aopkb.org/aop-xml}weight-of-evidence').find('{http://www.aopkb.org/aop-xml}emperical-support-linkage').text:
						overlapdict[KER.get('id')].append(item)
						ode[KER.get('id')].append(genedict[key][0])
						if not genedict[key][0] in ensembllist:
							ensembllist.append(genedict[key][0])
	else:
		KRd['nowoee'] += 1
	if not KER.find('{http://www.aopkb.org/aop-xml}weight-of-evidence').find('{http://www.aopkb.org/aop-xml}uncertainties-or-inconsistencies').text == None:
		KRd['haswoeu'] += 1
	else:
		KRd['nowoeu'] += 1
	if KER.find('{http://www.aopkb.org/aop-xml}taxonomic-applicability').text == None or KER.find('{http://www.aopkb.org/aop-xml}taxonomic-applicability').text == '\n    ':
		KRd['notaxapp'] += 1
	else:
		KRd['hastaxapp'] += 1
	if KER.find('{http://www.aopkb.org/aop-xml}evidence-supporting-taxonomic-applicability').text == None or KER.find('{http://www.aopkb.org/aop-xml}evidence-supporting-taxonomic-applicability').text == '\n    ':
		KRd['notaxappev'] += 1
	else:
		KRd['hastaxappev'] += 1
writtenoutput.write('\nKEY EVENT RELATIONSHIPS')
writtenoutput.write('The total amount of KERs:')
writtenoutput.write(str(amountofKERs))
writtenoutput.write('\n')
writtenoutput.write(str(KRd['nodescription'])+ 'no description\n')
writtenoutput.write(str(KRd['hasdescription'])+ 'has description\n')
writtenoutput.write(str(KRd['nowoeb']) + 'no WoE biological plausibility\n')
writtenoutput.write(str(KRd['haswoeb'])+ 'has WoE biological plausibility\n')
writtenoutput.write(str(KRd['nowoee']) + 'no WoE emperical-support-linkage\n')
writtenoutput.write(str(KRd['haswoee'])+ 'has WoE emperical-support-linkage\n')
writtenoutput.write(str(KRd['notaxapp']) +'no taxonomic-applicability\n')
writtenoutput.write(str(KRd['hastaxapp'])+ 'has taxonomic-applicability')
writtenoutput.write('\n Total amount of texts for text mapping:')
writtenoutput.write(n)
writtenoutput.write('\n Total amount of descriptions for text mapping:')
writtenoutput.write(ndes)
writtenoutput.write('\n Total amount of biological plausibility for text mapping:')
writtenoutput.write(nbiopl)
writtenoutput.write('\n Total amount of emperical support for text mapping:')
writtenoutput.write(nemps)

g = open('HGNCmappedgenesonKERs.txt', 'w')
for item in overlapdict:
	g.write(item )
	for gene in overlapdict[item]:
		g.write( '\t'+gene)
	g.write('\n')
g.close()
g = open('HGNCmappedgenesonKERse.txt', 'w')
for item in ode:
	g.write(item )
	for a in item[1:]:
		g.write('\t')
		g.write(str(ode[item]))
	g.write('\n')
g.close()

g = open('SPARQLqueryENSIDsafterHGNCmapKER.txt', 'w')
g.write('PREFIX wp:      <http://vocabularies.wikipathways.org/wp#>'+'\n'+'PREFIX rdfs:    <http://www.w3.org/2000/01/rdf-schema#>'+'\n'+'PREFIX dcterms: <http://purl.org/dc/terms/>'+'\n'+'SELECT DISTINCT ?pathway str(?ensId) as ?geneProduct'+'\n'+'WHERE {'+'\n'+'    ?geneProduct a wp:GeneProduct . '+'\n'+'    #?geneProduct rdfs:label ?label .'+'\n'+'    ?geneProduct wp:bdbEnsembl ?ensId .'+'\n'+'    ?geneProduct dcterms:isPartOf ?pathway .'+'\n'+'    ?pathway a wp:Pathway .'+'\n'+'    FILTER (')
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
strd['nodescription'] = 0
strd['hasdescription'] = 0
strd['nochem'] = 0
strd['haschem'] = 0
strd['noexp']=0
strd['hasexp'] =0

amountofStressors = 0
for st in root.findall('{http://www.aopkb.org/aop-xml}stressor'):
	amountofStressors += 1
	if st.find('{http://www.aopkb.org/aop-xml}description').text == None:
		strd['nodescription'] += 1
	else:
		strd['hasdescription'] += 1
	if st.find('{http://www.aopkb.org/aop-xml}chemicals') == None:
		strd['nochem'] += 1
	else:
		strd['haschem'] += 1
	if st.find('{http://www.aopkb.org/aop-xml}exposure-characterization').text == None:
		strd['noexp'] += 1
	else:
		strd['hasexp'] += 1
		
writtenoutput.write('\n STRESSORS\n')
writtenoutput.write('Total amount:'+ str(amountofStressors)+'\n')
writtenoutput.write(str(strd['nodescription'])+ 'no description\n')
writtenoutput.write(str(strd['hasdescription'])+ 'has description\n')
writtenoutput.write(str(strd['nochem'])+ 'no chem\n')
writtenoutput.write(str(strd['haschem'])+ 'has chem\n')
writtenoutput.write(str(strd['noexp'])+ 'no exp\n')
writtenoutput.write(str(strd['hasexp'])+ 'has exp\n')

total = 0
no = 0
for chemical in root.findall('{http://www.aopkb.org/aop-xml}chemical'):
	total+= 1
	if chemical.find('{http://www.aopkb.org/aop-xml}jchem-inchi-key').text == None:
		no += 1

a=total-no

percenthasinchi = a*100/total
writtenoutput.write(str(percenthasinchi)+ '% has inchi\n')







