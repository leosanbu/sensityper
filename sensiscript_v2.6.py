import os
import argparse as arg
from argparse import RawTextHelpFormatter

parser = arg.ArgumentParser(prog="sensiscript",
	formatter_class=RawTextHelpFormatter,
	description='Genomic antimicrobial susceptibility typing for Neisseria gonorrhoeae.\n',
	usage = '%(prog)s [options]')

parser = arg.ArgumentParser(description='sensiscript.py: genomic antimicrobial susceptibility typing for Neisseria gonorrhoeae', usage = '%(prog)s [options]')
parser.add_argument('-i', '--input_AMRtable', help='ARIBA output table containing genotypic AMR mechanisms', required=True)
parser.add_argument('-a', '--antibiotics', help='List and order of antibiotics to check separated by commas (options: ceftriaxone, ciprofloxacin, azithromycin, tetracycline, penicillin, spectinomycin, zoliflodacin) (default: ceftriaxone,azithromycin,ciprofloxacin,tetracycline,penicillin,spectinomycin,zoliflodacin)', required=False, default='ceftriaxone,azithromycin,ciprofloxacin,tetracycline,penicillin,spectinomycin,zoliflodacin')
parser.add_argument('-d', '--database', help='Path to the sensiscript.db file (default: sensiscript.db is in the same directory as the main script)', required=False, default='sensitype.db')
parser.add_argument('-p', '--pena', help='Path to the sensiscript.penA.db file (default: sensiscript.penA.db is in the same directory as the main script)', required=False, default='sensitype.penA.db')
parser.add_argument('-o', '--outfile', help='Print results to <outfile>. If not provided, results will be printed on screen.', required=False)
parser.add_argument('--suppress-html', help='Suppress automatic HTML output generation (used internally by pipeline mode)', action='store_true', required=False, default=False)
arg = parser.parse_args()

## Functions ##

def get_arguments(arg):
	args = {
		'intable': arg.input_AMRtable if arg.input_AMRtable else False,
		'antibiotics': arg.antibiotics.rstrip().split(',') if arg.antibiotics else ['ceftriaxone','azithromycin','ciprofloxacin','tetracycline','penicillin', 'spectinomycin', 'zoliflodacin'],
		'database': arg.database if arg.database else 'sensitype.db',
		'pena': arg.pena if arg.pena else 'sensitype.penA.db',
		'outfile': arg.outfile if arg.outfile else False,
		'suppress_html': arg.suppress_html if hasattr(arg, 'suppress_html') else False
	}
	return args

def read_databases(database, pena):
	amrdict = {}
	amrdictr = {}
	abxdict = {}
	mutnamedb = {}
	with open (database, 'r') as db:
		for line in db:
			if not line.startswith('antibiotic'): #skip header line
				linesplit = line.rstrip().split('\t')
				mutnamedb[linesplit[2]] = linesplit[3]
				if linesplit[0] not in abxdict:
					abxdict[linesplit[0]] = [linesplit[1]]
				else:
					if linesplit[1] not in abxdict[linesplit[0]]:
						abxdict[linesplit[0]].append(linesplit[1])
				if linesplit[1] not in amrdict:
					amrdict[linesplit[1]] = [linesplit[2]]
				else:
					if linesplit[2] not in amrdict[linesplit[1]]:
						amrdict[linesplit[1]].append(linesplit[2])
				amrdictr[linesplit[2]] = [linesplit[1]]
	penA_mosaic_vec = []
	with open (pena, 'r') as pdb:
		for line in pdb:
			linesplit = line.rstrip().split('\t')
			if linesplit[2] == 'yes':
				penA_mosaic_vec.append('penA.'+linesplit[0])
	return [amrdict, amrdictr, abxdict, penA_mosaic_vec, mutnamedb]

def initialize_abx(antibiotics, abxdict):
	recommended_treatment = {}
	found_mechanisms = {}
	wildtype_alleles = {}
	selected_determinants = []
	for i in antibiotics:
		recommended_treatment[i] = []
		found_mechanisms[i] = []
		wildtype_alleles[i] = []
		for x in abxdict[i]:
			if x not in selected_determinants:
				selected_determinants.append(x)
	return [recommended_treatment, found_mechanisms, wildtype_alleles, selected_determinants]

def write_header(antibiotics, outfile):
	abxlist =[] 
	for x in antibiotics:
		ab = x+'_NWT\t'+x+'_WT'
		abxlist.append(ab)
	#antibioticsR = [x+'_R' for x in antibiotics]
	#antibioticsS = [x+'_S' for x in antibiotics]
	outheader = 'isolate'+'\t'+'treatment recommendation'+'\t'+'\t'.join(abxlist)
	if outfile:
		outfilehandle = open(outfile, 'w+')
		outfilehandle.write(outheader+'\n')
	else:
		outfilehandle = False
		print(outheader)
	return outfilehandle

def find_columns_ariba(sel, header, linesplit, amrdict):
	sel_columns = []
	for s in amrdict[sel]:
		if s in header:
			sel_columns.append(s)
	extract_col = []
	for s in sel_columns:
		col_i = header.index(s)
		col_result = linesplit[col_i]
		extract_col.append(col_result)
		if s == '23S.23S.2045G' or s == '23S.23S.2597T':
			if s+'.%' in header:
				extract_col = get_mutated_23Sreads(s, header, linesplit, extract_col)
				sel_columns.append(s+'.%')
	return [extract_col, sel_columns]

def get_mutated_23Sreads(sel, header, linesplit, extract_col):
	col_i_per = header.index(sel+'.%')
	col_result_per = linesplit[col_i_per]
	extract_col.append(col_result_per)
	return extract_col

def cro_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic):
	rec_treat = True # Only DO NOT recommend when: penA.A501P or penA.A311+V316
	mech = []
	wt = []
	target_sites = abxdict[antibiotic]
	checkA311 = 0
	checkV316 = 0
	for sites in target_sites:
		for det in line_results:
			if sites in det:
				if sites == 'penA.A501' and 'yes' in line_results[det]:
					rec_treat = False
					mech.append(det)
				else:
					if 'no' in line_results[det]:
						wt.append(sites+'_WT')
				if sites == 'penA.A311' and 'yes' in line_results[det]:
					checkA311 = 1
					mech.append(det)
				else:
					if 'no' in line_results[det]:
						wt.append(sites+'_WT')
				if sites == 'penA.V316' and 'yes' in line_results[det]:
					checkV316 = 1
					mech.append(det)
				else:
					if 'no' in line_results[det]:
						wt.append(sites+'_WT')
	if checkA311 == 1 and checkV316 == 1:
		rec_treat = False
	if 'penA.mosaic' in line_results: #just reported for extra information but does not exclude ceftriaxone
		if line_results[mutnamedb['penA.ref_seq']] in penA_mosaic_vec:
			mosaic_nb = line_results[mutnamedb['penA.ref_seq']]
			mech.append(mosaic_nb)
	filtered_wt = exclude_mutated_codon(wt, mech)
	return [rec_treat, mech, set(filtered_wt)]

def cip_tet_spt_zol_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic):
	rec_treat = True
	mech = []
	wt = []
	target_sites = abxdict[antibiotic]
	for site in target_sites:
		for det in amrdict[site]:
			if mutnamedb[det] in line_results:
				if 'yes' in line_results[mutnamedb[det]]:
					mech.append(mutnamedb[det])
					rec_treat = False
				else:
					if 'no' in line_results[mutnamedb[det]]:
						if '16S.' in site:
							wt.append(site.rstrip(site[-1])+'_WT')
						else:
							if 'tetM' in site:
								wt.append('tetM.not_present')
							else:
								wt.append(site+'_WT')
			else:
				if site in line_results: #e.g. when no particular mutations for gyrB.D429 or gyrB.K450
					wt.append(site+'_WT')
	novel_mech = include_novel_mutations(line_results, abxdict, antibiotic)
	if len(novel_mech)>0:
		for n in novel_mech:
			mech.append(n)
		if rec_treat:
			rec_treat = 'novel'
	filtered_wt = exclude_mutated_codon(wt, mech)
	return [rec_treat, mech, set(filtered_wt)]

def azm_treatment(line_results, amrdict, amrdictr, abxdict, mutnamedb):
	rec_treat = True
	mech = []
	wt = []
	twentysmutation = False
	mtrDmosaic = False
	mtrCdisrupted = False
	target_sites = abxdict['azithromycin']
	for site in target_sites:
		for det in amrdict[site]:
			if mutnamedb[det] in line_results:
				if det == '23S.23S.2045G': 
					if 'yes' in line_results[mutnamedb[det]] or 'het' in line_results[mutnamedb[det]]:
						twentysmutation = True
						if det+'.%' in line_results:
							twentyspercent = line_results[det+'.%']
						else:
							twentyspercent = '99.9'
						mech.append(mutnamedb[det]+'['+twentyspercent+'%]')
					elif 'no' in line_results[mutnamedb[det]]:
							wt.append(site.rstrip(site[-1])+'_WT')
				elif det == '23S.23S.2597T': 
					if 'yes' in line_results[mutnamedb[det]] or 'het' in line_results[mutnamedb[det]]:
						twentysmutation = True
						if det+'.%' in line_results:
							twentyspercent = line_results[det+'.%']
						else:
							twentyspercent = '99.9'
						mech.append(mutnamedb[det]+'['+twentyspercent+'%]')
					elif 'no' in line_results[mutnamedb[det]]:
							wt.append(site.rstrip(site[-1])+'_WT')
				elif det == 'mtrD.ref_seq':
					if 'mosaic' in line_results[mutnamedb[det]]:
						mtrDmosaic = True
						mech.append(line_results[mutnamedb[det]])
					else:
						wt.append('mtrD.WT')
				elif det == 'mtrC.assembled':
					if 'interrupted' in line_results[mutnamedb[det]]: # mtrC reverses susceptibility for the mtrD mosaic only
						mtrCdisrupted = True
						mech.append(mutnamedb[det])
					else:
						wt.append('mtrC.WT')
	if twentysmutation:
		rec_treat = False  # 23S mutation always causes resistance
	elif mtrDmosaic:
		if mtrCdisrupted:
			rec_treat = True  # mtrC disruption reverts mtrD mosaic resistance
		else:
			rec_treat = False  # mtrD mosaic causes resistance
	novel_mech = include_novel_mutations(line_results, abxdict, "azithromycin")
	if len(novel_mech)>0:
		for n in novel_mech:
			mech.append(n)
		if rec_treat:
			rec_treat = 'novel'
	return [rec_treat, mech, wt]

def pen_treatment(line_results, amrdict, amrdictr, abxdict, penA_mosaic_vec, mutnamedb):
	rec_treat = True
	mech = []
	wt = []
	for x in abxdict['penicillin'][:-1]:
		for det in amrdict[x]:
			if mutnamedb[det] in line_results:
				if 'yes' in line_results[mutnamedb[det]]:
					if mutnamedb[det] not in mech:
						mech.append(mutnamedb[det])
					rec_treat = False
				else:
					if 'no' in line_results[mutnamedb[det]]:
						if 'bla' in mutnamedb[det]:
							wt.append('blaTEM.not_present')
						else:
							wt.append(x+'_WT')
	if 'penA.mosaic' in line_results:
		if line_results[mutnamedb['penA.ref_seq']] in penA_mosaic_vec:
			mosaic_nb = line_results[mutnamedb['penA.ref_seq']]
			mech.append(mosaic_nb)
			rec_treat = False
	novel_mech = include_novel_mutations(line_results, abxdict, "penicillin")
	if len(novel_mech)>0:
		for n in novel_mech:
			mech.append(n)
		if rec_treat:
			rec_treat = 'novel'
	filtered_wt = exclude_mutated_codon(wt, mech) # to exclude absent blaTEM from WT list
	return [rec_treat, mech, set(filtered_wt)]

def check_treatment(antibiotics, amrdict, recommended_treatment, found_mechanisms, wildtype_alleles, line_results, mutnamedb):
	for i in antibiotics:	
		if i == 'ceftriaxone':
			cro_check = cro_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic="ceftriaxone")
			recommended_treatment[i] = cro_check[0]
			found_mechanisms[i] = cro_check[1]
			wildtype_alleles[i] = cro_check[2]
		elif i == 'ciprofloxacin':
			cip_check = cip_tet_spt_zol_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic="ciprofloxacin")
			recommended_treatment[i] = cip_check[0]
			found_mechanisms[i] = cip_check[1]
			wildtype_alleles[i] = cip_check[2]
		elif i == 'azithromycin':
			azm_check = azm_treatment(line_results, amrdict, amrdictr, abxdict, mutnamedb)
			recommended_treatment[i] = azm_check[0]
			found_mechanisms[i] = azm_check[1]
			wildtype_alleles[i] = azm_check[2]
		elif i == 'tetracycline':
			tet_check = cip_tet_spt_zol_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic="tetracycline")
			recommended_treatment[i] = tet_check[0]
			found_mechanisms[i] = tet_check[1]
			wildtype_alleles[i] = tet_check[2]
		elif i == 'penicillin':
			pen_check = pen_treatment(line_results, amrdict, amrdictr, abxdict, penA_mosaic_vec, mutnamedb)
			recommended_treatment[i] = pen_check[0]
			found_mechanisms[i] = pen_check[1]
			wildtype_alleles[i] = pen_check[2]
		elif i == 'spectinomycin':
			spt_check = cip_tet_spt_zol_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic="spectinomycin")
			recommended_treatment[i] = spt_check[0]
			found_mechanisms[i] = spt_check[1]
			wildtype_alleles[i] = spt_check[2]
		elif i == 'zoliflodacin':
			zol_check = cip_tet_spt_zol_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic="zoliflodacin")
			recommended_treatment[i] = zol_check[0]
			found_mechanisms[i] = zol_check[1]
			wildtype_alleles[i] = zol_check[2]
	return [recommended_treatment, found_mechanisms, wildtype_alleles]

def call_wildtype(isolate, line_results, selected_determinants):
	with open(isolate, 'r') as report:
		headrep = report.readline().rstrip().split('\t')
		for r in report:
			rsplit = r.rstrip().split('\t')
			gene = rsplit[6]
			known_var = rsplit[13]
			known_var_change = rsplit[16]
			has_known_var = rsplit[17]
			ref_allele = rsplit[22]
			ctg_allele = rsplit[25]
			smtls_nt = rsplit[27]
			smtls_cov = rsplit[28]
			description = rsplit[29]
			mut_name = gene+'.'+known_var_change
			if known_var_change != '.':
				if not 'ins' in known_var_change: 
					codon_in_description = description.rsplit(':')[4].rsplit('.')[1] #get codon as in database, e.g. penA.502 must be penA.501
					codon_in_description = codon_in_description[:-1]
					for s in selected_determinants:
						if gene in s and codon_in_description in s:
							if known_var == '1':
								if has_known_var == '0': # check if known mutation is NOT present
									if len(ctg_allele)==3: # codon
										mut_name_cut = s+translate_codon(ctg_allele)
									else: # nucleotide
										mut_name_cut = mut_name
									if mut_name_cut in line_results:
										if ref_allele == ctg_allele:
											if line_results[mut_name_cut] != 'yes':
												if not '_WT' in line_results[mut_name_cut]:
													line_results[mut_name_cut] = line_results[mut_name_cut]+'_WT' # include base/codon in contig
									else:
										if ref_allele == ctg_allele:
											line_results[s] = 'no_WT'
			else: ## novel alleles
				ref_ctg_change = rsplit[18]
				if ref_ctg_change != '.': # discard mtrD/mtrC lines, which are not point changes
					novel_mut = gene+'.'+ref_ctg_change
					novel_mut_cut = gene+'.'+ref_ctg_change[:-1]
					for s in selected_determinants:
						if novel_mut_cut in s: #unknown change in known codon
							if gene == '23S' or gene == '16S':
								novel_allele_name = novel_mut_cut+'_'+ctg_allele
							else:
								aa = translate_codon(ctg_allele)
								novel_allele_name = novel_mut_cut+'_'+aa
							line_results[novel_allele_name] = 'novel'
							if gene == '23S': # check number of mutated copies
								smtls_split = smtls_nt.split(',')
								cov_split = smtls_cov.split(',')
								if len(smtls_split) == 1:
									if smtls_split[0] == ctg_allele:
										line_results[novel_allele_name] = line_results[novel_allele_name]+'[100.0%]'
								elif len(smtls_split) == 2:
									total = int(cov_split[0])+int(cov_split[1])
									nt_index = smtls_split.index(ctg_allele)
									percent = round(int(cov_split[nt_index])*100/total,1)
									line_results[novel_allele_name] = line_results[novel_allele_name]+'['+str(percent)+'%]'
	if 'penA.insD345' in line_results: # if penA.insD345 NOT present will NOT be in the report, so add here
		if line_results['penA.insD345'] == 'NA':
			line_results['penA.insD345'] = 'no_WT'
	else:
		line_results['penA.insD345'] = 'no_WT'
	return line_results

def exclude_mutated_codon(wt, mech):
	exclude_wt = []
	for w in wt:
		for m in mech:
			w2 = w.replace('_WT', '')
			if w2 in m: # mutated codon, exclude from WT list
				exclude_wt.append(w)
	filtered_wt = []
	for w in wt:
		if w not in exclude_wt:
			filtered_wt.append(w)
	return filtered_wt

def include_novel_mutations(line_results, abxdict, antibiotic): 
	include_novel = []
	for i in line_results:
		if line_results[i] == 'novel': #check whether the novel mutation is in a known gene of the target antibiotic
			tmpmut = i.split('_')[0]
			checkmut = False
			for x in abxdict[antibiotic]:
				if tmpmut in x:
					checkmut = True
					include_novel.append(i)
		elif 'novel[' in line_results[i]: #check whether the novel mutation is in a known gene of the target antibiotic
			tmpmut = i.split('_')[0]
			checkmut = False
			for x in abxdict[antibiotic]:
				if tmpmut in x:
					checkmut = True
			include_novel.append(i+line_results[i].replace('novel', ''))
	return include_novel

def translate_codon(codon):
	table = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
		'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
	}
	aa = table[codon]
	return aa

def create_treatment_line(isolate, recommended_treatment, found_mechanisms, wildtype_alleles, antibiotics):
	treatment_prediction = []
	mechanisms = []
	#wt_mechanisms = []
	for i in antibiotics:
		if recommended_treatment[i] == True:
			treatment_prediction.append(i)
		elif recommended_treatment[i] == 'novel':
			treatment_prediction.append(i+'(WARN:novel_mutation)')
		mechanisms.append('/'.join(sorted(found_mechanisms[i])))
		mechanisms.append('/'.join(sorted(wildtype_alleles[i])))
		#wt_mechanisms.append('/'.join(sorted(wildtype_alleles[i])))
		treatment_prediction_line = ','.join(treatment_prediction)+'\t'+'\t'.join(mechanisms)#+'\t'+'\t'.join(wt_mechanisms)
		if len(treatment_prediction)<2:
			treatment_prediction_line = '(UND) '+treatment_prediction_line
	short_isolate = os.path.basename(isolate.replace('_ARIBA/report_complete.tsv', ''))
	if outfile:
		outfilehandle.write(short_isolate+'\t'+treatment_prediction_line+'\n')
	else: 
		print(short_isolate+'\t'+treatment_prediction_line)

def checkAMR_and_predict(intable, amrdict, amrdictr, abxdict, penA_mosaic_vec, mutnamedb, recommended_treatment, found_mechanisms, wildtype_alleles, selected_determinants):
	with open(intable, 'r') as amrtable: #'ariba_summary.csv'
		header = amrtable.readline().rstrip().split(',')
		for line in amrtable:
			linesplit = line.rstrip().split(',')
			isolate = linesplit[0]
			line_results = {}
			for sel in selected_determinants:
				extract_col, sel_columns = find_columns_ariba(sel, header, linesplit, amrdict)
				for count, item in enumerate(sel_columns):
					if item in mutnamedb:
						line_results[mutnamedb[item]] = extract_col[count]
						if '.%' in item:
							line_results[item] = extract_col[1]
			line_results2 = call_wildtype(isolate, line_results, selected_determinants) # Check wildtypes
			process_antibiotics = check_treatment(antibiotics, amrdict, recommended_treatment, found_mechanisms, wildtype_alleles, line_results2, mutnamedb) # Check treatment
			recommended_treatment = process_antibiotics[0]
			found_mechanisms = process_antibiotics[1]
			wildtype_alleles = process_antibiotics[2]
			create_treatment_line(isolate, recommended_treatment, found_mechanisms, wildtype_alleles, antibiotics) # create treatment prediction line (results)


##########
## Main ##
##########


if __name__ == '__main__':

	# Get arguments #
	args = get_arguments(arg)
	intable = args['intable']
	antibiotics = args['antibiotics']
	database = args['database']
	pena = args['pena']
	outfile = args['outfile']

	# Read database of abx and amr mechanisms #
	amrdict, amrdictr, abxdict, penA_mosaic_vec, mutnamedb = read_databases(database, pena)

	# Initialize output directories in the antibiotic order specified by the user #
	recommended_treatment, found_mechanisms, wildtype_alleles, selected_determinants = initialize_abx(antibiotics, abxdict)

	# Write output header to screen or outfile #
	outfilehandle = write_header(antibiotics, outfile)

	# Check AMR mechanisms and predict treatment for each isolate #
	checkAMR_and_predict(intable, amrdict, amrdictr, abxdict, penA_mosaic_vec, mutnamedb, recommended_treatment, found_mechanisms, wildtype_alleles, selected_determinants)

	# Close output file so all buffered data is flushed to disk before HTML reads it #
	if outfilehandle:
		outfilehandle.close()

	# Generate HTML output (unless suppressed) #
	suppress_html = args.get('suppress_html', False)
	if outfile and not suppress_html:
		try:
			from html_generator import generate_resistance_profile_html
			html_path = outfile.replace('.tsv', '.html')
			generate_resistance_profile_html(
				tsv_path=outfile,
				output_html_path=html_path,
				antibiotics=antibiotics
			)
		except Exception as e:
			print(f"Warning: Could not generate HTML output: {e}")
