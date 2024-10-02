import argparse as arg
from argparse import RawTextHelpFormatter

parser = arg.ArgumentParser(prog="sensiscript",
	formatter_class=RawTextHelpFormatter,
	description='Genomic antimicrobial susceptibility typing for Neisseria gonorrhoeae.\n',
	usage = '%(prog)s [options]')

parser = arg.ArgumentParser(description='sensiscript.py: genomic antimicrobial susceptibility typing for Neisseria gonorrhoeae', usage = '%(prog)s [options]')
parser.add_argument('-i', '--input_AMRtable', help='ARIBA output table containing genotypic AMR mechanisms', required=True)
parser.add_argument('-a', '--antibiotics', help='List and order of antibiotics to check separated by commas (options: ceftriaxone, ciprofloxacin, azithromycin, tetracycline, penicillin, spectinomycin, zoliflodacin) (default: ceftriaxone,ciprofloxacin,azithromycin,tetracycline,penicillin)', required=False, default='ceftriaxone,ciprofloxacin,azithromycin,tetracycline,penicillin')
parser.add_argument('-d', '--database', help='Path to the sensiscript.db file (default: sensiscript.db is in the same directory as the main script)', required=False, default='sensiscript.db')
parser.add_argument('-p', '--pena', help='Path to the sensiscript.penA.db file (default: sensiscript.penA.db is in the same directory as the main script)', required=False, default='sensiscript.penA.db')
parser.add_argument('-o', '--outfile', help='Print results to <outfile>. If not provided, results will be printed on screen.', required=False)
arg = parser.parse_args()

## Functions ##

def get_arguments(arg):
	args = {
		'intable': arg.input_AMRtable if arg.input_AMRtable else False,
		'antibiotics': arg.antibiotics.rstrip().split(',') if arg.antibiotics else ['ceftriaxone','ciprofloxacin','azithromycin','tetracycline','penicillin'],
		'database': arg.database if arg.database else 'sensiscript.db',
		'pena': arg.pena if arg.pena else 'sensiscript.penA.db',
		'outfile': arg.outfile if arg.outfile else False
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

def initialize_abx(antibiotics):
	recommended_treatment = {}
	found_mechanisms = {}
	selected_determinants = []
	for i in antibiotics:
		recommended_treatment[i] = []
		found_mechanisms[i] = []
		for x in abxdict[i]:
			if x not in selected_determinants:
				selected_determinants.append(x)
	return [recommended_treatment, found_mechanisms, selected_determinants]

def write_header(antibiotics, outfile):
	outheader = 'isolate'+'\t'+'treatment recommendation'+'\t'+'\t'.join(antibiotics)
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
	return [extract_col, sel_columns]

def cro_zol_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic):
	rec_treat = True # Only DO NOT recommend when: penA.A501P or any mutation in penA.T483 or penA.A311
	mech = []
	target_sites = abxdict[antibiotic]
	for site in target_sites:
		for det in amrdict[site]:
			if mutnamedb[det] in line_results:
				if line_results[mutnamedb[det]] == 'yes':
					rec_treat = False # if ANY of the three
					mech.append(mutnamedb[det])
	if antibiotic == 'ceftriaxone' and 'penA.mosaic' in line_results: #just reported for extra information but does not exclude ceftriaxone
		if line_results[mutnamedb['penA.ref_seq']] in penA_mosaic_vec:
			mosaic_nb = line_results[mutnamedb['penA.ref_seq']]
			mech.append(mosaic_nb)
	return [rec_treat, mech]

def cip_tet_spt_treatment(line_results, amrdict, amrdictr, abxdict, mutnamedb, antibiotic):
	rec_treat = True
	mech = []
	target_sites = abxdict[antibiotic]
	for site in target_sites:
		for det in amrdict[site]:
			if mutnamedb[det] in line_results:
				if line_results[mutnamedb[det]] == 'yes':
					mech.append(mutnamedb[det])
					rec_treat = False
	return [rec_treat, mech]

def azm_treatment(line_results, amrdict, amrdictr, abxdict, mutnamedb):
	rec_treat = True
	mech = []
	twentysmutation = False
	mtrDmosaic = False
	mtrCdisrupted = False
	target_sites = abxdict['azithromycin']
	for site in target_sites:
		for det in amrdict[site]:
			if mutnamedb[det] in line_results:
				if det == '23S.23S.2045G': 
					if line_results[mutnamedb[det]] == 'yes':
						twentysmutation = True
						mech.append(mutnamedb[det])
				if det == '23S.23S.2597T': 
					if line_results[mutnamedb[det]] == 'yes':
						twentysmutation = True
						mech.append(mutnamedb[det])
				if det == 'mtrD.ref_seq':
					if 'mosaic' in line_results[mutnamedb[det]]:
						mtrDmosaic = True
						mech.append(mutnamedb[det])
				if det == 'mtrC.assembled':
					if 'interrupted' in line_results[mutnamedb[det]]: # mtrC reverses susceptibility for the mtrD mosaic only
						mtrCdisrupted = True
						mech.append(mutnamedb[det])
		if twentysmutation:
			rec_treat = False
		if mtrDmosaic:
			rec_treat = False
			if mtrCdisrupted:
				rec_treat = True
	return [rec_treat, mech]

def pen_treatment(line_results, amrdict, amrdictr, penA_mosaic_vec, mutnamedb):
	rec_treat = True
	mech = []
	for x in abxdict['penicillin'][:-1]:
		for det in amrdict[x]:
			if mutnamedb[det] in line_results:
				if line_results[mutnamedb[det]] == 'yes':
					if mutnamedb[det] not in mech:
						mech.append(mutnamedb[det])
					rec_treat = False
	if 'penA.mosaic' in line_results:
		if line_results[mutnamedb['penA.ref_seq']] in penA_mosaic_vec:
			mosaic_nb = line_results[mutnamedb['penA.ref_seq']]
			mech.append(mosaic_nb)
			rec_treat = False
	return [rec_treat, mech]

def check_treatment(antibiotics, amrdict, recommended_treatment, found_mechanisms, line_results, mutnamedb):
	for i in antibiotics:	
		if i == 'ceftriaxone':
			cro_check = cro_zol_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic="ceftriaxone")
			recommended_treatment[i] = cro_check[0]
			found_mechanisms[i] = cro_check[1]
		elif i == 'ciprofloxacin':
			cip_check = cip_tet_spt_treatment(line_results, amrdict, amrdictr, abxdict, mutnamedb, antibiotic="ciprofloxacin")
			recommended_treatment[i] = cip_check[0]
			found_mechanisms[i] = cip_check[1]
		elif i == 'azithromycin':
			azm_check = azm_treatment(line_results, amrdict, amrdictr, abxdict, mutnamedb)
			recommended_treatment[i] = azm_check[0]
			found_mechanisms[i] = azm_check[1]
		elif i == 'tetracycline':
			tet_check = cip_tet_spt_treatment(line_results, amrdict, amrdictr, abxdict, mutnamedb, antibiotic="tetracycline")
			recommended_treatment[i] = tet_check[0]
			found_mechanisms[i] = tet_check[1]
		elif i == 'penicillin':
			pen_check = pen_treatment(line_results, amrdict, amrdictr, penA_mosaic_vec, mutnamedb)
			recommended_treatment[i] = pen_check[0]
			found_mechanisms[i] = pen_check[1]
		elif i == 'spectinomycin':
			spt_check = cip_tet_spt_treatment(line_results, amrdict, amrdictr, abxdict, mutnamedb, antibiotic="spectinomycin")
			recommended_treatment[i] = spt_check[0]
			found_mechanisms[i] = spt_check[1]
		elif i == 'zoliflodacin':
			zol_check = cro_zol_treatment(line_results, amrdict, abxdict, mutnamedb, antibiotic="zoliflodacin")
			recommended_treatment[i] = zol_check[0]
			found_mechanisms[i] = zol_check[1]
	return [recommended_treatment, found_mechanisms]

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
	recommended_treatment, found_mechanisms, selected_determinants = initialize_abx(antibiotics)

	# Write output header to screen or outfile #
	outfilehandle = write_header(antibiotics, outfile)

	# Check AMR mechanisms and predict treatment for each isolate #
	final_results = []
	with open(intable, 'r') as amrtable: #'ARIBA_summary_sensityping2.csv'
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
					else:
						line_results[item] = extract_col[count]
			process_antibiotics = check_treatment(antibiotics, amrdict, recommended_treatment, found_mechanisms, line_results, mutnamedb)
			recommended_treatment = process_antibiotics[0]
			found_mechanisms = process_antibiotics[1]
			# create treatment prediction line (results)
			treatment_prediction = []
			mechanisms = []
			for i in recommended_treatment:
				if recommended_treatment[i]:
					treatment_prediction.append(i)
				mechanisms.append('/'.join(found_mechanisms[i]))
			treatment_prediction_line = ','.join(treatment_prediction)+'\t'+'\t'.join(mechanisms)
			if len(treatment_prediction)<2:
				treatment_prediction_line = '(UND) '+treatment_prediction_line
			if outfile:
				outfilehandle.write(isolate+'\t'+treatment_prediction_line+'\n')
			else:
				print(isolate+'\t'+treatment_prediction_line)

