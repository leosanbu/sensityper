#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sensitype (CLI orchestrator) — v0.6.7-py36

What this does
--------------
- Runs ARIBA in batch over input directories (with a retry loop per sample)
- Runs the Sensitype rules script (formerly sensiscript_v2.5.py) to generate per-antibiotic predictions
- Runs Sensitreat (this file) to select recommended_1, recommended_2 (no "second-line" wording)
- Provides a 'pipeline' mode to chain modules
- Regimen-category flags are written to treatment_output.tsv.
  Exactly one is "YES" (others "NO"). Priority is controlled by --sensitreat_order.

Breaking changes vs earlier internal versions
---------------------------------------------
- Command names:   sensityping -> sensitype,   treatment_prediction -> sensitreat
- Output headers:  1st_line_treatment / 2nd_line_treatment -> recommended_1 / recommended_2
- DB file names:   sensiscript.db -> sensitype.db,   sensiscript.penA.db -> sensitype.penA.db
- v0.6.6:          Path resolution precedence added (CLI > ENV > script-dir fallback),
                   with runtime config printouts; this build is Python 3.6+ compatible.
- v0.6.7:          Sensitreat comment strings updated to "RECOMMENDATION" labels; added
                   azithromycin monotherapy regimen category (only selectable if explicitly
                   included in --sensitreat_order). Antibiotic abbreviations removed from CLI
                   help/examples and output headers now use full regimen labels.

Important usage note (v0.6.7)
-----------------------------
- Regimen categories in --sensitreat_order MUST be written using full names:
    - ceftriaxone+azithromycin
    - ceftriaxone
    - azithromycin+spectinomycin
    - azithromycin
    - ciprofloxacin
    - spectinomycin
    - zoliflodacin
- Abbreviations (e.g., CRO, AZM, CIP, SPC, ZOL) are no longer accepted in sensitreat_order.
"""

import os
import subprocess
import argparse
import csv
import sys
from pathlib import Path
from typing import Optional  # 3.6+ compatible unions

__version__ = "0.6.7-py36"

# -------------------------------------------------------------------
# Path resolution helpers: CLI > ENV > script-dir fallback
# -------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent

def _env(name: str, default=None):
    v = os.environ.get(name)
    return v if v not in (None, "", "None") else default

def resolve_path(cli_value: Optional[str],
                 env_var: str,
                 default: Optional[Path],
                 must_exist: bool = True,
                 is_dir: bool = False) -> Optional[str]:
    """
    Resolution order:
      1) CLI value
      2) ENV variable (env_var)
      3) default Path (usually under SCRIPT_DIR)
    Returns a string path (or None if must_exist=False and nothing provided).
    """
    candidate = None  # type: Optional[Path]

    if cli_value:
        candidate = Path(cli_value)
    else:
        env_val = _env(env_var)
        if env_val:
            candidate = Path(env_val)
        elif default is not None:
            candidate = default

    if candidate is None:
        if must_exist:
            raise FileNotFoundError(
                "Path for {env} not provided. Pass a CLI flag or set {env}, or place the default next to the script."
                .format(env=env_var)
            )
        return None

    if must_exist:
        if is_dir and not candidate.is_dir():
            raise FileNotFoundError("{env}: directory not found: {p}".format(env=env_var, p=candidate))
        if not is_dir and not candidate.exists():
            raise FileNotFoundError("{env}: file not found: {p}".format(env=env_var, p=candidate))

    return str(candidate)

def config_print(label: str, value: Optional[str]):
    print("[CONFIG] {label:<22} = {val}".format(
        label=label,
        val=value if value is not None else "(not set)"
    ))

# -------------------------------------------------------------------
# Default fallbacks (script-dir)
# -------------------------------------------------------------------
# We keep helper scripts next to this file; databases in resources/
DEFAULT_ARIBA_BATCH = SCRIPT_DIR / "ariba_batch_v0.2.py"
DEFAULT_RULES_PATH  = SCRIPT_DIR / "sensiscript_v2.6.py"
DEFAULT_DB_MAIN     = SCRIPT_DIR / "resources" / "sensitype.db"
DEFAULT_DB_PENA     = SCRIPT_DIR / "resources" / "sensitype.penA.db"
DEFAULT_ARIBA_DB    = SCRIPT_DIR / "resources" / "ariba_db"   # put your bundled ARIBA DB here

# Will be set at runtime (ENV fallback or defaults above)
ARIBA_BATCH_PATH = None  # type: Optional[str]
SENSITYPE_RULES_PATH = None  # type: Optional[str]

# -------------------------------------------------------------------
# ARIBA batch + summary
# -------------------------------------------------------------------
def run_ariba_batch(input_dirs: str, output_dir: str, db_path: str, threads: int) -> None:
    """
    Calls external ariba_batch.py to run ARIBA across dirs.
    Then retries sample(s) missing report_complete.tsv up to 3 times with direct 'ariba run'.
    Finally writes filenames.txt and runs 'ariba summary' -> ariba_summary.*.
    """
    # 1) run ariba_batch.py
    cmd = [
        sys.executable, ARIBA_BATCH_PATH,
        '-d', input_dirs,
        '-o', output_dir,
        '--db_path', db_path,
        '-t', str(threads)
    ]
    print("[ARIBA batch] {c}".format(c=" ".join(cmd)))
    subprocess.run(cmd, check=True)

    print("\nChecking for failed samples to retry...\n")
    forward_suffixes = ['_1.fastq.gz', '_R1.fastq.gz', '_R1_001.fastq.gz']
    reverse_suffixes = ['_2.fastq.gz', '_R2.fastq.gz', '_R2_001.fastq.gz']
    dirs = input_dirs.split(',')

    for d in dirs:
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if not any(f.endswith(suf) for suf in forward_suffixes):
                continue

            # detect pair suffix
            sample_name = None
            rev_suffix = None
            for suf in forward_suffixes:
                if f.endswith(suf):
                    sample_name = f.rsplit(suf, 1)[0]
                    rev_suffix = reverse_suffixes[forward_suffixes.index(suf)]
                    break
            if sample_name is None:
                continue

            forward = os.path.join(d, f)
            reverse = os.path.join(d, sample_name + rev_suffix)
            outdir = os.path.join(output_dir, sample_name + '_ARIBA')
            report_complete = os.path.join(outdir, 'report_complete.tsv')

            if os.path.exists(report_complete):
                continue

            print("\nreport_complete.tsv missing for {s}. Retrying...".format(s=sample_name))

            for attempt in range(1, 4):
                print("Attempt {a}/3 for {s}".format(a=attempt, s=sample_name))
                if os.path.exists(outdir):
                    subprocess.run(['rm', '-rf', outdir], check=False)
                cmd = ['ariba', 'run', '--threads', str(threads), db_path, forward, reverse, outdir]
                subprocess.run(cmd, check=False)

                report_file = os.path.join(outdir, 'report.tsv')
                if os.path.exists(report_file):
                    with open(report_file, 'r') as infile, open(report_complete, 'w') as outfile:
                        for line in infile:
                            if 'D147_T148insT' in line:
                                pattern = r"0\t\.\tp\t\.\t0\tD147_T148insT"
                                replacement = "1\tSNP\tp\tD147_T148insT\t1\tD147_T148insT"
                                modified_line = re.sub(pattern, replacement, line)
                            elif 'R146_D147insR' in line:
                                pattern = r"0\t\.\tp\t\.\t0\tR146_D147insR"
                                replacement = "1\tSNP\tp\tR146_D147insR\t1\tR146_D147insR"
                                modified_line = re.sub(pattern, replacement, line)
                            else:
                                modified_line = line
                            outfile.write(modified_line)

                if os.path.exists(report_complete):
                    print("Success: {s}".format(s=sample_name))
                    break
                else:
                    print("Still failed: {s}".format(s=sample_name))
            else:
                print("Gave up after 3 failed attempts: {s}".format(s=sample_name))

    # 2) Generate filenames.txt for the ARIBA summary including absolute paths
    # Absolute paths are needed for sensitype to access the individual report_complete.tsv files in order to check wildtypes
    with open('filenames.txt', 'w') as f:
        report_files = [os.path.join(os.path.abspath(output_dir), d, 'report_complete.tsv') for d in os.listdir(output_dir) if d.endswith('_ARIBA')]
        for count, report_file in enumerate(report_files):
            f.write(report_file+'\n')

    # 3) ariba summary
    print("\nRunning ARIBA summary...")
    subprocess.run([
        'ariba', 'summary', 'ariba_summary',
        '-f', 'filenames.txt',
        '--cluster_cols', 'assembled,ref_seq,pct_id',
        '--col_filter', 'n', '--row_filter', 'n',
        '--no_tree', '--v_groups', '--known_variants'
    ], check=True)

def run_ariba(parsed_args_or_dict) -> None:
    """
    Supports both direct subcommand 'ariba' and pipeline mode (dict).
    """
    # Resolve helper scripts first (ENV or fallback to script dir)
    global ARIBA_BATCH_PATH, SENSITYPE_RULES_PATH
    ARIBA_BATCH_PATH = _env("SENSITYPE_ARIBA_BATCH") or str(DEFAULT_ARIBA_BATCH)
    SENSITYPE_RULES_PATH = _env("SENSITYPE_RULES_PATH") or str(DEFAULT_RULES_PATH)

    if isinstance(parsed_args_or_dict, dict):  # pipeline mode
        input_dirs = parsed_args_or_dict['--input_dirs']
        output_dir = parsed_args_or_dict['--output_dir']
        db_path = resolve_path(parsed_args_or_dict.get('--db_path'),
                               env_var="SENSITYPE_ARIBA_DB",
                               default=DEFAULT_ARIBA_DB,
                               must_exist=True, is_dir=True)
        threads = int(parsed_args_or_dict.get('--threads', 1) or 1)
    else:  # direct subcommand
        input_dirs = parsed_args_or_dict.input_dirs
        output_dir = parsed_args_or_dict.output_dir
        db_path = resolve_path(parsed_args_or_dict.db_path,
                               env_var="SENSITYPE_ARIBA_DB",
                               default=DEFAULT_ARIBA_DB,
                               must_exist=True, is_dir=True)
        threads = parsed_args_or_dict.threads

    # Config printout
    print("\n=== Effective configuration (ARIBA) ===")
    config_print("ARIBA_BATCH_PATH", ARIBA_BATCH_PATH)
    config_print("SENSITYPE_RULES_PATH", SENSITYPE_RULES_PATH)
    config_print("ARIBA_DB", db_path)
    print("======================================\n")

    # Ensure output dir exists
    create_output_dir(output_dir)
    run_ariba_batch(input_dirs, output_dir, db_path, threads)

# -------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------
def create_output_dir(output_dir: str) -> None:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("Created output directory: {d}".format(d=output_dir))
    else:
        print("Output directory already exists: {d}".format(d=output_dir))

def rename_files(directory: str, preview: bool = False) -> None:
    """
    Standardize FASTQ names to <sample>_R1.fastq.gz / <sample>_R2.fastq.gz.
    Recognizes common suffix variants.
    """
    fwd_suffixes = ['_1.fastq.gz', '_R1.fastq.gz', '_R1_001.fastq.gz']
    rev_suffixes = ['_2.fastq.gz', '_R2.fastq.gz', '_R2_001.fastq.gz']

    for fname in sorted(os.listdir(directory)):
        full = os.path.join(directory, fname)
        if not os.path.isfile(full):
            continue

        for i, fwd in enumerate(fwd_suffixes):
            if fname.endswith(fwd):
                sample = fname[:-len(fwd)]
                rev_guess = sample + rev_suffixes[i]
                fwd_new = sample + "_R1.fastq.gz"
                rev_new = sample + "_R2.fastq.gz"

                src_fwd = full
                src_rev = os.path.join(directory, rev_guess)
                dst_fwd = os.path.join(directory, fwd_new)
                dst_rev = os.path.join(directory, rev_new)

                if not os.path.exists(src_rev):
                    # reverse may already be normalized; try alternatives
                    for ralt in rev_suffixes:
                        alt = os.path.join(directory, sample + ralt)
                        if os.path.exists(alt):
                            src_rev = alt
                            break

                actions = []
                if os.path.exists(src_fwd) and src_fwd != dst_fwd:
                    actions.append(("mv", src_fwd, dst_fwd))
                if os.path.exists(src_rev) and src_rev != dst_rev:
                    actions.append(("mv", src_rev, dst_rev))

                if actions:
                    if preview:
                        for a in actions:
                            print("[PREVIEW] {op} {src} -> {dst}".format(op=a[0], src=a[1], dst=a[2]))
                    else:
                        for a in actions:
                            os.rename(a[1], a[2])
                            print("Renamed {src} -> {dst}".format(src=a[1], dst=a[2]))
                break

# -------------------------------------------------------------------
# Run Sensitype rules (external script)
# -------------------------------------------------------------------
def run_sensiscript(input_AMRtable: str,
                    database: str,
                    pena: str,
                    antibiotics: list,
                    outfile: str,
                    suppress_html: bool = False) -> None:
    """
    Calls the rules script (default: sensiscript_v2.5.py) to generate per-antibiotic predictions.
    DB names default to sensitype.db / sensitype.penA.db per v0.6.0 migration.
    """
    cmd = [
        sys.executable, SENSITYPE_RULES_PATH,
        '-i', input_AMRtable,
        '-d', database,
        '-p', pena,
        '-a', ','.join(antibiotics),
        '-o', outfile
    ]
    if suppress_html:
        cmd.append('--suppress-html')
    print("[Sensitype rules] {c}".format(c=" ".join(cmd)))
    subprocess.run(cmd, check=True)

def run_sensitype(parsed_args_or_dict) -> None:
    """
    New name for the 'sensityping' module.
    """
    # Resolve helper scripts first (ENV or fallback to script dir)
    global ARIBA_BATCH_PATH, SENSITYPE_RULES_PATH
    ARIBA_BATCH_PATH = _env("SENSITYPE_ARIBA_BATCH") or str(DEFAULT_ARIBA_BATCH)
    SENSITYPE_RULES_PATH = _env("SENSITYPE_RULES_PATH") or str(DEFAULT_RULES_PATH)

    if isinstance(parsed_args_or_dict, dict):  # pipeline mode
        is_pipeline_mode = True
        input_AMRtable = parsed_args_or_dict.get('--input_AMRtable', os.path.abspath('ariba_summary.csv'))

        sensiscript_db = resolve_path(parsed_args_or_dict.get('--sensiscript_db'),
                                      env_var="SENSITYPE_DB",
                                      default=DEFAULT_DB_MAIN,
                                      must_exist=True, is_dir=False)

        sensiscript_pena = resolve_path(parsed_args_or_dict.get('--sensiscript_pena'),
                                        env_var="SENSITYPE_PENA_DB",
                                        default=DEFAULT_DB_PENA,
                                        must_exist=True, is_dir=False)

        antibiotics = parsed_args_or_dict.get(
            '--sensiscript_antibiotics',
            'ceftriaxone,azithromycin,ciprofloxacin,tetracycline,penicillin,zoliflodacin'
        ).split(',')
        sensiscript_outfile = parsed_args_or_dict['--sensiscript_outfile']
    else:
        is_pipeline_mode = False
        input_AMRtable = parsed_args_or_dict.input_AMRtable

        sensiscript_db = resolve_path(parsed_args_or_dict.sensiscript_db,
                                      env_var="SENSITYPE_DB",
                                      default=DEFAULT_DB_MAIN,
                                      must_exist=True, is_dir=False)

        sensiscript_pena = resolve_path(parsed_args_or_dict.sensiscript_pena,
                                        env_var="SENSITYPE_PENA_DB",
                                        default=DEFAULT_DB_PENA,
                                        must_exist=True, is_dir=False)

        antibiotics = parsed_args_or_dict.sensiscript_antibiotics.split(',')
        sensiscript_outfile = parsed_args_or_dict.sensiscript_outfile

    # Config printout
    print("\n=== Effective configuration (Sensitype rules) ===")
    config_print("SENSITYPE_RULES_PATH", SENSITYPE_RULES_PATH)
    config_print("sensitype.db", sensiscript_db)
    config_print("sensitype.penA.db", sensiscript_pena)
    print("===============================================\n")

    # Suppress HTML in pipeline mode to avoid duplicate HTML files
    run_sensiscript(input_AMRtable, sensiscript_db, sensiscript_pena, antibiotics, sensiscript_outfile,
                   suppress_html=is_pipeline_mode)

# -------------------------------------------------------------------
# Sensitreat (recommendations)
# -------------------------------------------------------------------
def process_output(input_file: str,
                   available_antibiotics: list,
                   sensitreat_order: list,
                   alert_output: str,
                   treatment_output: str) -> None:

    def get_pred(antibiotic, recommended_antibiotics):#row, *keys):
        return_pred = 'yes'
        if antibiotic not in recommended_antibiotics:
            return_pred = 'no'
        return return_pred
        #for key in keys:
            #print(key)
            #print(row)
            #if key in row:
                #val = (row[key] or "").strip() #DICE NO PARA mtrC.disruped only
                #print(val)
                #if 'ceftriaxone' in key.lower():
                #    if 'A311V' in val.upper():
                #    return 'no'
                #return 'no' if val else 'yes'
        #return 'yes'

    def yesno(b: bool) -> str:
        return 'YES' if b else 'NO'

    def canon_regimen(label: str) -> str:
        """
        Canonicalise regimen-category labels (v0.6.7: full names only).

        Accepted examples:
          - ceftriaxone+azithromycin
          - azithromycin+spectinomycin
          - ceftriaxone
          - azithromycin
          - ciprofloxacin
          - spectinomycin
          - zoliflodacin
        """
        s = label.strip().lower().replace(' ', '')
        if s in ('ceftriaxone+azithromycin', 'azithromycin+ceftriaxone'):
            return 'ceftriaxone+azithromycin'
        if s in ('azithromycin+spectinomycin', 'spectinomycin+azithromycin'):
            return 'azithromycin+spectinomycin'
        if s in ('ceftriaxone', 'azithromycin', 'ciprofloxacin', 'spectinomycin', 'zoliflodacin'):
            return s
        return s

    with open(input_file, 'r') as infile, \
         open(alert_output, 'w', newline='') as alert_file, \
         open(treatment_output, 'w', newline='') as treated_file:

        reader = csv.DictReader(infile, delimiter='\t')
        headers = reader.fieldnames or []
        print("\nColumn headers in input file:\n{h}\n".format(h=headers))

        alert_writer = csv.writer(alert_file, delimiter='\t')
        treated_writer = csv.writer(treated_file, delimiter='\t')

        treated_writer.writerow([
            'isolate', 'recommended_1', 'recommended_2',
            'Predicted Profile', 'Recommended Treatment', 'Comment',
            'ceftriaxone+azithromycin',
            'ceftriaxone',
            'azithromycin',
            'azithromycin+spectinomycin',
            'ciprofloxacin',
            'spectinomycin',
            'zoliflodacin',
            'chosen_regimen'
        ])
        alert_writer.writerow(['Alert'] + headers)

        avail = set(a.strip().lower() for a in available_antibiotics)

        # NOTE: azithromycin monotherapy is supported as a regimen category, but is NOT
        # included in the default order. It is only selectable if the user explicitly
        # adds "azithromycin" to --sensitreat_order.
        canon_order = [canon_regimen(x) for x in sensitreat_order] or [
            'ceftriaxone+azithromycin',
            'ceftriaxone',
            'azithromycin+spectinomycin',
            'ciprofloxacin',
            'spectinomycin',
            'zoliflodacin'
        ]
        #print(sensitreat_order)
        for row in reader:
            isolate = row.get('isolate', '')

            # From rules TSV: list of antibiotics that are recommended (susceptible options)
            recommended_antibiotics = [
                x.strip().lower()
                for x in (row.get('treatment recommendation', '') or '').split(',')
                if x.strip()
            ]
            #print(recommended_antibiotics)
            rec_avail = [abx for abx in recommended_antibiotics if abx in avail]
            rec_avail_set = set(rec_avail)
            #print(rec_avail_set)
            # Maintain two simple "top" options among available recommendations, following the order
            # defined by available_antibiotics list as given (not the regimen order).
            treatment_recommendation = []
            #print(available_antibiotics)
            for abx in available_antibiotics:
                abx_l = abx.strip().lower()
                if abx_l in rec_avail_set and abx_l not in treatment_recommendation:
                    treatment_recommendation.append(abx_l)
                if len(treatment_recommendation) >= 2:
                    break
            #print(treatment_recommendation)

            rec1 = treatment_recommendation[0] if len(treatment_recommendation) > 0 else 'None'
            rec2 = treatment_recommendation[1] if len(treatment_recommendation) > 1 else 'None'
            #print(rec1, rec2)

            ceftriaxone_call = get_pred('ceftriaxone', recommended_antibiotics)#, 'CRO_predicted', 'ceftriaxone_predicted')
            azithromycin_call = get_pred('azithromycin', recommended_antibiotics)#row, 'azithromycin_NWT')#, 'AZM_predicted', 'azithromycin_predicted')
            ciprofloxacin_call = get_pred('ciprofloxacin', recommended_antibiotics)#row, 'ciprofloxacin_NWT')#, 'CIP_predicted', 'ciprofloxacin_predicted')
            spectinomycin_call = get_pred('spectinomycin', recommended_antibiotics)#row, 'spectinomycin_NWT')#, 'SPC_predicted', 'spectinomycin_predicted')
            zoliflodacin_call = get_pred('zoliflodacin', recommended_antibiotics)#row, 'zoliflodacin_NWT')#, 'ZOL_predicted', 'zoliflodacin_predicted')

            profile = "ceftriaxone={c}, azithromycin={a}, ciprofloxacin={i}".format(
                c=ceftriaxone_call.upper(),
                a=azithromycin_call.upper(),
                i=ciprofloxacin_call.upper()
            )
            if any(h.lower().startswith('spectinomycin') or h.lower().startswith('spc') for h in headers):
                profile += ", spectinomycin={s}".format(s=spectinomycin_call.upper())

            def has(drug: str) -> bool:
                return drug.lower() in rec_avail_set

            available_regimens = set()
            if has('ceftriaxone') and has('azithromycin'):
                available_regimens.add('ceftriaxone+azithromycin')
            if has('ceftriaxone'):
                available_regimens.add('ceftriaxone')
            if has('azithromycin'):
                available_regimens.add('azithromycin')
            if has('azithromycin') and has('spectinomycin'):
                available_regimens.add('azithromycin+spectinomycin')
            if has('ciprofloxacin'):
                available_regimens.add('ciprofloxacin')
            if has('spectinomycin'):
                available_regimens.add('spectinomycin')
            if has('zoliflodacin'):
                available_regimens.add('zoliflodacin')
            #print(available_regimens)

            pick = 'None'
            for opt in canon_order:
                if opt in available_regimens:
                    pick = opt
                    break
            #print(canon_order)
            #print(pick)

            treatment = ''
            comment = ''
            if pick == 'ceftriaxone+azithromycin':
                treatment = 'Ceftriaxone 1 g IM + Azithromycin 2 g orally'
                comment = 'Acceptable combination therapy (RECOMMENDATION 2)'
            elif pick == 'ceftriaxone':
                treatment = 'Ceftriaxone 1 g IM'
                comment = 'Guideline-based monotherapy (RECOMMENDATION 1)'
            elif pick == 'azithromycin+spectinomycin':
                treatment = 'Spectinomycin 2 g IM + Azithromycin 2 g orally'
                comment = 'Alternative regimen (RECOMMENDATION 3)'
            elif pick == 'azithromycin':
                treatment = 'Azithromycin 2 g orally (single dose)'
                comment = 'Avoid monotherapy due to rapid macrolide resistance selection'
            elif pick == 'ciprofloxacin':
                treatment = 'Ciprofloxacin 500 mg orally'
                comment = 'Use only when susceptibility confirmed (e.g., gyrA S91 wild-type)'
            elif pick == 'spectinomycin':
                treatment = 'Spectinomycin 2 g IM'
                comment = 'Lower cure rates in oropharyngeal infection; avoid for pharyngeal disease when possible (RECOMMENDATION 1)'
            elif pick == 'zoliflodacin':
                treatment = 'Zoliflodacin 3 g orally (single dose)'
                comment = 'Investigational oral option; phase 3 non-inferior to ceftriaxone+azithromycin for uncomplicated urogenital infection'
            else:
                treatment = 'XDR_isolate — manual follow-up'
                comment = 'Flag for review'

            flags = {
                'ceftriaxone+azithromycin': yesno(pick == 'ceftriaxone+azithromycin'),
                'ceftriaxone':             yesno(pick == 'ceftriaxone'),
                'azithromycin':            yesno(pick == 'azithromycin'),
                'azithromycin+spectinomycin': yesno(pick == 'azithromycin+spectinomycin'),
                'ciprofloxacin':           yesno(pick == 'ciprofloxacin'),
                'spectinomycin':           yesno(pick == 'spectinomycin'),
                'zoliflodacin':            yesno(pick == 'zoliflodacin'),
            }

            print("{iso} → ceftriaxone={cro}, azithromycin={azm}, ciprofloxacin={cip}, spectinomycin={spc}  | regimen={pick}".format(
                iso=isolate,
                cro=ceftriaxone_call,
                azm=azithromycin_call,
                cip=ciprofloxacin_call,
                spc=spectinomycin_call,
                pick=pick
            ))

            # write row
            treated_writer.writerow([ 
                isolate,
                rec1 if rec1 != 'None' else 'None',
                rec2 if rec2 != 'None' else 'None',
                profile,
                treatment, comment,
                flags['ceftriaxone+azithromycin'],
                flags['ceftriaxone'],
                flags['azithromycin'],
                flags['azithromycin+spectinomycin'],
                flags['ciprofloxacin'],
                flags['spectinomycin'],
                flags['zoliflodacin'],
                pick
            ])

            if pick == 'None' or (ceftriaxone_call == 'no' and azithromycin_call == 'no' and ciprofloxacin_call == 'no'):
                alert_state = 'XDR' if (ceftriaxone_call == 'no' and azithromycin_call == 'no' and ciprofloxacin_call == 'no') else 'None'
                rowout = [alert_state] + [row.get(h, '') for h in headers]
                alert_writer.writerow(rowout)

def run_sensitreat(parsed_args_or_dict) -> None:
    if isinstance(parsed_args_or_dict, dict):
        input_file = parsed_args_or_dict.get('--input_file') or parsed_args_or_dict.get('--sensiscript_outfile')
        available_antibiotics = parsed_args_or_dict.get(
            '--available_antibiotics',
            'ceftriaxone,azithromycin,ciprofloxacin,tetracycline,penicillin,spectinomycin,zoliflodacin'
        ).split(',')
        sensitreat_order = parsed_args_or_dict.get(
            '--sensitreat_order',
            'ceftriaxone+azithromycin,ceftriaxone,azithromycin+spectinomycin,ciprofloxacin,spectinomycin,zoliflodacin'
        )
        if sensitreat_order:
            sensitreat_order = sensitreat_order.split(',')
        else:
            sensitreat_order = 'ceftriaxone+azithromycin,ceftriaxone,azithromycin+spectinomycin,ciprofloxacin,spectinomycin,zoliflodacin'.split(',')
        alert_output = parsed_args_or_dict.get('--alert_output', 'alert_output.tsv')
        treatment_output = parsed_args_or_dict.get('--treatment_output', 'treatment_output.tsv')
    else:
        input_file = parsed_args_or_dict.input_file
        available_antibiotics = parsed_args_or_dict.available_antibiotics.split(',')
        sensitreat_order = parsed_args_or_dict.sensitreat_order.split(',')
        alert_output = parsed_args_or_dict.alert_output
        treatment_output = parsed_args_or_dict.treatment_output

    available_antibiotics = [a.strip() for a in available_antibiotics if a.strip()]
    sensitreat_order = [a.strip() for a in sensitreat_order if a.strip()]

    process_output(input_file, available_antibiotics, sensitreat_order, alert_output, treatment_output)

    # Generate combined HTML with tabbed interface
    try:
        from html_generator import generate_combined_tabbed_html
        combined_html = treatment_output.replace('.tsv', '.html')
        generate_combined_tabbed_html(
            sensiscript_tsv_path=input_file,
            treatment_tsv_path=treatment_output,
            output_html_path=combined_html,
            antibiotics=available_antibiotics,
            alert_tsv_path=alert_output
        )
    except Exception as e:
        print(f"Warning: Could not generate combined HTML output: {e}")

# -------------------------------------------------------------------
# Pipeline plumbing
# -------------------------------------------------------------------
def parse_arguments_to_dict(arguments: str) -> dict:
    """
    Convert a string of --flag value pairs to a dict.
    Example:
      "--input_file <FILE> --treatment_output <FILE>"
    -> {"--input_file":"<FILE>", "--treatment_output":"<FILE>"}
    """
    parts = arguments.split()
    parsed = {}
    i = 0
    while i < len(parts):
        if parts[i].startswith('--'):
            if i + 1 < len(parts) and not parts[i + 1].startswith('--'):
                parsed[parts[i]] = parts[i + 1]
                i += 2
            else:
                parsed[parts[i]] = None
                i += 1
        else:
            i += 1
    return parsed

def run_pipeline(selected_modules: list, other_arguments: str) -> None:
    arguments_dict = parse_arguments_to_dict(other_arguments)
    for name in selected_modules:
        if name in modules:
            print("Executing module: {m}".format(m=name))
            modules[name](arguments_dict)
        else:
            print("Module '{m}' not recognized.".format(m=name))

# -------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------
modules = {
    "ariba": run_ariba,
    "sensitype": run_sensitype,
    "sensitreat": run_sensitreat,
}

def main():
    parser = argparse.ArgumentParser(
        description="Neisseria gonorrhoeae Sensitype Pipeline (v{v})".format(v=__version__),
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(help="Available modules", dest='command')

    # ariba
    ariba_parser = subparsers.add_parser('ariba', help="Run ARIBA on input directories")
    ariba_parser.add_argument(
        '--input_dirs', required=True,
        help=("Comma-separated directories containing FASTQs.\n"
              "  e.g., --input_dirs <RUN_DIR1>,<RUN_DIR2>")
    )
    ariba_parser.add_argument(
        '--output_dir', required=True,
        help=("Output directory for per-sample ARIBA runs.\n"
              "  e.g., --output_dir <OUTPUT_DIR>")
    )
    ariba_parser.add_argument(
        '--db_path', required=False,
        help=("Path to ARIBA database directory.\n"
              "CLI > $SENSITYPE_ARIBA_DB > ./resources/ariba_db  (fallback)\n"
              "  e.g., --db_path <ARIBA_DB_DIR>")
    )
    ariba_parser.add_argument(
        '--threads', type=int, default=1,
        help=("Threads for ARIBA.\n"
              "  e.g., --threads 8   (default: 1)")
    )
    ariba_parser.set_defaults(func=run_ariba)

    # sensitype
    sensitype_parser = subparsers.add_parser('sensitype', help="Run Sensitype rules (susceptibility prediction)")
    sensitype_parser.add_argument(
        '--input_AMRtable', required=True,
        help=("ARIBA summary CSV/TSV (e.g., ariba_summary.csv/tsv).\n"
              "  e.g., --input_AMRtable <ARIBA_SUMMARY.(csv|tsv)>")
    )
    sensitype_parser.add_argument(
        '--sensiscript_outfile', required=True,
        help=("Output TSV from rules script.\n"
              "  e.g., --sensiscript_outfile <RULES_OUTFILE.tsv>")
    )
    sensitype_parser.add_argument(
        '--sensiscript_db', required=False,
        help=("Path to sensitype rules DB.\n"
              "CLI > $SENSITYPE_DB > ./sensitype.db")
    )
    sensitype_parser.add_argument(
        '--sensiscript_pena', required=False,
        help=("Path to penA DB.\n"
              "CLI > $SENSITYPE_PENA_DB > ./sensitype.penA.db")
    )
    sensitype_parser.add_argument(
        '--sensiscript_antibiotics',
        default='ceftriaxone,azithromycin,ciprofloxacin,tetracycline,penicillin,zoliflodacin',
        help=("Comma-separated antibiotics (evaluation order for rules stage).\n"
              "  e.g., --sensiscript_antibiotics ceftriaxone,ciprofloxacin,azithromycin,spectinomycin\n"
              "  (default: ceftriaxone,azithromycin,ciprofloxacin,tetracycline,penicillin,zoliflodacin)")
    )
    sensitype_parser.set_defaults(func=run_sensitype)

    # sensitreat
    treat_parser = subparsers.add_parser('sensitreat', help="Assign treatment recommendations")
    treat_parser.add_argument(
        '--input_file', required=True,
        help=("Sensitype rules output TSV (the file with mutation/mechanism columns).\n"
              "  e.g., --input_file <RULES_OUTFILE.tsv>")
    )
    treat_parser.add_argument(
        '--available_antibiotics',
        default='ceftriaxone,azithromycin,ciprofloxacin,tetracycline,penicillin,spectinomycin,zoliflodacin',
        help=("Comma-separated list of available antibiotics in your setting.\n"
              "Used to filter which recommendations can be chosen and to pick recommended_1/2.\n"
              "  e.g., --available_antibiotics ceftriaxone,ciprofloxacin,azithromycin,spectinomycin\n"
              "  (default: ceftriaxone,azithromycin,ciprofloxacin,tetracycline,penicillin,spectinomycin,zoliflodacin)")
    )
    treat_parser.add_argument(
        '--sensitreat_order',
        default='ceftriaxone+azithromycin,ceftriaxone,azithromycin+spectinomycin,ciprofloxacin,spectinomycin,zoliflodacin',
        help=("Priority order for regimen categories (paper categories).\n"
              "Accepts combos and singles (case-insensitive) using full names only, e.g.:\n"
              "  --sensitreat_order ceftriaxone+azithromycin,ceftriaxone,azithromycin+spectinomycin,ciprofloxacin,spectinomycin,zoliflodacin\n"
              "Canonical categories: ceftriaxone+azithromycin, ceftriaxone, azithromycin, azithromycin+spectinomycin,\n"
              "                     ciprofloxacin, spectinomycin, zoliflodacin\n"
              "  (default: ceftriaxone+azithromycin,ceftriaxone,azithromycin+spectinomycin,ciprofloxacin,spectinomycin,zoliflodacin)\n"
              "NOTE: azithromycin monotherapy is supported but NOT in the default order; include 'azithromycin' explicitly to enable it.")
    )
    treat_parser.add_argument(
        '--alert_output', default='alert_output.tsv',
        help=("Output TSV for flagged isolates (XDR/unassigned).\n"
              "  e.g., --alert_output <ALERTS_OUT.tsv>  (default: alert_output.tsv)")
    )
    treat_parser.add_argument(
        '--treatment_output', default='treatment_output.tsv',
        help=("Output TSV with recommended_1/2 and regimen flags.\n"
              "  e.g., --treatment_output <TREATMENT_OUT.tsv>  (default: treatment_output.tsv)")
    )
    treat_parser.set_defaults(func=run_sensitreat)

    # pipeline
    pipe_parser = subparsers.add_parser('pipeline', help="Run multiple modules in order")
    pipe_parser.add_argument(
        '--modules', required=True,
        help=("Comma-separated list: ariba,sensitype,sensitreat\n"
              "  e.g., --modules ariba,sensitype,sensitreat")
    )
    pipe_parser.add_argument(
        '--other_arguments', required=True,
        help=('Quoted string of flags for modules. Example:\n'
              '  --other_arguments "--input_dirs <RUN_DIRS> '
              '--output_dir <OUTPUT_DIR> '
              '--db_path <ARIBA_DB_DIR> '
              '--sensiscript_outfile <RULES_OUTFILE.tsv> '
              '--available_antibiotics ceftriaxone,ciprofloxacin,azithromycin,spectinomycin '
              '--sensitreat_order ceftriaxone+azithromycin,ceftriaxone,azithromycin+spectinomycin,ciprofloxacin,spectinomycin,zoliflodacin '
              '--treatment_output <TREATMENT_OUT.tsv>"')
    )

    # rename
    rename_parser = subparsers.add_parser('rename', help="Standardize fastq file names")
    rename_parser.add_argument(
        '--directories', required=True,
        help=("Comma-separated directories to normalize FASTQ names in.\n"
              "  e.g., --directories <RUN_DIR1>,<RUN_DIR2>")
    )
    rename_parser.add_argument(
        '--pre', action='store_true',
        help=("Preview only (no changes).\n"
              "  e.g., --pre")
    )

    args = parser.parse_args()

    if args.command == 'pipeline':
        selected_modules = args.modules.split(',')
        other_arguments = args.other_arguments
        run_pipeline(selected_modules, other_arguments)
    elif args.command == 'rename':
        for directory in args.directories.split(','):
            rename_files(directory.strip(), preview=args.pre)
    elif hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
