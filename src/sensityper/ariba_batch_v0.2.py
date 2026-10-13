import os
import re
import argparse
import subprocess

def run_ariba(input_dirs, output_dir, db_path, threads):
    # Split the input directories by comma
    dirs = input_dirs.split(',')

    # Define possible forward and reverse suffix patterns
    forward_suffixes = ['_1.fastq.gz', '_R1.fastq.gz', '_R1_001.fastq.gz']
    reverse_suffixes = ['_2.fastq.gz', '_R2.fastq.gz', '_R2_001.fastq.gz']

    # Loop through each directory
    for d in dirs:
        # Get all the forward read files in the directory based on possible suffixes
        forward_files = [f for f in os.listdir(d) for suffix in forward_suffixes if f.endswith(suffix)]

        for f in forward_files:
            # Extract the sample name based on the detected suffix
            for suffix in forward_suffixes:
                if f.endswith(suffix):
                    sample_name = f.rsplit(suffix, 1)[0]
                    reverse_suffix = reverse_suffixes[forward_suffixes.index(suffix)]
                    break

            # Construct the reverse read file name
            reverse_file = os.path.join(d, sample_name + reverse_suffix)

            # Construct each output folder inside the main output directory
            output_folder = os.path.join(output_dir, sample_name + '_ARIBA')

            # Check if the output folder already exists
            if os.path.exists(output_folder):
                print(f"Skipping {output_folder} as it already exists.")
                #continue

            # Run the ARIBA command with the specified number of threads
            cmd = [
                'ariba', 'run', '--verbose', '--threads', str(threads), db_path,
                os.path.join(d, f), reverse_file, output_folder
            ]
            subprocess.run(cmd)

            # Modify the report.tsv file to include the penA.insD345 if detected
            report_file = os.path.join(output_folder, 'report.tsv')
            new_report_file = os.path.join(output_folder, 'report_complete.tsv')
            with open(report_file, 'r') as infile, open(new_report_file, 'w') as outfile:
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

    # Generate filenames.txt for the ARIBA summary including absolute paths
    #with open('filenames.txt', 'w') as f:
    #    report_files = [os.path.join(output_dir, d, 'report_complete.tsv') for d in os.listdir(output_dir) if d.endswith('_ARIBA')]
    #    #report_names = [d.replace('_ARIBA', '') for d in os.listdir(output_dir) if d.endswith('_ARIBA')]
    #    for count, report_file in enumerate(report_files):
    #        #print(report_file)
    #        f.write(report_file+'\n')

    # Execute ARIBA summary
    #cmd = [
    #    'ariba', 'summary', 'ariba_summary', '-f', 'filenames.txt', '--cluster_cols',
    #    'assembled,ref_seq,pct_id', '--col_filter', 'n',
    #    '--row_filter', 'n', '--no_tree', '--v_groups', '--known_variants'
    #]
    #subprocess.run(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ARIBA on multiple directories.")
    parser.add_argument('-d', '--input_dirs', required=True, help="Comma-separated list of input directories containing FASTQ files.")
    parser.add_argument('-o', '--output_dir', required=True, help="Output directory for ARIBA results.")
    parser.add_argument('--db_path', required=True, help="Path to the ARIBA database.")
    parser.add_argument('-t', '--threads', type=int, default=1, help="Number of threads for ARIBA run.")

    args = parser.parse_args()

    # Construct the output folder
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        create_folder = ['mkdir', output_dir]
        subprocess.run(create_folder)
    output_dir = os.path.abspath(output_dir)
    print(output_dir)

    run_ariba(args.input_dirs, output_dir, args.db_path, args.threads)

