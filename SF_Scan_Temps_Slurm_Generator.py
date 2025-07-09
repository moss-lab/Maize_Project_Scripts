import os
import pathlib
from pathlib import Path
import glob
import sys

def shell_build_start(filename, job, email='masoneis@iastate.edu', time=5, nodes=1, mem=250, tasks=98, notify='ALL'):
    # Build shell scripts with Pronto/HPC settings and modules
    with open(filename, 'w', newline='\n') as writefile:
        writefile.writelines('#!/bin/bash -l\n')
        writefile.writelines('#SBATCH --partition=nova\n')  # Partition to submit to
        writefile.writelines(f'#SBATCH --time={time}-00:00:00\n')  # Time limit for this job in days
        writefile.writelines(f'#SBATCH --nodes={nodes}\n')  # Nodes to be used for this job during runtime
        if mem != 0:
            writefile.writelines(f'#SBATCH --mem={mem}G\n')  # Optional memory allocation
        writefile.writelines(f'#SBATCH --ntasks-per-node={tasks}\n')
        writefile.writelines(f'#SBATCH --job-name={job}\n')  # Name of this job in work queue
        writefile.writelines(f'#SBATCH --mail-user={email}\n')  # Email to send notifications to
        writefile.writelines(f'#SBATCH --mail-type={notify}\n')  # Email notification type (BEGIN, END, FAIL, ALL)
        writefile.writelines('\n')
        # Add Pronto modules
        writefile.writelines('#SBATCH --export=NONE\n\n')
        writefile.writelines('module purge\n')
        writefile.writelines('module load micromamba\n')
        writefile.writelines('eval "$(micromamba shell hook --shell=bash)"\n')
        writefile.writelines('micromamba activate /lustre/hdd/LAS/wmoss-lab/programs/envs/ScanFold2\n')
        writefile.writelines('wait;\n\n')

# Reads all fasta files and outputs a shell script for sf2
count = 1
current_size = 0
fasta_size = 0
max_size = 100  # Number of runs per shell script
max_fasta_size = 50000000  # Files over 8GB are run separately to aid with HPC out of memory errors
fasta_directory = Path.cwd()
sf2_com = 'python /lustre/hdd/LAS/wmoss-lab/masoneis/testing/ScanFold-Scan_Temperatures.py'
print("Gate 1 Success")
print(fasta_directory)


# for fasta_filepath in fasta_directory.glob('Zm*/*.fasta'):
#    print(fasta_filepath)
#    # Exit if any motifs with known errors are present
#    print(f'Error with ENSG {fasta_filepath}, check error and move on')
#    sys.exit()


for fasta_filepath in fasta_directory.glob('Zm*/*.fasta'):
    print("Gate 2 Success")
    print(f'fasta file: {fasta_filepath.name}')
    if current_size == 0:
        # Build initial shell scripts
        shell_build_start(f'delta_sf_scan_writefile_{count}.sh', f'delta_sf_scan_premRNA_{count}')
    with open(f'delta_sf_scan_writefile_{count}.sh', 'a', newline='\n') as writefile:
        writefile.writelines(f'{sf2_com} -i {fasta_filepath} &\n')  #--folder_name {pathlib.Path(fasta_filepath).stem}
        # Find matching sequence, structure, database files
        fasta_size += int(os.path.getsize(fasta_filepath))
        print(f'DB size: {fasta_size}\tMax size: {max_fasta_size}')
        # Check size and separate out large database files
        if fasta_size > max_fasta_size:
            writefile.writelines('wait;\n')
            current_size = 0
            count += 1
            fasta_size = 0
        elif current_size >= max_size:
            # Close out shell scripts when full, increase shell count and reset size count
            writefile.writelines('wait;\n')
            current_size = 0
            count += 1
            fasta_size = 0
        else:
            current_size += 1
with open(f'delta_sf_scan_writefile_{count}.sh', 'a', newline='\n') as writefile:
    writefile.writelines('wait;\n')