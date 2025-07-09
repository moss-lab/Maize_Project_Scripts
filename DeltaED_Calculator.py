# Written by Mason Eisenhauer - Moss Lab              masoneis@iastate.edu
# This script will gather the Ensemble Diversity (ED) values from TWO ScanFold output files to calculate the Delta ED and output 
# into a wig file format for easier viewing in IGV. The script will also output the standard deviation of the given values.
# The script can reverse the wig file if "--r" is called. 

# USAGE: python DeltaED_Wig_Calculator.py <filename(HIGHERTEMP)> <filename(LOWERTEMP)> <gene_name> <chromosome_number(FORMAT: "chr#"")> <start_coord(OPTIONAL)> <--r(OPTIONAL)>
# EXAMPLE: python DeltaED_Calculator.py 42CTemp.out 28CTemp.out ZmHsf04 chr5 129840392 --r
# EXAMPLE: python DeltaED_Calculator.py 42CTemp.out 28CTemp.out ZmHsf12 chr10 --r
# Example: python DeltaED_Calculator.py 42CTemp.out 28CTemp.out ZmHsf31 chr2

import sys
import math
import os

file1 = sys.argv[1] # Retrieve higher temp file path from system argument
file2 = sys.argv[2] # Retrieve lower temp file path from system argument
outname = sys.argv[3] # Retrieve the name of the gene to use as part of the file name
chrom = sys.argv[4] # Retrieve the chromosome number in format chr#

# Retrieve the starting position of the gene (Defaults to 1 if not provided)
try:
    start_coord = sys.argv[5]
except:
    start_coord = 1

# Retrieve the flag to run reversal of the wig file
if sys.argv[5] == "--r" or sys.argv[6] == "--r":
    runReverseScript = True

with open(file1, 'r') as f1, open(file2, 'r') as f2:    # Open both files in read mode

    next(f1)    # Skip the first line
    next(f2)    # Skip the first line

    ed1 = []
    ed2 = []

    # Gather the 7th value from the files and convert to a float
    for line in f1:
        ed1.append(float(line.strip().split('\t')[6]))
    for line in f2:
        ed2.append(float(line.strip().split('\t')[6])) 

    deltaED = [] 

    # Append all of the DeltaED values to deltaED
    for i in range(len(ed1)):
        deltaED.append(ed1[i] - ed2[i])  

output_file = f"{outname}_DeltaED.wig" # Name the output file

with open(output_file, 'w') as f:   

    # Write the header of the wig file
    f.write("fixedStep chrom=" + chrom + " start=" + str(start_coord) + " step=1 span=1\n") 

    # Write all of the DeltaED values into the wig file
    for j in deltaED:
        f.write(f"{j}\n") 

# Calculate the standard deviation of the deltaED values
mean_ED = sum(deltaED) / len(deltaED)
squared_deviations = [math.pow(value - mean_ED, 2) for value in deltaED]
variance = sum(squared_deviations) / len(squared_deviations)
standard_deviation = math.sqrt(variance)

std_file = f"{outname}_std.out"

# Write the standard deviation in a new out file
with open(std_file, 'w') as f_std:
    f_std.write(f"Standard deviation:\n{standard_deviation}\n")

# Only run if the reverse strand flag is called
if runReverseScript == True:

    #Close all files in preparation for reverse script
    f1.close
    f2.close
    f.close
    f_std.close

    # Reverse the wig coordinates to be in line for reverse strands (Pulled from Rev_Strand_Wig_Convert.py by Warren with minimal modification)
    with open(output_file, 'r') as wig, open("Reversed_" + outname +"_DeltaED.wig", 'w') as reverse:
        reverse.write(str(f"fixedStep chrom={chrom} start={start_coord} step=1 span=1\n"))
        for line in reversed(open(output_file).readlines()[1:]):
            reverse.write(line)
    
    # Close and delete the forward strand output file (Only called when reverse strand is made)
    wig.close
    os.remove(output_file)