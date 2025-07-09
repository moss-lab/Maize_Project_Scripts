# IGV_File_Preparation

'''
By: Mason W. Eisenhauer, Moss Lab, Iowa State University
masoneis@iastate.edu

# Some code adapted from scripts by Warren B. Rouse

This program will pull all of the required filepaths for each item that is meant to be loaded into IGV and use them
to create IGV ready files that will line up on a genome wide scale. This will allow easier usage of this data.

Usage: Place into the folder that includes all of the directories of each gene and run through the terminal.

New directories will be created inside of each gene's directory called ''IGV_Prepared'' that will include the converted
files from usage of this program. 
'''

import os
import glob

# Create or obtain directories within each folder in the CWD // DONE
def open_define_filepaths():

    scanfold_data_path = []                         # Create a list that will hold all directory paths of the unprepared .wig and .bp files
    target_folder_path = []                         # Create a list that will hold all target path locations for prepared .wig and .bp files
    fasta_path = []                                 # Create a list that will hold all file paths for the .fasta files of each gene
    
    i = 0                                           # Hold the value for each iteration of the following code 

    for item in os.listdir(cwd):                                                                # Do the following for every item in the current working directory
        if os.path.isdir(item):                                                                 # If the item found is a directory, do the following (Should only find gene directories)
            #print(f'Initial Directory: {item}')

            temp_path = os.path.join(cwd, item)                                                 # Build a temporary path 
            #print(f'ScanfoldPath:{temp_path}')
            for item1 in glob.glob(os.path.join(temp_path, '*')):                               # For each item in the gene folder, do the following
                #print(f'Secondary Directory/Items: {item1}')
                if os.path.isdir(item1):                                                        # If the item found is a directory, do the following (Should only find scanfold output data directories)
                    #print(f'InnerItemRecognized: {item1}')

                    scanfold_data_path.append(os.path.join(temp_path, item1, "igv_files"))      # Append "igv_files" to the directory name, this will be the location of the files that need prepared
            
            target_folder_path.append(os.path.join(cwd, item, "IGV_Prepared"))                  # Build the target folder directory
            os.makedirs(target_folder_path[i])                                                  # Create the directories for the target folders
            fasta_path.append(temp_path)                                                        # Append the temporary path to the fasta_path list

            #print(f'item:   {item}')
            print(f'Scanfold Path:    {scanfold_data_path[i]}')
            #print(f'Item1:    {item1}')
            print(f'Target Path:    {target_folder_path[i]}')

            i += 1                                                                              # Update the iteration
    print("Filepaths found and created successfully")
    return scanfold_data_path, target_folder_path, fasta_path           # Returns a list of directories to be used in later processes

# Determine all of the coordinates and strand directions for each path and appends them to a list for later processing
def determine_coords_and_direction(fasta_path):

    start_coord = []                                # Create a list to hold all of the start coordinates for each gene
    end_coord = []                                  # Create a list to hold all of the end coordinates for each gene
    direction = []                                  # Create a list to hold all of the gene strand directions
    chromosome = []                                 # Create a list to hold all of the chromosomes that each gene is one

    i = 0                                           # Hold the value for each iteration of the following code

    for each in range(len(fasta_path)):                                                         # For every gene, do the following
        for item in glob.glob(os.path.join(fasta_path[i], '*.fasta')):                          # For the current iteration .fasta file, do the following
            #print(f'Fasta_Path I: {fasta_path[i]}')
            with open(item, 'r') as fasta:                                                      # Open the .fasta file and do the following
                firstline = fasta.readline().strip()                                            # Strip the first line of the fasta
                #print(firstline)                            
                parts = firstline.split(':')                                                    # Split the first line into separate parts
                #print(parts)
                num_parts = parts[1].split('-')                                                 # Split on of the parts into further parts
                #print(num_parts)
                if len(num_parts) == 3:                                                         # Determine the direction of the strand and append
                    direction.append('-')                   
                else:                                       
                    direction.append('+')
                chrom_fragments = parts[0].split('_')                                           # With another part, separate into further parts

                if len(chrom_fragments) == 2:                                                   # Determine the name of the chromosome or scaffold and append
                    chromosome.append(chrom_fragments[1])
                else: 
                    chrom_stitch = chrom_fragments[1] + '_' + chrom_fragments[2]
                    chromosome.append(chrom_stitch)

                start_coord.append(num_parts[0])                                                # Append the start coordinates
                end_coord.append(num_parts[1].split('+')[0])                                    # Append the end coordinates
                #print(start_coord[i])
                #print(end_coord[i])
                #print(direction[i])
                #print(chromosome[i])
                
                i += 1                                                                          # Update the iteration
    print("Coordinates, directions, and chromosomes found and stored successfully")
    return start_coord, end_coord, direction, chromosome            # Return the start coordinates, end coordinates, strand direction, and chromosome/scaffold

# Create all of the IGV prepared files in the target location using previous data
def create_igv_prepared_files(sf_data_path, target_path, start, end, direction, chromosome):
    
    i = 0                                           # Hold the value for each iteration of the following code

    for each in range(len(sf_data_path)):                                                       # For each file directory, do the following
        if direction[i] == '+':                                                                 # If positive strand, do the following

            #print("Forward Strand Initiated")
            for fwd_wig_file in glob.glob(os.path.join(sf_data_path[i], '*.wig')):                  # For each wig file, convert and output
                fwd_strand_wig_conversion(fwd_wig_file, target_path[i], start[i], chromosome[i])
            for fwd_bp_file in glob.glob(os.path.join(sf_data_path[i], '*.bp')):                    # For each bp file, convert and output
                fwd_strand_bp_conversion(fwd_bp_file, target_path[i], start[i], chromosome[i])

        else:                                                                                   # If negative strand, do the following

            #print("Negative Strand Initiated")
            for neg_wig_file in glob.glob(os.path.join(sf_data_path[i], '*.wig')):                  # For each wig file, convert and output
                neg_strand_wig_conversion(neg_wig_file,target_path[i], start[i], chromosome[i])
            for neg_bp_file in glob.glob(os.path.join(sf_data_path[i], '*.bp')):                    # For each bp file, convert and output
                neg_strand_bp_conversion(neg_bp_file, target_path[i], end[i], chromosome[i])

        i += 1                                                                                      # Update the iteration
    print("IGV files created and stored successfully")
    
# Create the forward strand wig file in the target path location. // Code adapted from Fwd_Strand_Wig_Convert.py by Warren B. Rouse
def fwd_strand_wig_conversion(wig, path, start, chrom):

    filename = wig                      #wig file
    basename = os.path.basename(wig)    # Take the basename of the file
    path += "/Genomic_" + basename      # Append the basename of the file path to fully create the target file path

    with open(filename , 'r') as wig, open(path, "a") as genomic:                       #Open the input wig file in read mode and open the new wig file in write mode
        genomic.write(str(f"fixedStep chrom={chrom} start={start} step=1 span=1\n"))    #Write the header using the user input coordinates and chromosome number
        for line in open(filename).readlines()[1:]:                                     #Read the metrics from the input file
            genomic.write(line)                                                         #Write the metrics to the new file
    #print("fwd_strand_wig completed successfully")

# Create the forward strand bp file in the target path location // Code adapted from Transcript_to_Genomic_Coord.py by Warren B. Rouse
def fwd_strand_bp_conversion(bp, path, start, chrom):

    bpi = []				#Defining the column in the bp file that will be turned into a list
    bpj = []                #Defining the column in the bp file that will be turned into a list
    bpi2 = []               #Defining the column in the bp file that will be turned into a list
    bpj2 = []               #Defining the column in the bp file that will be turned into a list
    color = []              #Defining the column in the bp file that will be turned into a list

    basename = os.path.basename(bp)     # Take the basename of the file
    path += "/Genomic_" + basename      # Append the basename to the file path to fully create the target file path

    with open(bp , 'r') as bpfile, open(path, "w") as genomic:						                        #Open the input bp file in read mode and open the new bp file in write mode
        color_head = bpfile.readlines()[0:7]										                        #Read the first 7 lines of the input bp file that do not need modified
        for line in color_head:
            genomic.write(line)														                        #Write these first 7 lines to the new bp file
    with open(bp , 'r') as bpfile, open(path, "a") as genomic:						                        #Open the input bp file in read mode again and open the new bp file with the added 7 lines in write mode
        lines = bpfile.readlines()[7:]												                        #Read the lines of the input starting at line 8, which is the first line that needs modified
        for line in lines:
            data = line.split('\t')													                        #Create new variable called data that contains all parts of the input file broken up at the tab
            bpi = (int(start)+int(data[1])-int(1))									                        #Add the input coordinate to the second peice of data, bp i, and subtract 1 to give new genomic coordinate as the starting postition
            bpi2= (int(start)+int(data[2])-int(1))
            bpj = (int(start)+int(data[3])-int(1))									                        #Add the input coordinate to the fourth peice of data, bp j, and subtract 1 to give new genomic coordinate as the ending postition
            bpj2 = (int(start)+int(data[4])-int(1))
            color = data[5]															                        #Color value that does not need modified
            genomic.write(str(f"{chrom}\t{str(bpi)}\t{str(bpi2)}\t{str(bpj)}\t{str(bpj2)}\t{color}"))		#Use an f string literal to print all new coordinates and color column under the first 7 lines
    #print("fwd_strand_bp completed successfully")

# Create the reverse strand wig file in the target path location // Code adapted from Rev_Strand_Wig_Convert.py by Warren B. Rouse
def neg_strand_wig_conversion(wig, path, start, chrom):

    filename = wig                              #wig file
    basename = os.path.basename(wig)            # Take the basename of the file
    path += "/Reversed_Genomic_" + basename     # Append the basename to the file path to fully create the target file path

    with open(filename , 'r') as wig, open(path, "a") as genomic:                       #Open the input wig file in read mode and open the new wig file in write mode
        genomic.write(str(f"fixedStep chrom={chrom} start={start} step=1 span=1\n"))    #Write the header using the user input coordinates and chromosome number
        for line in reversed(open(filename).readlines()[1:]):                           #Reverse the order of the metrics 
            genomic.write(line)                                                         #Write the reversed order metrics to the new file
    #print("neg_strand_wig completed successfully")

# Create the reverse strand bp file in the target path location // Code adapted from Rev_Transcript_to_Genomic_Coord.py by Warren B. Rouse
def neg_strand_bp_conversion(bp, path, end, chrom):

    bpi = []				#Defining the column in the bp file that will be turned into a list
    bpj = []				#Defining the column in the bp file that will be turned into a list
    bpi2 = []				#Defining the column in the bp file that will be turned into a list
    bpj2 = []				#Defining the column in the bp file that will be turned into a list
    color = []				#Defining the column in the bp file that will be turned into a list

    basename = os.path.basename(bp)         # Take the basename of the file
    path += "/Reverse_Genomic_" + basename  # Append the basename to the file path to fully create the target file path

    with open(bp , 'r') as bpfile, open(path, "w") as genomic:						                        #Open the input bp file in read mode and open the new bp file in write mode
        color_head = bpfile.readlines()[0:7]															    #Read the first 7 lines of the input bp file that do not need modified
        for line in color_head:
            genomic.write(line)																			    #Write these first 7 lines to the new bp file
    with open(bp , 'r') as bpfile, open(path, "a") as genomic:						                        #Open the input bp file in read mode again and open the new bp file with the added 7 lines in write mode
        lines = bpfile.readlines()[7:]																	    #Read the lines of the input starting at line 8, which is the first line that needs modified
        for line in lines:
            data = line.split('\t')																		    #Create new variable called data that contains all parts of the input file broken up at the tab
            bpi = (int(end)-int(data[1])+int(1))													        #Add the input coordinate to the second peice of data, bp i, and subtract 1 to give new genomic coordinate as the starting postition
            bpi2= (int(end)-int(data[2])+int(1))
            bpj = (int(end)-int(data[3])+int(1))													        #Add the input coordinate to the fourth peice of data, bp j, and subtract 1 to give new genomic coordinate as the ending postition
            bpj2 = (int(end)-int(data[4])+int(1))
            color = data[5]																				    #Color value that does not need modified
            genomic.write(str(f"{chrom}\t{str(bpi)}\t{str(bpi2)}\t{str(bpj)}\t{str(bpj2)}\t{color}"))		#Use an f string literal to print all new coordinates and color column under the first 7 lines

# Obtain the current working directory
cwd = os.getcwd()
#print(f'CWD:{cwd}')


#Get all filepaths prepared
scanfold_data_path, target_folder_path, fasta_path = open_define_filepaths()    
#print(scanfold_data_path[1])
#print(target_folder_path[1])
#print(fasta_path[1])

# Determine the coordinates, directions, and chromosomes of each gene
start_list, end_list, direction_list, chromosome_list = determine_coords_and_direction(fasta_path)                         
#print(s)
#print(e)
#print(d)

# Create the output files in the target directory
create_igv_prepared_files(scanfold_data_path, target_folder_path, start_list, end_list, direction_list, chromosome_list)   
print("Completed Successfully")