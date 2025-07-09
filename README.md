# Maize_Project_Scripts
Extra scripts created to analyze and visualize data from the Zea Mays database

IMPORTANT:  all scripts: Modifications will be needed based on file setup, environment setup, and systems by which runs are completed. These serve as files to be modified for individual use.


SCRIPTS:
DeltaED_Calculator.py - Can be used individually to generate DeltaED files to be loaded into IGV from ED files fed into the script. This script is meant for individual genes as a tool for visualization, as ScanFold-Scan_Temperatures.py will be a more complete analysis.

ED_Stats_Generator.py - Generates files based on the lowest, highest, and average DeltaED across multiple genes (or the entire database). The content of these files can be easily converted into excel format to create graphics based on that. Reads from the output given by ScanFold-Scan_Temperatures.py file.

IGV_File_Preparation.py - Scans the entire set of genes ran by ScanFold2.0, and pulls IGV_Prepared files, converts them to the specification needed, and outputs them within the folder. This is useful for having files ready and correct for visualization within IGV.

SF_Scan_Temps_Slurm_Generator.py - Creates the SLURM job files needed to run on servers. Creates multiple job files such that one job doesn't have too many individual tasks, conserving resources and allowing for jobs to be split across time as needed.

ScanFold-Scan_Temperatures.py - Modified ScanFold-Scan.py script that runs the program twice, exempting z-score and p-value calculations. Outputs a text file with data on MFE values and ED values for each temperature and for the changes in temperatures. This file will also include the predicted dot-bracket notation structures. 
