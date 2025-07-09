import os
import glob

#USAGE: Run in the parent Maize directory. Output will be a text file

def obtain_file_directories(root_dir): #Get all of the file paths for each .txt file in the maize folders

    for item in os.listdir(root_dir):   #For each item in the main directory...
        
        if os.path.isdir(item):             #If it is a subdirectory...
            #print(os.path.isdir(item))
            temp_path = os.path.join(cwd, item) #Join the subdirectory to the path

            for item1 in glob.glob(os.path.join(temp_path, '*.txt')):   #For each item with a .txt extension
                #print(f'Item: {item1}')
                #print(os.path.join(temp_path, '*.tsv'))
                file_directories.append(os.path.join(temp_path, item1))     #Add item to the data path
                
# def get_gene_name():    #Get all of the gene names from the file paths... For later use in sorting
    
#     i = 0
    
#     for item in file_directories:               #Get the names based on the structure of the directory
#         temp = file_directories[i].split('/')
#         gene_names.append(temp[7])
#         i+=1

def get_dED_stats():    #Get the average, highest, and lowest EDs for each gene
   
    highest = float(0)  #Hold the highest dED
    lowest = float(0)   #Hold the lowest dED
    #j = 0
    for file in file_directories:   #For each file ...
        
        f = open(file, "r") #Open the file in read mode
        i = 0               #Hold number of scores added together
        total_dED = 0       #Hold the total ED (added together)
        gene_names.append(f.readline().strip().split('\t')[13])
        for line in f.readlines()[1:]:  #For each line (skipping the header line) ...

            line_set = line.strip().split('\t') #Split the line with each tab
            dED = float(line_set[7])            #Define column 8 as the dED scores
            # print(dED)
            total_dED += dED    #Add the dED to the running total
            i += 1              #Increment the running total by 1
            all_dED.append(dED)
            if (dED > highest): #Take care of the highest and lowest scores for each gene
                highest = dED
            if (dED < lowest):
                lowest = dED
        
        f.close()

        # f = open(gene_names[j]+".txt", 'w')   #Create files for each gene including every dED value
        # j += 1
        # f.write("dED\n")
        # for item in all_dED:
        #     f.write(str(item) + "\n")
        # f.close()                             #End creating files ^^^

        avg_dED.append(total_dED / i)
        highest_dED.append(highest)
        lowest_dED.append(lowest)

def create_tsv_of_values():
    f = open("higest_dED_stats.tsv", 'w')
    f.write("gene_name\tdED\n")
    i = 0
    while (i < len(gene_names)):
        f.write(gene_names[i] + "\t" + str(highest_dED[i]) +"\n")
        i+=1

    f = open("lowest_dED_stats.tsv", 'w')
    f.write("gene_name\tdED\n")
    i = 0
    while (i < len(gene_names)):
        f.write(gene_names[i] + "\t" + str(lowest_dED[i]) +"\n")
        i+=1

    f = open("average_dED_stats.tsv", 'w')
    f.write("gene_name\tdED\n")
    i = 0
    while (i < len(gene_names)):
        f.write(gene_names[i] + "\t" + str(avg_dED[i]) +"\n")
        i+=1



if __name__ == "__main__":      #MAIN SCRIPT

    cwd = os.getcwd()   #Get CWD
    #print(cwd)

    file_directories = []       #Hold all of the file directories for later use
    obtain_file_directories(cwd)    #Get all the file paths

    gene_names = []             #Hold all gene names for later use (Same order as file directories)
    # get_gene_name()                 #Get all the gene names

    #print(gene_names)
    highest_dED = []            #Hold all of the highest dEDs in a file (Same order as file directories)
    lowest_dED = []             #Hold all of the lowest dEDs in a file (Same order as file directories)
    avg_dED = []                #Hold all of the average dEDs in a file (Same order as file directories)
    all_dED = []
    #print(file_directories)

    get_dED_stats()             #Get all the average, highest and lowest dEDs
    # print(highest_dED)
    # print(lowest_dED)
    # print(avg_dED)

    create_tsv_of_values()