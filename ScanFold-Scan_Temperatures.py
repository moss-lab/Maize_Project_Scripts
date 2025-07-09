# UPDATED TO RUN SCANFOLD_SCAN WITHOUT CALCULATING z-scores or p-values AT 28C and 42C AND PACKAGING THE RESULTS INTO A .TXT FILE #

"""

ScanFold-Scan
Contact: Ryan Andrews - randrews@iastate.edu

This program takes a fasta input file and uses a scanning window approach to
calculate thermodynamic z-scores for individual windows.

Usage:
$ python3.6 ScanFold-Scan.py filename [options]

"""

import sys
import argparse
import string
import re
import numpy as np
sys.path.append('/usr/local/lib/python3.6/site-packages')
import RNA
import random
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from Bio import SeqIO

#### Parsing arguments ####
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--filename', type=str,
                    help='input filename')
parser.add_argument('-s', type=int, default=1,
                    help='step size')
parser.add_argument('-w', type=int, default=120,
                    help='window size')
parser.add_argument('-r', type=int, default=50,
                    help='randomizations')
parser.add_argument('-t', type=int, default=37,
                    help='Folding temperature')
parser.add_argument('-type', type=str, default='mono',
                    help='randomization type')
parser.add_argument('-p', '--print_to_screen', action='store_true',
                    help='print to screen option (default off)')
parser.add_argument('--print_random', type=str, default='off',
                    help='print to screen option (default off)')
parser.add_argument('-c', '--constraints', type=str,
                    help='optional | input constraint file')


args = parser.parse_args()
myfasta = args.filename
step_size = int(args.s)
window_size = int(args.w)
randomizations = int(args.r)
temperature = int(args.t)
type = str(args.type)
print_to_screen = args.print_to_screen
print_random = str(args.print_random)
constraints = args.constraints


#### Defining global variables ###############
temperature = 28
# w = open(myfasta+".forward.win_"+str(window_size)+".stp_"+str(step_size)+".rnd_"+str(randomizations)+".shfl_"+str(type)+"28C.txt", 'w')

##### Establish global "folding model" incorporting temperature (or other settings) #####
md = RNA.md()
md.temperature = int(temperature)

def multiprocessing(func, args,
                    workers):
    with ProcessPoolExecutor(workers) as ex:
        res = ex.map(func, args)
    return list(res)

#### Defining Dinucleotide function #####
# Taken from
# altschulEriksonDinuclShuffle.py
# P. Clote, Oct 2003
# NOTE: One cannot use function "count(s,word)" to count the number
# of occurrences of dinucleotide word in string s, since the built-in
# function counts only nonoverlapping words, presumably in a left to
# right fashion.
def computeCountAndLists(s):
  #WARNING: Use of function count(s,'UU') returns 1 on word UUU
  #since it apparently counts only nonoverlapping words UU
  #For this reason, we work with the indices.

  #Initialize lists and mono- and dinucleotide dictionaries
  List = {} #List is a dictionary of lists
  List['A'] = []; List['C'] = [];
  List['G'] = []; List['U'] = [];
  nuclList   = ["A","C","G","U"]
  s       = s.upper()
  s       = s.replace("T","U")
  nuclCnt    = {}  #empty dictionary
  dinuclCnt  = {}  #empty dictionary
  for x in nuclList:
    nuclCnt[x]=0
    dinuclCnt[x]={}
    for y in nuclList:
      dinuclCnt[x][y]=0

  #Compute count and lists
  nuclCnt[s[0]] = 1
  nuclTotal     = 1
  dinuclTotal   = 0
  for i in range(len(s)-1):
    x = s[i]; y = s[i+1]
    List[x].append( y )
    nuclCnt[y] += 1; nuclTotal  += 1
    dinuclCnt[x][y] += 1; dinuclTotal += 1
  assert (nuclTotal==len(s))
  assert (dinuclTotal==len(s)-1)
  return nuclCnt,dinuclCnt,List

def chooseEdge(x,dinuclCnt):
  numInList = 0
  for y in ['A','C','G','U']:
    numInList += dinuclCnt[x][y]
  z = random.random()
  denom=dinuclCnt[x]['A']+dinuclCnt[x]['C']+dinuclCnt[x]['G']+dinuclCnt[x]['U']
  numerator = dinuclCnt[x]['A']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['A'] -= 1
    return 'A'
  numerator += dinuclCnt[x]['C']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['C'] -= 1
    return 'C'
  numerator += dinuclCnt[x]['G']
  if z < float(numerator)/float(denom):
    dinuclCnt[x]['G'] -= 1
    return 'G'
  dinuclCnt[x]['U'] -= 1
  return 'U'

def connectedToLast(edgeList,nuclList,lastCh):
  D = {}
  for x in nuclList: D[x]=0
  for edge in edgeList:
    a = edge[0]; b = edge[1]
    if b==lastCh: D[a]=1
  for i in range(2):
    for edge in edgeList:
      a = edge[0]; b = edge[1]
      if D[b]==1: D[a]=1
  ok = 0
  for x in nuclList:
    if x!=lastCh and D[x]==0: return 0
  return 1

def eulerian(s):
  nuclCnt,dinuclCnt,List = computeCountAndLists(s)
  #compute nucleotides appearing in s
  nuclList = []
  for x in ["A","C","G","U"]:
    if x in s: nuclList.append(x)
  #compute numInList[x] = number of dinucleotides beginning with x
  numInList = {}
  for x in nuclList:
    numInList[x]=0
    for y in nuclList:
      numInList[x] += dinuclCnt[x][y]
  #create dinucleotide shuffle L
  firstCh = s[0]  #start with first letter of s
  lastCh  = s[-1]
  edgeList = []
  for x in nuclList:
    if x!= lastCh: edgeList.append( [x,chooseEdge(x,dinuclCnt)] )
  ok = connectedToLast(edgeList,nuclList,lastCh)
  return ok,edgeList,nuclList,lastCh

def shuffleEdgeList(L):
  n = len(L); barrier = n
  for i in range(n-1):
    z = int(random.random() * barrier)
    tmp = L[z]
    L[z]= L[barrier-1]
    L[barrier-1] = tmp
    barrier -= 1
  return L

def dinuclShuffle(s):
  ok = 0
  while not ok:
    ok,edgeList,nuclList,lastCh = eulerian(s)
  nuclCnt,dinuclCnt,List = computeCountAndLists(s)

  #remove last edges from each vertex list, shuffle, then add back
  #the removed edges at end of vertex lists.
  for [x,y] in edgeList: List[x].remove(y)
  for x in nuclList: shuffleEdgeList(List[x])
  for [x,y] in edgeList: List[x].append(y)

  #construct the eulerian path
  L = [s[0]]; prevCh = s[0]
  for i in range(len(s)-2):
    ch = List[prevCh][0]
    L.append( ch )
    del List[prevCh][0]
    prevCh = ch
  L.append(s[-1])
 # print(L)
  t = "".join(L)
  return t

#### Defining my functions #####
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])


###### Function to calculate ZScore on list of MFEs #################
def pscore_function(energy_list, randomizations):
    below_native = 0
    total_count = len(energy_list)
    native_mfe = float(energy_list[0])
    #scrambled_mean_mfe = np.mean(energy_list[1:randomizations])
    for MFE in energy_list:
        if float(MFE) < float(native_mfe):
            below_native += 1

    pscore = float(float(below_native) / float(total_count))

    return pscore;

###### Function to calculate ZScore on list of MFEs #################
def zscore_function(energy_list, randomizations):
    mean = np.mean(energy_list)
    sd = np.std(energy_list)
    native_mfe = energy_list[0]
    scrambled_mean_mfe = np.mean(energy_list[1:randomizations])
    #scrambled_sd = np.std(energy_list[1:randomizations])
    if sd != 0:
        zscore = (native_mfe - scrambled_mean_mfe)/sd
    if sd == 0:
        zscore = "#DIV/0!"
    return zscore;

def rna_folder(frag):
    (structure, MFE) = RNA.fold(str(frag))
    return MFE;

def randomizer(frag):
    result = ''.join(random.sample(frag,len(frag)))
    return result;

###### Function to calculate MFEs using RNAfold #################
def energies(seq_list):
    energy_list = []

    energy_list = multiprocessing(rna_folder, [sequence for sequence in seq_list], 12)
    # for sequence in seq_list:
    #     #fc = RNA.fold_compound(str(sequence))
    #     (structure, MFE) = RNA.fold(str(sequence)) # calculate and define variables for mfe and structure
    #     energy_list.append(MFE) # adds the native fragment to list

    return energy_list;

######Function to create X number of scrambled RNAs in list #################
#test
def scramble(text, randomizations, type):
    frag = str(text)
    frag_seqs = []
    if type == "di":
        for _ in range(randomizations):
            result = dinuclShuffle(frag)
            frag_seqs.append(result)
    elif type == "mono":
        frag_seqs = multiprocessing(randomizer, [frag for i in range(randomizations)], 12)

        # for _ in range(int(randomizations)):
        #     result = ''.join(random.sample(frag,len(frag)))
        #     frag_seqs.append(result)
    else:
        print("Shuffle type not properly designated; please input \"di\" or \"mono\"")

    return frag_seqs;

##################### Main Script #########################################
#w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
start_nucleotide_28 = []
end_nucleotide_28 = []
temp_28 = "28"
MFE_total_28 = []
ED_total_28 = []
frag_28 = []
structure_28 = []
centroid_28 = []





with open(myfasta, 'r') as forward_fasta:

    for cur_record in SeqIO.parse(forward_fasta, "fasta") :

            read_name = cur_record.name
            full_seq_length = len(cur_record.seq)
            print("Scanning sequence "+str(read_name)+"\nSequence Length: "+str(len(cur_record.seq))+"nt long.")

            #### this will change based on input fasta file header format #########
            #print(read.name)
            #fasta_header = read_name.split('|')
            #print(fasta_header)
            #gene_id = fasta_header[0]
            #transcript_id = fasta_header[1]
            #chromosome = "chr"+fasta_header[2]
            #gene_start = fasta_header[3]
            #gene_end = fasta_header[4]
            #strand = fasta_header[5]
        ##### Establish empty lists to capture calculated metrics per window ######
            zscore_total = []
            numerical_z = []
            pscore_total = []
            numerical_p = []
            MFE_total = []
            ED_total = []


            ##### Load constraints #####
            constraint_dict = {}
            constraint_list = []
            if args.constraints != None:
                print("Considering constraint input")
                constraint_file = open(args.constraints, "r")
                constraints = constraint_file.readlines()[2]

                #print(constraints)
                i = 1
                for nt in constraints:
                    #print(i, nt)
                    constraint_list.append(nt)
                    constraint_dict[i] = nt
                    i += 1


                if len(constraint_list)-1 == len(cur_record.seq):
                    print("Constraint list is "+str(len(constraint_list)-1)+"nt long.")
                else:
                    print("Constraint list is "+str(len(constraint_list)-1)+"nt long.")
                    raise ValueError("Error detected. Sequence and Constraints must be same length.")


            #gff3file = open(read.name+'.gff3', 'w')
            #pscore_wig = open(read.name+'.pscore.wig', 'w')
            #zscore_wig = open(read.name+".zscore.wig", 'w')
            #ED_wig = open(read.name+".ED.wig", 'w')
            #MFE_wig = open(read.name+".MFE.wig", 'w')
            #print(read.name, read.sequence)
            length = len(cur_record.seq)

            ##### This will ignore any sequences in fasta file which are smaller than window size #####
            if length < int(window_size):
                continue
            seq = cur_record.seq

            ##### Write the header to the output file #####
            # w.write("i\tj\tTemperature\tNative_dG\tZ-score\tP-score\tEnsembleDiversity\tSequence\tStructure\tCentroid\t"+read_name+"\n")

    ##### Main routine using defined functions: ##########################################


            i = 0 #Important to establish i=0 here...
            while i == 0 or i <= (length - window_size):
                start_nucleotide = i + 1 # This will just define the start nucleotide coordinate value
                frag = seq[i:i+int(window_size)] # This breaks up sequence into fragments
                #print(frag)
                #print(str(len(frag)))
                start_nucleotide = i + 1
                end_nucleotide = i + window_size
                if -1 == 0:
                    print("Magic")
                # if 'N' in frag:
                #     w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str(frag)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\n")
                #     i += step_size #this ensures that the next iteration increases by "step size" length
                else:
                    #print(start_nucleotide)
                    #print(end_nucleotide)
                    frag = frag.transcribe()

                    ##### Dealing with completely ambigious (All "N") sequence #####
                    if frag == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
                        MFE = int(0.0)
                        zscore = "#DIV/0"
                        ED = int(0.0)
                        pscore = int(0.0)
                        structure = "........................................................................................................................"
                        centroid = "........................................................................................................................"
                    else:
                        fc = RNA.fold_compound(str(frag), md) #creates "Fold Compound" object using folding model
                        fc.pf() # performs partition function calculations
                        frag_q = (RNA.pf_fold(str(frag)), md) # calculate partition function "fold" of fragment
                        (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                        MFE = round(MFE, 2)
                        MFE_total.append(MFE)
                        (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                        ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                        ED_total.append(ED)            #print(structure)
                        #fmfe = fc.pbacktrack()
                        #print(str(fmfe))
                        if constraints == None:
                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            fc.pf()# performs partition function calculations
                            #frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment

                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)
                            #fmfe = fc.pbacktrack()
                            #print(str(fmfe))
                        elif constraints != None:
                            window_constraint_list = constraint_list[start_nucleotide-1:end_nucleotide]
                            window_constraints = ''.join(window_constraint_list)
                            #print(frag)
                            #print(len(window_constraints))
                            fc.hc_add_from_db(window_constraints)
                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            fc.pf()# performs partition function calculations
                            #frag_q = fc.pf_fold() # calculate partition function "fold" of fragment
                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)
                        
                        # seqlist = [] # creates the list we will be filling with sequence fragments
                        # seqlist.append(frag) # adds the native fragment to list
                        # scrambled_sequences = scramble(frag, randomizations, type)
                        # seqlist.extend(scrambled_sequences)
                        # energy_list = energies(seqlist)
                        # if print_random == "on":
                        #     print(energy_list)


                        # try:
                        #     zscore = round(zscore_function(energy_list, randomizations), 2)
                        # except:
                        #     zscore = zscore_function(energy_list, randomizations)
                        # zscore_total.append(zscore)

                        # #print(zscore)
                        # pscore = round(pscore_function(energy_list, randomizations), 2)
                        # #print(pscore)
                        # pscore_total.append(pscore)



                    if print_to_screen == True:
                        if constraints != None:
                            print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(ED)+"\n"+str(frag)+"\n"+str(window_constraints)+"\n"+str(structure)+"\n"+str(centroid)+"\n")
                        else:
                            print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(ED)+"\n"+str(frag)+"\n"+str(structure)+"\n"+str(centroid)+"\n")
                    
                    
                    start_nucleotide_28.append(start_nucleotide)
                    end_nucleotide_28.append(end_nucleotide)
                    MFE_total_28.append(MFE)
                    ED_total_28.append(ED)
                    frag_28.append(frag)
                    structure_28.append(structure)
                    centroid_28.append(centroid)


                    #w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
                    #gff3file.write()
                    #pscore_wig.write()
                    #zscore_wig.write()
                    #ED_wig.write()
                    #MFE_wig.write()

                    i += step_size #this ensures that the next iteration increases by "step size" length

            #print(len(zscore_total))



            # for z in zscore_total:
            #     try:
            #         numerical_z.append(float(z))
            #     except ValueError:
            #         continue
            # #print(len(numerical_z))

            # #print(len(pscore_total))
            # for p in pscore_total:
            #     try:
            #         numerical_p.append(float(p))
            #     except ValueError:
            #         continue
            # window_count = len(zscore_total)



            #print(len(numerical_p))
            #print(window_count)
            #print(type(window_count))
            #print(step_size)
            #print(type(step_size))
            #print(length)
            #print(type(length))
            #coverage = round((float(window_count)*float(step_size))/float(length), 2)
            #print(coverage)
            #print(len(MFE_total))
            #print(len(ED_total))

            # mean_pscore = round(np.mean(numerical_p), 2)
            # mean_zscore = round(np.mean(numerical_z), 2)
            
            #mean_MFE = round(np.mean(MFE_total), 2)
            #mean_ED = round(np.mean(ED_total), 2)
            #w.write("---\t---\t---\t---\t---\t---\t---\t---\tSummary:\tLength\tMeanMFE\tMeanZ\tMeanPscore\tMeanED\n---\t---\t---\t---\t---\t---\t---\t---\t---\t"+str(length)+"\t"+str(mean_MFE)+"\t"+str(mean_zscore)+"\t"+str(mean_pscore)+"\t"+str(mean_ED)+"\n\n")
            #s.write(str(read_name)+"\t"+str(length)+"\t"+str(mean_MFE)+"\t"+str(mean_zscore)+"\t"+str(mean_pscore)+"\t"+str(mean_ED)+"\n")

# for i in range(len(start_nucleotide_28)):

#     # start_nucleotide_28.append(start_nucleotide)
#     # end_nucleotide_28.append(end_nucleotide)
#     # MFE_total_28.append(MFE)
#     # ED_total_28.append(ED)
#     # frag_28.append(frag)
#     # structure_28.append(structure)
#     # centroid_28.append(centroid)

#     w.write(str(start_nucleotide_28[i])+"\t"+str(end_nucleotide_28[i])+"\t"+temp_28+"\t"+str(MFE_total_28[i])+"\t"+str(ED_total_28[i])+"\t"+str(frag_28[i])+"\t"+str(structure_28[i])+"\t"+str(centroid_28[i])+"\n")

# w.close


print("DONE WITH 28C")



MFE_total_42 = []
ED_total_42 = []
structure_42 = []
centroid_42 = []

temperature = 42
##### Establish global "folding model" incorporting temperature (or other settings) #####
md = RNA.md()
md.temperature = int(temperature)

# w = open(myfasta+".forward.win_"+str(window_size)+".stp_"+str(step_size)+".rnd_"+str(randomizations)+".shfl_"+str(type)+"42C.txt", 'w')
##################### Main Script #########################################
with open(myfasta, 'r') as forward_fasta:

    for cur_record in SeqIO.parse(forward_fasta, "fasta") :

            read_name = cur_record.name
            full_seq_length = len(cur_record.seq)
            print("Scanning sequence "+str(read_name)+"\nSequence Length: "+str(len(cur_record.seq))+"nt long.")

            #### this will change based on input fasta file header format #########
            #print(read.name)
            #fasta_header = read_name.split('|')
            #print(fasta_header)
            #gene_id = fasta_header[0]
            #transcript_id = fasta_header[1]
            #chromosome = "chr"+fasta_header[2]
            #gene_start = fasta_header[3]
            #gene_end = fasta_header[4]
            #strand = fasta_header[5]
        ##### Establish empty lists to capture calculated metrics per window ######
            zscore_total = []
            numerical_z = []
            pscore_total = []
            numerical_p = []
            MFE_total = []
            ED_total = []


            ##### Load constraints #####
            constraint_dict = {}
            constraint_list = []
            if args.constraints != None:
                print("Considering constraint input")
                constraint_file = open(args.constraints, "r")
                constraints = constraint_file.readlines()[2]

                #print(constraints)
                i = 1
                for nt in constraints:
                    #print(i, nt)
                    constraint_list.append(nt)
                    constraint_dict[i] = nt
                    i += 1


                if len(constraint_list)-1 == len(cur_record.seq):
                    print("Constraint list is "+str(len(constraint_list)-1)+"nt long.")
                else:
                    print("Constraint list is "+str(len(constraint_list)-1)+"nt long.")
                    raise ValueError("Error detected. Sequence and Constraints must be same length.")


            #gff3file = open(read.name+'.gff3', 'w')
            #pscore_wig = open(read.name+'.pscore.wig', 'w')
            #zscore_wig = open(read.name+".zscore.wig", 'w')
            #ED_wig = open(read.name+".ED.wig", 'w')
            #MFE_wig = open(read.name+".MFE.wig", 'w')
            #print(read.name, read.sequence)
            length = len(cur_record.seq)

            ##### This will ignore any sequences in fasta file which are smaller than window size #####
            if length < int(window_size):
                continue
            seq = cur_record.seq

            ##### Write the header to the output file #####
            # w.write("i\tj\tTemperature\tNative_dG\tZ-score\tP-score\tEnsembleDiversity\tSequence\tStructure\tCentroid\t"+read_name+"\n")

    ##### Main routine using defined functions: ##########################################


            i = 0 #Important to establish i=0 here...
            while i == 0 or i <= (length - window_size):
                start_nucleotide = i + 1 # This will just define the start nucleotide coordinate value
                frag = seq[i:i+int(window_size)] # This breaks up sequence into fragments
                #print(frag)
                #print(str(len(frag)))
                start_nucleotide = i + 1
                end_nucleotide = i + window_size
                if -1 == 0:
                    print("Magic")
                # if 'N' in frag:
                #     w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\t"+str(frag)+"\t"+str("Not Available - N in fragment")+"\t"+str("Not Available - N in fragment")+"\n")
                #     i += step_size #this ensures that the next iteration increases by "step size" length
                else:
                    #print(start_nucleotide)
                    #print(end_nucleotide)
                    frag = frag.transcribe()

                    ##### Dealing with completely ambigious (All "N") sequence #####
                    if frag == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
                        MFE = int(0.0)
                        zscore = "#DIV/0"
                        ED = int(0.0)
                        pscore = int(0.0)
                        structure = "........................................................................................................................"
                        centroid = "........................................................................................................................"
                    else:
                        fc = RNA.fold_compound(str(frag), md) #creates "Fold Compound" object using folding model
                        fc.pf() # performs partition function calculations
                        frag_q = (RNA.pf_fold(str(frag)), md) # calculate partition function "fold" of fragment
                        (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                        MFE = round(MFE, 2)
                        MFE_total.append(MFE)
                        (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                        ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                        ED_total.append(ED)            #print(structure)
                        #fmfe = fc.pbacktrack()
                        #print(str(fmfe))
                        if constraints == None:
                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            fc.pf()# performs partition function calculations
                            #frag_q = (RNA.pf_fold(str(frag))) # calculate partition function "fold" of fragment

                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)
                            #fmfe = fc.pbacktrack()
                            #print(str(fmfe))
                        elif constraints != None:
                            window_constraint_list = constraint_list[start_nucleotide-1:end_nucleotide]
                            window_constraints = ''.join(window_constraint_list)
                            #print(frag)
                            #print(len(window_constraints))
                            fc.hc_add_from_db(window_constraints)
                            (structure, MFE) = fc.mfe() # calculate and define variables for mfe and structure
                            fc.pf()# performs partition function calculations
                            #frag_q = fc.pf_fold() # calculate partition function "fold" of fragment
                            MFE = round(MFE, 2)
                            MFE_total.append(MFE)
                            (centroid, distance) = fc.centroid() # calculate and define variables for centroid
                            ED = round(fc.mean_bp_distance(), 2) # this caclulates ED based on last calculated partition funciton
                            ED_total.append(ED)            #print(structure)
                        
                        # seqlist = [] # creates the list we will be filling with sequence fragments
                        # seqlist.append(frag) # adds the native fragment to list
                        # scrambled_sequences = scramble(frag, randomizations, type)
                        # seqlist.extend(scrambled_sequences)
                        # energy_list = energies(seqlist)
                        # if print_random == "on":
                        #     print(energy_list)


                        # try:
                        #     zscore = round(zscore_function(energy_list, randomizations), 2)
                        # except:
                        #     zscore = zscore_function(energy_list, randomizations)
                        # zscore_total.append(zscore)

                        # #print(zscore)
                        # pscore = round(pscore_function(energy_list, randomizations), 2)
                        # #print(pscore)
                        # pscore_total.append(pscore)



                    if print_to_screen == True:
                        if constraints != None:
                            print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(ED)+"\n"+str(frag)+"\n"+str(window_constraints)+"\n"+str(structure)+"\n"+str(centroid)+"\n")
                        else:
                            print(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(ED)+"\n"+str(frag)+"\n"+str(structure)+"\n"+str(centroid)+"\n")
                    
                    

                    MFE_total_42.append(MFE)
                    ED_total_42.append(ED)
                    structure_42.append(structure)
                    centroid_42.append(centroid)
                    
                    #w.write(str(start_nucleotide)+"\t"+str(end_nucleotide)+"\t"+str(temperature)+"\t"+str(MFE)+"\t"+str(ED)+"\t"+str(frag)+"\t"+str(structure)+"\t"+str(centroid)+"\n")
                    #gff3file.write()
                    #pscore_wig.write()
                    #zscore_wig.write()
                    #ED_wig.write()
                    #MFE_wig.write()

                    i += step_size #this ensures that the next iteration increases by "step size" length

            #print(len(zscore_total))



            # for z in zscore_total:
            #     try:
            #         numerical_z.append(float(z))
            #     except ValueError:
            #         continue
            # #print(len(numerical_z))

            # #print(len(pscore_total))
            # for p in pscore_total:
            #     try:
            #         numerical_p.append(float(p))
            #     except ValueError:
            #         continue
            # window_count = len(zscore_total)



            #print(len(numerical_p))
            #print(window_count)
            #print(type(window_count))
            #print(step_size)
            #print(type(step_size))
            #print(length)
            #print(type(length))
            #coverage = round((float(window_count)*float(step_size))/float(length), 2)
            #print(coverage)
            #print(len(MFE_total))
            #print(len(ED_total))

            # mean_pscore = round(np.mean(numerical_p), 2)
            # mean_zscore = round(np.mean(numerical_z), 2)
            
            mean_MFE = round(np.mean(MFE_total), 2)
            mean_ED = round(np.mean(ED_total), 2)
            #w.write("---\t---\t---\t---\t---\t---\t---\t---\tSummary:\tLength\tMeanMFE\tMeanZ\tMeanPscore\tMeanED\n---\t---\t---\t---\t---\t---\t---\t---\t---\t"+str(length)+"\t"+str(mean_MFE)+"\t"+str(mean_zscore)+"\t"+str(mean_pscore)+"\t"+str(mean_ED)+"\n\n")
            #s.write(str(read_name)+"\t"+str(length)+"\t"+str(mean_MFE)+"\t"+str(mean_zscore)+"\t"+str(mean_pscore)+"\t"+str(mean_ED)+"\n")


# for i in range(len(start_nucleotide_42)):

#     # start_nucleotide_28.append(start_nucleotide)
#     # end_nucleotide_28.append(end_nucleotide)
#     # MFE_total_28.append(MFE)
#     # ED_total_28.append(ED)
#     # frag_28.append(frag)
#     # structure_28.append(structure)
#     # centroid_28.append(centroid)

#     w.write(str(start_nucleotide_42[i])+"\t"+str(end_nucleotide_42[i])+"\t"+temp_42+"\t"+str(MFE_total_42[i])+"\t"+str(ED_total_42[i])+"\t"+str(frag_42[i])+"\t"+str(structure_42[i])+"\t"+str(centroid_42[i])+"\n")

# w.close

w = open(myfasta+".forward.win_"+str(window_size)+".stp_"+str(step_size)+".rnd_"+str(randomizations)+".shfl_"+str(type)+"42C-28CDelta.txt", 'w')

##### Write the header to the output file #####
w.write("i\tj\tNative_dG_28\tNative_dG_42\tDelta_Native_dG_42-28\tEnsembleDiversity_28\tEnsembleDiversity_42\tDelta_EnsembleDiversity_42-28\tSequence\tMFE_Structure_28\tMFE_Structure_42\tCentroid_28\tCentroid_42\t"+read_name+"\n")
for i in range(len(start_nucleotide_28)):
   w.write(str(start_nucleotide_28[i])+"\t"+str(end_nucleotide_28[i])+"\t"+str(MFE_total_28[i])+"\t"+str(MFE_total_42[i])+"\t"+str(MFE_total_42[i]-MFE_total_28[i])+"\t"+str(ED_total_28[i])+"\t"+str(ED_total_42[i])+"\t"+str(ED_total_42[i]-ED_total_28[i])+"\t"+str(frag_28[i])+"\t"+str(structure_28[i])+"\t"+str(structure_28[i])+"\t"+str(centroid_28[i])+"\t"+str(centroid_42[i])+"\n")

w.close()