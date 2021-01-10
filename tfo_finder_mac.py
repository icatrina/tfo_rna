#!/anaconda/bin/python3
import csv
import sys
import pandas as pd
import os
from Bio.SeqUtils import MeltingTemp as mt
import subprocess
import re
from Bio.Blast import NCBIXML
print("\n"*5)
print('TFOfinder  Copyright (C) 2017  Irina E. Catrina' + '\n'+
'This program comes with ABSOLUTELY NO WARRANTY;'  + '\n' +
'This is free software, and you are welcome to redistribute it' + '\n'+
'under certain conditions; for details please read the GNU_GPL.txt file included with PinMol' + '\n')
undscr = "->"*30
print(undscr)
print("\n"+"WARNING: The input sscount file and any previous output files will be overwritten!  Save them in a "+"\n"+"different location than the current file, or rename them to "+"\n"+"ensure they are not misused (e.g. use probes from a different target)."+"\n")
print(undscr)

def readSScountFile(filename): #sscount file for molecular beacon design
    userpath = os.path.dirname(filename) #use the path of input to save all files
    with open(filename) as infile: #check if the input file is txt and if Us are present not Ts
        if 'T' in infile.read():
            print('Wrong type of file (not txt) OR you may have used DNA parameters to fold your target RNA!')
            sys.exit('Try again using RNA parameters!')
    baselines = []
    with open(filename) as infile: #clean-up the file remove whitespaces from left side
        for line in infile:
            line = line.lstrip(" ")
            newlines = re.sub(' +', ' ', line)
            if not line:
                continue
            baselines.append(newlines.upper())
    so = int(baselines[0]) #number of suboptimal structures
    return (so, userpath, baselines)

def seqTarget(f): #sequence of target & sscount for each probe as fraction (1 for fully single stranded)
    bl = [[],[],[]]
    reader = csv.reader(f)
    for row in reader:
        for col in range(3):
            bl[col].append(row[col])
    sscount = bl[1]
    position = bl[0]
    max_base = bl[0][-1]
    bases = bl[2]
    seq = ''.join(bases)
    size = len(seq)
    return(sscount, position, max_base, bases, seq, size)

def itersplit_into_x_chunks(argum, size, chunksize): #split sequence in chunks of probe size
    for pos in range(0, size-chunksize+1):
        yield argum[pos:pos+chunksize]

def reverse_complement(seq): #generate RNA complement
    return seq.translate(basecomplement)[::-1]

def probeLength(probe): #input desired probe length; limited to range[6,11)]
    if probe <6 or probe >13:
        print('The value you entered is incorrect!')
        sys.exit('Try again!')
    else:
        probe = probe
    return (probe)

def seqProbes(mb_seq, mb_size, mb_sscount, probe):
    result = list(itersplit_into_x_chunks(mb_seq, mb_size, probe))
    basesl = []
    for i in result:
        i = reverse_complement(i)
        basesl.append(i)

    Tml = []
    for i in basesl:
        Tmx = mt.Tm_NN(i, dnac1 = 50000, dnac2 = 50000, Na = 100, nn_table = mt.RNA_NN1, saltcorr = 1)
        Tml.append(int(Tmx))
    result_bases = list(itersplit_into_x_chunks(mb_bases, mb_size, probe)) #list of lists of each base for each probe
    #base number as j and list of these numbers as jl, list of percent of Gs and Cs as perl
    j = 0
    perl = []
    jl = []
    for i in result_bases:
        j += 1
        aas = i.count('A')
        gs = i.count('G')
        per = int((aas+gs)/probe*100)
        perl.append(per)
        jl.append(j)
    size2=len(mb_sscount)
    result2 = list(itersplit_into_x_chunks(mb_sscount, size2, probe))
    sumsl = []
    for i in result2:
        i = list(map(int, i))
        sums = sum(i)/(probe*mb_so)
        sumsl.append(sums)
    return (jl, perl, sumsl, basesl, Tml) #put together all data as indicated in header

def regionTarget(tg_start, tg_end): #if only a region of the target needs to be considered
    tg_diff = tg_end - tg_start
    assert(tg_start>0 and tg_end>0), "The base numbers should be positive and larger then zero!"
    if tg_diff < 0: #consider only probes within the requested region of target RNA, but check that the program can finish
        print("The number of end base cannot be smaller than the number of the start base!")
        sys.exit('Try again!')
    elif tg_diff < probe:
        print("You have to enter a region with a size larger than the probe size!")
        sys.exit('Try again!')
    elif tg_start >=1 and tg_end <= int(mb_max_base):
        df = pd.read_csv(mb_userpath+'/all_probes.csv')
        tgstart = tg_start - 1
        tgend = tg_end - probe + 2
        slice2 = df[tgstart:tgend]
        for row in slice2:
            slice2 = slice2.sort_values(by='sscount', ascending=[True]) #sort ascending by sscount = larger sscount more accessible target region
        return(slice2)

    elif no_pb > 1 and no_pb <= row_no and no_pb <= 50:
        no_pb = no_pb
        input1 = open(mb_userpath+'/probes_sortedby5.csv', 'r').read().split('\n')
        outputData = input1[:no_pb+1]
        output = open(mb_userpath +'/DG_probes.csv', 'w')
        output.write('\n'.join(outputData))
        output.close()
    return(no_pb)

def stemDesign(): #design the stem of the final probes
    i = 0
    with open(mb_userpath+'/mb_picks.csv') as ff:
        tw =[[],[]]
        reader = csv.reader(ff)
        for row in reader:
            for col in range(2):
                tw[col].append(row[col])
        bs_ps =tw[0]
        fseq =tw[1]
        seq_slc = []

        for item in fseq:
            i += 1
            seql = list(item)
            aseq = (str(i) + " TFO sequence at base number " + bs_ps[i-1] + ' is:  '+ item)
            print(aseq + '\n')
            with open (mb_userpath+'/Final_TFOs.csv', 'a') as outputf:
                outputf.write(aseq+'\n')
            with open(mb_userpath+'/Seq'+str(i)+'.seq', 'w') as fiseq:
                fiseq.write(';'+'\n'+ str(i) + ' at base # ' + bs_ps[i-1] + ' TFO' + '\n' + item.strip() +'1')
        return(i)

if __name__ == "__main__":
    filename = input('Enter a file name: ') # request ss-count file path and name
    mb_so, mb_userpath, mb_baselines = readSScountFile(filename)

    while True: #probe length?
        try:
            probe=int(input("Enter the length of probe; a number between 6 and 13: "))
        except:
            print('You must type a number between 6 and 13, try again:')
            continue

        else:
            probe = probeLength(probe)
            break

    with open (filename, 'w') as outfile:
        for line in mb_baselines:
            outfile.write(line)

    with open(filename, 'r') as infile, open(mb_userpath+'/sscounttxt_tocsv.csv', 'w') as csv_file:
        next(infile)
        reader = csv.reader(infile, delimiter = ' ')
        writer = csv.writer(csv_file, delimiter = ",", lineterminator = '\n')
        writer.writerows(reader)


    with open(mb_userpath+'/sscounttxt_tocsv.csv', 'r') as f:
        mb_sscount, mb_position, mb_max_base, mb_bases, mb_seq, mb_size = seqTarget(f)

    os.remove(mb_userpath+'/sscounttxt_tocsv.csv')

    basecomplement = str.maketrans({'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'})
    mb_jl, mb_perl, mb_sumsl, mb_basesl, mb_Tml = seqProbes(mb_seq, mb_size, mb_sscount, probe)

    with open(mb_userpath+'/all_probes.csv', 'w') as csv_file:
        writer = csv.writer(csv_file, lineterminator='\n')
        writer.writerow(["Base number", "%GA", "sscount", "Probe sequence", "Tm"])
        rows = zip(mb_jl,mb_perl,mb_sumsl,mb_basesl,mb_Tml)
        for row in rows:
            writer.writerow(row)

    tg_start = int(input("If a specific region within the target is needed, please enter the number of start base, or 1: "))
    tg_end = int(input("  and the number of end base or max number of bases " + str(mb_max_base) + ": "))
    mb_slice2 = regionTarget(tg_start, tg_end)
    mb_slice2.to_csv(mb_userpath+'/all_probes_sorted_ss.csv', index=False)

    os.remove(mb_userpath+'/all_probes.csv')

    with open(mb_userpath+'/all_probes_sorted_ss.csv','r') as fin, open(mb_userpath+'/GA_probes.csv', 'w') as fout:
        next(fin)
        reader = csv.reader(fin)
        writer = csv.writer(fout, lineterminator='\n')
        writer.writerow(["Base number", "%GA", "sscount", "Probe sequence", "Tm"])
        for row in reader:
            if int(row[1]) > 70 and int(row[1]) < 101:
                writer.writerow(row)


    #sort by sscount?
    flistsort = pd.read_csv(mb_userpath+'/GA_probes.csv', sep=',', usecols=[0,1,2,3,4])

    #new file with only sequences of probes for calculating free energies using oligoscreen
    flistsort2 = pd.read_csv(mb_userpath+'/GA_probes.csv', usecols=[3], skiprows=1)
    flistsort2.to_csv(mb_userpath+'/oligoscreen_input.lis', index=False)

    with open (mb_userpath+'/GA_probes.csv', 'r') as GA_probes, open(mb_userpath+'/all_probes_sorted_ss.csv', 'r') as all_probes:
        no_GAprobes = sum(1 for row in GA_probes)
        no_probes = sum(1 for row in all_probes)

    subprocess.check_output(["oligoscreen", mb_userpath+'/oligoscreen_input.lis', mb_userpath+'/oligoscreen_output.csv'])
    read_oligosc = pd.read_csv(mb_userpath+'/oligoscreen_output.csv', delimiter = '\t', usecols=[1,2,3])

    os.remove(mb_userpath+'/oligoscreen_input.lis')
    os.remove(mb_userpath+'/oligoscreen_output.csv')

    #keep only probes that meet the energy requirements and sort them
    data_comb = pd.concat([flistsort, read_oligosc], axis=1)
    data_filter = data_comb[(data_comb.DGbimolecular > -7.5) & (data_comb.DGunimolecular > -2.5)]
    for row in data_filter:
        data_sort = data_filter.sort_values(['sscount','DGunimolecular', 'DGbimolecular', '%GA', 'DGduplex'], ascending=[True, False, False, True, True]) #sort descending by sscount = larger sscount more accessible target region
        data_sort.to_csv(mb_userpath+'/probes_sortedby5.csv', index=False)

    #determine the total number of probes that meet the eg criteria for the selected target (region or full)
    with open(mb_userpath+'/probes_sortedby5.csv', 'r') as flistsort3:
        reader = csv.reader(flistsort3)
        row_no = sum(1 for row in reader)-1
        print("Maximum number of possible probes is: "+str(row_no)+"\n")
 
    read_sscnt = pd.read_csv(mb_userpath+'/all_probes_sorted_ss.csv', delimiter = ',')
    no_ss = (read_sscnt["sscount"] < 0.1).sum()
    #os.remove(mb_userpath+'/usetxtfile.txt')

    with open (mb_userpath+'/Final_TFOs.csv', 'a') as add_output:
        add_output.write('Results for ' + '\"'+filename +'\"'+ ' using ' + str(probe) + ' as probe length, and for a target region between ' +
        str(tg_start) +' and ' + str(tg_end) + ' nucleotides:  ' + '\n' +
        '\t'+'1. Total number of possible probes =  ' + str(no_probes-1)+ '\n'+
        '\t'+'2. Number of probes that meet GA > 70% and energetic criteria =  ' + str(row_no) + '\n'
        '\t'+'3. Number of probes that have an ss-count fraction smaller than 0.5 =  ' + str(no_ss)+ '\n')
    print('Results for ' + '\"'+filename +'\"'+ ' using ' + str(probe) + ' as probe length, and for a target region between ' +
        str(tg_start) +' and ' + str(tg_end) + ' nucleotides:  ' + '\n')
    print('\t'+'1. Total number of possible probes =  ' + str(no_probes-1)+ '\n')
    print('\t'+'2. Number of probes that have a GA content > 70% =  '+ str(no_GAprobes-1)+ '\n' )
    print('\t'+'3. Number of probes that have an ss-count fraction smaller than 0.5 =  ' + str(no_ss)+ '\n')

    print("\n"+"This information can be also be found in the file Final_TFOs.csv"+"\n")
    print("\n"+"To select the best TFO, check GA_probes.csv file for the smallest sscount and largest GA%, note that the GA% is the purine content in the target not in the probe sequence."+"\n")
    #print("\n"+"Finally, to avoid false positives you need to check for mostly single stranded sequences complementary to the selected probe(s) in other mRNAs expressed in the tissue of interest."+"\n")
