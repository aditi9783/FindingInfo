#!/usr/bin/python

import getopt, sys
import math
import numpy as np

fname = None
w = 10 # window size (on both sides of a position, thus net window is 2w+1) for computing linguistic complexity and shannon's entropy

# start of readUserGenomeSeq #######
def readUserGenomeSeq(fname):
    fh = open(fname, 'r')
    fh.readline() # remove fasta header
    seq = ""
    for line in fh:
        seq = seq+line.rstrip("\n")
    fh.close()
    return seq
# end of readUserGenomeSeq #######

# start of addbase ############
def addbase(mers, alphabet):
    newmers = []
    for m in mers:
        for a in alphabet:
            newmers.append(m+a)
    return newmers
# end of addbase ############

# start of getallmers #######
def getallmers(i): # get all possible oligomers of size i
    alphabet = ["A", "T", "C", "G"]
    i_mers = alphabet[:] # all mers start with the single base
    n = i-1 # i is len of mers, and we already added 1 base (alphabet), thus n is the remaining len of mer left
    while n != 0:
        i_mers = addbase(i_mers, alphabet)
        n -= 1
    return i_mers 
# end of getallmers #######

# start of getseqmers #########
def getseqmers(i, seqarr): # get all mers of size i present in the seq
    seqmerdict = {} # key: unique mer present in seq, val: num of times the mer occured in seq
    for idx in range(0, len(seqarr)-i+1):
        mer = "".join(seqarr[idx:idx+i])
        if mer in seqmerdict: # this mer has been seen before
            seqmerdict[mer] += 1
        else:
            seqmerdict[mer] = 1 # initialize count for this mer
    return seqmerdict # return dict of unique mers seen in the seq
# end of getseqmers #########

# start of getallmersize ########
def getallmersize(i):
    alphabet = ["A", "T", "C", "G"]
    asize = len(alphabet)
    n = i-1
    while n != 0:
        asize = asize * asize
        n -= 1
    return asize
# end of getallmersize #########

# start of computeLC #############
def computeLC( seqarr ): # compute Trifonov linguistic complexity (LC) score for the seq
    seql = len(seqarr)
    LC = 1.0

    for i in range(1,seql): # consider all mers from 1 to seql-1
        seqimerdict = getseqmers(i, seqarr) # get all mers of size i that are present in seq
        actualcount = sum(seqimerdict.values()) # total count of all mers found in this seq
        expmercount = 1 # expected count for each mer
        diffcount = 0 # difference in total actual mers found and total expected mer count
        allimers = 1 # this is a dummy val for num of all possible mers. This number is 4294967296 for i =5. Thus, dont need to compute allmersize for i > 5
        if i <= 5:
            allimers = getallmersize(i) # get num of all possible mers of size i
            if allimers < seql+1: # num possible mers is less than seqlen, thus each mer can occur multiple times. Calculate uniform dist val.
                expmercount = max(1, int(float(seql)/allimers)) # expected count for each mer assuming uniform occurences for all mers
            # determine what is the total expected count
            totalexp = expmercount * allimers
            if totalexp < actualcount: # if there are actual mer counts that are missed because exp counts can only be integers, determine this diff
                diffcount = actualcount - totalexp # 
        excess = 0 - diffcount
        for sm in seqimerdict:
            excess += seqimerdict[sm] - expmercount # 1 here is the min expected count. If a mer occurs more than that => excess
        LC = LC * (float(actualcount)-excess)/actualcount
        #print i, allimers, seqimerdict, actualcount-excess, actualcount
    return LC
# end of computeLC #############

# start of shannonE ############
def shannonE( seqarr ):
    basecount = {"A":1, "T":1, "C":1, "G":1} # pseudocount of 1 so that natural log is not 0 if a base is not encountered
    for base in seqarr: # get base counts
        basecount[base] += 1

    total = float(sum(basecount.values()))
    basefrac = [basecount[nt]/total for nt in basecount]
    SE = 0.0 # shannon's entropy
    for bf in basefrac:
        SE += bf * math.log(bf,4) # this is log 4, thus entropy is in range 0-1 for alphabet of size 4. 
    SE = -1.0 * SE
    return SE
# end of shannonE ############

# start of usage #############
def usage(): # print usage message
    print ("\nThis is the help message. Type -h to print this message.")
    print ("FindingInfo: calculating linguistic complexity and Shannon's entropy for overlapping windows along a genome.")
    print ("Dependencies: numpy and math packages.")
    print ("User options:")
    print ("\t-f: full path of the genome sequence in fasta format (required).")
    print ("\t-w: the window size for computing linguistic complexity (LC) and Shannon's entropy (H). Default is 10, i.e., LC and H are computed for windows of length 21: from 10 bases upstream to 10 bases downstream of a given genomic position.\n")
# end of usage ###############

# start of main ##############
def main():
    global fname, w
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'hw:f:', ['help','window=','file='])
    except getopt.GetoptError as err:
        # print help information and exit:
        print (str(err))  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    if len(opts) == 0:
        usage()
        sys.exit()
    for o, a in opts:
        if o in ('-f', '--file'): # read genome fasta file supplied by user
            fname = a
        elif o in ('-h', '--help'):
            usage()
            sys.exit()
        elif o in ('-w', '--window'):
            w = int(a)
        else:
            assert False, "unhandled option"
    if fname == None:
        print ("Please provide a sequence file in fasta format.")
        usage()
        sys.exit()
# end of main ################

if __name__ == '__main__':
    main()
    refseq = readUserGenomeSeq(fname)
    refseqarr = list(refseq)
    print ("Genome length: ", len(refseqarr))
    refseqLC = [] # LC for all windows
    refseqSE = [] # shannon entropy for all windows
    for i in range(0,len(refseqarr)-w+1): # from beginning of seq till the pos where the window of size w will exist
        wstart = max(0,i-w) # window start
        wend = min(len(refseqarr),i+w+1) # window end
        seqregion = refseqarr[wstart:wend] # get seq region from i-w to i+w
        refseqLC.append(computeLC(seqregion))
        refseqSE.append(shannonE(seqregion))
    LCarr = np.array(refseqLC)
    SEarr = np.array(refseqSE)
    avgLC = np.mean(LCarr) 
    avgSE = np.mean(SEarr) 
    stdLC = np.std(LCarr) 
    stdSE = np.std(SEarr)
    factorLC = 2*stdLC # how much far from the mean should I consider
    factorSE = 2*stdSE

    foutLC = open("flankingregion_LCscore_w"+str(w)+".txt", 'w') # LC score for all windows
    foutSE = open("flankingregion_SEscore_w"+str(w)+".txt", 'w') # SE score for all windows
    foutlowcomp = open("flankingregion_lowcomplexity_w"+str(w)+".txt", 'w') # regions of low complexity: both LC and SE score are below 2*SD
    for i in range(len(refseqLC)):
        foutLC.write(str(i+1)+"\t"+str(refseqLC[i])+"\n")
        foutSE.write(str(i+1)+"\t"+str(refseqSE[i])+"\n")
        if refseqLC[i] < avgLC-factorLC and refseqSE[i] < avgSE-factorSE: # if LC and SE score are 2*SD below their respective means => Low complexity
            foutlowcomp.write(str(i+1)+"\t"+str(refseqLC[i])+"\t"+str(refseqSE[i])+"\n")
    foutLC.write("Mean:"+str(avgLC)+"\tSD:"+str(stdLC)+"\n")
    foutSE.write("Mean:"+str(avgSE)+"\tSD:"+str(stdSE)+"\n")

    foutLC.close()
    foutSE.close()
    foutlowcomp.close()

