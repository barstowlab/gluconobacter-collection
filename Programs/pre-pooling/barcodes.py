#!/sw/bin/python3.5



# ------------------------------------------------------------------------------------------------ #
def ShannonEntropyDNA(sequence, bases):

    import math

    entropy = 0

    for x in bases:
        p_x = float(sequence.count(x))/len(sequence)
        if p_x > 0:
            entropy += -1.0*p_x*math.log(p_x, 4)

    return entropy
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def hamdist(str1, str2):
    # Count the # of differences between equal length strings str1 and str2
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
# kosudoku-barcode
# Created by Buz Barstow 2016-03-12
# Last modified by Sean Medin 2019-08-20
# Generates new barcodes for Knockout Sudoku
# ------------------------------------------------------------------------------------------------ #


import sys
sys.path.insert(0, './')
from utilities.input import get_input
from utilities.utils import ensure_dir
import itertools
from copy import deepcopy
import operator
import pdb
from Bio.Seq import Seq
from math import ceil,floor,sqrt


# ------------------------------------------------------------------------------------------------ #
# Import the input parameter file

inputParameters = ['rowPools', 'colPools', 'nPlates', 'nExtraBarcodes', 'barcodeSeq1', \
                   'barcodeSeq2', 'barcodeLength', 'outputFileName']

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/pre-pool-files/'
file_name = file_intro + 'barcodes.inp'

inputParameterValues = get_input(file_name, inputParameters)

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Parse the input data
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
nPlates = int(inputParameterValues['nPlates'])
nExtraBarcodes = int(inputParameterValues['nExtraBarcodes'])
barcodeLength = int(inputParameterValues['barcodeLength'])
barcodeSeq1 = str(inputParameterValues['barcodeSeq1'])
barcodeSeq2 = str(inputParameterValues['barcodeSeq2'])
outputFileName = file_intro + str(inputParameterValues['outputFileName'])
# ------------------------------------------------------------------------------------------------ #

ensure_dir(outputFileName)

# ------------------------------------------------------------------------------------------------ #
# Generate all possible arbitrary bp sequences of the desired length and calculate their 
# Shannon entropy
bases = 'ACTG'
seqArray = []

for num in itertools.product(bases, repeat=barcodeLength):
    seq = ''.join(num)
    seqArray.append([seq, ShannonEntropyDNA(seq, bases)])
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Calculate how many barcodes we'll need to generate
plateCols = ceil(sqrt(nPlates))
plateRows = floor(sqrt(nPlates))

while plateRows*plateCols < nPlates:
    plateRows += 1

i = 0
plateRowPools = []
while i < plateRows:
    plateRowPools.append('PR' + '{:02d}'.format(i + 1))
    i += 1

i = 0
plateColPools = []
while i < plateCols:
    plateColPools.append('PC' + '{:02d}'.format(i + 1))
    i += 1

extraPools = []
while i < nExtraBarcodes:
    extraPools.append('X' + '{:01d}'.format(i + 1))
    i += 1

pools = rowPools + colPools + plateRowPools + plateColPools + extraPools

nBarcodes = len(pools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Generate a list of barcode sequences that have the highest entropy and are maximally different
# from one another.

seqArraySorted = sorted(deepcopy(seqArray), key=operator.itemgetter(1), reverse=False)

mostInterestingSeqs = []
min_entropy = 5
min_ham = barcodeLength


i = 0
while i < nBarcodes:

    mostInterestingSeq = [deepcopy(seqArraySorted.pop())[0]]
    mostInterestingSeqs.append(mostInterestingSeq)

    # Go through the remaining elements in the sequence array and calculate their Hamming distances
    # to all of the interesting sequences.
    j = 0
    while j < len(seqArraySorted):
        seqArraySorted[j].append(hamdist(seqArraySorted[j][0], mostInterestingSeqs[-1][0]))
        seqArraySorted[j][2:] = sorted(seqArraySorted[j][2:]) # sorts hamming distance because only the closest distance between any two sequences matter
        j += 1

    # 	pdb.set_trace()

    # Sort the newly appended sequence array and sort it by its Shannon entropy and Hamming
    # distances to all other sequences. To do this, I'm going to define a function on the fly.

    j = 0
    testKeys = [1]
    while j < len(mostInterestingSeqs):
        testKeys.append(j+2)
        j += 1


    def sortTuple(item):
        returnList = []
        for key in testKeys:
            returnList.append(item[key])
        returnTuple = tuple(returnList)
        return returnTuple


    seqArraySorted = sorted(seqArraySorted, key=lambda x: sortTuple(x), reverse=False)
    min_entropy = min(min_entropy,seqArraySorted[-1][1])
    min_ham = min(min_ham, seqArraySorted[-1][2])


    i += 1
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Output the barcodes

outputHandle = open(outputFileName, 'w')

outputStr = 'Pool,Barcode,Reverse Complement\n'
outputHandle.write(outputStr)



i = 0
while i < len(pools):

    outputStr = str(pools[i]) + ',' + mostInterestingSeqs[i][0] + ',' \
                + str(Seq(mostInterestingSeqs[i][0]).reverse_complement()) + '\n'

    outputHandle.write(outputStr)

    i += 1

outputHandle.close()
# ------------------------------------------------------------------------------------------------ #

# prints min entropy and min hamming distance used
print('Minimum Entropy Used:' + str(min_entropy))
print('Minimum hamming distance used:' + str(min_ham))

