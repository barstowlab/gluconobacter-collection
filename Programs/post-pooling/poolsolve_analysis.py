# generates graphs to help figure out how to optimize poolsolve

import sys
sys.path.insert(0, './')
from operator import itemgetter

from utilities.input import get_input
from utilities.utils import ensure_dir
from utilities.grid import ImportSudokuGridLayout, FlattenSudokuGridLookupDict, \
    AnnotateSudokuGridWithFeautureNames, WriteSudokuGridSummaryTable
from utilities.pool import ImportPoolPresenceTable, SolvePoolPresenceTable
from utilities.xml import ExportSudokuGridToXML
from utilities.feature import ImportFeatureArrayFromGenBank, UpdateFeatureArrayWithSudokuWells, \
    CalculateFeatureArrayTaxonomy
from utilities.poolanalysis import ReadBayesianInferenceParameters
from utilities.poolsolve_analysis_helper import *

import pdb
import os
import warnings
import csv
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------- #
warnings.simplefilter("ignore")
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]

file_intro = '../../projects/' + proj_name + '/post-pool-files/'
file_name = file_intro + 'poolsolve_analysis.inp'


inputParameters = [ \
    'rowPools', 'colPools', 'controlPools', 'sudokuGridLayout', 'genBankFileName', \
    'poolPresenceTableFileName', \
    'readCountThreshold', 'maxGapForCoordGrouping', 'voigtScoreThreshold', 'maxTotalCoords', \
    'maxSinglePoolCoordNumber', \
    'sudokuGridXMLFile', 'sudokuGridSummaryFileName', \
    'bayesianParameterFileName', 'barcodeFile', 'matingFile']

inputParameterValues = get_input(file_name, inputParameters)

# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
controlPools = inputParameterValues['controlPools'].split(',')
poolPresenceTableFileName = file_intro + inputParameterValues['poolPresenceTableFileName']
genBankFileName = file_intro + inputParameterValues['genBankFileName']
sudokuGridLayout = file_intro + inputParameterValues['sudokuGridLayout']
maxGapForCoordGrouping = int(inputParameterValues['maxGapForCoordGrouping'])
readCountThreshold = int(inputParameterValues['readCountThreshold'])
voigtScoreThreshold = float(inputParameterValues['voigtScoreThreshold'])
maxTotalCoords = int(inputParameterValues['maxTotalCoords'])
maxSinglePoolCoordNumber = int(inputParameterValues['maxSinglePoolCoordNumber'])
sudokuGridXMLFile = file_intro + inputParameterValues['sudokuGridXMLFile']
sudokuGridSummaryFileName = file_intro + inputParameterValues['sudokuGridSummaryFileName']
bayesianParameterFileName = file_intro + inputParameterValues['bayesianParameterFileName']
barcodeFile = file_intro + inputParameterValues['barcodeFile']
matingFile = file_intro + inputParameterValues['matingFile']
# ----------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table

[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, \
                                                                  colPools)

numToMating = {}
with open(matingFile, newline='', encoding="utf-8-sig") as csvfile:
    mating_reader = csv.reader(csvfile, delimiter=',')
    mating_reader.__next__()
    for row in mating_reader:
        numToMating[row[0]] = row[1]


with open(sudokuGridLayout, newline='', encoding="utf-8-sig") as csvfile:
    grid_reader = csv.reader(csvfile, delimiter=',')
    cols = grid_reader.__next__()
    for row in grid_reader:
        pr = row[0]
        for i in range(1, len(row)):
            if row[i].isdigit():
                sudokuGridLookupDict[pr][cols[i]].mating = numToMating[row[i]]
            else:
                sudokuGridLookupDict[pr][cols[i]].mating = -1

if controlPools[0].lower() == 'none':
    sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools
    controlPools = None
else:
    sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools \
                        + controlPools




poolPresenceTable = ImportPoolPresenceTable(poolPresenceTableFileName, sudokuPoolColumns)

# reads barcode file
barcodeMap = {}
with open(barcodeFile, newline='', encoding="utf-8-sig") as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in csvreader:
        if row[0] in sudokuPoolColumns:
            barcodeMap[row[0]] = row[1]

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Read in the Bayesian inference parameters

[logReadNumberRatioHistogramFitDict, logReadNumberRatioHistogramIntegralDict] = \
    ReadBayesianInferenceParameters(bayesianParameterFileName)

# counts number of times each pool combo appears
all_combos = {}
individual_counts = {}

for i in range(1, len(sudokuPoolColumns)):
    individual_counts[sudokuPoolColumns[i]] = 0
    for j in range(i + 1, len(sudokuPoolColumns)):
        itype = sudokuPoolColumns[i]
        jtype = sudokuPoolColumns[j]
        all_combos[(itype, jtype)] = 0

for row in range(poolPresenceTable.size):
    for i in range(1, len(sudokuPoolColumns)):
        if poolPresenceTable[row][i] >= readCountThreshold:
            individual_counts[sudokuPoolColumns[i]] += 1
        for j in range(i + 1, len(sudokuPoolColumns)):
            if poolPresenceTable[row][i] >= readCountThreshold and poolPresenceTable[row][j] >= readCountThreshold:
                tup = (sudokuPoolColumns[i], sudokuPoolColumns[j])
                all_combos[tup] += 1

examine_mutations(all_combos, barcodeMap)

examine_mutliple_pickings(poolPresenceTable, sudokuGridLookupDict, prPools, pcPools, sudokuPoolColumns, readCountThreshold)

examine_contamination(all_combos, individual_counts)

# print(len(rowPools))
# print(len(colPools))
# print(len(prPools))
# print(len(pcPools))
#
# # finds distribution of number of plate rows, plate columns, columns, and rows for each entry (looks for contamination)
# dup_rows = np.zeros(len(rowPools) + 1)
# dup_cols = np.zeros(len(colPools) + 1)
# dup_prows = np.zeros(len(prPools) + 1)
# dup_pcols = np.zeros(len(pcPools) + 1)
# for row in range(poolPresenceTable.size):
#     num_rows = 0
#     num_cols = 0
#     num_prows = 0
#     num_pcols = 0
#     for i in range(1, len(poolPresenceTable.dtype.names)):
#         type = poolPresenceTable.dtype.names[i]
#         if poolPresenceTable[row][i] > readCountThreshold:
#             if type in rowPools:
#                 num_rows += 1
#             elif type in colPools:
#                 num_cols += 1
#             elif type in prPools:
#                 num_prows += 1
#             else:
#                 num_pcols += 1
#     dup_rows[num_rows] += 1
#     dup_cols[num_cols] += 1
#     dup_prows[num_prows] += 1
#     dup_pcols[num_pcols] += 1
#
# plt.figure()
# plt.bar(list(range(len(dup_rows))), dup_rows)
# plt.title('Rows duplicate distribution')
# print(sum(dup_rows[2:]))
#
# plt.figure()
# plt.bar(list(range(len(dup_prows))), dup_prows)
# plt.title('Plate rows duplicate distribution')
# print(sum(dup_prows[2:]))
#
# plt.figure()
# plt.bar(list(range(len(dup_cols))), dup_cols)
# plt.title('Cols duplicate distribution')
# print(sum(dup_cols[2:]))
#
# plt.figure()
# plt.bar(list(range(len(dup_pcols))), dup_pcols)
# plt.title('Plate cols duplicate distribution')
# print(sum(dup_pcols[2:]))
#
#
#
# # finds distribution of overlaps between indices (treats same type and different type separately)
# # will look for spills and barcode mutations
# expected_combos = {}
# unexpected_combos = {}
# for i in range(1, len(sudokuPoolColumns)):
#     for j in range(i + 1, len(sudokuPoolColumns)):
#         itype = sudokuPoolColumns[i]
#         jtype = sudokuPoolColumns[j]
#         if (itype in rowPools and jtype in rowPools) or (itype in colPools and jtype in colPools) or (itype in prPools and jtype in prPools) \
#             or (itype in pcPools and jtype in pcPools):
#             unexpected_combos[(i,j)] = 0
#         else:
#             expected_combos[(i,j)] = 0
#
# for row in range(poolPresenceTable.size):
#     for i in range(1, len(sudokuPoolColumns)):
#         for j in range(i + 1, len(sudokuPoolColumns)):
#             if poolPresenceTable[row][i] >= readCountThreshold and poolPresenceTable[row][j] >= readCountThreshold:
#                 if (i,j) in unexpected_combos:
#                     unexpected_combos[(i,j)] += 1
#                 else:
#                     expected_combos[(i,j)] += 1
#
# unexpected_combo_counts = []
# expected_combo_counts = []
# for combo in unexpected_combos:
#     unexpected_combo_counts.append(unexpected_combos[combo])
#
# for combo in expected_combos:
#     expected_combo_counts.append(expected_combos[combo])
#
# plt.figure()
# plt.hist(unexpected_combo_counts)
# plt.title('Combos that should not happen')
#
# plt.figure()
# plt.hist(expected_combo_counts)
# plt.title('Combos that should happen about 400 times')
#
# # calculates hamming distance between barcodes and compares overlap frequency for each case
# hamming_combos = {}
# for i in range(1, len(sudokuPoolColumns)):
#     for j in range(i + 1, len(sudokuPoolColumns)):
#         code1 = barcodeMap[i]
#         code2 = barcodeMap[j]
#         ham = 0
#         for c1, c2 in zip(code1, code2):
#             if c1 != c2:
#                 ham += 1
#         hamming_combos[(i,j)] = ham
#
# unexpected_hammings_tots = np.zeros(9)
# expected_hammings_tots = np.zeros(9)
#
# unexpected_hammings_count = np.zeros(9)
# expected_hammings_count = np.zeros(9)
#
# for combo in unexpected_combos:
#     ham = hamming_combos[combo]
#     unexpected_hammings_tots[ham] += unexpected_combos[combo]
#     unexpected_hammings_count[ham] += 1
#
# for combo in expected_combos:
#     ham = hamming_combos[combo]
#     expected_hammings_tots[ham] += expected_combos[combo]
#     expected_hammings_count[ham] += 1
#
# plt.figure()
# plt.bar(list(range(9)), unexpected_hammings_tots / unexpected_hammings_count)
# plt.title('Unexpected Combos')
#
# plt.figure()
# plt.bar(list(range(9)), expected_hammings_tots / expected_hammings_count)
# plt.title('Expected Combos')
#
# # prints list of most commonly appearing expected and unexpected
# unexpected_tups = []
# expected_tups = []
#
# for combo in expected_combos:
#     expected_tups.append(((sudokuPoolColumns[combo[0]], sudokuPoolColumns[combo[1]]), expected_combos[combo]))
#
# for combo in unexpected_combos:
#     unexpected_tups.append(((sudokuPoolColumns[combo[0]], sudokuPoolColumns[combo[1]]), unexpected_combos[combo]))
#
# expected_tups.sort(key=itemgetter(1), reverse=True)
# print(expected_tups)
#
# unexpected_tups.sort(key=itemgetter(1), reverse=True)
# print(unexpected_tups)
#
# # finds how often each entry occurs
#
# plt.ioff()
# plt.show()
