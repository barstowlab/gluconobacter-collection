#!/sw/bin/python3.5


# ----------------------------------------------------------------------------------------------- #
# condense.py
# Created by Buz Barstow 2016-02-05
# Modified by Sean Medin 2019-09-25

# Code needed to produce re-array instructions for a non-redundant set of mutants from the
# progenitor collection that strike the best balance between isolatability and likelihood of 
# function disruption. 
# ----------------------------------------------------------------------------------------------- #

import sys
sys.path.insert(0, './')

from utilities.grid import FlattenSudokuGridLookupDict, AnnotateSudokuGridWithFeautureNames
from utilities.input import get_input
from utilities.utils import ensure_dir
from utilities.xml import ExportSudokuGridToXML, ImportSudokuGridFromXML

from utilities.feature import ImportFeatureArrayFromGenBank, UpdateFeatureArrayWithSudokuWells, \
    CalculateFeatureArrayTaxonomy

from utilities.condense import PickBestSudokuWellToDisruptFeature, \
    FindFeaturesWithDisruptionsInProgenitorButNotInQC, AssignProgenitorWellsToCondensedCollection, \
    WriteSudokuGridCondensationInstructions, SortBestSudokuWellsIntoSinglyAndMultiplyOccupiedSets, \
    rows96, columns96, AssignPlateGridCoordinatesToCondensedCollectionPlates


# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [ \
    'sudokuGridXMLFile', 'referenceGenomeGenBankFileName', 'rearrayInstructionsFileName', \
    'rearrayFillOrder', 'rearrayFillPattern']

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/post-pool-files/'
file_name = file_intro + 'condense.inp'

inputParameterValues = get_input(file_name, inputParameters)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Parse the input data
sudokuGridXMLFile = file_intro + inputParameterValues['sudokuGridXMLFile']
referenceGenomeGenBankFileName = file_intro + inputParameterValues['referenceGenomeGenBankFileName']
rearrayInstructionsFileName = file_intro + inputParameterValues['rearrayInstructionsFileName']
rearrayFillOrder = inputParameterValues['rearrayFillOrder']
rearrayFillPattern = inputParameterValues['rearrayFillPattern']
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
ensure_dir(rearrayInstructionsFileName)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Import the progenitor collection catalog
sudokuGridLookupDict, prPools, pcPools, rowPools, colPools = \
    ImportSudokuGridFromXML(sudokuGridXMLFile)

# Flatten the progenitor collection catalog
sudokuWellArray = FlattenSudokuGridLookupDict(sudokuGridLookupDict, prPools, pcPools, rowPools, \
                                              colPools)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Import a feature array
featureArray = ImportFeatureArrayFromGenBank(referenceGenomeGenBankFileName)

# Make a gene feature array
geneFeatureArray = []
for feature in featureArray:
    if feature.featureType == 'CDS':
        if 'locus_tag' not in feature.tagDict and 'db_xref' in feature.tagDict:
            feature.tagDict['locus_tag'] = feature.tagDict['db_xref']
        geneFeatureArray.append(feature)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Update the gene feature array with data from the progenitor collection
UpdateFeatureArrayWithSudokuWells(geneFeatureArray, sudokuWellArray)
# CalculateFeatureArrayTaxonomy(geneFeatureArray)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Calculate the features in the genome that have representatives in the progenitor only, and in the
# progenitor and QC collection
[featuresWithDisruptionsInQCAndProgenitor, featuresWithDisruptionsInProgenitorOnly] = \
    FindFeaturesWithDisruptionsInProgenitorButNotInQC(geneFeatureArray)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Pick the best sudoku well objects from the features that only appear in the progenitor collection
featuresAndBestWellArray = PickBestSudokuWellToDisruptFeature( \
    featuresWithDisruptionsInProgenitorOnly, sudokuGridLookupDict)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Separate out the representatives into ones that are singly occupied and ones that are multiply
# occupied.

[singlyOccupiedRepresenativeWells, multiplyOccupiedRepresentativeWells] = \
    SortBestSudokuWellsIntoSinglyAndMultiplyOccupiedSets(featuresAndBestWellArray, \
                                                         rearrayPickOrder='columns')

# Assign wells for rearray of mutants from singly occupied wells
[condensedWellArrayForSingles, currentPlateNumber, currentRowNumber, currentColNumber] = \
    AssignProgenitorWellsToCondensedCollection(singlyOccupiedRepresenativeWells, 1, 1, 1, rows96, \
                                               columns96, fillOrder=rearrayFillOrder, fillPattern=rearrayFillPattern, fillLastPlateWithBlanks=True)

# Assign wells for the colony purification of mutants from multiply occupied wells
[condensedWellArrayForMultiples, currentPlateNumber, currentRowNumber, currentColNumber] = \
    AssignProgenitorWellsToCondensedCollection(multiplyOccupiedRepresentativeWells, currentPlateNumber, \
                                               currentRowNumber, currentColNumber, rows96, \
                                               columns96, fillOrder=rearrayFillOrder, fillPattern=rearrayFillPattern, fillLastPlateWithBlanks=True)

condensedWellArray = condensedWellArrayForSingles + condensedWellArrayForMultiples
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Assign the condensed well array to a new Sudoku grid for re-pooling
condensedWellArray = AssignPlateGridCoordinatesToCondensedCollectionPlates(condensedWellArray)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Write out rearraying instructions
# ----------------------------------------------------------------------------------------------- #
WriteSudokuGridCondensationInstructions(rearrayInstructionsFileName, condensedWellArray)