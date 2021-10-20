#!/sw/bin/python3.5

# ----------------------------------------------------------------------------------------------- #
# kosudoku-isitinthere
# Created by Buz Barstow 2015-6-14
# Updated by Buz Barstow 2016-05-05
# Last updated by Alexa Schmitz 2020-06-03
# ----------------------------------------------------------------------------------------------- #

import sys
sys.path.insert(0, './')
import pdb
import os

from utilities.input import get_input
from utilities.utils import ensure_dir
from utilities.pool import ImportPoolPresenceTable, SolvePoolPresenceTable
from utilities.xml import ExportSudokuGridToXML, ImportSudokuGridFromXML

from utilities.feature import ImportFeatureArrayFromGenBank, UpdateFeatureArrayWithSudokuWells, \
CalculateFeatureArrayTaxonomy

from utilities.grid import ImportSudokuGridLayout, WriteSudokuGridSummaryTable, \
AnnotateSudokuGridWithFeautureNames, FlattenSudokuGridLookupDict

from utilities.condense import ImportPredictedCondensedCollectionCatalog

from utilities.isitinthere import UpdateCondensedCollectionCatalogWithProgenitorCollection, \
ExpandOutUniqueCoords, \
ReducePoolPresenceTable, \
SolvePoolPresenceTableForColonyPurifiedCollection, \
WriteCondensedCollectionSummaryTable, \
AnnotateCondensedCollectionSudokuGridWithFeautureNames, \
UpdateCondensedCollectionCatalogWithInternalAddressFormats, \
UpdateSudokuGridLookupDictToSudokuColonyPurifiedWell, \
AssignLikelyContentsToSudokuGridLookupDictAndCheckPredictions, \
UpdateLikelyReadAlignmentCoords

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [\
'rowPools', 'colPools', 'poolPresenceTableFileName', \
'sudokuGridLayout', 'outputLog', 'maxGapForCoordGrouping', \
'readCountThreshold', 'maxTotalCoords', 'maxSinglePoolCoordNumber', \
'pgsudokuGridXMLFile', 'referenceGenomeGenBankFileName', \
'predictedCondensedCollectionCatalogFileName', \
'condSudokuGridSummaryFileName', 'condSudokuGridXMLFile']

argv = sys.argv
if len(argv) != 2:
	print("Error: IndexSummarize takes 1 argument, the input file name")
	sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/post-pool-files/'
file_name = file_intro + 'isitinthere.inp'

inputParameterValues = get_input(file_name, inputParameters)
# ----------------------------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
# Parse the input data
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
sudokuGridLayout = file_intro + inputParameterValues['sudokuGridLayout']
maxGapForCoordGrouping = int(inputParameterValues['maxGapForCoordGrouping'])
readCountThreshold = int(inputParameterValues['readCountThreshold'])
maxTotalCoords = int(inputParameterValues['maxTotalCoords'])
maxSinglePoolCoordNumber = int(inputParameterValues['maxSinglePoolCoordNumber'])

pgsudokuGridXMLFile = file_intro + inputParameterValues['pgsudokuGridXMLFile']
poolPresenceTableFileName = file_intro + inputParameterValues['poolPresenceTableFileName']
referenceGenomeGenBankFileName = file_intro + inputParameterValues['referenceGenomeGenBankFileName']

predictedCondensedCollectionCatalogFileName = \
file_intro + inputParameterValues['predictedCondensedCollectionCatalogFileName']

condSudokuGridSummaryFileName = file_intro + inputParameterValues['condSudokuGridSummaryFileName']
condSudokuGridXMLFile = file_intro + inputParameterValues['condSudokuGridXMLFile']
# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
# Initialize the output files
ensure_dir(condSudokuGridSummaryFileName)
ensure_dir(condSudokuGridXMLFile)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table

[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, \
colPools)

sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools

poolPresenceTable = ImportPoolPresenceTable(poolPresenceTableFileName, sudokuPoolColumns)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Read the detailed colony purification data file in. This contains the progenitor collection 
# locations for the colony purified bugs

predictedCondensedCollectionCatalog = \
ImportPredictedCondensedCollectionCatalog(predictedCondensedCollectionCatalogFileName)

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the progenitor collection Sudoku Grid

pSudokuGridLookupDict, pgprPools, pgpcPools, pgrowPools, pgcolPools = \
ImportSudokuGridFromXML(pgsudokuGridXMLFile)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Make a list of transposons that we added into the colony purified collection

flatPSudokuGrid = FlattenSudokuGridLookupDict(pSudokuGridLookupDict, pgprPools, pgpcPools, \
pgrowPools, pgcolPools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Go through the progenitor sudoku grid and add data from it to the colonyPurificationDict
updatedCondensedCollectionCatalog = \
UpdateCondensedCollectionCatalogWithProgenitorCollection(predictedCondensedCollectionCatalog, \
pSudokuGridLookupDict, pgprPools, pgpcPools, pgrowPools, pgcolPools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Expand out the unique coords put in so that there is some wiggle room for alignment error. 
# Expand it out by about 4 bp either side.
uniqueCoordsPutInExpanded = ExpandOutUniqueCoords(updatedCondensedCollectionCatalog)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Reduce down the pool presence table
reducedPoolPresenceTable = ReducePoolPresenceTable(poolPresenceTable, uniqueCoordsPutInExpanded)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Update plate name addresses in updatedColonyPurificationDict with internal names #SEANCOMMENT: Redundant?
updatedCondensedCollectionCatalog = \
UpdateCondensedCollectionCatalogWithInternalAddressFormats(updatedCondensedCollectionCatalog, \
sudokuGridLookupDict, prPools, pcPools)			
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Solve the pool presence table
# This is a much simpler solution than that used for the progenitor collection, so I'm using
# a special function for this case. 
SolvePoolPresenceTableForColonyPurifiedCollection(sudokuGridLookupDict, reducedPoolPresenceTable, \
rowPools, colPools, prPools, pcPools, readCountThreshold)

# Update the sudoku well objects in the sudokugridlookupdict to a sub-class (colony purified well)
# that has a little extra data in it. 
UpdateSudokuGridLookupDictToSudokuColonyPurifiedWell(sudokuGridLookupDict, rowPools, colPools,\
prPools, pcPools)

# Check the sudoku grid lookup dict against the predicted catalog and assign likely contents
AssignLikelyContentsToSudokuGridLookupDictAndCheckPredictions(sudokuGridLookupDict, \
updatedCondensedCollectionCatalog)

flattenedSudokuWellList = FlattenSudokuGridLookupDict(sudokuGridLookupDict, prPools, pcPools, \
rowPools, colPools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Update the flattened sudoku well list to make a list of sudoku genomic coord objects from the 
# simplified likely read alignment coords

UpdateLikelyReadAlignmentCoords(flattenedSudokuWellList)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
featureArray = ImportFeatureArrayFromGenBank(referenceGenomeGenBankFileName)

geneFeatureArray = []
for feature in featureArray:
	if feature.featureType == 'CDS':
		geneFeatureArray.append(feature)

# UpdateFeatureArrayWithSudokuWells(geneFeatureArray, flattenedSudokuGrid)
# featureArrayTaxonomyDict = CalculateFeatureArrayTaxonomy(geneFeatureArray)

AnnotateCondensedCollectionSudokuGridWithFeautureNames(flattenedSudokuWellList, geneFeatureArray)
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
# Write out the colony purified collection summary
WriteCondensedCollectionSummaryTable(flattenedSudokuWellList, condSudokuGridSummaryFileName, \
cpCollectionAddressDictKey='colonyPurified', progenitorCollectionAddressDictKey='progenitor')
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Write out the new sudoku grid as xml
ExportSudokuGridToXML(sudokuGridLookupDict, prPools, pcPools, rowPools, colPools, \
condSudokuGridXMLFile, useSudokuColonyPurifiedWells=True)
# ------------------------------------------------------------------------------------------------ #
