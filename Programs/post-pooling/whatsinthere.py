#!/sw/bin/python3.5

# ----------------------------------------------------------------------------------------------- #
# whatsinthere
# Created by Sean Medin 2021-03-20 (based on isitinthere)
# Updated by Sean Medin 2021-03-20
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
ReverseReducePoolPresenceTable, \
WriteCondensedCollectionSummaryTable, \
UpdateLikelyReadAlignmentCoords, ReviewPastData, DocumentUnexpectedFinds,WriteUnexpectedFindsOutputTable

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [\
'rowPools', 'colPools', 'poolPresenceTableFileName', \
'sudokuGridLayout', 'outputLog', 'maxGapForCoordGrouping', \
'readCountThreshold', 'maxTotalCoords', 'maxSinglePoolCoordNumber', \
'referenceGenomeGenBankFileName', 'previousResults', 'newResults']

argv = sys.argv
if len(argv) != 2:
	print("Error: IndexSummarize takes 1 argument, the input file name")
	sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/post-pool-files/'
file_name = file_intro + 'whatsinthere.inp'

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

poolPresenceTableFileName = file_intro + inputParameterValues['poolPresenceTableFileName']
referenceGenomeGenBankFileName = file_intro + inputParameterValues['referenceGenomeGenBankFileName']


previousResultsFileName = file_intro + inputParameterValues['previousResults']

newResultsFileName = file_intro + inputParameterValues['newResults']

# ----------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
# Initialize the output files
ensure_dir(newResultsFileName)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table

[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, \
colPools)

sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools

poolPresenceTable = ImportPoolPresenceTable(poolPresenceTableFileName, sudokuPoolColumns)

# collects relevant information from isitinthere data
well_plates, coords_already_found, plate_mapping = ReviewPastData(previousResultsFileName)

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Reduce down the pool presence table
reducedPoolPresenceTable = ReverseReducePoolPresenceTable(poolPresenceTable, coords_already_found)
# ------------------------------------------------------------------------------------------------ #

featureArray = ImportFeatureArrayFromGenBank(referenceGenomeGenBankFileName)

geneFeatureArray = []
for feature in featureArray:
	if feature.featureType == 'CDS':
		geneFeatureArray.append(feature)

# documents unexpected finds
unexpectedFinds = DocumentUnexpectedFinds(sudokuGridLookupDict, reducedPoolPresenceTable, \
rowPools, colPools, prPools, pcPools, readCountThreshold,well_plates, featureArray)


# ------------------------------------------------------------------------------------------------ #
# Write out the unexpected finds summary
WriteUnexpectedFindsOutputTable(newResultsFileName, unexpectedFinds, plate_mapping)
# ------------------------------------------------------------------------------------------------ #

