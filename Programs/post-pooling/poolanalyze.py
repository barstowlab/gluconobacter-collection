#!/sw/bin/python3.5

# ----------------------------------------------------------------------------------------------- #
# poolanalyze
# Created by Buz Barstow 2015-6-14
# Last updated by Sean Medin 2019-08-23

# All code necessary to analyze the pool presence table in preparation for calculating Bayesian
# inference parameters.
# ----------------------------------------------------------------------------------------------- #

import sys
sys.path.insert(0, './')
import warnings
from matplotlib.pyplot import ion, ioff, show

from utilities.input import get_input
from utilities.utils import ensure_dir
from utilities.grid import ImportSudokuGridLayout
from utilities.pool import ImportPoolPresenceTable
from utilities.poolanalysis import CalculateReadCountHistograms, \
CalculateReadCountHistogramsForSingleAddressLines, MakeAndPlotReadCountHistogram, \
PlotEffectOfReadThresholdOnPoolPresenceTableTaxonomy

from numpy import log

# ----------------------------------------------------------------------------------------------- #
warnings.simplefilter("ignore")
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Turn the matplotlib interactive features on to allow multiple plots to be shown from the command
# line
ion()
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Import the input parameter file

inputParameters = [\
'rowPools', 'colPools', 'controlPools', \
'poolPresenceTableFileName', 'sudokuGridLayout', \
'poolPresenceTaxonomyTableFileName', \
'readCountThresholdForReadCountRatioHistogramCalculation', \
'startReadCountThresholdForTrialSolution', 'endReadCountThresholdForTrialSolution', 'readCountThresholdIncrement']

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/post-pool-files/'
file_name = file_intro + 'poolanalyze.inp'
inputParameterValues = get_input(file_name, inputParameters)
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Parse the input data

controlPools = inputParameterValues['controlPools'].split(',')
rowPools = inputParameterValues['rowPools'].split(',')
colPools = inputParameterValues['colPools'].split(',')
poolPresenceTableFileName = file_intro + inputParameterValues['poolPresenceTableFileName']
sudokuGridLayout = file_intro + inputParameterValues['sudokuGridLayout']

readCountThresholdForReadCountRatioHistogramCalculation = \
int(inputParameterValues['readCountThresholdForReadCountRatioHistogramCalculation'])

startReadCountThresholdForTrialSolution = \
int(inputParameterValues['startReadCountThresholdForTrialSolution'])

endReadCountThresholdForTrialSolution = \
int(inputParameterValues['endReadCountThresholdForTrialSolution'])

readCountThresholdIncrement = int(inputParameterValues['readCountThresholdIncrement'])

poolPresenceTaxonomyTableFileName = file_intro + inputParameterValues['poolPresenceTaxonomyTableFileName']
# ----------------------------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------------------------- #
# Initialize the output files
ensure_dir(poolPresenceTaxonomyTableFileName)
# ----------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------ #
# Get the sudoku grid plate layout and the plate row and plate column pools from the sudoku
# grid definition file
[sudokuGridLookupDict, prPools, pcPools] = ImportSudokuGridLayout(sudokuGridLayout, rowPools, \
colPools)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the sudoku pool columns
if controlPools[0].lower() == 'none':
    sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools
    controlPools = None
else:
    sudokuPoolColumns = ['readAlignmentCoord'] + rowPools + colPools + prPools + pcPools \
    + controlPools
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the pool presence table
poolPresenceTable = ImportPoolPresenceTable(poolPresenceTableFileName, sudokuPoolColumns)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Make and Plot histograms of read counts for all lines
[totalReadCountPerLineArray, totalReadCountPerLineNoControlsArray, \
totalRowPoolCountArray, totalColPoolCountArray, totalPRPoolCountArray, totalPCPoolCountArray, \
rowPoolCountArray, colPoolCountArray, prPoolCountArray, pcPoolCountArray] = \
CalculateReadCountHistograms(\
poolPresenceTable, rowPools, colPools, prPools, pcPools, controlPools, \
readCountThresholdForReadCountRatioHistogramCalculation)


MakeAndPlotReadCountHistogram(rowPoolCountArray, plotCumulative=False, \
plotHistogram=True, \
title="Distribution of Row Pool Read Count", \
xlabel="Row Pool Read Count Per Line", ylabel="Lines")

MakeAndPlotReadCountHistogram(colPoolCountArray, plotCumulative=False, \
plotHistogram=True, \
title="Distribution of Col Pool Read Count", \
xlabel="Col Pool Read Count Per Line", ylabel="Lines")

MakeAndPlotReadCountHistogram(prPoolCountArray, plotCumulative=False, \
plotHistogram=True, \
title="Distribution of PR Pool Read Count", \
xlabel="PR Pool Read Count Per Line", ylabel="Lines")

MakeAndPlotReadCountHistogram(pcPoolCountArray, plotCumulative=False, \
plotHistogram=True, \
title="Distribution of PC Pool Read Count", \
xlabel="PC Pool Read Count Per Line", ylabel="Lines")
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Make histograms of read counts for single address lines
[saRowPoolCountArray, saColPoolCountArray, saPrPoolCountArray, saPcPoolCountArray] = \
CalculateReadCountHistogramsForSingleAddressLines(poolPresenceTable, rowPools, colPools, \
prPools, pcPools, readCountThresholdForReadCountRatioHistogramCalculation)

logSaRowPoolCountArray = log(saRowPoolCountArray)
logSaColPoolCountArray = log(saColPoolCountArray)
logSaPrPoolCountArray = log(saPrPoolCountArray)
logSaPcPoolCountArray = log(saPcPoolCountArray)

MakeAndPlotReadCountHistogram(logSaRowPoolCountArray, plotCumulative=False, \
plotHistogram=True, \
title="Distribution of Log Row Pool Read Count for Single Address Lines", \
xlabel="Log Row Pool Read Count Per Line", ylabel="Lines", binWidth=0.2)

MakeAndPlotReadCountHistogram(logSaColPoolCountArray, plotCumulative=False, \
plotHistogram=True, \
title="Distribution of Log Col Pool Read Count for Single Address Lines", \
xlabel="Log Col Pool Read Count Per Line", ylabel="Lines", binWidth=0.2)

MakeAndPlotReadCountHistogram(logSaPrPoolCountArray, plotCumulative=False, \
plotHistogram=True, \
title="Distribution of Log PR Pool Read Count for Single Address Lines", \
xlabel="Log PR Pool Read Count Per Line", ylabel="Lines", binWidth=0.2)

MakeAndPlotReadCountHistogram(logSaPcPoolCountArray, plotCumulative=False, \
plotHistogram=True, \
title="Distribution of Log PC Pool Read Count for Single Address Lines", \
xlabel="Log PC Pool Read Count Per Line", ylabel="Lines", binWidth=0.2)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the effect of read count threshold on the pool presence table taxonomy.
# This step is really important in deciding upon the read count threshold for inclusion of 
# coordinates into the Knockout Sudoku solution of the pool presence table.  
PlotEffectOfReadThresholdOnPoolPresenceTableTaxonomy(poolPresenceTable, \
rowPools, colPools, prPools, pcPools, controlPools, startReadCountThresholdForTrialSolution, \
endReadCountThresholdForTrialSolution, readCountThresholdIncrement, poolPresenceTaxonomyTableFileName)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Show the plots generated in the script
# Turn the interactive mode off to allow plots to stay on screen. 
ioff()
show()
# ------------------------------------------------------------------------------------------------ #