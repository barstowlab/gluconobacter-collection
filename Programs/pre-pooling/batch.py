#!/sw/bin/python3.5

# ------------------------------------------------------------------------------------------------ #
# kosudoku-batch
# Created by Buz Barstow 2016-05-06
# ------------------------------------------------------------------------------------------------ #

from utilities.input import get_input
from utilities.batch import CalculateNPick, CalculateNPoolDays, CalculateNPool
import sys

# ------------------------------------------------------------------------------------------------ #
# Import the input parameter file

inputParameters = ['Nplate', 'tLiquid', 'tSolid', 'tSaturation', 'rPick', 'tRest', \
'mutantsPerPlate']

argv = sys.argv
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../../projects/' + proj_name + '/pre-pool-files/'
file_name = file_intro + 'batch.inp'

inputParameterValues = get_input(file_name, inputParameters)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Parse the input data
Nplate = float(inputParameterValues['Nplate'])
tLiquid = float(inputParameterValues['tLiquid'])
tSolid = float(inputParameterValues['tSolid'])
tSaturation = float(inputParameterValues['tSaturation'])
rPick = float(inputParameterValues['rPick'])
tRest = float(inputParameterValues['tRest'])
mutantsPerPlate = int(inputParameterValues['mutantsPerPlate'])

# Calculate the number of mutants that we can pick, and the number of days we can do it over
Npick, pickDays, tPickPerDay = CalculateNPick(tLiquid, tSolid, rPick=rPick, tIncubate=tSaturation)


# Calculate the number of days that we can pool over
NPoolDays = CalculateNPoolDays(tLiquid/24, tRest/24, pickDays, tSaturation)

# Calculate the number of mutants that we can pool
Npool = CalculateNPool(60, tLiquid/24, NPoolDays, mutantsPerPlate)

# Calculate the number of mutants in a batch
Nbatch = min([Npool, Nplate, Npick])

print("Maximum number of mutants that can be picked, Npick = " + str(Npick) + " colonies\n")
print("Maximum number of days over which we can pick, pickDays = " + str(pickDays) + " days\n")
print("Picking time per day: " + '{:0.1f}'.format(tPickPerDay) + " hours\n")
print("Maximum number of mutants that can be pooled, Npool = " + str(Npool) + " colonies\n")
print("Number of mutants per batch that can be pooled, Nbatch = " + str(Nbatch) + " colonies\n")

# ------------------------------------------------------------------------------------------------ #
