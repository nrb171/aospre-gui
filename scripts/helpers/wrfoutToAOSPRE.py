### Convert WRF output to AOSPRE compatible format
# usage: python wrfoutToAOSPRE.py -in /path/to/wrfout/files/ -out /path/to/output/files/ -start
# -in and -out can be the same directory
# if you don't know the start time, it can be found with: `ncdump -h /path/to/wrfout.nc | grep :SIMULATION_START_DATE`

### Load modules
import datetime as dt
import argparse
import os

### parse command line arguments
parser = argparse.ArgumentParser(description='Rename WRF output files to AOSPRE compatible format')
parser.add_argument('-in' , '--inputDir', help='Directory containing WRF output files')
parser.add_argument('-out','--outputDir', help='Directory to save AOSPRE compatible files')
parser.add_argument('-o', '--options', help='Options for renaming files ln or mv (default: ln)', default='ln')
parser.add_argument('-start', '--startTime', 
                    help='Start time for the WRF simulation in format YYYY-MM-DD_HH:MM:SS. If not provided, will use the first file found in the directory.',
                    default=None)
args = parser.parse_args()



wrfoutDir = args.inputDir
outputDir = args.outputDir
option = args.options
startString = args.startTime

if option not in ['ln','mv']:
    print('symlink')
    option = 'ln'

if option == 'ln':
    option = 'ln -s'


### find files in wrfoutDir, find start time, 
# find files
wrfFiles = os.listdir(wrfoutDir)
for i in range(len(wrfFiles)):
    wrfFiles[i].replace('/',':')
mask = ['wrfout' in file for file in wrfFiles]
# print(mask)
wrfFiles = [wrfFiles[i] for i in range(len(wrfFiles)) if mask[i]]
wrfFiles.sort()



# print(wrfFiles)

# find start time
if startString is None:
    raise(ValueError('No start time provided . Please provide a start time with -start.'))
else:
    startDt = dt.datetime.strptime(startString, '%Y-%m-%d_%H:%M:%S')

### convert to seconds since start time for each file and rename
checked = 0
for file in wrfFiles:
    if 'wrfout' not in file:
        continue
    fileString = file.split('_')
    prefix = fileString[0]+'_'+fileString[1]+'_'
    timeString = fileString[2]+'_'+fileString[3]
    # check if file is wrfout file
    fileTime = dt.datetime.strptime(timeString, '%Y-%m-%d_%H:%M:%S')
    
    # calculate time difference
    timeDelta = fileTime - startDt
    seconds = timeDelta.total_seconds()

    # rename file
    newFile = prefix+'{:06.0f}'.format(seconds)+'.nc'
    # command = 'mv '+wrfoutDir+'/'+file+' '+outputDir+'/'+newFile
    command = option+' '+wrfoutDir+'/'+file+' '+outputDir+'/'+newFile

    print(command)
    if checked == 0:
        if input('is this correct? (y/n): ') != 'y':
            break
        checked = 1

    os.system(command)
    



