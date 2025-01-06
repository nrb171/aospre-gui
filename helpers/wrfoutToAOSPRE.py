### Convert WRF output to AOSPRE compatible format
# usage: python wrfoutToAOSPRE.py -in /path/to/wrfout/files/ -out /path/to/output/files/
# -in and -out can be the same directory

### Load modules
import datetime as dt
import sys
import argparse
import os

### parse command line arguments
parser = argparse.ArgumentParser(description='Rename WRF output files to AOSPRE compatible format')
parser.add_argument('-in' , '--inputDir', help='Directory containing WRF output files')
parser.add_argument('-out','--outputDir', help='Directory to save AOSPRE compatible files')
args = parser.parse_args()

wrfoutDir = args.inputDir
outputDir = args.outputDir

### find files in wrfoutDir, find start time, 
# find files
wrfFiles = os.listdir(wrfoutDir)
for i in range(len(wrfFiles)):
    wrfFiles[i].replace('/',':')
mask = ['wrfout' in file for file in wrfFiles]
print(mask)
wrfFiles = [wrfFiles[i] for i in range(len(wrfFiles)) if mask[i]]
wrfFiles.sort()



print(wrfFiles)

# find start time
startString = wrfFiles[0].split('_')
timeString = startString[2]+'_'+startString[3]
startDt = dt.datetime.strptime(timeString, '%Y-%m-%d_%H:%M:%S')

### convert to seconds since start time for each file and rename
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
    command = 'mv '+wrfoutDir+'/'+file+' '+outputDir+'/'+newFile
    print(command)
    os.system(command)
    



