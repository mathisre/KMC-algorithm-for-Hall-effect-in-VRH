# BEFORE USAGE: Define masterData to be folder containing data folders
# BEFORE USAGE: Define combinedDataLoc to be folder where plottable data files are to be placed

# Usage: python3 analyze_directory.py directoryToBeAnalyzed nameOfOutputFile
from __future__ import print_function

from sys import argv, version_info
from os.path import isfile
from os import listdir, chdir
from sys import argv
import time
from statistics import mean, stdev
from math import atan, atan2
from sampler import sampler
from IPython import embed

get_input = input
#Backwards compatability with python 2.7
if version_info[:2] <= (2, 7):
    get_input = raw_input

try:
	if argv[1] == "curr":
		location = ""
		relData  = ""
		saveName = argv[2]
		masterData = ""
		combinedDataLoc = ""
		dataLocation = ""
	else:
		location = argv[1]
		saveName = argv[2]
		relData = "/" + location + "/"
		masterData = "/home/mathisre/Master/data"
		dataLocation = masterData + relData
		chdir(str(dataLocation))
		combinedDataLoc = masterData + "/combinedData/"
except (IndexError, ValueError):
	message = ("Put folder name and save filename in cmd. For example python3 dir_measureConductivity.py" 
		  "kj_HT1 kj_HT1 will combine data in masterData/kj_HT1/ and save as kj_HT1_conductivity.dat. \nTo analyze and save to" 
		  "current directory, use python3 dir_measureConductivity.py curr saveName. This will also save "
		  "file to current directory.")
	print(message)
	exit()

saveFilename = saveName + ".dat"

loc = sampler.extractLocalizationLength(location)
fill = 0.5

Q = -1; L = 100; 

Hz = []; Ex = []; temp = []; variedL = []; Lname = [];
fillList = [];  fillName = []
files = listdir(".")
if "noTri" in files[0]:
	shift = 1
	filestruct = "_jumps_noTri_Hz_"
else:
	shift = 0
	filestruct = "_jumps_Hz_"


# MAKE THE PROGRAM ABLE TO USE VARIABLE n AND L without having L in filename for all the files
for file in files:
	name = file.split("_")

	if "jump" in file and "Hz" in file and "Ex" in file and "T"  in file:
		#ifL in nsame, thena dd L, if not add 100

		if "L" in file:

			if int(name[0]) == 0:
				variedL.append(int(name[2]))
				fillList.append(float(0))
				Lname.append(True) 
				fillName.append(False)
				Hz.append(float(name[5+shift]))
				Ex.append(float(name[7+shift]))
				temp.append(float(name[9+shift]))
		elif "fill" in file:

			if int(name[0]) == 0:
				fillList.append(float(name[2]))
				variedL.append(L)
				fillName.append(True)
				Lname.append(False)
				Hz.append(float(name[5+shift]))
				Ex.append(float(name[7+shift]))
				temp.append(float(name[9+shift]))

		else:

			if int(name[0]) == 0:
				variedL.append(L)
				fillList.append(float(0))
				Lname.append(False)
				fillName.append(False)
				Hz.append(float(name[3+shift]))
				Ex.append(float(name[5+shift]))
				temp.append(float(name[7+shift]))

nList = []
for E,T,H,fill in zip(Ex,temp,Hz,fillList):
	n_i = 0
	if fillName[0] == False:
		for file in files:
			if "Ex_%.5f" % E in file and "T_%.5f" % T in file and "Hz_%.5f" % H in file:
				name = file.split("_")	

				if int(name[0])+1 > n_i:
					n_i = int(name[0]) + 1
	else:		
		for file in files:
			if "Ex_%.5f" % E in file and "T_%.5f" % T in file and "Hz_%.5f" % H in file and "fill_%.5f" % fill in file:
				name = file.split("_")	
				if int(name[0])+1 > n_i:
					n_i = int(name[0]) + 1

	nList.append(n_i)

fileParameterMatrix = (temp, Hz, Ex, Lname, fillName, variedL, nList,fillList)
temp, Hz, Ex, Lname, fillName, variedL, nList, fillList = sampler.sortMatrixByArray((fileParameterMatrix), 0)

varyL = False
varyN = False
varyFill = False
if variedL != [variedL[0] for _ in range(len(variedL))]:
	varyL = True
if nList != [nList[0] for _ in range(len(nList))]:
	varyN = True
if fillList != [fillList[0] for _ in range(len(fillList))]:
	varyFill = True

if not varyFill:
	fillList = [fill for _ in range(len(nList))]


print("\nTemperatures:")
[print("{:5.4f}, ".format(temp[i]), end='') for i in range(len(temp)-1)]
print( "{:5.4f}".format(temp[len(temp)-1]))

print("Magnetic fields:")
[print("{:5.4f}, ".format(Hz[i]), end='') for i in range(len(Hz)-1)]
print( "{:5.4f}".format(Hz[len(Hz)-1]))

print("Electric fields:")
[print("{:5.4f}, ".format(Ex[i]), end='') for i in range(len(Ex)-1)]
print( "{:5.4f}".format(Ex[len(Ex)-1]))

if varyN:
	print("Number of samples:")
	[print("{:3d}, ".format(nList[i]), end='') for i in range(len(nList)-1)]
	print( "{:3d}".format(nList[len(nList)-1]))
else:
	print("\nn =", nList[0])

if varyL:
	print("Lengths:")
	[print("{:6d}, ".format(variedL[i]), end='') for i in range(len(variedL)-1)]
	print( "{:6d}".format(variedL[len(variedL)-1]))


if varyFill:
	print("Fill:")
	[print("{:.4f}, ".format(fillList[i]), end='') for i in range(len(fillList)-1)]
	print( "{:.4f}".format(fillList[len(fillList)-1]))


if varyL: 
	print("\nA = {:4.2f}. Is that ok?".format(loc))
else:
	print("\nL = {:3d} and A = {:4.2f}. Is that ok?".format(L,loc))

print("Filename:", saveFilename)


variedL, L, loc = sampler.everythingOK(variedL, L, loc)


print("A = ", loc)
print("L = ", L)
time.sleep(1)
sampler = sampler(Q, loc)

print(variedL)

# Check if all files are present
for E,T,H,L,n,sizeInFilename, fillInFilename, fill in zip(Ex, temp, Hz, variedL, nList, Lname, fillName, fillList): # Check if all files are present before starting analysis
	for i in range(n):
		LName = "_L_" + str(L) if sizeInFilename else ""		
		FName = "_fill_%.4f" % (fill) if fillInFilename else ""
		filename = dataLocation + str(i) + LName + FName + filestruct + "{0:.5f}_Ex_{1:.5f}_T_{2:.5f}_run_0.dat".format(H,E,T)
		if isfile(filename) != True:
			print("The archives are incomplete!")
			print(filename + " is missing.")
			print("Exiting program.")
			exit()
# Analyze files
for E,T,H,L,n,sizeInFilename, fillInFilename, fill in zip(Ex, temp, Hz, variedL, nList, Lname, fillName, fillList): # analyze files, one combination of T,E,H,L at a time
	sampler.initETHLnConfig(E,T,H,L,n,fill)
	for i in range(n):		
		LName = "_L_" + str(L) if sizeInFilename else ""	
		FName = "_fill_%.4f" % (fill) if fillInFilename else ""
		spesfile =  str(i) + LName + FName + filestruct + "{0:.5f}_Ex_{1:.5f}_T_{2:.5f}_run_0.dat".format(H,E,T)
		filename = dataLocation + spesfile
		print(spesfile)
		sampler.openAndFitData(dataLocation, spesfile)
		sampler.sampleFirstLine(filename)
	sampler.storeToLists()

sampler.writeToFile(combinedDataLoc + saveFilename)
sampler.finalize()


