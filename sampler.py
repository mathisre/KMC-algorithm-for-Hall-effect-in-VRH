# Read dx jump data and uses linear fit to measure conductivity for variations in temperature and electric field
# Saves data in files of constant E and constant T
# Put data location file (inside data/) and 
from __future__ import print_function

from sys import version_info
import numpy as np
from scipy.optimize import curve_fit
from os.path import isfile
from statistics import mean, stdev, median, StatisticsError
from math import atan
from IPython import embed


get_input = input
#Backwards compatability with python 2.7
if version_info[:2] <= (2, 7):
    get_input = raw_input

def linear(x, a, b):
	return a*x + b

class sampler:
	def __init__(self, Q, loc, dataCutPercent = 0.50):
		self.Q = Q; self.loc = loc;
		self.dataCutPercent = dataCutPercent

		self.cond_x = []; self.cond_y = []; self.err_x = []; self.err_y = []
		self.E_list = []; self.T_list = []; self.H_list = []; self.L_list = []
		self.jumpLength = []; self.jumpLengthInter = []; self.err_jl = []; self.err_jl_i = []
		self.dxlist = []; self.dylist = []; self.err_dx = []; self.err_dy = []
		self.dxlist_i = []; self.dylist_i = []; self.err_dx_i = []; self.err_dy_i = []
		self.acc_2_list = []; self.acc_3_list = []; self.err_acc_2 = []; self.err_acc_3 = []
		self.driftAngle = []; self.driftAngle_err = []
		self.meanArea = []; self.meanArea_err = []
		self.meanWeightedArea = []; self.meanWeightedArea_err = []
		self.meanXCurrWeightedJL = []; self.meanXCurrWeightedJL_err = []
		self.meanYCurrWeightedJL = []; self.meanYCurrWeightedJL_err = []
		self.dataLength = []; self.differentSizeCombinations = []
		self.fill_list = []

		self.dataSizeFlag = False
		self.printOutlierFunctionInformation = False


	def initETHLnConfig(self, E, T, H, L, n, fill):
		# Initialize lists for data for this config
		self.manyCond_x = []; self.manyCond_y= []
		self.manyJL = []; self.manyJL_i = []; self.manydx = []; self.manydx_i = []
		self.manydy = []; self.manydy_i = []; self.manyAcc_2 = []; self.manyAcc_3 = []
		self.manyDriftRatio = []; self.manyWeightedArea = []; self.manyArea = []
		self.manyXCurrWeightedJL = []
		self.manyYCurrWeightedJL = []
		self.FirstLineOk = True

		self.E = E
		self.T = T
		self.H = H
		self.L = L
		self.n = n
		self.fill = fill

	def openAndFitData(self, dataLocation, spesfile):
		filename = dataLocation + spesfile
		data = np.loadtxt(filename, skiprows=1, unpack=True)		
		time = data[0]
		dxData = data[4]
		dyData = data[5]
		data_size = len(dxData)
		discardData = int(round(self.dataCutPercent*data_size))
	
		dxParam, dxPcov = curve_fit(linear, time[discardData:], self.Q*dxData[discardData:]/(self.L*self.L*1))
		dyParam, dyPcov = curve_fit(linear, time[discardData:], self.Q*dyData[discardData:]/(self.L*self.L*1))	
		
		self.manyCond_x.append(dxParam[0])
		self.manyCond_y.append(dyParam[0])

		self.manyDriftRatio.append(dyData[-1]/dxData[-1])

		self.dataLength.append(data_size)
		if self.dataLength[-1] != self.dataLength[0]:
			self.dataSizeFlag = True
			if [self.E,self.T,self.H] not in self.differentSizeCombinations:
				fileNumber = int(spesfile.split("_")[0])
				self.differentSizeCombinations.append([self.E,self.T,self.H, fileNumber])

	def sampleFirstLine(self, filename):
		file = open(filename, 'r')
		firstLine = file.readline()		
		firstLine = firstLine.split()

		try: # Read first line of file, contains jump length etc information
			if firstLine[9].split("=")[0] != "<jl>": 
				self.FirstLineOk = False

			if float(firstLine[9].split("=")[1]) != 0: # If file actually contains information			
				self.manyJL.append(float(firstLine[9].split("=")[1]))
				self.manydx.append(float(firstLine[10].split("=")[1]))
				self.manydy.append(float(firstLine[11].split("=")[1]))
				self.manyJL_i.append(float(firstLine[12].split("=")[1]))
				self.manydx_i.append(float(firstLine[13].split("=")[1]))
				self.manydy_i.append(float(firstLine[14].split("=")[1]))
				self.manyAcc_2.append(float(firstLine[16].split("=")[1]))
				self.manyAcc_3.append(float(firstLine[18].split("=")[1]))

			
		except (IndexError, ValueError):
			pass

		file.close()


	def storeToLists(self): 
		# Save means and errors into lists
		if not self.FirstLineOk:
			print("First line in previous format. Unable to read")

		n2 = len([i for i in self.manyAcc_3 if i!= 0]) #Find number of non-zero elements in list
		if n2 == 0:
			n2 = 1
			currentJumpLength   = 0
			jumpLengthError     = 0
			currentJumpLength_i = 0
			jumpLengthError_i   = 0
	
			currentDx    = 0
			currentDx_i  = 0
			currentDy    = 0
			currentDy_i  = 0
	
			errorCurrentDx   = 0
			errorCurrentDx_i = 0
			errorCurrentDy   = 0
			errorCurrentDy_i = 0
	
			acc_2 = 0
			acc_3 = 0
			error_acc_2 = 0
			error_acc_3 = 0

			currentArea             = 0
			currentWeightedArea     = 0
			currentArea_err         = 0
			currentWeightedArea_err = 0

			currentXCurrWeightedJL     = 0
			currentXCurrWeightedJL_err = 0
			currentYCurrWeightedJL     = 0
			currentYCurrWeightedJL_err = 0

		else:			
			currentJumpLength   = mean(self.manyJL)
			jumpLengthError     = stdev(self.manyJL)
			currentJumpLength_i = mean(self.manyJL_i)
			jumpLengthError_i   = stdev(self.manyJL_i)
	
			currentDx   = mean(self.manydx)
			currentDx_i = mean(self.manydx_i)
			currentDy   = mean(self.manydy)
			currentDy_i = mean(self.manydy_i)
	
			errorCurrentDx   = stdev(self.manydx)
			errorCurrentDx_i = stdev(self.manydx_i)
			errorCurrentDy   = stdev(self.manydy)
			errorCurrentDy_i = stdev(self.manydy_i)
	
			acc_3 = mean(self.manyAcc_3)
			acc_2 = mean(self.manyAcc_2)
			error_acc_2 = stdev(self.manyAcc_2)
			error_acc_3 = stdev(self.manyAcc_3)

			# These datas can fluctuate a lot, this process removes outliers
			currentArea, currentArea_err 				 	   = self.remove_outliers_return_mean_error(self.manyArea,"manyArea")
			currentWeightedArea, currentWeightedArea_err 	   = self.remove_outliers_return_mean_error(self.manyWeightedArea,"manyWeightedArea")
			currentXCurrWeightedJL, currentXCurrWeightedJL_err = self.remove_outliers_return_mean_error(self.manyXCurrWeightedJL,"manyXCurrWeightedJL") 
			currentYCurrWeightedJL, currentYCurrWeightedJL_err = self.remove_outliers_return_mean_error(self.manyYCurrWeightedJL,"manyYCurrWeightedJL")
		

		meanCond_x = mean(self.manyCond_x)
		errCond_x  = stdev(self.manyCond_x)
		meanCond_y = mean(self.manyCond_y)
		errCond_y  = stdev(self.manyCond_y)
	
		self.err_x.append(errCond_x/np.sqrt(self.n))
		self.err_y.append(errCond_y/np.sqrt(self.n))
		
		self.cond_x.append(meanCond_x)
		self.cond_y.append(meanCond_y)
	
		self.jumpLength.append(currentJumpLength)
		self.jumpLengthInter.append(currentJumpLength_i)
	
		self.err_jl.append(jumpLengthError/np.sqrt(n2))
		self.err_jl_i.append(jumpLengthError_i/np.sqrt(n2))
	
		self.dxlist.append(currentDx)
		self.dxlist_i.append(currentDx_i)
		self.dylist.append(currentDy)
		self.dylist_i.append(currentDy_i)
	
		self.err_dx.append(errorCurrentDx/np.sqrt(n2))
		self.err_dx_i.append(errorCurrentDx_i/np.sqrt(n2))
		self.err_dy.append(errorCurrentDy/np.sqrt(n2))
		self.err_dy_i.append(errorCurrentDy_i/np.sqrt(n2))
	
		self.acc_2_list.append(acc_2)
		self.acc_3_list.append(acc_3)
		self.err_acc_2.append(error_acc_2/np.sqrt(n2))
		self.err_acc_3.append(error_acc_3/np.sqrt(n2))
	
		self.meanDriftAngle = atan(mean(self.manyDriftRatio))
		self.errorDriftAngle = atan(stdev(self.manyDriftRatio)/np.sqrt(n2))
	
		self.driftAngle.append(self.meanDriftAngle)
		self.driftAngle_err.append(self.errorDriftAngle)

		self.meanArea.append(currentArea)
		self.meanWeightedArea.append(currentWeightedArea)
		self.meanArea_err.append(currentArea_err)
		self.meanWeightedArea_err.append(currentWeightedArea_err)

		self.meanXCurrWeightedJL.append(currentXCurrWeightedJL)
		self.meanXCurrWeightedJL_err.append(currentXCurrWeightedJL_err)
		self.meanYCurrWeightedJL.append(currentYCurrWeightedJL)
		self.meanYCurrWeightedJL_err.append(currentYCurrWeightedJL_err)
	
		self.E_list.append(self.E)
		self.H_list.append(self.H)
		self.T_list.append(self.T)
		self.L_list.append(self.L)
		self.fill_list.append(self.fill)

	def writeToFile(self, filename):

		E_arr = np.asarray(self.E_list)
		T_arr = np.asarray(self.T_list)
		H_arr = np.asarray(self.H_list)
		L_arr = np.asarray(self.L_list)
		fill_arr = np.asarray(self.fill_list)

		cond_x_arr = np.asarray(self.cond_x)
		err_x_arr  = np.asarray(self.err_x)
		cond_y_arr = np.asarray(self.cond_y)
		err_y_arr  = np.asarray(self.err_y)		
		
		dxlist_arr     = np.asarray(self.dxlist)
		dylist_arr     = np.asarray(self.dylist)
		err_dx_arr     = np.asarray(self.err_dx)
		err_dy_arr     = np.asarray(self.err_dy)
		dxlist_i_arr   = np.asarray(self.dxlist_i) 
		dylist_i_arr   = np.asarray(self.dxlist_i)
		err_dx_i_arr   = np.asarray(self.err_dx_i) 
		err_dy_i_arr   = np.asarray(self.err_dy_i)
		acc_2_list_arr = np.asarray(self.acc_2_list) 
		acc_3_list_arr = np.asarray(self.acc_3_list)
		err_acc_2_arr  = np.asarray(self.err_acc_2) 
		err_acc_3_arr  = np.asarray(self.err_acc_3)
		
		err_jl_arr   = np.asarray(self.err_jl)
		err_jl_i_arr = np.asarray(self.err_jl_i)
		jl_arr       = np.asarray(self.jumpLength)
		jl_i_arr     = np.asarray(self.jumpLengthInter)

		driftAngle_arr 	   = np.asarray(self.driftAngle)    *180/np.pi
		driftAngle_err_arr = np.asarray(self.driftAngle_err)*180/np.pi


		meanArea_arr 		     = np.asarray(self.meanArea)
		meanWeightedArea_arr     = np.asarray(self.meanWeightedArea)
		meanArea_err_arr 	     = np.asarray(self.meanArea_err)
		meanWeightedArea_err_arr = np.asarray(self.meanWeightedArea_err)

		meanXCurrWeightedJL_arr     = np.asarray(self.meanXCurrWeightedJL)
		meanXCurrWeightedJL_err_arr = np.asarray(self.meanXCurrWeightedJL_err)
		meanYCurrWeightedJL_arr     = np.asarray(self.meanYCurrWeightedJL)
		meanYCurrWeightedJL_err_arr = np.asarray(self.meanYCurrWeightedJL_err)


		loc_arr      = np.asarray([self.loc]*len(self.E_list))	

		headertext = "E-field, temp, H-field, current-x, error-x,, current-y, error-y, jl, loc, jl, error-jl, jl inter, error-jl-inter, \
		dx, error-dx, dy, error-dy, dx_I, error-dx_I, dy_I, error-dy_I, acc-2, error-acc-2, acc-3, error-acc-3, driftAngle, driftAngle_err, \
		 meanArea, meanWeightedArea, meanArea_err, meanWeightedArea_err \
		 meanXCurrWeightedJL, meanXCurrWeightedJL_err, meanYCurrWeightedJL_arr, meanYCurrWeightedJL_err, loc_arr, L_arr, fill_arr"

		np.savetxt(filename, np.transpose([E_arr, T_arr, H_arr, cond_x_arr, err_x_arr, cond_y_arr, err_y_arr, \
			jl_arr, err_jl_arr, jl_i_arr, err_jl_i_arr, dxlist_arr, err_dx_arr, dylist_arr, err_dy_arr, dxlist_i_arr, \
			err_dx_i_arr, dylist_i_arr, err_dy_i_arr, acc_2_list_arr, err_acc_2_arr, acc_3_list_arr, err_acc_3_arr, driftAngle_arr, driftAngle_err_arr, \
			meanArea_arr, meanArea_err_arr, meanWeightedArea_arr, meanWeightedArea_err_arr, \
			meanXCurrWeightedJL_arr, meanXCurrWeightedJL_err_arr, meanYCurrWeightedJL_arr, meanYCurrWeightedJL_err_arr, loc_arr, L_arr, fill_arr]),fmt = "%1.6e", header=headertext)

	@staticmethod
	def extractLocalizationLength(location):
		loc = 0
		try:
			if "A" in location:	
				number = location[location.index("A")+1:]
			
				if "_" in number:
					number_split = number.split("_")
					number = int(number_split[0]) + int(number_split[1])*0.1
				else:
					number = int(number)
				loc = number
			else: 
				loc = 3
		except ValueError:
			loc = get_input("Put loc manually: ")
			while type(loc) != float:
				try: 
					loc = float(loc)
				except ValueError:			
					loc = get_input("Put loc manually: ")
		return loc

	@staticmethod
	def sortMatrixByArray(matrix, index):
		transMatrix = np.array(matrix).transpose()
		sortedMatrix = np.array(sorted(transMatrix, key=lambda a:a[index])).transpose()
		listWithLists = []
		for quantity in sortedMatrix:
			if quantity[0] < 1:
				quantityList = [float(x) for x in quantity]
			if quantity[0] >= 1:
				quantityList = [int(x) for x in quantity]
			listWithLists.append(quantityList)
		return tuple(listWithLists)

	@staticmethod
	def everythingOK(variedL, L, loc):      
	        ok = get_input("Filename ok? (y/n) " )
	        while ok != "y":
	                while ok != "n" and "loc" not in ok and "L=" not in ok and "embed" not in ok:
	                        print(ok)
	                        ok = get_input("Filename ok? (y/n) " )
	                variedL, L, loc = readInput(ok, variedL, L, loc)
	                ok = get_input("Filename ok? (y/n) " )
	        return variedL, L, loc

	@staticmethod
	def readInput(ok, variedL, L, loc):
		if ok == "y":
			return variedL, L, loc
		elif ok == "n":
			exit()
		elif "loc" in ok:
		    loc = float(ok.split("=")[1])
		    return variedL, L, loc
		elif "L" in ok:
		    L = int(ok.split("=")[1])
		    variedL = [L for _ in range(len(variedL))]
		    return variedL, L, loc
		elif "embed" in ok:
		    embed()
		    return variedL, L, loc



	def finalize(self):
		if self.dataSizeFlag:
			st = ""
			for ele in self.differentSizeCombinations:
				st += "(%d, %.3f, %.3f, %.3f)" % (ele[3], ele[0], ele[1], ele[2])
			print("\n Data files (n, E, T, H) = {:s} were of different lengths than the first file. Is everything as it should?\n".format(st))

	def remove_outliers_return_mean_error(self, data, name):

		try:
			data_mean = mean(data)
			data_stdev = stdev(data)

			if self.printOutlierFunctionInformation == True:
				print("\n" + name + ":")
				print("Before")	
				print(np.asarray(data))
				print("Mean, stdev: %4.4f, %4.4f" % (data_mean, stdev(data)))
	
			data = [x for x in data if (x >= 0.0)]
			data_mean = mean(data)
			data_stdev = stdev(data)
	
			if data_stdev > 0.2:
				# if data_mean > 0.01:
	
				for i in range(1):
						
					data_stdev = stdev(data)
					data_mean = mean(data)
			
					data = [x for x in data if (x > data_mean - 2.5 * data_stdev)]		
					data = [x for x in data if (x < data_mean + 2.5 * data_stdev)]
	
				data_stdev = stdev(data)
				if data_stdev > 0.2:
					data_mean = mean(data)

					data = [x for x in data if (x > data_mean - 2.5 * data_stdev)]		
					data = [x for x in data if (x < data_mean + 2.5 * data_stdev)]
		
							# data_mean = mean(data)
							# data_stdev = stdev(data)
				
				data_stdev = stdev(data)
		
			
			data_mean = mean(data)
			data_mean_stdev = stdev(data)/np.sqrt(len(data))

			if self.printOutlierFunctionInformation == True:
				print("After")
				print(np.asarray(data))
				print("Mean, stdev: %4.4f, %4.4f" % (data_mean, stdev(data)))
		except StatisticsError:
			data_mean = 0
			data_mean_stdev = 0
		return data_mean, data_mean_stdev
		



