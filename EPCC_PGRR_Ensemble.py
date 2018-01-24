############################################################################
# Project: The Lenard effect of preciptation at the RUAO,
# Title: Ensemble processing of the PG, Time and Rain Rate data,
# Author: James Gilmore,
# Email: james.gilmore@pgr.reading.ac.uk.
# Version: 1.2.1
# Date: 12/08/16
# Status: Operational (Basic)
# Change: Added in support for calling from outside modules. i.e. the function PGRR_Ensembler now states a unique file location in PGRR_Loc
############################################################################

#Initialising the python script
from __future__ import absolute_import, division, print_function
from scipy import stats, interpolate
from array import array
import sys,csv,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from PIL import Image
import time as systime

sys.path.insert(0, '/home/th863480/PhD/Global_Functions')
#from externals import PGRR_Plots #Cyclical import issue!!!
from lowess import lowess
from Gilly_Utilities import nan_removal_dd

def PGRR_Ensembler(bincount=None, loop=None, selectcase=None, method=None, case=None, PGRR_Loc=None, pertcount=None, dataimport=True, year=None, month=None, time=None, rainrate=None, pg=None, Storage_Path='../../', Plots_Path='Plots/Ensemble/', Data_Path='Processed_Data/Ensemble/'):
		
	#Import the data gathered by PGRainRate.py
	#if dataimport == True: year, month, time, rainrate, pg = np.genfromtxt('../../Processed_Data/Rain_Days/PG_Data_' + str(PGRR_Loc) + '.csv', dtype=float, delimiter=',', unpack=True, usecols=(0,1,2,3,4))
	if dataimport == True: year, month, time, rainrate, pg = np.genfromtxt('/glusterfs/phd/users/th863480/WC2_Surface_Electrification/Processed_Data/Radar/PG_LWC_Median_Comparison.csv', dtype=float, delimiter=',', unpack=True, usecols=(0,1,2,3,4))

	############################################################################
	#Initialise
	"""Convert Data types to numpy arrays"""
	month = np.asarray(month) if month is not None else None
	time = np.asarray(time) if time is not None else None
	rainrate = np.array(rainrate, dtype=float) if rainrate is not None else None
	pg = np.array(pg, dtype=float) if pg is not None else None
	year = np.asarray(year) if year is not None else None
	
	"""Remove zero values from processed data. Used "year" as criteria as theres a chance of real
	zero values cropping up in other columns."""
	Month = month.copy()[year.copy()!=0] if month is not None else np.nan
	Time = time.copy()[year.copy()!=0] if time is not None else np.nan
	Rainrate = rainrate.copy()[year.copy()!=0] if year is not None else rainrate.copy()
	PG = pg.copy()[year.copy()!=0] if year is not None else pg.copy()
	Year = year.copy()[year.copy()!=0] if year is not None else np.nan
	
	"""Define arrays for statistical analysis"""
	slope = np.zeros(int(loop))
	intercept = np.zeros(int(loop))
	r_value = np.zeros(int(loop))
	p_value = np.zeros(int(loop))
	std_err =  np.zeros(int(loop))
	lowessval = np.zeros([int(loop),bincount+30*loop])

	"""Sort the data by the date (YEAR, YEARDAY)"""
	PGRR = np.asarray(zip(Rainrate, PG))
	PGRRsort = PGRR[np.lexsort((PGRR[:, 1], PGRR[:, 0]))]
	#print("mean = ", np.mean(PGRRsort[:, 1]), " median = ", np.median(PGRRsort[:, 1]))
	PGsort = PGRRsort[:,1]
	RRsort = PGRRsort[:,0]
	
	Rainrate, PG = nan_removal_dd(np.array([rainrate,pg]), unpack=True)
	
	"""Removes Zeros from Data"""
	PGRRsort = PGRRsort[PGRRsort[:,0] != 0]
	PGRR = PGRRsort

	############################################################################
	"""Loop around multiple bin sizes (if that was what we selected)"""
	
	
	for k in xrange(loop):
		"""Initalise the matrices and vectors for each loop"""
		RainRateBin = np.zeros((bincount+30*k)-1)
		RainRateBinLimit = np.zeros(bincount+30*(k))
		TimeTipBin = np.zeros(bincount+30*(k))
		PGTipBin = np.zeros(bincount+30*(k))
		TotalBin = np.zeros(bincount+30*(k))
		PGTipBinMedian = np.zeros([bincount+30*(k),len(Year)]) if year is not None else None
		#PGTipBinMedianConst = np.zeros([bincount+30*k, len(PGRR)/(bincount+30*k)]) probs dont need it now (for case ==3)
		PGTipPosition = np.zeros(bincount+30*(k))
		PGTipBinMedianFinal = np.zeros(bincount+30*(k))
		PGTipBinMedianFinal_se = np.zeros(bincount+30*(k))
		eps = sys.float_info.epsilon
		PGTipBinValues = np.zeros([len(PGRR)/(bincount+30*k)+1,bincount])

		"""Define the Rain Rate for each bin with the centred values determined as well."""
		for i in range(bincount+30*(k)):
			RainRateBinLimit[i] = i*25/(bincount+30*(k))
		for i in range((bincount+30*(k))-1):
			RainRateBin[i] = 0.5*(RainRateBinLimit[i+1]-RainRateBinLimit[i])
			
		if selectcase == 1:
			############################################################################
			"""Define the mean (ensemble) PG and Tip Times for the statistically significant data.
			Equal Bin Spacing, Variable Bin Counts"""
			for j in range(len(Year)):
				#print("PG[j]", PG[j])
				for i in range(1,bincount+30*(k)):
					if (Rainrate[j] < i*5/(bincount+30*(k)) and Rainrate[j] > (i-1)*5/(bincount+30*(k))):
						PGTipBin[i] += PG[j]
						TimeTipBin[i] += Time[j]
						TotalBin[i] += 1
			PGTipBinned = PGTipBin.copy()/(TotalBin.copy()+eps)
			PGTipBinned_se = 0 #######NEEDS FIXING (MAKE PGTipBin a container for all PG and not +=)
			TimeTipBinned = TimeTipBin.copy()/(TotalBin.copy()+eps)

			#Removes NaN values
			PGTipBinned = [0 if np.isnan(x) else x for x in PGTipBinned]
			#PGTipBinned_se = [0 if np.isnan(x) else x for x in PGTipBinned_se]
			TimeTipBinned = [0 if np.isnan(x) else x for x in TimeTipBinned]
		
			#Select values for plotting
			yvalue = np.asarray(PGTipBinned)
			yvalue_se = np.asarray(PGTipBinned_se)
			amethod = "Static_Mean"
		
			############################################################################

		elif selectcase == 2:	
			############################################################################
			"""Define the median PG and Tip Times for the statistically significant data.
			Equal Bin Spacing, Variable Bin Counts."""

			for j in range(len(Year)):
				for i in range(1,bincount+30*(k)):
					if (Rainrate[j] < i*25/(bincount+30*(k)) and Rainrate[j] > (i-1)*25/(bincount+30*(k))):
						PGTipBinMedian[i,PGTipPosition[i]] = PG[j]
						PGTipPosition[i]+=1
			for i in range(bincount+30*(k)):
				PGTipBinMedianFinal[i] = np.median(PGTipBinMedian[i,:].copy()[PGTipBinMedian[i,:].copy()!=0])
				PGTipBinMedianFinal_se[i] = np.std(PGTipBinMedian[i,:].copy()[PGTipBinMedian[i,:].copy()!=0])/np.sqrt(PGTipPosition[i])
				PGTipBinMedianFinal[np.isnan(PGTipBinMedianFinal)] = 0
				PGTipBinMedianFinal_se[np.isnan(PGTipBinMedianFinal_se)] = 0
			
			#Select values for plotting
			yvalue = np.asarray(PGTipBinMedianFinal)
			yvalue_se = np.asarray(PGTipBinMedianFinal_se)
			amethod = "Static_Median"
			print("PGTipBinMedian", PGTipBinMedianFinal)

			############################################################################
		
		elif selectcase == 3:
			############################################################################
			"""Define the mean (ensemble) PG and Tip Times for the statistically significant data.
			Variable Bin Spacing, Equal Bin Counts."""
			
			PGTipBinned = np.zeros(bincount+30*(k))
			PGTipBinned_se = np.zeros(bincount+30*(k))
			PGTipBinned[0] = np.mean(PGRRsort[0:(len(PGRR)/(bincount+30*k)), 1].copy())
			PGTipBinned_se[0] = np.std(PGRRsort[0:(len(PGRR)/(bincount+30*k)), 1].copy())/np.sqrt(len(PGRR)/(bincount+30*k))
			for i in range(1,(bincount+30*(k))):
				PGTipBinned[i] = np.mean(PGRRsort[(len(PGRR)/(bincount+30*k)*i):(len(PGRR)/(bincount+30*k)*(i+1)), 1].copy())
				PGTipBinned_se[i] =  np.std(PGRRsort[(len(PGRR)/(bincount+30*k)*i):(len(PGRR)/(bincount+30*k)*(i+1)), 1].copy())/np.sqrt(len(PGRR)/(bincount+30*k))
				PGTipBinned[np.isnan(PGTipBinned)] = 0
				PGTipBinned_se[np.isnan(PGTipBinned_se)] = 0
				
			#Select values for plotting
			yvalue = np.asarray(PGTipBinned)
			yvalue_se = np.asarray(PGTipBinned_se)
			amethod = "Dynamic_Mean"
			#print("PGTipBinned", PGTipBinned)

			#Define the Rain Rate for each bin with the centred values determined as well.
			RainRateBinLimit[0] = 0.5*PGRRsort[(len(PGRR)/(bincount+30*k)), 0]
			for i in range(1,bincount+30*(k)):				
				RainRateBinLimit[i] = 0.5*(PGRRsort[(len(PGRR)/(bincount+30*k)*i), 0]-PGRRsort[(len(PGRR)/(bincount+30*k)*(i-1)), 0])+PGRRsort[(len(PGRR)/(bincount+30*k)*(i-1)), 0]
			print(RainRateBinLimit)
			############################################################################

		elif selectcase == 4:	
			############################################################################
			"""Define the median PG and Tip Times for the statistically significant data.
			Variable Bin Spacing, Equal Bin Counts."""
			
			PGTipBinMedianFinal[0] = np.median(PGRRsort[0:(len(PGRR)/(bincount+30*k)), 1].copy())
			PGTipBinMedianFinal_se[0] = np.std(PGRRsort[0:(len(PGRR)/(bincount+30*k)), 1].copy())/np.sqrt(len(PGRR)/(bincount+30*k))
			for i in range(1,(bincount+30*(k))):
				PGTipBinMedianFinal[i] = np.median(PGRRsort[np.floor(len(PGRR)/(bincount+30*k)*i):np.floor(len(PGRR)/(bincount+30*k)*(i+1)), 1].copy())
				PGTipBinMedianFinal_se[i] = np.std(PGRRsort[(len(PGRR)/(bincount+30*k)*i):(len(PGRR)/(bincount+30*k)*(i+1)), 1].copy())/(np.sqrt(len(PGRR)/(bincount+30*k))+eps)
				#PGTipBinMedianFinal = PGTipBinMedianFinal[~np.isnan(PGTipBinMedianFinal)]
				#PGTipBinMedianFinal_se = PGTipBinMedianFinal_se[~np.isnan(PGTipBinMedianFinal_se)]
				
				#Puts all the values found for each bin into a matrix
				try:
					PGTipBinValues[:, i] = PGRRsort[np.floor(len(PGRR)/(bincount+30*k)*i):np.floor(len(PGRR)/(bincount+30*k)*(i+1)), 1].copy()
				except:
					PGTipBinValues[:, i] = np.append(PGRRsort[np.floor(len(PGRR)/(bincount+30*k)*i):np.floor(len(PGRR)/(bincount+30*k)*(i+1)), 1].copy(),0)
			
			############################################################################
			#Select values for plotting
			yvalue = np.asarray(PGTipBinMedianFinal)
			yvalue_se = np.asarray(PGTipBinMedianFinal_se)
			amethod = "Dynamic_Median"
			#print("PGTipBinMedian", np.sort(PGTipBinMedianFinal))

			#Define the Rain Rate for each bin with the centred values determined as well.
			#print(PGRRsort)
			#print('RainRateBinLimit', RainRateBinLimit.shape, np.array(PGRRsort).shape)
			RainRateBinLimit[0] = PGRRsort[(len(PGRR)/(bincount+30*k)), 0]
			for i in range(1,bincount+30*(k)):				
				RainRateBinLimit[i] = 0.5*(PGRRsort[(len(PGRR)/(bincount+30*k)*i), 0]-PGRRsort[(len(PGRR)/(bincount+30*k)*(i-1)), 0])+PGRRsort[(len(PGRR)/(bincount+30*k)*(i-1)), 0]

			############################################################################
		
		elif selectcase == 5:
			############################################################################
			"""Define the median PG and Tip Times for the statistically significant data
			Variable Bin Spacing, Equal Bin Counts, Perturbed over a larger area"""
			
			#Need to redefine these variables to expand to the pertcount size rather than bincount size. The reason its pertcount+1 because we have an extra bin created due to this method.
			PGTipBinMedianFinal = np.zeros((pertcount+1)+30*(k))
			PGTipBinMedianFinal_se = np.zeros((pertcount+1)+30*(k))
			RainRateBinLimit = np.zeros((pertcount+1)+30*(k))
			
			#Define the jump size for each pertub. i.e. by how many events do we change the lower and upper boundaries in the next virtual bin
			dynlength = (len(PGRRsort)-len(PGRRsort)/(bincount+30*(k)))/pertcount
	
			# for i in xrange(pertcount+1):
				# PGTipBinMedianFinal[i] = np.median(PGRRsort[np.floor(0+dynlength*i):np.floor(len(PGRRsort)/bincount+dynlength*i),1].copy())
				# PGTipBinMedianFinal_se[i] = np.std(PGRRsort[np.floor(0+dynlength*i):np.floor(len(PGRRsort)/bincount+dynlength*i),1].copy())/np.sqrt((pertcount+1)+30*(k)+eps)
			
			for i in xrange(pertcount+1):
				PGTipBinMedianFinal[i] 		= np.median(PG[np.floor(0+dynlength*i):np.floor(len(PGRRsort)/bincount+dynlength*i)].copy())
				PGTipBinMedianFinal_se[i] 	= np.std(PG[np.floor(0+dynlength*i):np.floor(len(PGRRsort)/bincount+dynlength*i)].copy())/np.sqrt((pertcount+1)+30*(k)+eps)
			
			
			#print(PGTipBinMedianFinal)
			
			############################################################################
			#Select values for plotting
			yvalue = np.asarray(PGTipBinMedianFinal)
			yvalue_se = np.asarray(PGTipBinMedianFinal_se)
			amethod = "Virtual_Median"

			#Define the Rain Rate for each bin with the centred values determined as well.
			for i in range(pertcount+1):				
				RainRateBinLimit[i] = 0.5*(PGRRsort[np.floor(0+dynlength*i), 0]+PGRRsort[np.floor(len(PGRRsort)/bincount+dynlength*i)-1, 0])

			#print(RainRateBinLimit)	
				
			############################################################################
		
		else:
			sys.exit("Please select either the Mean (1) or Median (2) case.")
		
		#print("Bin Counts", PGTipPosition)
				
		#Calculation of the linear regression model along with statistical parameters.
		slope[k], intercept[k], r_value[k], p_value[k], std_err[k] = stats.linregress(RainRateBinLimit, yvalue)
		#print("RainRateBinLimit", RainRateBinLimit)
		#print("yvalue", yvalue)
		
		#try:
		#	for m in xrange(len(lowess(RainRateBinLimit+eps, yvalue+eps, 1/2))):
		#		lowessval[k,m] = lowess(RainRateBinLimit+eps, yvalue+eps, 1/2)[m]	
		#except:	
		lowessval[k,:] = 0
	
	#print("Method", method)
	
	try:
		PGEnsembleData = zip(RainRateBinLimit, yvalue, yvalue_se)
	except:
		PGEnsembleData = zip(RainRateBinLimit, yvalue, np.zeros(len(yvalue)))
		
	# with open(Storage_Path + Data_Path + 'PGEnsembleData_Normal_' + str(PGRR_Loc) + str(amethod) + str(bincount)+ "-" + str(pertcount) + str(method) + ".csv", "wb") as output:
			# writer = csv.writer(output, lineterminator='\n')
			# writer.writerows(PGEnsembleData)
	
	# PGEnsembleLowess = zip(RainRateBinLimit, lowessval[0,:])
	# with open(Storage_Path + Data_Path + 'PGEnsembleLowess_Normal_' + str(PGRR_Loc) + str(amethod) + str(bincount) + "-" + str(pertcount) + str(method) + ".csv", "wb") as output:
			# writer = csv.writer(output, lineterminator='\n')
			# writer.writerows(PGEnsembleLowess)

	# if selectcase == 4:
		# with open(Storage_Path + Data_Path + 'PGEnsembleBinValues_' + str(PGRR_Loc) + str(amethod) + str(bincount) + "-" + str(pertcount) + str(method) + ".csv", "wb") as output:
			# writer = csv.writer(output, lineterminator='\n')
			# writer.writerows(PGTipBinValues)		
			
	############################################################################
	"""Plot the ensemble PG against Rain Rate. See external.py for the source function."""
	
	# plot = PGRR_Plots()
	# try:
		# plot.PGRainSlim(np.max(RainRateBinLimit)+0.2, np.max(yvalue)+0.2, Storage_Path + Plots_Path, 'PGEnsemble_' + str(PGRR_Loc) + str(amethod) + str(bincount) + str(case) + str(pertcount) + str(method), "png", RainRateBinLimit, yvalue, os.path.abspath(os.path.dirname(__file__)))
	# except:
		# print("Plotting Error")
	
	#img = Image.open('../Plots/Ensemble/PGEnsemble_Normal_v26_' + str(amethod) + str(bincount) + str(case) + str(pertcount) + '.png')
	#img.show()
	
	return np.asarray(PGEnsembleData)
	
if __name__ == "__main__":

	y = "y"
	n = "n"
	############################################################################
	"""User input for the further processing of the PGRR data."""
	print("####################################################################")
	print("The Lenard effect of preciptation at the RUAO. Using the processed ")
	print("data collected from the PGRainRate.py script the average for each ")
	print("rain rate can be found.")
	print("####################################################################\n")
	selectcase = input("Please select the averaging method: Type '1' for Mean, Type '2' for Median: ")
	loop = str(input('Do you want to loop over many bins? y/n: '))
	if loop == "n":
		bincount = input("How many bins for the averaging would you like (recommended = 100): ")	
		loop = 1
	elif loop == "y":
		bincount = 30
		loop = 10
	if selectcase == 5:
		pertcount = input("How many times do you want to perturb over the data. The higher the number the better the resolution (recommended = 100):")
	else:
		pertcount=None

	############################################################################
	"""Launch the ensembler"""
	PGRR_Ensembler(bincount, loop, selectcase,None,None,'Chil_PG_LWC_Median', pertcount)
	