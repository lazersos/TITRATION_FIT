#!/usr/bin/env python3
"""
This class provides a quick fit routine for 
a titration curve using the method of 
equal areas.
Lazer Industries
'Makeing tomorrow today'
"""

# Libraries
import sys
import numpy as np

class TITRATION():
	"""Class which handles doing the titration fit

	"""
	def __init__(self):
		self.filename = None
		self.lverb = True
		self.xdata = None
		self.ydata = None
		self.fit_func = lambda x,a,b,c,d,e:0.5*a*(1.0-np.tanh((x-b)/c))+d*x+e
		self.coefs = None
		self.XLabel = 'Volume [mL]'
		self.YLabel = 'Potential [mV]'

	def titrationFunc(self,x,a,b,c,d,e):
		""" Test function for titration curve
		"""
		H = -np.power(10,x)
		H2 = H*H
		H3 = H2*H
		#top = H3 +      Ka * H2 - (Kw+Ka*Ca) * H - Ka*Kw
		#bot = H3 + (Ka+Cd) * H2 + (Ka*Cb-Kw) * H - Ka*Kw
		#Vb = -Va * top / bot
		top = H3 +      a  * H2 + b * H - c
		bot = H3 +      d  * H2 + e * H - c
		Vb  =  e * top / bot
		return Vb

	def setVerbosity(self,lverb):
		""" Sets the verbosity of the class

		This routine sets the verbosity of the class

		Parameters
		----------
		lverb : boolean
			Verbosity of the class
		"""
		self.lverb = lverb

	def setFileName(self,filename):
		""" Sets the Filename to work on

		This routine filename used for the fit data

		Parameters
		----------
		filename : string
			Filename to operate on
		"""
		if self.lverb: print('  DATA FILE : '+str(filename))
		self.filename = filename

	def setData(self,xdata=None,ydata=None):
		""" Sets the x and y data to fit

		This routine sets the x andy data for the fitting

		Parameters
		----------
		x : real
			X data in [mL]
		y : real
			Y data in [mV]
		"""
		if (type(xdata) is type(None)) or (type(ydata) is type(None)): 
			if self.lverb: print('!!!!! Issue with XData or YData in setData !!!!!')
			return
		self.xdata = xdata
		self.ydata = ydata
		if self.lverb:
			print('  XDATA = ['+str(min(xdata))+','+str(max(xdata))+'] : LEN = '+str(len(xdata)))
			print('  YDATA = ['+str(min(ydata))+','+str(max(ydata))+'] : LEN = '+str(len(ydata)))

	def readFile(self):
		""" Reads the data in filename

		Reads the data from a CSV formated file
		"""
		if ".csv" not in self.filename:
			if self.lverb: print('!!!!!  '+str(self.filename)+ ' may not be a CSV file  !!!!!')
			return
		temparr = np.genfromtxt(self.filename, delimiter=',')
		self.setData(xdata=temparr[:,0],ydata=temparr[:,1])

	def fitData(self):
		""" Fit the data

		Fits the data using the function
		"""
		from scipy.optimize import curve_fit
		if (type(self.xdata) is type(None)) or (type(self.ydata) is type(None)): return
		bounds =([min(self.ydata),min(self.xdata),1E-6,-5,0],[2.*max(self.ydata),max(self.xdata),max(self.xdata)-min(self.xdata),5,max(self.ydata)])
		p0 = np.array([max(self.ydata),np.mean(self.xdata),1,-0.01,10])
		popt, pcov = curve_fit(self.fit_func, self.xdata, self.ydata,p0=p0,bounds=bounds)
		self.coefs = popt

	def makePlot(self):
		""" Plot the data

		Creates a plot of the data
		"""
		import matplotlib.pyplot as plt
		plt.plot(self.xdata, self.ydata, 'ko', label='Data')
		x2 = np.linspace(0,max(self.xdata)*1.2,1000)
		plt.plot(x2,self.fit_func(x2, *self.coefs), 'r-', label='Fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f, e=%5.3f' % tuple(self.coefs))
		plt.xlabel(self.XLabel)
		plt.ylabel(self.YLabel)
		plt.legend()
		plt.show()



if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description='''A tool/library to fit
		a titration curve to a provided dataset.  A user provided
		dataset in CSV format with the first column being 
		volume in mL and the second column being potential
		in mV.''')
	parser.add_argument("--file", dest="filename",
		help="Filename containting data (CSV format)", default = None)
	parser.add_argument("-q", "--quiet", dest="lquiet", action='store_true',
		help="Verbose screen output.", default = False)
	args = parser.parse_args()
	lverb = not args.lquiet
	if (lverb): print('-----  TITRATION FITTING PROGRAM -----')
	FIT = TITRATION()
	FIT.setVerbosity(lverb)
	FIT.setFileName(args.filename)
	FIT.readFile()
	FIT.fitData()
	FIT.makePlot()
	if (lverb): print('--------------------------------------')
	sys.exit(0)


