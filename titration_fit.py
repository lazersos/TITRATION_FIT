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
		# d/dx tahn = 1/cosh2 = 1-tanh2
		self.fit_func_slope = lambda x,a,b,c,d,e:-0.5*a*(1.0/(np.cosh((x-b)/c))**2)/c+d
		self.slope_func = lambda x,a,b:a*x+b
		self.slope_func_integrand = lambda x,a,b:0.5*a*x*x+b*x
		self.coefs = None
		self.x0_tp = None
		self.y0_tp = None
		self.slope_tp = None
		self.ycept_tp = None
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

	def calcTitrationPoint(self):
		""" Sets the titration point

		Returnst the titration point and slope from the dataset

		"""
		# Get the derivative information
		x = np.linspace(0,max(self.xdata)*1.2,1000)
		fprime = self.fit_func_slope(x,*self.coefs)
		self.slope_tp = min(fprime) # slope is negative so minimum value
		maxindex = list(fprime).index(self.slope_tp)
		self.x0_tp = x[maxindex]
		self.y0_tp = self.fit_func(self.x0_tp,*self.coefs)
		self.ycept_tp  = self.y0_tp - self.x0_tp*self.slope_tp

	def getPerpLine(self,m,x0,y0):
		""" Get the titration point

		Returnst the titration point and slope from the dataset

		Parameters
		----------
		m  : real
			Slope of the titration line
		x0 : real
			X location of the titration point
		y0 : real
			Y location of the titration point

		Returns
		-------
		m  : real
			Slope of the titration curve at (X0,Y0)
		b  : real
			Y-intercept of line passing through (X0,Y0) with slope m
		"""
		# Normal is -1/dfdx
		m = -1.0/m
		b = y0 - x0*m
		return m,b

	def calcTangentLines(self):
		""" Calculates the paramteters for the tangent lines
		"""
		fact = 2
		x = self.x0_tp-self.coefs[2]*fact
		y = self.fit_func(x,*self.coefs)
		m = self.fit_func_slope(x,*self.coefs)
		b = y - x*m
		self.x1_tan = x
		self.y1_tan = y
		self.slope1_tan = m
		self.ycept1_tan = b 
		x = self.x0_tp+self.coefs[2]*fact
		y = self.fit_func(x,*self.coefs)
		m = self.fit_func_slope(x,*self.coefs)
		b = y - x*m
		self.x2_tan = x
		self.y2_tan = y
		self.slope2_tan = m
		self.ycept2_tan = b

	def calcAreas(self):
		""" Calculates the area of the triangles
		"""
		import scipy.integrate as integrate
		# First calcualte the points where the curves meet
		[m2,b2] = self.getPerpLine(self.slope_tp,self.x0_tp,self.y0_tp)
		xr    = (b2-self.ycept1_tan)/(self.slope1_tan-m2)
		xl    = (b2-self.ycept2_tan)/(self.slope2_tan-m2)
		# Now calculate the areas
		result1 = integrate.quad(lambda x: self.slope_func(x,m2,b2), xl, self.x0_tp)
		result2 = integrate.quad(lambda x: self.slope_func(x,self.slope2_tan,self.ycept2_tan), xl, self.x0_tp)
		self.areal = result1[0]-result2[0]
		result1 = integrate.quad(lambda x: self.slope_func(x,m2,b2), self.x0_tp, xr)
		result2 = integrate.quad(lambda x: self.slope_func(x,self.slope1_tan,self.ycept1_tan), self.x0_tp, xr)
		self.arear = result2[0]-result1[0]
		# Now get the centerpoints of the triangles
		yl = self.slope2_tan*xl+self.ycept2_tan
		yr = self.slope1_tan*xr+self.ycept1_tan
		self.x0_al = xl + (self.x0_tp-xl)/3.0
		self.y0_al = yl - (self.y0_tp-yl)/3.0
		self.x0_ar = xr + (self.x0_tp-xr)/3.0
		self.y0_ar = yr - (self.y0_tp-yr)/3.0





	def makePlot(self):
		""" Plot the data

		Creates a plot of the data
		"""
		import matplotlib.pyplot as plt
		# Plot Dataset
		plt.plot(self.xdata, self.ydata, 'ko', label='Data')
		# Create an x array for plotting  and plot fit function
		x2 = np.linspace(0,max(self.xdata)*1.2,1000)
		plt.plot(x2,self.fit_func(x2, *self.coefs), 'r-', label='Fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f, e=%5.3f' % tuple(self.coefs))
		# Plot the titration point
		plt.plot(self.x0_tp, self.y0_tp, 'bo', label='Titration Point (%5.3f,%5.3f)' % (self.x0_tp,self.y0_tp))
		# Plot the tripple point normal line
		[m2,b2] = self.getPerpLine(self.slope_tp,self.x0_tp,self.y0_tp)
		plt.plot(x2,self.slope_func(x2,m2,b2),'c:')
		# Plot the tangent lines
		plt.plot(x2,self.slope1_tan*x2+self.ycept1_tan, 'b:')
		plt.plot(x2,self.slope2_tan*x2+self.ycept2_tan, 'b:')
		# Annotations
		plt.annotate('%5.3f [mV*mL]' % (self.areal),
            xy=(self.x0_al, self.y0_al), xycoords='data',
            xytext=(-25, 45), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.05),
            horizontalalignment='right', verticalalignment='bottom')
		plt.annotate('%5.3f [mV*mL]' % (self.arear),
            xy=(self.x0_ar, self.y0_ar), xycoords='data',
            xytext=(85, -45), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.05),
            horizontalalignment='right', verticalalignment='bottom')
		# Cleanup the plot
		plt.ylim(bottom=0,top=1.4*max(self.ydata))
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
	FIT.calcTitrationPoint()
	FIT.calcTangentLines()
	FIT.calcAreas()
	FIT.makePlot()
	if (lverb): print('--------------------------------------')
	sys.exit(0)


