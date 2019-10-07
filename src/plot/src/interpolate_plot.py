# ***************************************************************************
# *                                                                         *
# *          AFLOWpi - Central Michigan University University, 2017         *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2017 - Andrew Supka and Marco Fornari - AFLOW.ORG consortium
#
#  This file is part of AFLOWpi software.
#
#  AFLOWpi is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ***************************************************************************

import AFLOWpi
import numpy
import matplotlib
matplotlib.use('PDF')
from matplotlib import pylab
from matplotlib import pyplot
import os
import logging
import io
import glob
import re
import copy
import scipy


def interpolatePlot(calcs,variable1,variable2,zaxis='Energy',xaxisTitle=None, yaxisTitle=None,zaxisTitle=None,title=None,fileName='interpolatePlot.pdf',delta=False,text_min=False,vhline_min=False,circle_min=False,delta_min=True,rel_center=False,plot_color='jet',find_min=False):    
    '''
    Takes a list of calculations and plots the energy of the calculations as a function of two input variables
    the first value is the baseline for the energy value and the energy plotted is the difference between that
    energy and the other energies in the grid

    Arguments:
          calcs (dict): dictionary of dictionaries of calculations
          variable1 (str): a string of the variable in the calculations that you want as your x axis
          variable2 (str): a string of the variable in the calculations that you want to your y axis

    Keyword Arguments:
          title (str): Title of plot (default: None)
          zaxis (str): Choice out of the keywords in each calc to plot in the Z axis (default: Energy)
          xaxisTitle (str): title of xaxis (default: same as variable1)
          yaxisTitle (str): title of yaxis (default: same as variable2)
          zaxisTitle (str): title of zaxis (default: same as zaxis)
          fileName (str): Name (and path where default is directory where script is run from ) of the 
                          output file (default: 'interpolatePlot.pdf')
          delta (bool): Z-axis scale relative to its first value
          delta_min (bool): Z-axis scale relative to its lowest value
          text_min (bool): Display text of minimum value next to the minimum value if find_min=True
          vhline_min (bool): Display text of minimum value next to the minimum value if find_min=True
          circle_min (bool): Display text of minimum value next to the minimum value if find_min=True
          delta_min (bool): Display text of minimum value next to the minimum value if find_min=True
          rel_center (bool): Display text of minimum value next to the minimum value if find_min=True
          plot_color (str): string of the matplotlib colormap to be used 
          find_min (bool): Interpolate between points and find value of 
                           variable1 and variable2 for the minimum value of Z axis

    Returns:
          None

    '''

    if xaxisTitle is None:
	    xaxisTitle=variable1
    if yaxisTitle is None:
	    yaxisTitle=variable2
    X=[]
    Y=[]
    Z=[]
    try:
        if zaxis=='Energy':
            calcs=AFLOWpi.retr.grabEnergyOut(calcs)
    except:
	    pass

    for key,oneCalc, in list(calcs.items()):
        X.append(float(oneCalc[variable1]))
        Y.append(float(oneCalc[variable2]))
        Z.append(float(oneCalc[zaxis]))


    X=numpy.asarray(X)
    Y=numpy.asarray(Y)
    X=numpy.unique(X)
    Y=numpy.unique(Y)
    Z=numpy.array(Z)

    logging.debug(X)
    logging.debug(Y)
    logging.debug(Z)
    try:
        logging.debug(X.shape)
        logging.debug(Y.shape)
        logging.debug(Z.shape)
    except:
        logging.debug("noshape found")
    x_bkup=copy.deepcopy(X)
    y_bkup=copy.deepcopy(Y)

#    X=numpy.around(X,decimals=7)
#    Y=numpy.around(Y,decimals=7)

    xmin,xmax = numpy.amin(X),numpy.amax(X)
    ymin,ymax = numpy.amin(Y),numpy.amax(Y)


    if delta_min==True:
	    Z-=numpy.amin(Z)
    z_min, z_max = Z.min(), Z.max()

    Z=numpy.reshape(Z,(X.size,int(float(Z.size)/float(X.size))))

    if delta==True:
	    X, Y = numpy.meshgrid(X-X[0], Y-Y[0])
    else:
	    X, Y = numpy.meshgrid(X, Y)

    Z=numpy.ma.masked_values(Z,0.0).T

    fig = pyplot.figure(figsize=(10,6 ))
    interp='bicubic'
    ax = fig.add_subplot(111)


    if zaxis=='Energy':
        try:
                if find_min==True:
                    valEnergy = AFLOWpi.pseudo._getMinimization(calcs,fitVars=[variable1,variable2],return_energy=True)
                    plotVals=valEnergy[0]

                if rel_center==True:
                        pyplot.plot([(xmax+xmin)/2.0 ,plotVals[variable1]], [(ymax+ymin)/2.0 ,plotVals[variable2] ], 'k-')
                if text_min==True:
                        coordText='%s,\n%s'%(numpy.around(plotVals[variable1],decimals=6),numpy.around(plotVals[variable2],decimals=6))
                        pyplot.figtext(0.98*plotVals[variable1],0.98*plotVals[variable2],coordText,color='k')
                if circle_min==True:
                        circ = matplotlib.patches.Ellipse((plotVals[variable1],plotVals[variable2] ), width=0.01000000000*numpy.abs(xmax-xmin),height=0.01000000000*numpy.abs(ymax-ymin), color='k')
                        ax.add_patch(circ)
                        if rel_center==True:
                                circ = matplotlib.patches.Ellipse(((xmax+xmin)/2,(ymax+ymin)/2 ),  width=0.01000000000*numpy.abs(xmax-xmin),height=0.01000000000*numpy.abs(ymax-ymin), color='k')
                                ax.add_patch(circ)


        except Exception as e:
                AFLOWpi.run._fancy_error_log(e)    
    
    if int(matplotlib.__version__[0])<2:
        im = pylab.imshow(Z,cmap=plot_color,aspect='equal',vmin=z_min,vmax=z_max,
                          interpolation=interp,extent=[X.min(), X.max(), Y.min(), Y.max()],origin='lower',hold=True)
    else:
        im = pylab.imshow(Z,cmap=plot_color,aspect='equal',vmin=z_min,vmax=z_max,
                          interpolation=interp,extent=[X.min(), X.max(), Y.min(), Y.max()],origin='lower')

    cbar = pylab.colorbar()

    if zaxisTitle!=None:
	    cbar.set_label(zaxisTitle,size=18)
    
    ax.set_xlabel(xaxisTitle)
    ax.set_ylabel(yaxisTitle)

    if title!=None:
	    pyplot.title(title)

    pyplot.savefig(fileName,bbox_inches='tight')
    try:
	    AFLOWpi.retr._moveToSavedir(fileName)
    except Exception as e:
	    pass



def interpolatePlot1D(calcs,variable1,yaxis='Energy',xaxisTitle=None, yaxisTitle=None,title=None,fileName='interpolatePlot.pdf',delta=False,circle_min=False):    
    '''
    Takes a list of calculations and plots the energy of the calculations as a function of two input variables
    the first value is the baseline for the energy value and the energy plotted is the difference between that
    energy and the other energies in the grid

    Arguments:
          calcs (dict): dictionary of dictionaries of calculations
          variable1 (dict): a string of the variable in the calculations that you want as your x axis

    Keyword Arguments:
          yaxis (str): Choice out of the keywords in each calc to plot in the Z axis (default: Energy)
          xaxisTitle (str): title of xaxis (default: same as variable1)
          yaxisTitle (str): title of yaxis (default: same as yaxis)
          title (str): Title of plot (default: None)
          fileName (str): Name (and path where default is directory where script is run from ) of the 
                          output file (default: 'interpolatePlot.pdf')
          delta (bool): Z-axis scale relative to its first value
          circle_min (bool): Display text of minimum value next to the minimum value

    '''

    if xaxisTitle is None:
	    xaxisTitle=variable1
    if yaxisTitle is None:
	    yaxisTitle=yaxis
    X=[]
    Y=[]

    try:
	    calcs=AFLOWpi.retr.grabEnergyOut(calcs)
    except:
	    pass

    for key,oneCalc, in list(calcs.items()):
        X.append(oneCalc[variable1])
        Y.append(oneCalc[yaxis])
	
    xmin,xmax = numpy.amin(X),numpy.amax(X)
        
    Y=numpy.array(Y)
    Y=numpy.ma.masked_values(Y,0.0)
    origYMin=numpy.amin(Y)
    origYMax=numpy.amax(Y)
    Y-=numpy.amin(Y)

    fig = pyplot.figure()

    ax = fig.add_subplot(111)

    x_interp = numpy.linspace(xmin, xmax, num=100, endpoint=True)
    interp_func = scipy.interpolate.interp1d(X, Y, kind='cubic')

    try:
	    valEnergy = AFLOWpi.pseudo._getMinimization(calcs,fitVars=[variable1],return_energy=True)
	    plotVals=valEnergy[0]
	    foundMin=plotVals[yaxis]-origYMin
	    
	    minBound=numpy.amin(Y)
	    if minBound>foundMin:
		    minBound=foundMin

	    minBound-=numpy.abs((numpy.amax(Y)-numpy.amin(Y)))*0.02

	    pyplot.axis([xmin,xmax,minBound,numpy.amax(Y)*1.01])
	    pyplot.plot(X,Y,'o',x_interp,interp_func(x_interp),'-',[plotVals[variable1]],[foundMin],'x')

	    pyplot.legend(['data','interpolated','min'], loc='best')	    

    except Exception as e:
	    AFLOWpi.run._fancy_error_log(e)    
    
    
    ax.set_xlabel(xaxisTitle)
    ax.set_ylabel(yaxisTitle)

    if title!=None:
	    pyplot.title(title)

	    
    pyplot.savefig(fileName,bbox_inches='tight')
    try:
	    AFLOWpi.retr._moveToSavedir(fileName)
    except Exception as e:
	    AFLOWpi.run._fancy_error_log(e)

