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
import StringIO
import glob
import re


def grid_plot(calcs,xaxis,yaxis,zaxis='Energy',colorbarUnits=None,zaxis_title=None,plot_title=None,xAxisStr=None,yAxisStr=None,fileName='grid_plot.pdf',runlocal=True):

    calcs_aux=copy.deepcopy(calcs)
    if zaxis=='Energy':

        try:
            calcs = AFLOWpi.retr.grabEnergyOut(calcs)
            zaxis='aux_energy'            
            for ID,oneCalc in calcs.iteritems():
                calcs[ID][zaxis]=calcs[ID]['Energy']
        except:
            pass


        if zaxis_title==None:
            zaxis_title='Energy (Ry)'

    if len(calcs)>1 and type(calcs)==type([1,2,3]):
        if zaxis_title!=None:
            if type(zaxis_title)==type([1,2,3]):
                pass
            else:
                zaxis_title=[zaxis_title]

        calcs_aux=calcs[0]
        AFLOWpi.plot.__distortionEnergy(calcs[0],xaxis,yaxis,zaxis=zaxis,calcs=calcs[1:],colorbarUnits=colorbarUnits,titleArray=zaxis_title,plotTitle=plot_title,xAxisStr=xAxisStr,yAxisStr=yAxisStr,fileName=fileName,percentage=False)

        
    else:
        for ID,oneCalc in calcs.iteritems():
            calcs_aux[ID][zaxis]=0.0

        AFLOWpi.plot.__distortionEnergy(calcs_aux,xaxis,yaxis,zaxis=zaxis,calcs=[calcs],colorbarUnits=colorbarUnits,titleArray=[zaxis_title],plotTitle=plotTitle,xAxisStr=xAxisStr,yAxisStr=yAxisStr,fileName=fileName,percentage=False)

def distortionEnergy(calcs1,xaxis,yaxis,zaxis='Energy',calcs=None,colorbarUnits=None,titleArray=None,plotTitle=None,xAxisStr=None,yAxisStr=None,fileName='distortionEnergy.pdf',percentage=False,runlocal=True):

        if runlocal==True:
            AFLOWpi.plot.__distortionEnergy(calcs1,xaxis,yaxis,zaxis=zaxis,calcs=calcs,colorbarUnits=colorbarUnits,titleArray=titleArray,plotTitle=plotTitle,xAxisStr=xAxisStr,yAxisStr=yAxisStr,fileName=fileName,percentage=percentage)
        else:
            pass


        calc_array=[calcs.extend(calcs1)]
	command = """
	   calcArray=[]

	   labels= ['CUB','RHO','TET']
	   refDict={}
	   import os
	   curDir = os.path.abspath(os.curdir)
	   for label in labels:

		calcs =  AFLOWpi.prep.loadlogs('CRT', label,label,config='/Users/supka/fornari-research-dev/tests/AFLOWpi_tests.config')


	## append the calculations to a list for input them into the plotting function
		calcArray.append(calcs)

	   AFLOWpi.plot.__distortionEnergy(calcArray[0],'_AFLOWPI_A_','_AFLOWPI_B_',calcs=calcArray[1:],plotTitle='$ABO_{3}$ Distortion $\Delta$E',titleArray=[labels[1]+'$\Delta$E Ry',labels[2]+'$\Delta$E Ry'],xAxisStr='A',yAxisStr='B')
	"""
	AFLOWpi.prep.runAfterAllDone(calcArray,command)



def __distortionEnergy(calcs1,xaxis,yaxis,zaxis='Energy',calcs=None,colorbarUnits=None,titleArray=None,plotTitle=None,xAxisStr=None,yAxisStr=None,fileName='distortionEnergy.pdf',percentage=False):
    '''
    Plots the change in energy of a certain chemistry when the structure is changed

    Arguments:
          calcs1 (dict):  calculations of calculations with which to compare calcs2 and calcs3 with
          xaxis  (str):  variable in the calculations you want as yaxis
          yaxis  (str): variable in the calculations you want as yaxis

    Keyword Arguments:
          calcs (list): a set of calcs that has distorted parameters from calcs1 (default:None)
          title (str): title of plot (default:None)
          colorbarUnits (str): the units of the elements in the array.(default:same as xaxis)
          xaxisStr (str): The label for the horizontal axis of the plot.(default:same as yaxis)
          yaxisStr (str): The label for the vertical axis of the plot.(default:None)
          plotTitle (str) The title of the plot.(default:None)
          titleArray (list): an array for the labels for the colorbar for the sets (default:None)
          fileName (str): name (and path where default is directory where script is run from ) 
                          of the output file (default: 'distortionEnergy.pdf')
	
    Returns:
          None

    '''

    if len(calcs)>3:
	    print 'List of distortion calculations must not exceed 3'
	    return

    if xAxisStr==None:
	    xAxisStr=xaxis
    if yAxisStr==None:
	    yAxisStr=yaxis

    
    def grabMatrix(calcs):
	    calcs1 = AFLOWpi.retr.grabEnergyOut(calcs)
	    X=set()
	    Y=set()
	    Z=[]
	    for key,oneCalc, in calcs1.iteritems():
		    X.add(oneCalc[xaxis])
		    Y.add(oneCalc[yaxis])
                    val = oneCalc[zaxis]

		    Z.append(val)

	    Z=numpy.ma.masked_array(Z)    
	    Z=Z.reshape(len(Y),len(X))  
	    return Z


    def grabLabels(calcs,xLabel,yLabel):
	    xVals=[]
	    yVals=[]
	    counter=0
	    for ID,oneCalc in calcs.iteritems():
		    if oneCalc[xLabel] not in xVals:
			    xVals.append(oneCalc[xLabel])
		    if oneCalc[yLabel] not in yVals:
			    yVals.append(oneCalc[yLabel])
	    
	    return [xVals,yVals]

    energy1 = grabMatrix(calcs1)
    energyArrayList=[]
    printBool=0

    fig = pyplot.figure(figsize=(12,6))
    neg_flag=True    

    for item in range(len(calcs)):
        diff=grabMatrix(calcs[item])-energy1
        if percentage==True:
            diff = 100.0*diff/energy1


        try:
            for j in diff:
                for k in j:
                    if k > 0.1:
                        neg_flag=False
        except Exception,e:
            print e

        energyArrayList.append(diff)

    energyArrayListMask = copy.deepcopy(energyArrayList)
    z_max = 0
    z_min = 0

    for item in range(len(energyArrayList)):
	    if energyArrayList[item].max()>z_max:
		    z_max = energyArrayList[item].max()
	    if energyArrayList[item].min()<z_min:
		    z_min = energyArrayList[item].min()
    for item in range(len(energyArrayList)):
	    for item1 in range(len(energyArrayList)):
                

                energyArrayListMask[item] = numpy.ma.masked_where(energyArrayList[item]>energyArrayList[item1], energyArrayListMask[item])
    
    varLabels= grabLabels(calcs1,xaxis,yaxis)

    
    cbcolor = ['Greens_r','Blues_r','Reds_r']
    if neg_flag==False:
        cbcolor = ['Greens','Blues','Reds']

    for item in range(len(energyArrayListMask)):
	    parr=energyArrayListMask[item]

	    plot2 = pyplot.pcolor(parr, vmin=z_min, vmax=z_max,cmap=cbcolor[item])
	    parrMask= numpy.ma.getmaskarray(parr)
	    if printBool:
		    for y in range(parrMask.shape[0]):
			    for x in range(parrMask.shape[1]):
				    if not parrMask[y,x]:
					    pyplot.text(x + 0.5, y + 0.5, '%.4f' % parr[y, x],
							horizontalalignment='center',
							verticalalignment='center',
							)

	    try:
		    pyplot.colorbar().set_label(titleArray[item])
	    except:
		    pass


    pyplot.xlabel(xAxisStr)
    pyplot.ylabel(yAxisStr)
    pyplot.title(plotTitle)
    
    pyplot.xticks(numpy.arange(0.5, len(varLabels[0])), varLabels[0])
    pyplot.yticks(numpy.arange(0.5, len(varLabels[1])), varLabels[1])
    
    pyplot.savefig(fileName,bbox_inches='tight')

    pyplot.cla()
    pyplot.clf()
    pyplot.close()

    AFLOWpi.retr._moveToSavedir(fileName)
