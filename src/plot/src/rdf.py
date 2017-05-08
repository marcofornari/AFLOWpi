# ***************************************************************************
# *                                                                         *
# *          AFLOWpi - Central Michigan University University, 2017         *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2017 - Andrew Supka and Marco Fornari - AFLOW.ORG consortium
#
#  This file is part of AFLOW software.
#
#  AFLOW is free software: you can redistribute it and/or modify
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

import os
import datetime
import cPickle
import logging 
import matplotlib
import glob
matplotlib.use('PDF')
from matplotlib import pylab
import numpy 
import StringIO
import copy 
import re 
import AFLOWpi.prep
import AFLOWpi.retr
import matplotlib.patches
import matplotlib.image 
import pprint
import math
import collections
from matplotlib import pyplot
import scipy.interpolate



def radialPDF(calcs,atomNum,filterElement=None,runlocal=False,inpt=False,outp=True,title='',n_bins=30,y_range=[0,3],**kwargs):
	'''
	kwargs get passed to pyplot.hist
	'''
        AFLOWpi.retr.atomicDistances(calcs,runlocal=runlocal,inpt=inpt,outp=outp)
	returnDict = {}
	for ID,oneCalc in calcs.iteritems():
		if runlocal:
                    if outp==True:
			figObj = AFLOWpi.plot.__radialPDF(oneCalc,ID,atomNum,filterElement=filterElement,outp=True,title=title,**kwargs)
                        figObj = AFLOWpi.plot.__radialPDF(oneCalc,ID,atomNum,filterElement=filterElement,outp=False,title=title,**kwargs)                    

			oneCalc.update({'radialPDF_atom%s' % atomNum:figObj})
		else:
			kwargsStringList = ['%s = %s' % (key,value) for key,value in kwargs.iteritems() if key != 'runlocal' and type(value) != type('string')]
			kwargsStringList.extend(["%s = '%s'" % (key,value) for key,value in kwargs.iteritems() if key != 'runlocal' and type(value) == type('string')])

			kwargsString = ','.join(kwargsStringList)
			kwargsString+=',title="%s",' % title
			if outp==True:
                            AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__radialPDF(oneCalc,ID,outp=True,%s,%s)" % (atomNum,kwargsString))	
                        if outp==False:
                            AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__radialPDF(oneCalc,ID,outp=False,%s,%s)" % (atomNum,kwargsString))	

	return calcs



###################################################################################################
def __radialPDF(oneCalc,ID,atomNum,filterElement=None,title='',file_prefix='',file_postfix='',outp=True,y_range=None,**kwargs):
    """
    kwargs get passed onto the hist function inside
    """
    try:
        if outp==True:
            inOrOut='output'
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_output_dist.out' % ID),'r') as inFile:
                inFileString = inFile.read()
        else:
            inOrOut='input'
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_input_dist.out' % ID),'r') as inFile:
                inFileString = inFile.read()

    except Exception,e:
        print e
        print e
        print e
        logging.error('Could not find %s_dist.out in %s. Are you sure you ran AFLOWpi.retr.atomicDistances?' % (ID,oneCalc['_AFLOWPI_FOLDER_']))
        print 'Could not find %s_dist.out in %s. Are you sure you  ran AFLOWpi.retr.atomicDistances?' % (ID,oneCalc['_AFLOWPI_FOLDER_'])
	return
    fig = pyplot.figure()
    ax1 = pyplot.subplot(111)
    inFileStringSplit = inFileString.split('\n')

    
    splitLines = sorted([line.split()[:4] for line in inFileString.split('\n') if len(line.split()) < 7])

    lineList=[]
    for line in splitLines:
        try:
            if int(line[0]) == atomNum:
		    if filterElement!=None:
			    if filterElement==line[2].split('-')[1]:
				    lineList.append(float(line[3]))
		    else:
			    lineList.append(float(line[3]))			    
        except:
            pass


    bins=numpy.linspace(float(y_range[0]),float(y_range[1]),float(n_bins))

    try:
	    del kwargs['n_bins']
	    del kwargs['y_range']
	    del kwargs['filterElements']
    except KeyError:
	    pass

    
    n, bins, patches = pyplot.hist(lineList,bins=bins,**kwargs)

    if y_range!=None:
        try:
            axes = pyplot.gca()
            axes.set_ylim([y_range[0],y_range[1]])
        except Exception,e:
            print e
    newbin = []
    for entry in range(len(bins)):
	    try:
		    newbin.append((bins[entry+1]+bins[entry])/2)
	    except Exception,e:
		    pass

    if file_prefix!='':
        file_prefix+='_'

    if file_postfix!='':
        file_postfix='_'+file_postfix


    fileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%sradialPDF_%s%s_atom_%s_%s%s.pdf' % (file_prefix,AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atomNum,inOrOut,file_postfix))


    if title=='':
	    figtitle = '%s%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID) 
    else:
	    figtitle=title
    
    t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=14,horizontalalignment='center') #[x,y]
    newplot = pyplot.plot(newbin,n,linestyle=' ',marker='')
    pylab.grid(color='k', linestyle='--')
    pylab.xlabel('Radial Distance from Atom ($\AA$)')
    pylab.ylabel('No. of Atoms')
    plotObj=newplot

    pyplot.savefig(fileName)
   
    try:
	    AFLOWpi.retr._moveToSavedir(fileName)
    except Exception,e:
	    pass

    AFLOWpi.retr._moveToSavedir(fileName)
    pyplot.cla()
    pyplot.clf()
    pyplot.close()

    return lineList


###################################################################################################
