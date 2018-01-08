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
matplotlib.use('pdf')
from matplotlib import pylab
from matplotlib import pyplot
import os
import logging
import StringIO
import glob
import re
import cPickle
import matplotlib.lines as mlines


def band_topology(calcs,yLim=[-10,10],DOSPlot='',runlocal=False,postfix='',tight_banding=False):
	'''
	Generates electronic band structure plots for the calculations in the dictionary of dictionaries
	of calculations with the option to have a DOS or PDOS plot to accompany it.

	Arguments:
              calcs (dict): dictionary of dictionaries representing the set of calculations
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)
              DOSPlot (str): DOS or the PDOS plots next to the eletronic band structure plot (assuming you ran 
                             either ppDOS or ppPDOS) (default: NONE) 
              postfix (str): Postfix to the filename of the plot
              tight_banding (bool): Whether to treat the input data as from Quantum Espresso or WanT bands.x

        Returns:
              None

	'''

	for oneCalc in calcs.items():
		if runlocal:
			AFLOWpi.plot.__band_topology_Plot(oneCalc,DOSPlot,postfix=postfix,tight_banding=tight_banding)
		else:
			AFLOWpi.prep._addToBlock(oneCalc[1],oneCalc[0],'PLOT',"AFLOWpi.plot.__band_topology(oneCalc,ID,DOSPlot='%s',postfix='%s',tight_banding=%s)" % (DOSPlot,postfix,tight_banding))	


def __band_topology(oneCalc,ID,DOSPlot='',postfix='',tight_banding=False):
	'''
        Wrapper function for AFLOWpi.plot.__bandPlot 
        OBSOLETE. NEEDS REMOVAL


	Arguments:
              oneCalc (dict): Single calculation that is the value of the dictionary of dictionaries of calculations
              ID (str): ID of the calculation
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)
              DOSPlot (str): DOS or the PDOS plots next to the eletronic band structure plot (assuming you ran 
                             either ppDOS or ppPDOS) (default: NONE) 
              postfix (str): Postfix to the filename of the plot
              tight_banding (bool): Whether to treat the input data as from Quantum Espresso or WanT bands.x

        Returns:
              None

	'''



	ipol=0
	jpol=1
	spol=2

        oneCalc = (ID,oneCalc)

	for topo_type in ['BC','SBC']:


		if topo_type=='SBC':
			dem_files = glob.glob(oneCalc[1]['_AFLOWPI_FOLDER_']+'/Omegaj_*')
		else:
			dem_files = glob.glob(oneCalc[1]['_AFLOWPI_FOLDER_']+'/Omega_*')

		dir_dict = {'x':0,'y':1,'z':2}

		for dat_file in dem_files:

			dat_name = dat_file[:-4]
			dat_split_name = dat_name.split('_')

			ipol=dir_dict[dat_split_name[-1][0]]
			jpol=dir_dict[dat_split_name[-1][1]]
			if topo_type=='SBC':
			    spol=dir_dict[dat_split_name[-2]]
			else:
			    spol=2

			AFLOWpi.plot.__band_topology_Plot(oneCalc,ipol,jpol,spol=spol,postfix=postfix,
							  tight_banding=tight_banding,topo_type=topo_type)


def __getPath_Topo(oneCalc,ID):
    '''


    Arguments:


    Keyword Arguments:

        
    Returns:


    '''

    try:
        want_stdout_path = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/kpath_points.txt')[-1]
    except:
        want_stdout_path = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_bands_up.out'%ID)[-1]

    with open(want_stdout_path,"r") as ofo:
	    lines=ofo.readlines()

    output_path_string=""
    flag=False
    points_list=[]
    for l in lines:
	    if len(l.strip())==0:
		    flag=True
	    if flag==False:
		    lspl=l.split()
		    output_path_string+="0.0 0.0 0.0 %s ! %s\n"%(lspl[1],lspl[0])
	    else:
		    points_list.extend([float(x) for x in l.split()])

    points=numpy.reshape(numpy.asarray(points_list),(len(points_list)/3,3))
    

    r = numpy.diff(points,axis=0)

    dist=numpy.cumsum(numpy.sqrt(numpy.sum(r**2,axis=1)))
    dist = numpy.concatenate((numpy.array([0.0]),dist),axis=0)
    

    calcID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='PAO-TB',last=True)

    nspin=2
    try:
	    with open("bands_1.dat","r") as ofo:
		    by_band = numpy.array([map(float,x.split()) for x in ofo.readlines()]).T
	    ofs=""
	    for band in xrange(by_band.shape[0]):
		    for kpt in xrange(by_band.shape[1]):
			    ofs+="%s %s\n"%(dist[kpt],by_band[band,kpt])
		    if band!=by_band.shape[0]-1:
			    ofs+="\n"

	    filebands = os.path.join(oneCalc["_AFLOWPI_FOLDER_"],'%s_bands_paopy_down_cleaned.dat'%calcID)
	    with open(filebands,"w") as ofo:
		    ofo.write(ofs)
	    
    except:
	    nspin=1

    if nspin==2:
	    filebands = os.path.join(oneCalc["_AFLOWPI_FOLDER_"],'%s_bands_paopy_up_cleaned.dat'%calcID)
    else:
	    filebands = os.path.join(oneCalc["_AFLOWPI_FOLDER_"],'%s_bands_paopy_cleaned.dat'%calcID)

    with open("bands_0.dat","r") as ofo:
	    by_band = numpy.array([map(float,x.split()) for x in ofo.readlines()]).T



    try:
	    ofs=""
	    for band in xrange(1,by_band.shape[0]):
		    for kpt in xrange(by_band.shape[1]):
			    ofs+="%s %s\n"%(dist[kpt],by_band[band,kpt])
		    if band!=by_band.shape[0]-1:
			    ofs+="\n"	    

	    with open(filebands,"w") as ofo:
		    ofo.write(ofs)
    except Exception,e:
	    AFLOWpi.run._fancy_error_log(e)
	    raise SystemExit
	    pass
	    
    return  output_path_string



def _clean_topo_bands(filebands,Efermi_shift):
     	x = []
	y = []
	k_x = []
	k_y = []


	'''
        looks at the filebands file and reads in the columns for the energy of the
        band and the value of the k point to make the path. this looks at the output
        of the plot_bands.x script and when it finds a blank line in the data it 
        appends a list of k point values and energy values for each band into a list
        of the bands
	'''

	try:

		data = numpy.loadtxt(filebands)

		k_y = [data[:,1].tolist()]
		k_x = [data[:,0].tolist()]



	

	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)
		logging.warning("output from bands calculation not found. Are you sure you ran ppBands and it completed properly?")
		print "Are you sure you ran ppBands and it completed properly?"
		return
	'''
        Eliminating the gaps between the paths in band structure plots by looking 
        for gaps in the data in the k point values that are greater than some value
	'''
 	gapThreshold = 2.00


 	for band in range(len(k_x)):
 		for kpoint in range(2,len(k_x[band])):
			if k_x[band][kpoint] - k_x[band][kpoint-1] > (k_x[band][kpoint-1] - k_x[band][kpoint-2])*gapThreshold:
				if k_x[band][kpoint-1] - k_x[band][kpoint-2]!=0:
					difference = k_x[band][kpoint] - k_x[band][kpoint-1]
					higher_vals = k_x[band][kpoint:]
				else:
					difference=0
			
				k_x[band][kpoint:] = [x - difference for x in k_x[band][kpoint:]]
		
				




	return k_x,k_y

def __band_topology_Plot(oneCalc,ipol,jpol,spol=2,yLim=[-10,10],DOSPlot='',postfix='',tight_banding=False,topo_type='AHC'): 
	'''
        Function to take the data files generated by the sumpdos ppBands functions and plots 
        the electronic band structure and the projected density of states with energy shifted 
        relative to the Fermi Energy.


	Arguments:
              oneCalc (tuple): Single calculation dictionary and the ID of the calculation 
                               NEEDS TO BE CHANGED TO oneCalc,ID FORMAT SOON
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)
              DOSPlot (str): DOS or the PDOS plots next to the eletronic band structure plot (assuming you ran 
                             either ppDOS or ppPDOS) (default: NONE) 
              tight_banding (bool): Whether to treat the input data as from Quantum Espresso or WanT bands.x
              postfix (str): Postfix to the filename of the plot

        Returns:
              None

	'''
	

	extension='png'


	calcCopy = oneCalc
	calcID = oneCalc[0]
	oneCalc = oneCalc[1]
        if tight_banding==True:
            calcID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='PAO-TB',last=True)
        else:
            calcID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='bands',last=True)

        LSDA=False
        nspin = int(AFLOWpi.scfuj.chkSpinCalc(oneCalc,calcID))
        if nspin!=1:
            LSDA=True


	SOC=False




	if SOC: 
		LSDA=False
		nspin=1

        
	if DOSPlot != '' and DOSPlot != 'APDOS' and DOSPlot != 'DOS':
		print "Not a valid choice for DOSPlot. Valid options are:'APDOS','DOS'"
		return

	if postfix!='':
		postfix='_'+postfix

	if DOSPlot=='APDOS' or DOSPlot=='DOS':
		AFLOWpi.plot.__sumpdos(oneCalc,calcID,TB=tight_banding)	
	

	if tight_banding:
		bandSym = AFLOWpi.plot.__getPath_Topo(oneCalc,calcID)
	else:
		bandSym = AFLOWpi.retr._getPathFromFile(calcCopy)
	
	if bandSym==None:
		print 'ERRORRRR!'
		return
	try:
                if tight_banding:
		#	Efermi=AFLOWpi.retr._getEfermi(oneCalc,'%s_WanT_dos'%calcID,directID=True)

			Efermi=0.0
                else:
			Efermi=AFLOWpi.retr._getEfermi(oneCalc,calcID)


		if type(Efermi)!=type(0.5):
			Efermi=Efermi[0]

	except Exception,e:
            print e
            Efermi=0.0


	dir_dict = ['x','y','z']

	subdir=oneCalc['_AFLOWPI_FOLDER_']
	if topo_type=='SBC':
		fileplot = os.path.join(subdir,'SPIN_BERRY_CURVATURE_%s_%s%s_%s_%s%s.%s' % (dir_dict[spol],dir_dict[ipol],dir_dict[jpol],AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix,extension))
	if topo_type=='BC':
		fileplot = os.path.join(subdir,'BERRY_CURVATURE_%s%s_%s_%s%s.%s' % (dir_dict[ipol],dir_dict[jpol],AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix,extension))

	"""get the path to the subdirectory of the calc that you are making plots for"""


	#	try:
#			AFLOWpi.prep._clean_want_bands(oneCalc,calcID)
#		except:
#			return



	if topo_type=='SBC':
		filebands = os.path.join(subdir,'Omegaj_%s_%s%s.dat'%(dir_dict[spol],dir_dict[ipol],dir_dict[jpol]))
	if topo_type=='BC':
		filebands = os.path.join(subdir,'Omega_%s%s.dat'%(dir_dict[ipol],dir_dict[jpol]))
			
	if not os.path.exists(filebands):
		return

	Efermi_shift=Efermi






        
	

	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
	

       #to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=20)             #to set the font size

        
        width = 20
	height = 14
	pylab.figure(figsize=(width, height))#to adjust the figure size
	
     #to do the gaussian smoothing               
 
     #################################

	if nspin==2:
		k_x_up,k_y_up = AFLOWpi.plot._clean_topobands(filebands_up,Efermi_shift)
		k_x_dn,k_y_dn = AFLOWpi.plot._clean_topo_bands(filebands_dn,Efermi_shift)
		k_x,k_y=k_x_up,k_y_up
	else:
		k_x,k_y = AFLOWpi.plot._clean_topo_bands(filebands,Efermi_shift)


	a=k_x[0]
	b=k_x[0]

	k_y = numpy.array(k_y)
	
	pylab.xlim(min(k_x[0]),max(k_x[0])) 
	yLim = [numpy.amin(k_y)*1.1,numpy.amax(k_y)*1.1]
	pylab.ylim((yLim[0],yLim[1]))


	if DOSPlot != '':
		ax1=pylab.subplot(121)

	if DOSPlot == '':
		ax1=pylab.subplot(111)	


	'''
        Plot each band (k_x[i]),(k_y[i]) on the band structure plot from list of
        the values of energy and position from their respective list by k points
	'''
	if nspin==2:
		for i in range(len(k_x_up)):
			try:
				
				if tight_banding==True:
					pylab.plot((k_x_up[i]),(k_y_up[i]),'r',alpha=1.0,marker=".",linestyle=" ",label="$\uparrow$",linewidth=2)
					pylab.plot((k_x_dn[i]),(k_y_dn[i]),'k',alpha=1.0,marker=".",linestyle=" ",label="$\downarrow$",linewidth=2)
				else:
					pylab.plot((k_x_up[i]),(k_y_up[i]),'r',alpha=1.0,marker=".",linestyle=" ",label="$\uparrow$",linewidth=2)
					pylab.plot((k_x_dn[i]),(k_y_dn[i]),'k',alpha=1.0,marker=".",linestyle=" ",label="$\downarrow$",linewidth=2)
			except:
				pass
		handles, labels = ax1.get_legend_handles_labels()
#		ax1.legend(handles[-2:], labels[-2:],numpoints=1)
		up_legend_lab = mlines.Line2D([], [], color='red',label='$\uparrow$')
		dn_legend_lab = mlines.Line2D([], [], color='black',label='$\downarrow$')
		ax1.legend(handles=[up_legend_lab,dn_legend_lab])
	else:
		SOC=False



		if tight_banding==True:
			for i in range(len(k_x)):
				if topo_type=='BC':
					pylab.plot((k_x[i]),(k_y[i]),'k',marker=" ",linewidth=3.0,color='darkorange')
				if topo_type=='SBC':
#					k_ys =floatdos=AFLOWpi.plot.__smoothGauss(k_y[i])
#					k_xs =floatdos=AFLOWpi.plot.__smoothGauss(k_x[i])
					pylab.plot(k_x[i],k_y[i],'k',marker=" ",linewidth=3.0,color='darkred')
	#


	if topo_type=='BC':
		pylab.ylabel(r"$\Omega_{%s%s}(\vec{k}$) (a.u.)"%(dir_dict[ipol],dir_dict[jpol]),
			     fontsize=20)
	if topo_type=='SBC':
		pylab.ylabel(r"$\Omega^{%s}_{%s%s}(\vec{k}$) (a.u.)"%(dir_dict[spol],dir_dict[ipol],dir_dict[jpol]),
			     fontsize=20)


#	yLim = [numpy.amin,numpy.amax()]
#	print yLim

#	pylab.yticks([])

	
	labs = numpy.arange(yLim[0],yLim[1]+1,1).tolist()

	for i in range(len(labs)):
		if labs[i]%5!=0:
			labs[i]=''
		else:
			labs[i]=str(labs[i])
			

#	pylab.yticks(numpy.arange(yLim[0],yLim[1]+1,5),labs)
#	print numpy.arange(yLim[0],yLim[1]+1,5)

	'''
        takes in a list of k points that was used as pw.x input for the 'bands'
        calculation as a string. It puts parts of that string into lists and 
        manipulates them to display the symmetry point boundary lines on the 
        band structure plot
	'''
          
	bandSymSplit =  bandSym.split()
	HSPList = []
	HSPSymList = []
	buf = StringIO.StringIO(bandSym)



	for line in buf:
		splitLine = line.split()
		if len(splitLine)==2: # new style kpoint path input case
			HSPList.append(splitLine[1])
			specialPointName = splitLine[-1].rstrip()
			
 	                #renames gG to greek letter for capital gamma
			if specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gG':
				specialPointName = r"$\Gamma$"
			elif specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gS':
				specialPointName = r"$\Sigma$"
			elif specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gS1':
				specialPointName = r"$\Sigma_{1}$"

			elif len(specialPointName) != 1:
				specialPointName = "$"+specialPointName[0]+r'_{'+specialPointName[1]+'}$' #if there is a subscript it makes the print out on the plot have the number subscripted 
			else:
				specialPointName = "$"+specialPointName[0]+"$" #formats with internal math renderer so all the labels look the same
			HSPSymList.append(specialPointName)
			
		elif len(splitLine)==6: # old style kpoint path input case with kpoint names
			HSPList.append(splitLine[3])
			specialPointName = splitLine[-1].strip()

			if specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gG': #renames gG to greek letter for capital gamma
				specialPointName = r"$\Gamma$"
			elif specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gS':
				specialPointName = r"$\Sigma$"
			elif specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gS1':
				specialPointName = r"$\Sigma_{1}$"

			elif len(specialPointName) != 1:
				specialPointName = "$"+specialPointName[0]+r'_{'+specialPointName[1]+'}$' #if there is a subscript it makes the print out on the plot have the number subscripted 
			else:
				specialPointName = "$"+specialPointName[0]+"$"  #formats with internal math renderer so all the labels look the same
			HSPSymList.append(specialPointName)                     

		elif len(splitLine)==5: # old style kpoint path input case without kpoint names
			try:
				
				if HSPSymList[-1]!=splitLine[4]:
					HSPList.append(counter)
					HSPSymList.append(splitLine[4])
					counter=1

				else:
					counter+=1
			except Exception,e:
				print e
				counter=1
				HSPSymList.append(splitLine[4])


	'''
        takes the number of k points between each point in the k point paths and
        figures out if they are separate paths (where the end of one path and the 
        next begin have zero k points between them). it also takes the labels for 
        the k points that were in the bands calculation input and creates a list 
        of them and makes a special label for the path boundary e.g. X|Q. All of
        these symmetry lines and symmetry path symbols are put into their own list
        symIndex: for the symmetry line's index in the data set
        symPrint: for the symbols that represent the special symmetry path k points 
	'''

	symIndex = [0]
	totalX =0
	SymPrint = []

	for i in range(len(HSPList)-1):
             try:
		if i==0: # for the first k point in the first path
			SymPrint.append(HSPSymList[i])
		if int(HSPList[i]) == 0 and i<len(HSPList)-2: # for the end of a path (where the number of k points between one point and another is zero)
			continue
		elif int(HSPList[i+1]) == 0 and i!=len(HSPList)-2: # for the point that begins a new path where the end of the last path (which has zero k points from it to this k point)
			if not tight_banding:
				totalX +=(int(HSPList[i])+1)
			else:
				totalX +=(int(HSPList[i]))
			symIndex.append(totalX)
			mid = '|'
			pathBetweenString = HSPSymList[i+1]+mid+HSPSymList[i+2]
			SymPrint.append(pathBetweenString)
		elif int(HSPList[i+1]) != 0: # for kpoints that are not at the beginning or end of paths
			SymPrint.append(HSPSymList[i+1])
			totalX +=int(HSPList[i])
			symIndex.append(totalX)
		elif i==len(HSPList)-2: # for the end of the last path
			symIndex.append(totalX+int(HSPList[i]))
			SymPrint.append(HSPSymList[i+1])
                elif int(HSPList[i-1]) == 0 and int(HSPList[i]) == 0 and i!=len(HSPList)-2:
			logging.debug('can not find HSP in __bandPlot. This shouldnt be able to be tripped')
	     except:
		     pass
     #add symmetry lines to the band structure plot


	for sym in symIndex:
            try:

                pylab.axvline(a[sym], color = 'k',linewidth=2)
            except Exception,e:

                pylab.axvline(a[-1], color = 'k',linewidth=2)
		
                pass
     #Print path labels to band structure x-axis
	try:
		bars=[]
		for index in symIndex:
			try:
				bars.append(a[index] )
			except:
				bars.append(a[-1] )

		matplotlib.rcParams['xtick.major.pad'] = 8

		pylab.xticks(bars,SymPrint, fontsize = 24)
	except Exception,e:
                print e
		pylab.xticks([a[-1] for index in symIndex],SymPrint,fontsize=24)

	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Femi level line
	locs, labels = pylab.xticks()

	ax1.spines['bottom'].set_linewidth(1.5)
	ax1.spines['left'].set_linewidth(1.5)
	ax1.spines['right'].set_linewidth(1.5)
	ax1.spines['top'].set_linewidth(1.5)
	ax1.set_frame_on(True) #or False 

#	low_tick_bound=int(numpy.ceil(yLim[0]))
#	high_tick_bound=int(numpy.floor(yLim[1])+1.0)

#	ax1.yaxis.set_ticks(numpy.arange(low_tick_bound,high_tick_bound))

#	pylab.ylim(yLim[0],yLim[1])
	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line

     ##############################
     #to increase the linewidth of the axis, on both subplots
	if topo_type=='BC':
		description='Berry Curvature'
	if topo_type=='SBC':		
		description='Spin Berry Curvature'
#        if DOSPlot=='APDOS':
#            description+=' and Atom Projected DOS'
#        if DOSPlot=='DOS':
#            description+=' and DOS'
        
	'''gives the name of the compound in the calculation in the name of the file for the band structure plot'''
	figtitle = ''
        compoundNameLatex = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)
	figtitle = '%s: %s' % (description,compoundNameLatex) 
	ax1.set_title(figtitle,fontsize=24)
#	ax1.axes.xaxis.set_label_position('top')
#	t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=20,horizontalalignment='center') #[x,y]

	matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')

	try:
		AFLOWpi.retr._moveToSavedir(fileplot)
	except Exception,e:
		pass

        pyplot.cla()
	pyplot.clf()
	pyplot.close()

