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
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle


def _get_plot_ext_type():
	plot_ext_type = AFLOWpi.prep._ConfigSectionMap('plot','plot_file_type')
	if plot_ext_type == "":
		plot_ext_type = "pdf"	
	if plot_ext_type.lower() not in ["png","pdf"]:
		plot_ext_type = "pdf"	
	return plot_ext_type

def bands(calcs,yLim=[-10,10],DOSPlot='',runlocal=False,postfix='',tight_banding=False):
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
			__bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding)
		else:
			AFLOWpi.prep._addToBlock(oneCalc[1],oneCalc[0],'PLOT',"AFLOWpi.plot.__bands(oneCalc,ID,yLim=[%s,%s],DOSPlot='%s',postfix='%s',tight_banding=%s)" % (yLim[0],yLim[1],DOSPlot,postfix,tight_banding))	


def __bands(oneCalc,ID,yLim=[-10,10],DOSPlot='',postfix='',tight_banding=False,SBC=False,label_map={}):
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

        oneCalc = (ID,oneCalc)
	SOC=False
	try:
		nc = AFLOWpi.retr._splitInput(oneCalc[1]['_AFLOWPI_INPUT_'])['&system']['noncolin']
		if nc.lower()=='.true.':
			SOC=True
	except Exception,e: pass

	if tight_banding:
		SOC=False

	if SBC==True:
		__bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding,SBC=SBC,label_map=label_map)

		return

	if not SOC:
		__bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding,label_map=label_map)
	else:
		__bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding,spin_dir='X',label_map=label_map)
		__bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding,spin_dir='Y',label_map=label_map)
		__bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding,spin_dir='Z',label_map=label_map)




def __getPath_WanT(oneCalc,ID):
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



def _clean_bands_data_qe(filebands,Efermi_shift):
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

		with open(filebands,'r') as datafile:
			data =datafile.readlines()

		for line in data:
			if line:
				try:
					x_val = float(line.split()[0])
                                        #to shift all the y values with respect to the Fermi level
					y_val = float(line.split()[1])-Efermi_shift
					x.append(x_val)
					y.append(y_val)
				except Exception,e:
					pass
			if not line.split():

				k_x.append(x)
				k_y.append(y)
				x=[]
				y=[]


	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)
#		logging.warning("output from bands calculation not found. Are you sure you ran ppBands and it completed properly?")
#		print "Are you sure you ran ppBands and it completed properly?"
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

def __bandPlot(oneCalc,yLim=[-10,10],DOSPlot='',postfix='',tight_banding=False,spin_dir='',SBC=False,label_map={}): 
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
	

	matplotlib.rcParams['figure.dpi'] = 300
	matplotlib.rcParams['savefig.dpi'] = 300
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

        
	if DOSPlot != '' and DOSPlot != 'APDOS' and DOSPlot != 'DOS':
		print "Not a valid choice for DOSPlot. Valid options are:'APDOS','DOS'"
		return

	if postfix!='':
		postfix='_'+postfix

	if DOSPlot=='APDOS' or DOSPlot=='DOS':
		AFLOWpi.plot.__sumpdos(oneCalc,calcID,TB=tight_banding)	
	

	if tight_banding:
		bandSym = AFLOWpi.plot.__getPath_WanT(oneCalc,calcID)
	else:
		bandSym = AFLOWpi.retr._getPathFromFile(calcCopy)
	
	if bandSym==None:
		print 'ERRORRRR!'
		return
	try:
                if tight_banding:
		#	Efermi=AFLOWpi.retr._getEfermi(oneCalc,'%s_WanT_dos'%calcID,directID=True)
			Efermi=AFLOWpi.retr._getEfermi(oneCalc,calcID+"_TB",directID=True)
#			Efermi=0.0
                else:
			Efermi=AFLOWpi.retr._getEfermi(oneCalc,calcID)


		if type(Efermi)!=type(0.5):
			Efermi=Efermi[0]

	except Exception,e:            
            Efermi=0.0
	
	if spin_dir!='':
		postfix+='_'+spin_dir
	subdir=oneCalc['_AFLOWPI_FOLDER_']

	petype = AFLOWpi.plot._get_plot_ext_type()
       	fileplot = os.path.join(subdir,'BANDS_%s_%s%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix,petype))

	"""get the path to the subdirectory of the calc that you are making plots for"""

	if tight_banding==True:
	#	try:
#			AFLOWpi.prep._clean_want_bands(oneCalc,calcID)
#		except:
#			return
		filebands = os.path.join(subdir,'%s_bands_paopy_cleaned.dat'%calcID)
		if not os.path.exists(filebands):

			filebands_up = os.path.join(subdir,'%s_bands_paopy_up_cleaned.dat'%calcID)
			filebands_dn = os.path.join(subdir,'%s_bands_paopy_down_cleaned.dat'%calcID)
                Efermi_shift=Efermi

	else:
		filebands = os.path.join(subdir,'%s_bands.xmgr'%calcID)
		if not os.path.exists(filebands):
			filebands_up = os.path.join(subdir,'%s_up_bands.xmgr'%calcID)
			filebands_dn = os.path.join(subdir,'%s_dn_bands.xmgr'%calcID)

                Efermi_shift=Efermi



        dos_ID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='bands',last=True)
	filedos = os.path.join(subdir,'%s_dos.dat'%dos_ID)

	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
	

       #to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=20)             #to set the font size
	
#	pylab.style.use('dark_background')
        width = 20
	height = 14
	fig = pylab.figure(figsize=(width, height))#to adjust the figure size
	
     #to do the gaussian smoothing               
 
     #################################
	if nspin==2:
		k_x_up,k_y_up = AFLOWpi.plot._clean_bands_data_qe(filebands_up,Efermi_shift)
		k_x_dn,k_y_dn = AFLOWpi.plot._clean_bands_data_qe(filebands_dn,Efermi_shift)
		k_x,k_y=k_x_up,k_y_up
	else:
		k_x,k_y = AFLOWpi.plot._clean_bands_data_qe(filebands,Efermi_shift)




	a=k_x[1]   # a set of k point values for one band for axis scaling purposes
       	b=k_y[1]




	if DOSPlot == '':
		ax1=pylab.subplot(111)	
		print 'Plotting electronic band structure of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
	else:
		ax1=pylab.subplot(121)	

	'''
        Plot each band (k_x[i]),(k_y[i]) on the band structure plot from list of
        the values of energy and position from their respective list by k points
	'''
	if nspin==2:
		for i in range(len(k_x_up)):
			try:
				
				if tight_banding==True:
					pylab.plot((k_x_up[i]),(k_y_up[i]),'r',alpha=1.0,marker=".",
						   linestyle=" ",label="$\uparrow$",linewidth=2)
				else:
					pylab.plot((k_x_up[i]),(k_y_up[i]),'r',alpha=1.0,marker=".",
						   linestyle=" ",label="$\uparrow$",linewidth=2)
			except:
				pass
		for i in range(len(k_x_up)):
			try:
				
				if tight_banding==True:
					pylab.plot((k_x_dn[i]),(k_y_dn[i]),'k',alpha=1.0,marker=".",
						   linestyle=" ",label="$\downarrow$",linewidth=2)
				else:
					pylab.plot((k_x_dn[i]),(k_y_dn[i]),'k',alpha=1.0,marker=".",
						   linestyle=" ",label="$\downarrow$",linewidth=2)
			except:
				pass
		handles, labels = ax1.get_legend_handles_labels()
#		ax1.legend(handles[-2:], labels[-2:],numpoints=1)
		up_legend_lab = mlines.Line2D([], [], color='red',label='$\uparrow$')
		dn_legend_lab = mlines.Line2D([], [], color='black',label='$\downarrow$')
		ax1.legend(handles=[up_legend_lab,dn_legend_lab],loc=1)
	else:
		SOC=False
		try:
			filebands_xspin = os.path.join(subdir,'%s_bands.xspin'%calcID)
			zspin_kpath,xspin = AFLOWpi.plot._clean_bands_data_qe(filebands_xspin,0.0)
			filebands_yspin = os.path.join(subdir,'%s_bands.yspin'%calcID)
			zspin_kpath,yspin = AFLOWpi.plot._clean_bands_data_qe(filebands_yspin,0.0)
			filebands_zspin = os.path.join(subdir,'%s_bands.zspin'%calcID)
			zspin_kpath,zspin = AFLOWpi.plot._clean_bands_data_qe(filebands_zspin,0.0)

			xspin=numpy.ravel(numpy.array(xspin),order="C")
			yspin=numpy.ravel(numpy.array(yspin),order="C")
			zspin=numpy.ravel(numpy.array(zspin),order="C")
		
			if spin_dir=='':
				rgb_arr=numpy.array([zspin,yspin,xspin]).T
			if spin_dir=='X':
				rgb_arr=xspin
			if spin_dir=='Y':
				rgb_arr=yspin
			if spin_dir=='Z':
				rgb_arr=zspin


			SOC=True
		except Exception,e:  pass

		if tight_banding==True:
#			if negative_SBC!=True:
#				cmap = pylab.cm.copper
#				my_Oranges = cmap(numpy.arange(cmap.N))
#				my_Oranges[:,-1] = numpy.linspace(0, 1, cmap.N)
#				my_Oranges = ListedColormap(my_Oranges)
#				my_Oranges = pylab.cm.gnuplot
#				my_Oranges = pylab.cm.autumn_r
#			else:
#				cmap = pylab.cm.copper
#				my_Oranges = cmap(numpy.arange(cmap.N))
#				my_Oranges[:,-1] = numpy.linspace(0, 1, cmap.N)
#				my_Oranges = ListedColormap(my_Oranges[::-1])
#				my_Oranges = pylab.cm.gnuplot_r
#				my_Oranges = 
#

			Ojn_path = os.path.join(subdir,'Omegajn_z_xy.dat')
			BC=True
			if BC:
				Ojn_path = os.path.join(subdir,'Omegajn_z_xy.dat')


			if SBC==True and os.path.exists(Ojn_path):
				cmap_pos = pylab.cm.get_cmap('YlOrRd')
				cmap_neg = pylab.cm.get_cmap('YlGnBu_r')

				import time
				st = time.time()
				#invert those colormaps
				cmap_neg = 1.0-cmap_neg(numpy.arange(0,255,1))[:,:-1]
				cmap_neg = matplotlib.colors.LinearSegmentedColormap.from_list('YlGnBu_r_inv',
												cmap_neg , N=256)
				cmap_pos = 1.0-cmap_pos(numpy.arange(0,255,1))[:,:-1]
				cmap_pos = matplotlib.colors.LinearSegmentedColormap.from_list('YlOrRd_inv',  
												cmap_pos , N=256)

				print time.time()-st





				rgba_test = cmap_pos(0.0)
				#plot the black band structure then add color for SBC
#				for i in range(len(k_x)):
#					pylab.scatter((k_x[i]),(k_y[i]),c=rgba_test,s=3,edgecolors='none')

				#remove bad values
				Omega_znk=numpy.loadtxt(Ojn_path)
				inf_ind = numpy.where(numpy.isinf(Omega_znk))
				Omega_znk[inf_ind] = 1.e20
				nan_ind = numpy.where(numpy.isnan(Omega_znk))
				Omega_znk[nan_ind] = 1.e-20

				nan_ind = numpy.where(Omega_znk<-1.e6)
				Omega_znk[nan_ind] = 1.e-20
				nan_ind = numpy.where(Omega_znk>1.e6)
				Omega_znk[nan_ind] = 1.e-20


				#limit plot to only 'good' bands
				k_x = numpy.array(k_x)
				k_y = numpy.array(k_y)

				bnd = Omega_znk.shape[0]
				nk  = Omega_znk.shape[1]

				k_x = numpy.ascontiguousarray(k_x[:bnd])
				k_y = numpy.ascontiguousarray(k_y[:bnd])

				#copy arrays
				zspinr = Omega_znk
				rgb_arr=zspinr


                                # ########################################################
				# rgb_arr_r = numpy.ravel(rgb_arr)
				# k_x_r = numpy.ravel(k_x)
				# k_y_r = numpy.ravel(k_y)
				# order = numpy.argsort(numpy.abs(rgb_arr_r))
				# k_x_r      = k_x_r[order]  
				# k_y_r      = k_y_r[order]    
				# rgb_arr_r  = rgb_arr_r[order]
                                # ########################################################

                                ########################################################
                                # ########################################################
				# for i in xrange(rgb_arr.shape[0]):
				# 	order = numpy.argsort(numpy.abs(rgb_arr[i]))
				# 	k_x[i]      = k_x[i][order]  
				# 	k_y[i]      = k_y[i][order]    
				# 	rgb_arr[i]  = rgb_arr[i][order]
                                # ########################################################
				# for i in xrange(rgb_arr.shape[1]):
  				# 	order = numpy.argsort(rgb_arr[:,i])
				# 	k_x[:,i]      = k_x[:,i][order]  
				# 	k_y[:,i]      = k_y[:,i][order]    
				# 	rgb_arr[:,i]  = rgb_arr[:,i][order]
                                # ########################################################


				#mask positive SBC values 
				k_ymn = numpy.ma.masked_where(zspinr>=0.0,k_y)
				k_xmn = numpy.ma.masked_where(zspinr>=0.0,k_x)
				rgb_arrmn = numpy.ma.masked_where(zspinr>0.0,rgb_arr)

				rgb_arrmn[numpy.where(rgb_arrmn>-1.0)]=-1.0

				print numpy.amin(rgb_arrmn)
				print numpy.amax(rgb_arrmn)

				#sort from least negative plotted first to most negative plotted last
                                ########################################################
				# for i in xrange(rgb_arrmn.shape[1]):
  				# 	neg_order = numpy.argsort(rgb_arrmn[:,i])
				# 	k_xmn[:,i]      = k_xmn[:,i][neg_order]  
				# 	k_ymn[:,i]      = k_ymn[:,i][neg_order]    
				# 	rgb_arrmn[:,i]  = rgb_arrmn[:,i][neg_order]
                                # #######################################################
				# for i in xrange(rgb_arrmn.shape[0]):
				# 	neg_order = numpy.argsort(numpy.abs(rgb_arrmn[i]))
				# 	k_xmn[i]      = k_xmn[i][neg_order]  
				# 	k_ymn[i]      = k_ymn[i][neg_order]    
				# 	rgb_arrmn[i]  = rgb_arrmn[i][neg_order]
                                # # ########################################################



				#do log base 10 for SBC values
				rgb_arrmn= -1.0*numpy.log10(numpy.abs(rgb_arrmn))

				#normalize the colormap for negative values
				normalize_neg = matplotlib.colors.Normalize(vmin=numpy.ceil(numpy.amin(rgb_arrmn)), 
									    vmax=0.0)
			        ###################################################################################
				#mask negative SBC values 
				k_ymp = numpy.ma.masked_where(zspinr<0.0,k_y)
				k_xmp = numpy.ma.masked_where(zspinr<0.0,k_x)
				rgb_arrmp = numpy.ma.masked_where(zspinr<0.0,rgb_arr)

				rgb_arrmp[numpy.where(rgb_arrmp<1.0)]=1.0


				print numpy.amin(rgb_arrmp)
				print numpy.amax(rgb_arrmp)


				# #sort from least positive plotted first to most positive plotted last
                                # ########################################################
				# for i in xrange(rgb_arrmp.shape[1]):
				# 	pos_order = numpy.argsort(rgb_arrmp[:,i])
				# 	k_xmp[:,i]      = k_xmp[:,i][pos_order]  
				# 	k_ymp[:,i]      = k_ymp[:,i][pos_order]    
				# 	rgb_arrmp[:,i]  = rgb_arrmp[:,i][pos_order]
                                # # ########################################################
				# for i in xrange(rgb_arrmp.shape[0]):
				# 	pos_order = numpy.argsort(rgb_arrmp[i])
				# 	k_xmp[i]      = k_xmp[i][pos_order]  
				# 	k_ymp[i]      = k_ymp[i][pos_order]    
				# 	rgb_arrmp[i]  = rgb_arrmp[i][pos_order]
                                # ########################################################


				#do log base 10 for SBC values
				rgb_arrmp=numpy.log10(rgb_arrmp)


				#normalize the colormap for postiive values
				normalize_pos = matplotlib.colors.Normalize(vmin=0.0, 
									vmax=numpy.floor(numpy.amax(rgb_arrmp)))





				rgb_arr = numpy.zeros_like(rgb_arr)

				mask_neg = numpy.ma.getmask(k_xmn)
				mask_pos = numpy.ma.getmask(k_xmp)
#				k_x[mask_neg] = k_xmn[mask_neg]
#				k_y[mask_neg] = k_ymn[mask_neg]
				rgb_arr[mask_neg] = rgb_arrmn[mask_neg]
#				k_x[mask_pos] = k_xmp[mask_pos]
#				k_y[mask_pos] = k_ymp[mask_pos]
				rgb_arr[mask_pos] = rgb_arrmp[mask_pos]

				rgb_arr = numpy.ma.array(rgb_arr)

				mask = numpy.ma.MaskedArray.argsort(rgb_arr,axis=1,kind='mergesort')
				


				for i in xrange(k_xmn.shape[0]):
					k_xmn[i] = k_xmn[i,mask[i]]
					k_ymn[i] = k_ymn[i,mask[i]]
					rgb_arrmn[i] = rgb_arrmn[i,mask[i]]
					k_xmp[i] = k_xmp[i,mask[i]]
					k_ymp[i] = k_ymp[i,mask[i]]
					rgb_arrmp[i] = rgb_arrmp[i,mask[i]]

#				print numpy.prod(numpy.shape(rgb_arr))
#				print numpy.sum(mask_neg)+numpy.sum(mask_pos)


#				normalize_a = matplotlib.colors.Normalize(vmin=numpy.floor(numpy.amin(rgb_arr)), 
#									vmax=numpy.ceil(numpy.amax(rgb_arr)))
                                # #######################################################
				# for i in xrange(rgb_arrmn.shape[0]):
				# 	order = numpy.argsort(rgb_arr[i])
				# 	k_x[i]      = k_x[i][order]  
				# 	k_y[i]      = k_y[i][order]    
				# 	rgb_arr[i]  = rgb_arr[i][order]
                                # ########################################################


				# #plot in the order of eigenvalues for each k point on path
				# for i in xrange(k_x.shape[0]):
				# 	pylab.scatter(k_x[i],k_y[i],c=rgb_arr[i],s=3,cmap=pylab.cm.gist_rainbow,
				# 		      edgecolors='none',norm=normalize_a)



				# k_xmp     = numpy.ascontiguousarray(k_xmp.T)
				# k_ymp     = numpy.ascontiguousarray(k_ymp.T)
				# rgb_arrmp = numpy.ascontiguousarray(rgb_arrmp.T)

				# k_xmn     = numpy.ascontiguousarray(k_xmn.T)
				# k_ymn     = numpy.ascontiguousarray(k_ymn.T)
				# rgb_arrmn = numpy.ascontiguousarray(rgb_arrmn.T)


				order = numpy.argsort(k_xmp,axis=1)[0]
				print order.shape
				print k_xmp.shape
				print k_ymp.shape
				print rgb_arrmp.shape
				k_xmn      = k_xmn[:,order]  
				k_ymn      = k_ymn[:,order]    
				rgb_arrmn  = rgb_arrmn[:,order]

				k_xmp      = k_xmp[:,order]  
				k_ymp      = k_ymp[:,order]    
				rgb_arrmp  = rgb_arrmp[:,order]

				

				#plot in the order of eigenvalues for each k point on pat
				for i in xrange(k_xmp.shape[0]):
					pylab.scatter(k_xmp[i],k_ymp[i],c=rgb_arrmp[i],s=20,cmap=cmap_pos,
						      edgecolors='none',norm=normalize_pos,alpha=1.0,zorder=i+1,marker='_')
				for i in xrange(k_xmp.shape[0]):
					pylab.scatter(k_xmn[i],k_ymn[i],c=rgb_arrmn[i],s=20,cmap=cmap_neg,
					 	      edgecolors='none',norm=normalize_neg,alpha=1.0,zorder=i,marker='_') 
					


				fileplot = os.path.join(subdir,'SPIN_BERRY_CURV_%s_%s%s.png' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix))



######################################################################################################
			
			else:
				for i in range(len(k_x)):
					pylab.plot((k_x[i]),(k_y[i]),'k',marker=".",
						   linestyle=" ",linewidth=1.3)	
		else:
			if SOC:
				normalize = matplotlib.colors.Normalize(vmin=0.0, vmax=1.0)
				k_xr=numpy.ravel(numpy.array(k_x),order="C")
				k_yr=numpy.ravel(numpy.array(k_y),order="C")
				zspinr=numpy.ravel(numpy.array(zspin),order="C")

				order = numpy.argsort(numpy.abs(rgb_arr))
				rgb_arr+=0.5
				k_xr=k_xr[order]#[::-1]
				k_yr=k_yr[order]#[::-1]
				rgb_arr=rgb_arr[order]#[::-1]

				pylab.scatter(k_xr,k_yr,c=rgb_arr,s=20,cmap='coolwarm',
					      edgecolors='none',norm=normalize) 
			else:
				for i in range(len(k_x)):
					pylab.plot((k_x[i]),(k_y[i]),'k',marker=".",linestyle=" ",linewidth=1.3)			
	#


	pylab.ylabel('E(eV)')
	pylab.xlim(min(k_x[1]),max(k_x[1])) 
	pylab.ylim(yLim[0],yLim[1])    
#	pylab.yticks(numpy.arange(yLim[0],yLim[1]+1,2))

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

                pylab.axvline(a[sym], color = 'k',linewidth=2,zorder=10000)
            except Exception,e:
                pylab.axvline(a[-1], color = 'k',linewidth=2,zorder=10000)
		
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
              
		pylab.xticks([a[-1] for index in symIndex],SymPrint)

	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3,zorder=10000) #fermi
	locs, labels = pylab.xticks()

##########################################################################################################
#	if DOSPlot == '':
#		ax1.set_title(figtitle,fontsize=24)



	if DOSPlot == 'APDOS':
		print 'Plotting electronic band structure and projected DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure and projected DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
		petype = AFLOWpi.plot._get_plot_ext_type()
		fileplot = os.path.join(subdir,'BANDPDOS_%s_%s%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix,petype))	
		ax3=pylab.subplot(122)
		ax3.xaxis.set_ticks([])
		ax3.yaxis.set_ticks([])
		pylab.xlabel('arbitrary units',fontsize = 20)
		ax3.set_position([0.75,0.1,0.20,0.8]) 
		ax3.axes.xaxis.set_label_position('bottom')
		ax2=pylab.subplot(122)

		def getPlotData(sumpdosFile):
			try:
				with open(sumpdosFile,'rb') as dataFile:
					data = cPickle.load(dataFile)

			except Exception,e:
			
				with open(sumpdosFile,'r') as dataFile:
					data = dataFile.read()
					data = data.split('\n')
					for entry in range(len(data)):
						data[entry]=data[entry].split()


			en = []
			pdos = []
			ldos = []
			ldosDOWN = []
			for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
				try:
					#to shift all the y values with respect to the Fermi level
					en.append(float(data[i][0])-Efermi) 
					ldos.append(float(data[i][1]))
					pdos.append(float(data[i][2]))
					ldosDOWN.append(-1*float(data[i][2]))
					
				except Exception, e:
					pass

			endos=map(float,en)  #to convert the list x from float to numbers
			floatdos=map(float,ldos)
			floatdosDOWN=map(float,ldosDOWN)
			enshift = numpy.array(endos) #to treat the list b as an array?

			return enshift,floatdos,floatdosDOWN

		maxDOS=0
		minDOS=0
		atomName= []
		pDOSNoPath = []
		atomList=[]

		if os.path.exists(os.path.join(subdir,'%s_dos.dat'%dos_ID)):
			pDOSNoPath.append(os.path.join(subdir,'%s_dos.dat'%dos_ID))


		alab = AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])

		
		atomList = []
		for i in alab:
			if i not in atomList:
				atomList.append(i)



		for species in atomList:
			filePath = os.path.join(subdir,'%s_All.sumpdos' % species)
			if os.path.exists(filePath):
				pDOSNoPath.append(filePath)
		color = 'k'

		colors = ['b','g','r','m','c', 'y', 'k',"orange"]
		for fi in range(len(pDOSNoPath)):
			filename =  pDOSNoPath[fi].split('/')[-1]

			try:
				species = re.findall(r'(.+?)_.+',filename)[0]
			except:
				if filename=='%s_dos.dat'%dos_ID:
					species='TOTAL'
			
			'''gets the energy and the DOS for the orbital for the atom'''
			enshift, floatdos,floatdosDOWN = getPlotData(pDOSNoPath[fi])
			"""makes a sum of all the pdos for a total"""
#			floatdos=AFLOWpi.plot.__smoothGauss(floatdos)
#			floatdosDOWN=AFLOWpi.plot.__smoothGauss(floatdosDOWN)
#			enshift=AFLOWpi.plot.__smoothGauss(enshift)


			''' scales DOS to larges value of DOS in the given energy range and finds the largest DOS between the different orbitals'''
			try:
				if max([floatdos[k] for k in range(len(floatdos)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ]) > maxDOS:
					maxDOS = max([floatdos[k] for k in range(len(floatdos)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ])
			
				if min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ]) < minDOS:
					minDOS = min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ])
			except:
                            pass


			if species in label_map.keys():
				species_lab = label_map[species]
			else:
				species_lab = species

			if species=='TOTAL':	
				pylab.plot(floatdos,enshift,'-',label=species_lab,color='k',linewidth=2)       
			else:
				pylab.plot(floatdos,enshift,'-',label=species_lab,linewidth=2,color=colors[fi%len(colors)])			   
#				pyplot.fill(floatdos,enshift,color=colors[fi%len(colors)])			   



			if LSDA:
				if species=='TOTAL':
					pylab.plot(floatdosDOWN,enshift,'-',label=species_lab,color='k',linewidth=2) 
				else:
					pylab.plot(floatdosDOWN,enshift,'-',label=species_lab,linewidth=2,color=colors[fi%len(colors)])		      

		handles, labels = ax2.get_legend_handles_labels()

		if LSDA:
			ax2.legend(handles[::2], labels[::2],fontsize=14,loc=1)
			dosRange=max([nump5Ay.abs(minDOS),numpy.abs(maxDOS)])
			pylab.xlim(1.1*minDOS,1.1*maxDOS) # scales DOS to larges value of DOS in the given energy range 
			pylab.axvline(0.0, color = 'k', linewidth = 2.0) #line separating up and down spin
		else:
#			ax2.legend(handles[::-1], labels[::-1],fontsize=14,loc=1)
			ax2.legend(handles, labels,fontsize=14,loc=1)
			pylab.xlim(0,1.1*maxDOS) # scales DOS to larges value of DOS in the given energy range

#		pylab.yticks(numpy.arange(yLim[0],yLim[1]+1,2))

		ax2.spines['bottom'].set_linewidth(1.5)
		ax2.spines['left'].set_linewidth(1.5)
		ax2.spines['right'].set_linewidth(1.5)
		ax2.spines['top'].set_linewidth(1.5)

		ax2.set_yticklabels([])	     #to hide ticklabels
		ax1.set_position([0.07,0.1,0.67,0.8]) #[left,bottom,width,height] 
		ax2.set_position([0.75,0.1,0.20,0.8]) #other useful options for the frame! :D


		ax2.yaxis.set_ticks_position('left')
		pylab.xlabel('Density of States (States/eV)',fontsize=20)
		ax2.axes.xaxis.set_label_position('top')
		locs, labels = pylab.xticks()
		
		# for item in range(len(labels)):
		# 	if item == len(labels)/2:
		# 		labels[item]=
		# 	else:
		# 		labels[item]=''
	

#		ax3.set_xticklabels(labels,fontsize = 20)

##########################################################################################################
	     #to plot the DOS
	if DOSPlot == 'DOS':
		print 'Plotting electronic band structure and DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure and DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
		petype = AFLOWpi.plot._get_plot_ext_type()
		fileplot = os.path.join(subdir,'BANDDOS_%s%s%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),oneCalc['_AFLOWPI_PREFIX_'],postfix,petype))
		ax2=pylab.subplot(122)

		try:
			data = open(filedos,'r').readlines()
		except Exception:
#			logging.warning("output from dos calculation not found. Are you sure you ran ppDOS and it completed properly?")
#			print "Are you sure you ran ppDOS and it completed properly?"
			return

		en = []
		dos = []
		dosdw = []

		for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
			try:
				#to shift all the y values with respect to the Fermi level 
				en.append(float(data[i].split()[0])-Efermi) 
				dos.append(float(data[i].split()[1]))
				dosdw.append(-1*float(data[i].split()[2]))

			except Exception, e:
				pass

		endos=map(float,en)  #to convert the list x from float to numbers
		floatdos=map(float,dos)
		floatdosDOWN=map(float,dosdw)
		enshift = numpy.array(endos) #to treat the list b as an array?

#		floatdos=AFLOWpi.plot.__smoothGauss(floatdos)
#		floatdosDOWN=AFLOWpi.plot.__smoothGauss(floatdosDOWN)
#		enshift=AFLOWpi.plot.__smoothGauss(enshift)		

		ax2=pylab.subplot(122)

		dosMAX = max([dos[k] for k in range(len(dos)) if en[k] < yLim[1] and en[k] > yLim[0] ])
		pylab.plot(floatdos,enshift,'k') #to plot the smoothed data
		
		if LSDA:
			dosMIN = min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if en[k] < yLim[1] and en[k] > yLim[0] ])
			pylab.plot(floatdosDOWN,enshift,'k') #to plot the smoothed data
			pylab.xlim(1.1*dosMIN,1.1*dosMAX) # scales DOS to larges value of DOS in the given energy range
			pylab.axvline(0.0, color = 'k', linewidth = 1.3) #line separating up and down spin
		else:
			pylab.xlim(0,1.1*dosMAX) # scales DOS to larges value of DOS in the given energy range

		ax2.spines['bottom'].set_linewidth(1.5)
		ax2.spines['left'].set_linewidth(1.5)
		ax2.spines['right'].set_linewidth(1.5)
		ax2.spines['top'].set_linewidth(1.5)

		locs, labels = pylab.xticks()
		
		for item in range(len(labels)):
			if item == len(labels)/2:
				labels[item]='arbitrary units'
			else:
				labels[item]=''
		
		ax2.set_xticklabels(labels,fontsize = 14)
		ax1.set_position([0.07,0.1,0.67,0.8]) #[left,bottom,width,height]
		ax2.set_position([0.75,0.1,0.20,0.8])

	     #other useful options for the frame! :D
		ax2.yaxis.set_ticks([])
		ax2.yaxis.set_ticks_position('left')
		pylab.xlabel('Density of States (States/eV)')
		ax2.axes.xaxis.set_label_position('top')
###########################################################################################

	ax1.spines['bottom'].set_linewidth(1.5)
	ax1.spines['left'].set_linewidth(1.5)
	ax1.spines['right'].set_linewidth(1.5)
	ax1.spines['top'].set_linewidth(1.5)
	ax1.set_frame_on(True) #or False 

	low_tick_bound=int(numpy.ceil(yLim[0]))
	high_tick_bound=int(numpy.floor(yLim[1])+1.0)

#	ax1.yaxis.set_ticks(numpy.arange(low_tick_bound,high_tick_bound))

	pylab.ylim(yLim[0],yLim[1])
	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line

     ##############################
     #to increase the linewidth of the axis, on both subplots
        description='Electronic Band Structure'
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


	#this has to go at the end for now....




		   

	if SBC:
		compoundNameLatex = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)
		description='Electronic Band Structure and Spin Berry Curvature'
		if BC: description='Electronic Band Structure and Berry Curvature'
			
		figtitle = '%s: %s' % (description,compoundNameLatex) 
		ax1.set_title(figtitle,fontsize=24)


		neg_cax = fig.add_axes([0.98,0.2875,0.02,0.25])
		pos_cax = fig.add_axes([0.98,0.5375,0.02,0.25])


		
		cbarneg = matplotlib.colorbar.ColorbarBase(neg_cax, cmap=cmap_neg,
							norm=normalize_neg,
							extend="neither",
							spacing='uniform',
							orientation='vertical')








 		labs = cbarneg.ax.get_yticklabels()
		labs = [item.get_text() for item in labs]


		num_ticks= numpy.abs(normalize_neg.vmin)+1
		tm = numpy.linspace(0, 1,num_ticks+1,endpoint=True)
		cbarneg.ax.yaxis.set_ticks(tm)

 		labs = cbarneg.ax.get_yticklabels()
		labs = [item.get_text() for item in labs]
		
		for item in xrange(len(labs)):
			if item==len(labs)-1:
				labs[item] =  u'    $0$'			
			else:
				labs[item] =  u'$-10^{%s}$'%(len(labs)-int(item)-1)

		cbarneg.ax.set_yticklabels(labs)

		cbarneg.ax.tick_params('both',width=3,length=7)
		cbarneg.ax.tick_params(direction='in')



# 		raise SystemExit


# 			elif t_item==0.0:
# 				labs[item]=u'     $0$'
				
# 			else:
# 				labs[item]=''




#		cbarneg.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))	


		cbarneg.ax.tick_params(direction='in')
		cbarneg.ax.set_yticklabels(labs)
		cbarneg.ax.tick_params('both',width=2,length=30)

		cbar = matplotlib.colorbar.ColorbarBase(pos_cax, cmap=cmap_pos,
							norm=normalize_pos,
							extend='neither',

							spacing='uniform',
							orientation='vertical')





		cbar.ax.tick_params(labelsize=18)
		cbarneg.ax.tick_params(labelsize=18)

		num_ticks= numpy.abs(normalize_pos.vmax)+1
		tm = numpy.linspace(0, 1,num_ticks+1,endpoint=True)
		cbar.ax.yaxis.set_ticks(tm)

 		labs = cbar.ax.get_yticklabels()
		labs = [item.get_text() for item in labs]
		
		for item in xrange(len(labs)):
			if item==0:
				labs[item] =  u''			
			else:
				labs[item] =  u'   $10^{%s}$'%(int(item))

		cbar.ax.set_yticklabels(labs)


		cbar.ax.tick_params('both',width=2,length=29)
		cbarneg.ax.tick_params('both',width=2,length=29)
		cbar.ax.tick_params(direction='in')

		if BC:
			cbar.ax.set_title('  $\Omega_{xy}$\n',size=44 )
		else:
			cbar.ax.set_title('  $\Omega^{z}_{xy}$\n',size=44 )


		ax1.set_position([0.07,0.1,0.90,0.95]) #[left,bottom,width,height] 		


	#save fig
	matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')

	try:
		AFLOWpi.retr._moveToSavedir(fileplot)
	except Exception,e:
		pass

        pyplot.cla()
	pyplot.clf()
	pyplot.close()

