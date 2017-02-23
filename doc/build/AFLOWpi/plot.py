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

    if xaxisTitle==None:
	    xaxisTitle=variable1
    if yaxisTitle==None:
	    yaxisTitle=variable2
    X=[]
    Y=[]
    Z=[]
    try:
        if zaxis=='Energy':
            calcs=AFLOWpi.retr.grabEnergyOut(calcs)
    except:
	    pass

    for key,oneCalc, in calcs.iteritems():
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


        except Exception,e:
                AFLOWpi.run._fancy_error_log(e)    
    

    im = pylab.imshow(Z,cmap=plot_color,aspect='auto',vmin=z_min,vmax=z_max,interpolation=interp,extent=[X.min(), X.max(), Y.min(), Y.max()],origin='lower',hold=True)
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
    except Exception,e:
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

    if xaxisTitle==None:
	    xaxisTitle=variable1
    if yaxisTitle==None:
	    yaxisTitle=yaxis
    X=[]
    Y=[]

    try:
	    calcs=AFLOWpi.retr.grabEnergyOut(calcs)
    except:
	    pass

    for key,oneCalc, in calcs.iteritems():
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

    except Exception,e:
	    AFLOWpi.run._fancy_error_log(e)    
    
    
    ax.set_xlabel(xaxisTitle)
    ax.set_ylabel(yaxisTitle)

    if title!=None:
	    pyplot.title(title)

	    
    pyplot.savefig(fileName,bbox_inches='tight')
    try:
	    AFLOWpi.retr._moveToSavedir(fileName)
    except Exception,e:
	    AFLOWpi.run._fancy_error_log(e)


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








###############################################################################################################################
def transport_plots(calcs,runlocal=False,postfix='',x_range=None):


	if runlocal:
		for ID,oneCalc in calcs.iteritems():
                    AFLOWpi.plot.__transport_plot(oneCalc,ID,postfix=postfix,x_range=x_range)
	else:
		for ID,oneCalc in calcs.iteritems():
			AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__transport_plot(oneCalc,ID,postfix='%s',x_range=%s)" %(postfix,x_range))


def optical_plots(calcs,runlocal=False,postfix='',x_range=None):


	if runlocal:
		for ID,oneCalc in calcs.iteritems():
                    AFLOWpi.plot.__transport_plot(oneCalc,ID,postfix=postfix,epsilon=True,x_range=x_range)
	else:
		for ID,oneCalc in calcs.iteritems():
			AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__transport_plot(oneCalc,ID,postfix='%s',epsilon=True,x_range=%s)" %(postfix,x_range))




def __transport_plot(oneCalc,ID,nm=False,postfix='',epsilon=False,x_range=None):
	'''


	Arguments:


	Keyword Arguments:


	'''

        trans_plot_dict={}
        trans_plot_dict['ZT']            = {'pf':'ZetaT_',
                                            'ft':'ZT:',
                                            'lc':'ZT ',
                                            'xl':'$\mu$ (eV)',
                                            'yl':'ZT (au)',
                                            'fp':'ZETAT',

                                            }
        trans_plot_dict['cond']          = {'pf':'cond_',
                                            'ft':'$Conduction$:',
                                            'lc':r'{\sigma/\tau}',
                                            'xl':'$\mu$ (eV)',
                                            'yl':r'$\sigma/\tau$ $(10^{20}$ $m/s/\Omega)$ ',
                                            'fp':'CONDUCTION',

                                            }
        trans_plot_dict['seebeck']       = {'pf':'seebeck_',
                                            'ft':'$Seebeck$:',
                                            'lc':'Seebeck ',
                                            'xl':'$\mu$ (eV)',
                                            'yl':'S $(10^{-3} $ $V/K)$',
                                            'fp':'SEEBECK',
                                            }

        trans_plot_dict['sig_seebeck']   = {'pf':'sigma_seebeck_',
                                            'ft':'$\sigma S$:',
                                            'lc':'\sigma S ',
                                            'xl':'$\mu$ (eV)',
                                            'yl':'$\sigma S$ ($Vm/\Omega/K)$ ',
                                            'fp':'SIGMA_SEEBECK',

                                            }
        trans_plot_dict['kappa']         = {'pf':'kappa_',
                                            'ft':'$\kappa$:',
                                            'lc':'\kappa ',
                                            'xl':'$\mu$ (eV)',
                                            'yl':'$\kappa$ $(10^{17}$ $W/m/K)$',
                                            'fp':'KAPPA',

                                            }
        # trans_plot_dict['epsilon_i']     = {'pf':'epsilon_imag_',
        #                                     'ft':'$\epsilon_{imaginary}$:',
        #                                     'lc':'\epsilon_{i} ',
        #                                     'xl':'$\hbar\omega$ (eV)',
        #                                     'yl':'$\epsilon_{imag}$ (au)',
        #                                     'fp':'EPSILON_IMAG',

	#}
        trans_plot_dict['epsilon']       = {'pf':'epsilon_',
                                            'ft':'$\epsilon$:',
                                            'lc':'\epsilon ',
                                            'xl':'$\hbar\omega$ (eV)',
                                            'yl':'$\epsilon$ (au)',
                                            'fp':'EPSILON',

                                            }
        


####################################################################################################################
        ##needed to make sure the real and imaginary parts of epsilon have the same scaling
        min_val=0.0
        max_val=0.0


        ID_list =  AFLOWpi.prep._return_ID(oneCalc,ID,step_type='PAO-TB',last=True,straight=True)

####################################################################################################################
	if epsilon==True:
		type_list=['epsilon',]
	else:
		type_list=['kappa','seebeck','sig_seebeck','cond','ZT',]

        for Type in type_list:
            try:
                max_val=0.0
                min_val=0.0


                file_name=[]
                file_name_up=[]
                file_name_down=[]
                for ID_ent in ID_list:

                    search=oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_%s*.dat'%(ID_ent,trans_plot_dict[Type]['pf']) 

                    file_name.extend(glob.glob(search))

                    file_name_up.extend(glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_%s*up*.dat'%(ID_ent,trans_plot_dict[Type]['pf']) ))
                    file_name_down.extend(glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_%s*down*.dat'%(ID_ent,trans_plot_dict[Type]['pf']) ))



                spin_polarized=False
                if len(file_name_down)!=0:
                        spin_polarized=True



                if spin_polarized==False:
                        ############################################################################################
                        ##Non Spin Polarized Case
                        ############################################################################################
                        data_by_temp=collections.OrderedDict()
                        for i in range(len(file_name)):
                            temperature = file_name[i][:-4].split('_')[-1][:-1]
                            if epsilon==True:

				    if temperature=="ima":
					    temperature="-100"
				    elif temperature=="rea":
					    temperature="-200"
				    else:
					    continue
			    else:
				    temperature = file_name[i][:-4].split('_')[-1][:-1]

                            data_by_temp[temperature] = read_transport_datafile(file_name[i])               

                if spin_polarized==True:
                        #############################################################################################
                        ##Spin Polarized Case
                        #############################################################################################
                        data_by_temp_dn=collections.OrderedDict()
                        data_by_temp_up=collections.OrderedDict()

                        for i in range(len(file_name_down)):
		            temperature = file_name_down[i][:-4].split('_')[-1][:-1]
                            if epsilon==True:

				    if temperature=="ima":
					    temperature="-100"
				    elif temperature=="rea":
					    temperature="-200"
				    else:

					    continue

			    else:
				    temperature = file_name_down[i][:-4].split('_')[-1][:-1]

                            data_by_temp_dn[temperature] = read_transport_datafile(file_name_down[i]) 

                        for i in range(len(file_name_up)):
		            temperature = file_name_up[i][:-4].split('_')[-1][:-1]
                            if epsilon==True:

				    if temperature=="ima":
					    temperature="-100"
				    elif temperature=="rea":
					    temperature="-200"
				    else:

					    continue

			    else:
				    temperature = file_name_up[i][:-4].split('_')[-1][:-1]

                            data_by_temp_up[temperature] = read_transport_datafile(file_name_up[i]) 



                #############################################################################################

                #############################################################################################
                width = 20
                height = 14
                if spin_polarized==True:
                    height = 28
                pylab.figure(figsize=(width, height)) #to adjust the figure size
		matplotlib.rc("font", size=30)             #to set the font size
                markers=["o","s","8","D","H","1"]
                colors =['b','g','c','r','m','0.75','y','orange']
                lines  =['-','--',':',]              
		if epsilon==True:
			lines  =['-',]              
			colors =['k','r']

                if spin_polarized==False:
                        ax2=pylab.subplot(111)

			ax2.tick_params(axis='both', which='major', labelsize=21)

                        sorted_all =  sorted([[float(x[0]),x[1]] for x in data_by_temp.items()])

                        for temp_index in range(len(sorted_all)):              
		            if epsilon==True:
				    temp = int(sorted_all[temp_index][0])

				    if temp==-100:
					    label_text='$'+trans_plot_dict[Type]['lc']+'_{2}$'
				    elif temp==-200:
					    label_text='$'+trans_plot_dict[Type]['lc']+'_{1}$'

				    else:
					    continue


			    else:
				    label_text='$'+trans_plot_dict[Type]['lc']+'$ @ %sK'%int(sorted_all[temp_index][0])


                            ls_choice  = lines[temp_index%len(lines)]
                            color_choice  = colors[temp_index%len(colors)]
                            marker_choice = markers[temp_index%len(markers)]

                            x_vals=sorted_all[temp_index][1][0]
                            y_vals=sorted_all[temp_index][1][1]



			    if Type=='kappa':
				    y_vals=numpy.asarray(y_vals)/1.0e17
			    elif Type=='sig_seebeck':
				    y_vals=numpy.asarray(y_vals)/1.0e-9
			    elif Type=='cond':
				    y_vals=numpy.asarray(y_vals)/1.0e20
			    if x_range==None:
				    max_x=max(x_vals)
				    min_x=min(x_vals)
			    else:
				    min_x=x_range[0]
				    max_x=x_range[1]




                            for i in [y_vals]:
                                set_max=max([i[j] for j in range(len(i)) if (x_vals[j]>min_x and  x_vals[j]<max_x)])

                                if set_max>max_val:
                                    max_val=set_max
				
			    for i in [y_vals]:
			        set_min=min([i[j] for j in range(len(i)) if (x_vals[j]>min_x and  x_vals[j]<max_x)])

                                if set_min<min_val:
                                    min_val=set_min

			    x_vals,y_vals = AFLOWpi.plot.__smoothGauss(x_vals,degree=5),AFLOWpi.plot.__smoothGauss(y_vals,degree=5)
                            ax2.plot(x_vals,y_vals,label=label_text,color=color_choice,linestyle=ls_choice,linewidth=4)


                if spin_polarized==True:
                    ax1=pylab.subplot(211)


                    sorted_dn =  sorted([[float(x[0]),x[1]] for x in data_by_temp_dn.items()])
                    sorted_up =  sorted([[float(x[0]),x[1]] for x in data_by_temp_up.items()])

                    for temp_index in range(len(sorted_dn)):                    
		            if epsilon==True:
				    try:
					    if int(sorted_dn[temp_index][0])==-100:
						    label_text='$'+trans_plot_dict[Type]['lc']+'^{down}_{2}$'
					    if int(sorted_dn[temp_index][0])==-200:
						    label_text='$'+trans_plot_dict[Type]['lc']+'^{down}_{1}$'
				    except:
					    continue

			    else:
				    label_text='$'+trans_plot_dict[Type]['lc']+'^{down}$ @ %sK'%int(sorted_dn[temp_index][0])
                            color_choice  = colors[temp_index%len(colors)]
                            ls_choice  = lines[temp_index%len(lines)]
                            marker_choice = markers[temp_index%len(markers)]
                            x_vals=sorted_dn[temp_index][1][0]
                            y_vals=sorted_dn[temp_index][1][1]


			    if Type=='kappa':
				    y_vals=numpy.asarray(y_vals)/1.0e17
			    elif Type=='sig_seebeck':
				    y_vals=numpy.asarray(y_vals)/1.0e-3
			    elif Type=='cond':
				    y_vals=numpy.asarray(y_vals)/1.0e20

			    if x_range==None:
				    max_x=max(x_vals)
				    min_x=min(x_vals)
			    else:
				    min_x=x_range[0]
				    max_x=x_range[1]

			    if epsilon==True:
				    if min_x<=0.2:
					    min_x=0.2

                            for i in [y_vals]:
                                set_max=max([i[j] for j in range(len(i)) if (x_vals[j]>min_x and  x_vals[j]<max_x)])
                                if set_max>max_val:
                                    max_val=set_max
                            for i in [y_vals]:
                                set_min=min([i[j] for j in range(len(i)) if (x_vals[j]>min_x and  x_vals[j]<max_x)])
                                if set_min<min_val:
                                    min_val=set_min

			    x_vals,y_vals = AFLOWpi.plot.__smoothGauss(x_vals,degree=5),AFLOWpi.plot.__smoothGauss(y_vals,degree=5)
                            ax1.plot(x_vals,y_vals,label=label_text,color=color_choice,linestyle=ls_choice,linewidth=4)

			    pyplot.ylabel(trans_plot_dict[Type]['yl'],{'fontsize':30})

                            

                    pylab.axhline(0.0, color = 'k',linestyle='dashed', linewidth = 1.3)
#                    ax1.legend(loc=0,fontsize=26)
                    ax1.legend(loc=0)
		    if Type in ['ZT']:
			    ax1.yaxis.set_ticks([0],)


#		    pylab.xlabel(trans_plot_dict[Type]['xl'],{'fontsize':26})
                    for temp_index in range(len(sorted_up)):                    
		            if epsilon==True:
				    try:
					    if int(sorted_up[temp_index][0])==-100:
						    label_text='$'+trans_plot_dict[Type]['lc']+'^{up}_{2}$'
					    if int(sorted_up[temp_index][0])==-200:
						    label_text='$'+trans_plot_dict[Type]['lc']+'^{up}_{1}$'
				    except:
					    continue

			    else:
				    label_text='$'+trans_plot_dict[Type]['lc']+'^{up}$ @ %sK'%int(sorted_up[temp_index][0])

                            color_choice  = colors[temp_index%len(colors)]
                            marker_choice = markers[temp_index%len(markers)]
                            ls_choice  = lines[temp_index%len(lines)]
                            x_vals=sorted_up[temp_index][1][0]
                            y_vals=sorted_up[temp_index][1][1]


			    if x_range==None:
				    max_x=max(x_vals)
				    min_x=min(x_vals)                    
			    else:
				    min_x=x_range[0]
				    max_x=x_range[1]

			    if epsilon==True:
				    if min_x<=0.2:
					    min_x=0.2

			    if Type=='kappa':
				    y_vals=numpy.asarray(y_vals)/1.0e17
			    elif Type=='sigma_seebeck':
				    y_vals=numpy.asarray(y_vals)/1.0e-3
			    elif Type=='cond':
				    y_vals=numpy.asarray(y_vals)/1.0e20



                            for i in [y_vals]:
                                set_max=max([i[j] for j in range(len(i)) if (x_vals[j]>min_x and  x_vals[j]<max_x)])
                                if set_max>max_val:
                                    max_val=set_max
                            for i in [y_vals]:
                                set_min=min([i[j] for j in range(len(i)) if (x_vals[j]>min_x and  x_vals[j]<max_x)])
                                if set_min<min_val:
                                    min_val=set_min

			    ax2=pylab.subplot(212)
#[left,bottom,width,height]

#			    ax2.tick_params(axis='both', which='major', labelsize=2)


			    x_vals,y_vals = AFLOWpi.plot.__smoothGauss(x_vals,degree=5),AFLOWpi.plot.__smoothGauss(y_vals,degree=5)
                            ax2.plot(x_vals,y_vals,label=label_text,color=color_choice,linestyle=ls_choice,linewidth=4)

			    pyplot.ylabel(trans_plot_dict[Type]['yl'],{'fontsize':30})
			    ax1.set_position([0.04,0.465,0.96,0.435]) 
			    ax2.set_position([0.04,0.02,0.96,0.435]) 
#			    pyplot.xlabel(trans_plot_dict[Type]['xl'],{'fontsize':22})
#			    pyplot.ylabel(trans_plot_dict[Type]['yl'],{'fontsize':22})
#                            ax2.xaxis.set_label(trans_plot_dict[Type]['xl'],)
 #                           ax2.yaxis.set_label(trans_plot_dict[Type]['yl'],)

 #line separating up and down spni

		if x_range==None:
			max_x=max(x_vals)
			min_x=min(x_vals)
		else:
			min_x=x_range[0]
			max_x=x_range[1]

		if epsilon==True:
			if min_x<=0.2:
				min_x=0.2

			

                mult_max=1.05
                if max_val<0.0:
                    mult_max=0.95
                mult_min=1.05
                if min_val>0.0:
                    mult_min=0.95
		


                try:
                    try:
                        ax1.set_ylim([min_val*mult_min,max_val*mult_max])
                        ax1.axhline(0.0, color = 'k',linestyle='dashed', linewidth = 1.3) 

			if Type in ['ZT']:
				ax1.yaxis.set_ticks([0.0])

                        ax1.xaxis.set_ticks([])

			if epsilon==True:
				if min_x<=0.2:
					min_x=0.2
                        ax1.set_xlim([min_x,max_x])

                    except:
                        pass

                    pyplot.xlabel(trans_plot_dict[Type]['xl'],{'fontsize':30})
                    pyplot.ylabel(trans_plot_dict[Type]['yl'],{'fontsize':30})

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)
                    pass
                ax2.set_ylim([min_val*mult_min,max_val*mult_max])
#                ax2.legend(loc=0,fontsize=26)
                ax2.legend(loc=0)

		if epsilon==True:
			if min_x<=0.2:
				min_x=0.2
                ax2.set_xlim([min_x,max_x])
		if Type in ['ZT']:
			ax2.yaxis.set_ticks([0.0])

#                pyplot.xlabel(trans_plot_dict[Type]['xl'],{'fontsize':22})


                compoundName = AFLOWpi.retr._getStoicName(oneCalc,strip=True)
                compoundNameLatex = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)
#                ax2.legend(loc=0,fontsize=26)
                ax2.legend(loc=0)

		if Type in ['ZT']:
			ax2.yaxis.set_ticks([0])

#                ax2.axes.yaxis.set_ticks_position('right')
		try:
			ax2.yaxis.set_label(trans_plot_dict[Type]['yl'])
			ax2.axhline(0.0, color = 'k',linestyle='dashed', linewidth = 1.3)
		except Exception,e: 
			AFLOWpi.run._fancy_error_log(e)

			pass
		try:
			ax1.yaxis.set_label(trans_plot_dict[Type]['yl'])
			ax1.axhline(0.0, color = 'k',linestyle='dashed', linewidth = 1.3) #line separating up and down spni

		except Exception,e:
			pass

                pylab.axhline(0.0, color = 'k',linestyle='dashed', linewidth = 1.3) #line separating up and down spin

                figtitle = r'%s %s' % (trans_plot_dict[Type]['ft'],compoundNameLatex)

                t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=50,horizontalalignment='center') #[x,y] 

                if postfix!='':
                    postfix='_'+postfix

                fileplot = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_%s_%s%s.pdf'%(trans_plot_dict[Type]['fp'],compoundName,ID,postfix))
                matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')


            except Exception,e:
		    AFLOWpi.run._fancy_error_log(e)
		    raise SystemExit



###############################################################################################################################
def read_transport_datafile(ep_data_file,mult_x=1.0,mult_y=1.0):
	'''


	Arguments:


	Keyword Arguments:


        Returns:

	'''

	with open(ep_data_file,'r') as epsilon_data_file:
		ep_data_string=epsilon_data_file.read()

	lines=ep_data_string.split('\n')[2:]
	en_array=[]
	val_array=[]

	for i in lines:
		try:
			float_array=[]
			split_line = i.split()

			for x in split_line:
				float_array.append(float(x))
			en_array.append(float_array[0])
			val_array.append(sum([float_array[1],float_array[5],float_array[9]]))
		except:
			pass


		en_array=[mult_x*x for x in en_array]
		val_array=[mult_y*y for y in val_array]

		
	return en_array,val_array	




###############################################################################################################################

import AFLOWpi
import os
from matplotlib import pyplot
from matplotlib import pylab
import numpy
import matplotlib
import StringIO
import logging 




def __plot_gruneisen(oneCalc,ID,optical=False):
	"""
	Function to take the data files generated by the sumpdos ppBands functions and plots the electronic band structure and the projected density of states with energy shifted relative to the Fermi Energy.
	
	Arguments:
              oneCalc (dict): Single calculation that is the value of the dictionary of dictionaries of calculations
              ID (str): ID of the calculation	
	Keyword Arguments:
              postfix (str): Output filename postfix
              THz (bool): Plot the frequencies in THz or cm^-1

        Returns:
              None

	"""
        grun_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal')
        grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phGRUN.gp'%grun_ID)

        with open(grun_file_name,'r') as gpfo:
            grun_data_str = gpfo.read()
        postfix=''
        THz=False


	
	calcCopy = oneCalc
	calcID = ID

        calcID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='phonon',last=True)        


	if postfix!='':
		postfix='_'+postfix

	subdir=oneCalc['_AFLOWPI_FOLDER_']
       	fileplot = os.path.join(subdir,'GRUNEISEN_%s_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,postfix))

	"""get the path to the subdirectory of the calc that you are making plots for"""

#        filebands = os.path.join(subdir,'%s.phBAND.gp'%ID)
        filebands = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phGRUN.gp'%grun_ID)
        filedos = os.path.join(subdir,'%s.phDOS.gp'%ID)
        filedos = os.path.join(subdir,'%s.phdos'%ID)

	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
	

       #to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=14)             #to set the font size

        
        width = 20
	height = 8
	pylab.figure(figsize=(width, height))#to adjust the figure size
     

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
        #scale for THz or cm^-1
        if THz==True:
            sf=0.0299792458
        else:
            sf=1.0

	try:
                max_val=0.0
                min_val=0.0
		with open(filebands,'r') as datafile:
			data =datafile.readlines()


		for line in data:
                    try:
                        if line.split()[0]==0.0:
                            continue
                        x_val = float(line.split()[0])
                        
                        y_val = [float(z)*sf for z in line.split()[1:]]
                        
                        x.append(x_val)
                        y.append(y_val)
                    except Exception,e:
                        continue

                k_y = numpy.asarray(y).T.tolist()

                for entry in k_y:
                    k_x.append(x)


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
		
				
	a=k_x[1]   # a set of k point values for one band for axis scaling purposes
	b=k_y[1]



        ax1=pylab.subplot(111)	

	'''
        Plot each band (k_x[i]),(k_y[i]) on the band structure plot from list of
        the values of energy and position from their respective list by k points
	'''
        colors =['b','g','c','r','m','y','orange']

	branch_list=['TA',"TA'",'LA','Optical']

	for i in range(len(k_x)):
#	    k_y_prime = [k_y[i][y] for y in range(len(k_y[i])) if not y%10]
#            k_y_prime=k_y[i]

	    k_y_prime = [k_y[i][y] for y in range(len(k_y[i])) if not y%51]
	    k_x_prime = [k_x[i][x] for x in range(len(k_x[i])) if not x%51]
#            k_y_prime=k_y[i]
#            k_x_prime=k_x[i]
	    k_y_prime = numpy.ma.masked_equal(k_y_prime,0.0)


            new_min_val=numpy.amin(k_y_prime)
            new_max_val=numpy.amax(k_y_prime)
            color_choice=colors[i%len(colors)]
	    if optical==False:
		    if i<3:
			    if new_max_val>max_val:
				    max_val=new_max_val
			    if new_min_val<min_val:
				    min_val=new_min_val
	    else:
		    if new_max_val>max_val:
			    max_val=new_max_val
		    if new_min_val<min_val:
			    min_val=new_min_val
				    


#	    k_y_prime = numpy.ma.masked_outside(k_y_prime,-10.0,10.0)

	    if i<3:
 		    pylab.plot(k_x_prime,k_y_prime,color=color_choice,linestyle=' ',marker='^',label=branch_list[i])
	    if optical==True:
		    if i==3:
			    pylab.plot(k_x_prime,k_y_prime,color='k',linestyle=' ',marker='^',label=branch_list[-1])
		    else:
			    pylab.plot(k_x_prime,k_y_prime,color='k',linestyle=' ',marker='^')
		    
#            pylab.plot((k_x[i]),(k_y[i]),color=color_choice)


	    
	pylab.ylabel('$Gr\ddotuneisen$ $parameter$')


	pylab.xlim(min(k_x[1]),max(k_x[1])) 
        
	print min_val,max_val
        pylab.ylim(min_val,max_val) 

	'''
        takes in a list of k points that was used as pw.x input for the 'bands'
        calculation as a string. It puts parts of that string into lists and 
        manipulates them to display the symmetry point boundary lines on the 
        band structure plot
	'''
        ph_band_in = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_matdyn_phBand.in'%ID)
#        print ph_band_in
        with open(ph_band_in,'r') as ph_band_in_file:
            ph_band_in_file_string = ph_band_in_file.read()
        bandSym=ph_band_in_file_string.split('/')[-1]
 	bandSymSplit =  bandSym.split()

 	HSPList = []
 	HSPSymList = []
 	buf = StringIO.StringIO(bandSym)

	for line in buf:
		splitLine = line.split()
		if len(splitLine)==2: # new style kpoint path input case
			HSPList.append(splitLine[1])
			specialPointName = splitLine[0].rstrip()
			
 	                #renames gG to greek letter for capital gamma
			if specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gG':
				specialPointName = r"$\Gamma$"

			elif len(specialPointName) != 1:
				specialPointName = "$"+specialPointName[0]+r'_{'+specialPointName[1]+'}$' #if there is a subscript it makes the print out on the plot have the number subscripted 
			else:
				specialPointName = "$"+specialPointName[0]+"$" #formats with internal math renderer so all the labels look the same
			HSPSymList.append(specialPointName)
			
		elif len(splitLine)==7: # old style kpoint path input case with kpoint names
			HSPList.append(splitLine[3])
			specialPointName = splitLine[4][1:].rstrip()

			if specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gG': #renames gG to greek letter for capital gamma
				specialPointName = r"$\Gamma$"
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
		if i==0: # for the first k point in the first path
			SymPrint.append(HSPSymList[i])
		if int(HSPList[i]) == 0 and i<len(HSPList)-2: # for the end of a path (where the number of k points between one point and another is zero)
			continue
		elif int(HSPList[i+1]) == 0 and i!=len(HSPList)-2: # for the point that begins a new path where the end of the last path (which has zero k points from it to this k point)
			totalX +=(int(HSPList[i])+1)
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

     #add symmetry lines to the band structure plot
#        print symIndex

	for sym in symIndex:
		try:
			pylab.axvline(a[sym], color = 'k')
		except Exception,e:

                    pass
     #Print path labels to band structure x-axis
	try:
            labs=[]
            for index in range(len(symIndex)):
                if index!=0:
                    labs.append(a[symIndex[index]-1])
                else:
                    labs.append(0)

            pylab.xticks(labs,SymPrint)
	except Exception,e:
#            print 'asdfasdfasdfasdf'
            print e
            pass
    #        return
	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Femi level line
	locs, labels = pylab.xticks()

##########################################################################################################


        print 'Plotting q resolved Gruneisen Parameter of  %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
        logging.info('Plotting q resolved Gruneisen Parameter of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))

#        ax2=pylab.subplot(122)

        try:
                data = open(filedos,'r').readlines()
        except Exception:
                logging.warning("output from dos calculation not found. Are you sure you ran ppDOS and it completed properly?")
                print "Are you sure you ran ppDOS and it completed properly?"
                return
        

        freq_dos=[]
        dos=[]
        plot_dos_x=[]
        for i in range(len(data)):      #convert to floats
            dat = [float(x) for x in data[i].split()]
            freq_dos.append(dat[0]*sf)

            dos.append(dat[1])

        plot_dos_x=[]
        plot_dos_y=[]
        pre_sort=[]
	#smooth the phDOS
	dos = AFLOWpi.plot.__smoothGauss(dos)   
	freq_dos = AFLOWpi.plot.__smoothGauss(freq_dos)   

#        pylab.plot(dos,freq_dos,'k',linestyle='-') #to plot the smoothed data
        pylab.ylim(min_val,max_val) 

#        ax2.spines['bottom'].set_linewidth(1.5)
#        ax2.spines['left'].set_linewidth(1.5)
#        ax2.spines['right'].set_linewidth(1.5)
#        ax2.spines['top'].set_linewidth(1.5)
#        ax1.yaxis.set_ticks([])
#        ax1.xaxis.set_ticks([])
#        ax2.yaxis.set_ticks_position('left')
#        pylab.xlabel('Phonon Density of States (arb. units)')
#        pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line


#        ax1.set_position([0.07,0.1,0.67,0.8]) #[left,bottom,width,height]
#        ax2.set_position([0.75,0.1,0.20,0.8])

 	figtitle = '$Gr\ddotuneisen$ $Parameter$: %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)) 
 	t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=14,horizontalalignment='center') #[x,y]

	handles, labels = ax1.get_legend_handles_labels()
	ax1.legend(handles, labels,numpoints=1)

  	matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')

        pyplot.cla()
 	pyplot.clf()
 	pyplot.close()




def __gruneisen_of_omega(oneCalc,ID,projected=True):
    matplotlib.rc("font", size=24)             #to set the font size

    therm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal')

    therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER.gp'%therm_ID)
#    therm_data =numpy.loadtxt(therm_file_name,dtype=numpy.float64,)

    therm_data=[]

    with open(therm_file_name,"r") as fo:
        data = fo.read()

    data=data.split('\n')
    therm_data = numpy.asarray([map(float,x.split()) for x in data if len(x.strip())!=0])
    
    therm_data=numpy.asarray(therm_data)
    therm_data[:,0]=numpy.around(therm_data[:,0],decimals=0)
    therm_data[:,1]=numpy.around(therm_data[:,1],decimals=2)

#    therm_data = numpy.unique(b).view(therm_data.dtype).reshape(-1, therm_data.shape[1])
    
    therm_data  = numpy.vstack({tuple(row) for row in therm_data})
#   therm_data[:,0] = numpy.ma.masked_where(therm_data[:,0] <= 1.0, therm_data, copy=False)
    for i in range(len(therm_data)):
	    if therm_data[i][0]<1.1:
		    therm_data[i][1]=0.0
	    



    therm_data[:,1] = numpy.ma.masked_equal(therm_data[:,1],0.0)
#    therm_data = numpy.ma.masked_equal(therm_data,0.0)
    therm_data[:,1] = numpy.ma.masked_greater(therm_data[:,1],100.0)
#    therm_data[:,0] = numpy.ma.masked_less_equal(therm_data[:,0],1.0)
    print therm_data.shape
    width = 10.0
    height = 8.0
    pylab.figure(figsize=(width, height))#to adjust the figure size
    
    ax1=pylab.subplot(111)	
    pylab.ylabel('$\gamma^{2}$')
    pylab.xlabel('$\omega$ $(cm^{-1})$')
    pylab.plot(therm_data[:,0],therm_data[:,1],'k',linestyle=' ',marker='o',fillstyle='none')

    figtitle = '$Gr\ddotuneisen$ $Parameter:$ %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)) 
    t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=24,horizontalalignment='center') #[x,y]
    pylab.xlim([0.0,numpy.amax(therm_data[:,0])])
    fileplot = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'SCATTER_%s_%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,))
	
    matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')

    pyplot.cla()
    pyplot.clf()
    pyplot.close()

    

# def __plot_thermal_conductivity(oneCalc,ID):

#     therm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal')
#     therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_thermal_cond.dat'%therm_ID)

#     with open(therm_file_name,'r') as tcfo:
#         therm_data_str = tcfo.read()


import 	matplotlib.pyplot

def _plot_lattice_TC(oneCalc,ID,temp_range=[80.0,800.0]):

	fname = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],)

	therm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal')
	therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_thermal_cond.dat'%therm_ID)

	data=numpy.loadtxt(therm_file_name, dtype=numpy.float, comments='#',)
	matplotlib.pyplot.plot(data[:,0],data[:,1])

	data_cut = scipy.where(data[0]>temp_range[0] and data[0]<temp_range[1],data)
	max_y = numpy.amax(data_cut[:,1])*1.1
	matplotlib.pyplot.ylim([0.0,max_y])

	matplotlib.pyplot.xlim(temp_range)
	matplotlib.pyplot.xlabel('T (K)')

	figtitle = 'Lattice Thermal Conductivity: %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)) 
	t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=14,horizontalalignment='center') #[x,y]

	matplotlib.pyplot.ylabel('$\kappa_{lat}$ $(\frac{W}{m\cdotK})$')


	fig_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"LATTICE_TC_%s_%s.pdf"%(AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID))
	matplotlib.pyplot.savefig(fig_file_name,bbox_inches='tight')

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



def phonon(calcs,runlocal=False,postfix='',THz=True,color_accoustic=False,color_optical=False):
	if runlocal:
		for ID,oneCalc in calcs.iteritems():
			AFLOWpi.plot.__plot_phonon(oneCalc,ID,postfix=postfix,THz=THz,color_accoustic=color_accoustic,color_optical=color_optical)
	else:
		AFLOWpi.prep._addToAll(calcs,'PLOT',"AFLOWpi.plot.__plot_phonon(oneCalc,ID,postfix=%s,THz=%s,color_accoustic=%s,color_optical=%s)"%(repr(postfix),THz,color_accoustic,color_optical))



def __plot_phonon(oneCalc,ID,postfix='',THz=True,color_accoustic=False,color_optical=False,DOSPlot='APDOS'): 
	"""
	Function to take the data files generated by the sumpdos ppBands functions and plots the electronic band structure and the projected density of states with energy shifted relative to the Fermi Energy.
	
	Arguments:
              oneCalc (dict): Single calculation that is the value of the dictionary of dictionaries of calculations
              ID (str): ID of the calculation	
	Keyword Arguments:
              postfix (str): Output filename postfix
              THz (bool): Plot the frequencies in THz or cm^-1

        Returns:
              None

	"""
	
	calcCopy = oneCalc
	calcID = ID

        calcID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='phonon',last=True)        


	if postfix!='':
		postfix='_'+postfix

	subdir=oneCalc['_AFLOWPI_FOLDER_']
       	fileplot = os.path.join(subdir,'PHONON_%s_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,postfix))

	"""get the path to the subdirectory of the calc that you are making plots for"""

        filebands = os.path.join(subdir,'%s.phBAND.gp'%ID)




	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''	


       #to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=20)             #to set the font size

        
        width = 20
	height = 8
	pylab.figure(figsize=(width, height))#to adjust the figure size
     

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
        #scale for THz or cm^-1
        if THz==True:
            sf=0.0299792458
        else:
            sf=1.0

	try:
                max_val=0.0
                min_val=0.0
		with open(filebands,'r') as datafile:
			data =datafile.readlines()



		for line in data:

                    try:
                        x_val = float(line.split()[0])
                        
                        y_val = [float(z)*sf for z in line.split()[1:]]
                        
                        x.append(x_val)
                        y.append(y_val)
                    except Exception,e:
                        pass

                k_y = numpy.asarray(y).T.tolist()

                for entry in k_y:
                    k_x.append(x)


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
		
				
	a=k_x[1]   # a set of k point values for one band for axis scaling purposes
	b=k_y[1]



        ax1=pylab.subplot(121)	

	'''
        Plot each band (k_x[i]),(k_y[i]) on the band structure plot from list of
        the values of energy and position from their respective list by k points
	'''
        colors =['b','g','c','r','m','y','orange']

#	for j in range(len(k_y[0])):
#		if k_y[0][j]==0.0 and k_y[1][j]==0.0 and k_y[1][j]==0.0:
#			for i in range(len(k_y)):
#				k_y[i][j]=0.0
#	gamma_index=
#	for i in range(len(k_y[0])):
	
####################################################################################
#	gamma_index=[]
#	for sym in range(len(SymPrint)):
#		if SymPrint[sym]=='$\\Gamma$':
#			gamma_index.append(symIndex[sym])


####################################################################################

	#set frequency of optical branches at gamma to either one before or one after 
	#so that the LOTO splitting at gamma doesn't result in an outlier point

	gamma_index=[]
	by_k = numpy.asarray(k_y).T
	for i in range(len(by_k)):
		num_zero =  len(by_k[i])-numpy.count_nonzero(by_k[i])
		if num_zero>=3:
			gamma_index.append(i)

	for j in gamma_index:
		for i in range(len(k_y)):
			k_y[i][j]=0.0

#		if k_y[0][j]==0.0:
#			for i in range(len(k_y)):
##				if i not in [0,1,2]:
#				else:
#				k_y[i][j]=0.0
#					try:
#						k_y[i][j]=k_y[i][j-1]
#					except:	
#						k_y[i][j]=k_y[i][j+1]



	for i in range(len(k_x)):
            new_min_val=min(k_y[i])
            new_max_val=max(k_y[i])
            color_choice=colors[i%len(colors)]
	    if color_accoustic==False:
		    if i<3:
			    color_choice='k'
	    if color_optical==False:
		    if i>2:
			    color_choice='k'

            if new_max_val>max_val:
                max_val=new_max_val
            if new_min_val<min_val:
                min_val=new_min_val

#	    if i>2:


            pylab.plot((k_x[i]),(k_y[i]),color=color_choice,linestyle='None',marker='.')



        if THz==True:
            pylab.ylabel('Frequency (THz)')
        else:
            pylab.ylabel('Frequency (cm$^{-1}$)')
	pylab.xlim(min(k_x[1]),max(k_x[1])) 
        



        pylab.ylim(1.1*min_val,max_val*1.1) 

	'''
        takes in a list of k points that was used as pw.x input for the 'bands'
        calculation as a string. It puts parts of that string into lists and 
        manipulates them to display the symmetry point boundary lines on the 
        band structure plot
	'''
        ph_band_in = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_matdyn_phBand.in'%ID)

        with open(ph_band_in,'r') as ph_band_in_file:
            ph_band_in_file_string = ph_band_in_file.read()
        bandSym=ph_band_in_file_string.split('/')[-1]
 	bandSymSplit =  bandSym.split()

 	HSPList = []
 	HSPSymList = []
 	buf = StringIO.StringIO(bandSym)

	for line in buf:
		splitLine = line.split()
		if len(splitLine)==2: # new style kpoint path input case
			HSPList.append(splitLine[1])
			specialPointName = splitLine[0].rstrip()
			
 	                #renames gG to greek letter for capital gamma
			if specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gG':
				specialPointName = r"$\Gamma$"

			elif len(specialPointName) != 1:
				specialPointName = "$"+specialPointName[0]+r'_{'+specialPointName[1]+'}$' #if there is a subscript it makes the print out on the plot have the number subscripted 
			else:
				specialPointName = "$"+specialPointName[0]+"$" #formats with internal math renderer so all the labels look the same
			HSPSymList.append(specialPointName)
			
		elif len(splitLine)==6: # old style kpoint path input case with kpoint names
			HSPList.append(splitLine[3])
			specialPointName = splitLine[5].rstrip()

			if specialPointName == 'G' or specialPointName == 'g' or specialPointName == 'Gamma' or specialPointName == 'gG': #renames gG to greek letter for capital gamma
				specialPointName = r"$\Gamma$"
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
		if i==0: # for the first k point in the first path
			SymPrint.append(HSPSymList[i])
		if int(HSPList[i]) == 0 and i<len(HSPList)-2: # for the end of a path (where the number of k points between one point and another is zero)
			continue
		elif int(HSPList[i+1]) == 0 and i!=len(HSPList)-2: # for the point that begins a new path where the end of the last path (which has zero k points from it to this k point)
			totalX +=(int(HSPList[i])+1)
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


     #add symmetry lines to the band structure plot
	for sym in symIndex:
		try:
			pylab.axvline(a[sym], color = 'k')
		except Exception,e:
			pass
     #Print path labels to band structure x-axis
	try:
		pylab.xticks([a[index] for index in symIndex],SymPrint)
	except Exception,e:
		pass
		return
	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Femi level line
	locs, labels = pylab.xticks()



##########################################################################################################
	ax2=pylab.subplot(122)

	if DOSPlot!='APDOS':
		print 'Plotting Phonons and phDOS of  %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting Phonons and phDOS of  %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))

		ax2=pylab.subplot(122)
		ax2.set_color_cycle(['k','r','g','b','c', 'm', 'y',])

		filedos = os.path.join(subdir,'%s.phdos'%ID)
		try:
			data = open(filedos,'r').readlines()
		except Exception:
			logging.warning("output from dos calculation not found. Are you sure you ran ppDOS and it completed properly?")
			print "Are you sure you ran ppDOS and it completed properly?"
			return

		freq_dos=[]
		dos=[]
		plot_dos_x=[]
		for i in range(len(data)):      #convert to floats
		    dat = [float(x) for x in data[i].split()]
		    freq_dos.append(dat[0]*sf)

		    dos.append(dat[1])

		plot_dos_x=[]
		plot_dos_y=[]
		pre_sort=[]
		#smooth the phDOS
	#	dos = AFLOWpi.plot.__smoothGauss(dos)   
	#	freq_dos = AFLOWpi.plot.__smoothGauss(freq_dos)   

		pylab.plot(dos,freq_dos,'k',linestyle='-', linewidth=3) #to plot the smoothed data

##########################################################################################################

	if DOSPlot=='APDOS':
		print 'Plotting Phonons and atom projected phDOS of  %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting Phonons and atom projected phDOS of  %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))




		filedos = os.path.join(subdir,'%s.aphdos'%ID)
		try:
			data = open(filedos,'r').read()
			data=data.split('\n')
			labels=data[0].split()
			data=data[1:]

		except Exception:
			logging.warning("output from dos calculation not found. Are you sure you ran ppDOS and it completed properly?")
			print "Are you sure you ran ppDOS and it completed properly?"
			return

		freq_dos=[]
		dos=[]
		plot_dos_x=[]
		for i in range(len(data)):      #convert to floats
		    dat = [float(x) for x in data[i].split()]
		    freq_dos.append(dat)

#		    dos.append(dat[1:])
		freq_dos=numpy.asarray(freq_dos)
		freq_dos[:,0]*=sf


		plot_dos_x=[]
		plot_dos_y=[]
		pre_sort=[]
		#smooth the phDOS
	#	dos = AFLOWpi.plot.__smoothGauss(dos)   
	#	freq_dos = AFLOWpi.plot.__smoothGauss(freq_dos)   

		min_xval = numpy.amin(freq_dos[:,2:])
		max_xval = numpy.amax(freq_dos[:,2:])
		for spec in range(2,len(labels)):
			ax2.plot(freq_dos[:,spec],freq_dos[:,0],linestyle='-', linewidth=3,label=labels[spec]) #to plot the smoothed data
		ax2.set_xlim(0.0,max_xval*1.2) 
		ax2.legend(fontsize=14)

#		handles, labels = ax2.get_legend_handles_labels()
##########################################################################################################

        ax1.set_ylim(1.1*min_val,max_val*1.1) 
        ax2.set_ylim(1.1*min_val,max_val*1.1) 

        ax2.spines['bottom'].set_linewidth(1.5)
        ax2.spines['left'].set_linewidth(1.5)
        ax2.spines['right'].set_linewidth(1.5)
        ax2.spines['top'].set_linewidth(1.5)
        ax2.yaxis.set_ticks([])
        ax2.xaxis.set_ticks([])
        ax2.yaxis.set_ticks_position('left')
        pylab.xlabel('Density of States (arb. units)')
        pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line



        ax1.set_position([0.07,0.1,0.69,0.8]) #[left,bottom,width,height]
        ax2.set_position([0.77,0.1,0.23,0.8])

 	figtitle = 'Phonon Dispersion and DOS: %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)) 
 	t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=24,horizontalalignment='center') #[x,y]

  	matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')

        pyplot.cla()
 	pyplot.clf()
 	pyplot.close()

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
import cPickle
import matplotlib.lines as mlines


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


def __bands(oneCalc,ID,yLim=[-10,10],DOSPlot='',postfix='',tight_banding=False):
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
        __bandPlot(oneCalc,yLim,DOSPlot,postfix=postfix,tight_banding=tight_banding)





def __getPath_WanT(oneCalc,ID):
    '''


    Arguments:


    Keyword Arguments:

        
    Returns:


    '''

    try:
        want_stdout_path = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_bands.out'%ID)[-1]
    except:
        want_stdout_path = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_bands_up.out'%ID)[-1]

    with open(want_stdout_path,'r') as in_file_obj:
        in_string = in_file_obj.read()

    path = AFLOWpi.retr._getPath(0.01,oneCalc,ID=ID)

    plot_bool=[]
    path_name=[]
    path_split = [x for x in  path.split('\n')[1:] if len(x.strip())]

    for i in path_split:
        path_name.append(i.split()[-1])
        if  int(i.split()[3]):

            plot_bool.append(True)
        else:
            plot_bool.append(False)

    gg = re.findall("\s*Number of kpts in each segment\n((?:.*:\W+(?:\d*)\n)*)",in_string)
    print gg

    num = [int(x) for x in re.findall('line\s*\d+:\s*(\d+)\s*\n',in_string)]
    print re.findall('line\s*\d+:\s*(\d+)\s*\n',in_string)

    total = 0
    include=[]

    output_path_string = ''
    print plot_bool
    for i in range(len(num)):

        total+=num[i]+1
	try:
		if plot_bool[i]:

		    if i==0:
			output_path_string+='%s %s %s %s ! %s\n' %(0.0,0.0,0.0,num[i],path_name[i])
		    else:
			output_path_string+='%s %s %s %s ! %s\n' %(0.0,0.0,0.0,num[i]+1,path_name[i],)
		else:
                    if i!=len(num):
			    output_path_string+='%s %s %s %s ! %s\n' %(0.0,0.0,0.0,0,path_name[i])
		for j in range(num[i]):
		    include.append(plot_bool[i])
#		else:
#			include.append(False)
		include.append(True)
		
	except Exception,e:
		print e

    output_path_string+='%s %s %s %s ! %s' %(0.0,0.0,0.0,0,path_name[-1])+'\n' 



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

def __bandPlot(oneCalc,yLim=[-10,10],DOSPlot='',postfix='',tight_banding=False): 
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
			Efermi=AFLOWpi.retr._getEfermi(oneCalc,'%s_WanT_dos'%calcID,directID=True)

			
                else:
			Efermi=AFLOWpi.retr._getEfermi(oneCalc,calcID)
			print 'EFERMI BANDS NO TB',Efermi

		if type(Efermi)!=type(0.5):
			Efermi=Efermi[0]

	except Exception,e:
            print e
            Efermi=0.0

	subdir=oneCalc['_AFLOWPI_FOLDER_']
       	fileplot = os.path.join(subdir,'BANDS_%s_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix))

	"""get the path to the subdirectory of the calc that you are making plots for"""

	if tight_banding==True:
		try:
			AFLOWpi.prep._clean_want_bands(oneCalc,calcID)
		except:
			return
		filebands = os.path.join(subdir,'%s_bands_want_cleaned.dat'%calcID)
		if not os.path.exists(filebands):

			filebands_up = os.path.join(subdir,'%s_bands_want_up_cleaned.dat'%calcID)
			filebands_dn = os.path.join(subdir,'%s_bands_want_down_cleaned.dat'%calcID)
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

        
        width = 20
	height = 14
	pylab.figure(figsize=(width, height))#to adjust the figure size
	
     #to do the gaussian smoothing               
 
     #################################
	if nspin==2:
#		AFLOWpi.prep._clean_want_bands(oneCalc,ID)
		k_x_up,k_y_up = AFLOWpi.plot._clean_bands_data_qe(filebands_up,Efermi_shift)
		k_x_dn,k_y_dn = AFLOWpi.plot._clean_bands_data_qe(filebands_dn,Efermi_shift)
		k_x,k_y=k_x_up,k_y_up
	else:
		k_x,k_y = AFLOWpi.plot._clean_bands_data_qe(filebands,Efermi_shift)

	a=k_x[1]   # a set of k point values for one band for axis scaling purposes
       	b=k_y[1]



	if DOSPlot != '':
		ax1=pylab.subplot(121)

	if DOSPlot == '':
		ax1=pylab.subplot(111)	
		print 'Plotting electronic band structure of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))

	'''
        Plot each band (k_x[i]),(k_y[i]) on the band structure plot from list of
        the values of energy and position from their respective list by k points
	'''
	if nspin==2:
		for i in range(len(k_x_up)):
			try:
				
				if tight_banding==True:
					pylab.plot((k_x_up[i]),(k_y_up[i]),'r',alpha=0.5,marker=".",linestyle=" ",label="$\uparrow$",linewidth=2)
					pylab.plot((k_x_dn[i]),(k_y_dn[i]),'k',alpha=0.5,marker=".",linestyle=" ",label="$\downarrow$",linewidth=2)
				else:
					pylab.plot((k_x_up[i]),(k_y_up[i]),'r',alpha=0.5,marker=".",linestyle=" ",label="$\uparrow$",linewidth=2)
					pylab.plot((k_x_dn[i]),(k_y_dn[i]),'k',alpha=0.5,marker=".",linestyle=" ",label="$\downarrow$",linewidth=2)
			except:
				pass
		handles, labels = ax1.get_legend_handles_labels()
#		ax1.legend(handles[-2:], labels[-2:],numpoints=1)
		up_legend_lab = mlines.Line2D([], [], color='red',label='$\uparrow$')
		dn_legend_lab = mlines.Line2D([], [], color='black',label='$\downarrow$')
		ax1.legend(handles=[up_legend_lab,dn_legend_lab])
	else:
		if tight_banding==True:
			for i in range(len(k_x)):
				pylab.plot((k_x[i]),(k_y[i]),'k',marker=".",linestyle=" ",linewidth=1.3)	
#				pylab.plot((k_x[i]),(k_y[i]),'k',linewidth=1.3)	
		else:
			for i in range(len(k_x)):
				pylab.plot((k_x[i]),(k_y[i]),'k',marker=".",linestyle=" ",linewidth=1.3)			
	#


	pylab.ylabel('E(eV)')
	pylab.xlim(min(k_x[1]),max(k_x[1])) 
	pylab.ylim(yLim[0],yLim[1])    
	pylab.yticks(numpy.arange(yLim[0],yLim[1]+1,2))

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
			totalX +=(int(HSPList[i])+1)
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
                print a[sym]
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
		pylab.xticks(bars,SymPrint, fontsize = 20)
	except Exception,e:
                print e
		pylab.xticks([a[-1] for index in symIndex],SymPrint)

	pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Femi level line
	locs, labels = pylab.xticks()

##########################################################################################################


	if DOSPlot == 'APDOS':
		print 'Plotting electronic band structure and projected DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure and projected DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))

		fileplot = os.path.join(subdir,'BANDPDOS_%s_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),calcID,postfix))	
		ax2=pylab.subplot(122)

		def getPlotData(sumpdosFile):
			try:
				with open(sumpdosFile,'rb') as dataFile:
					data = cPickle.load(dataFile)

			except Exception,e:
				print e
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


		atomList = list(set(AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])))

		for species in atomList:
			filePath = os.path.join(subdir,'%s_All.sumpdos' % species)
			if os.path.exists(filePath):
				pDOSNoPath.append(filePath)
		color = 'k'
		ax2.set_color_cycle(['r','g','b','c', 'm', 'y', 'k'])
		if LSDA:
			ax2.set_color_cycle(['r','r','g','g','b','b','c','c','m','m', 'y','y','k','k'])
		
		for filepath in pDOSNoPath:
			filename =  filepath.split('/')[-1]

			try:
				species = re.findall(r'(.+?)_.+',filename)[0]
			except:
				if filename=='%s_dos.dat'%dos_ID:
					species='TOTAL'
			
			'''gets the energy and the DOS for the orbital for the atom'''
			enshift, floatdos,floatdosDOWN = getPlotData(filepath)
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
			if species=='TOTAL':	
				pylab.plot(floatdos,enshift,'-',label=species,color='k',linewidth=2)		
			else:
				pylab.plot(floatdos,enshift,'-',label=species,linewidth=2)						

			if LSDA:
				if species=='TOTAL':
					pylab.plot(floatdosDOWN,enshift,'-',label=species,color='k',linewidth=2)		
				else:
					pylab.plot(floatdosDOWN,enshift,'-',label=species,linewidth=2)					   

		handles, labels = ax2.get_legend_handles_labels()

		if LSDA:
			ax2.legend(handles[::-2], labels[::-2],fontsize=14)
			dosRange=max([minDOS,maxDOS])
			pylab.xlim(-1.1*dosRange,1.1*dosRange) # scales DOS to larges value of DOS in the given energy range 
			pylab.axvline(0.0, color = 'k', linewidth = 1.3) #line separating up and down spin
		else:
			ax2.legend(handles[::-1], labels[::-1],fontsize=14)
			pylab.xlim(0,1.1*maxDOS) # scales DOS to larges value of DOS in the given energy range

		pylab.yticks(numpy.arange(yLim[0],yLim[1]+1,2))

		ax2.spines['bottom'].set_linewidth(1.5)
		ax2.spines['left'].set_linewidth(1.5)
		ax2.spines['right'].set_linewidth(1.5)
		ax2.spines['top'].set_linewidth(1.5)

		ax2.set_yticklabels([])	     #to hide ticklabels
		ax1.set_position([0.07,0.1,0.67,0.8]) #[left,bottom,width,height] 
		ax2.set_position([0.75,0.1,0.20,0.8]) #other useful options for the frame! :D

		ax2.yaxis.set_ticks([])
		ax2.yaxis.set_ticks_position('left')
		pylab.xlabel('Density of States (States/eV)')
		ax2.axes.xaxis.set_label_position('top')
		locs, labels = pylab.xticks()
		
		for item in range(len(labels)):
			if item == len(labels)/2:
				labels[item]='arbitrary units'
			else:
				labels[item]=''
		
		ax2.set_xticklabels(labels,fontsize = 20)

##########################################################################################################
	     #to plot the DOS
	if DOSPlot == 'DOS':
		print 'Plotting electronic band structure and DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))
		logging.info('Plotting electronic band structure and DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
		fileplot = os.path.join(subdir,'BANDDOS_%s%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),oneCalc['_AFLOWPI_PREFIX_'],postfix))
		ax2=pylab.subplot(122)

		try:
			data = open(filedos,'r').readlines()
		except Exception:
			logging.warning("output from dos calculation not found. Are you sure you ran ppDOS and it completed properly?")
			print "Are you sure you ran ppDOS and it completed properly?"
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

	ax1.yaxis.set_ticks(numpy.arange(low_tick_bound,high_tick_bound))

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
	ax1.set_title(figtitle)
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
import cPickle
import collections
import math

def __smoothGauss(list,strippedXs=False,degree=10):  

	'''


	Arguments:


	Keyword Arguments:

        
        Returns:


	'''

	window=degree*2-1  
	weight=numpy.array([1.0]*window)  
	weightGauss=[]  
	for i in range(window):  
		i=i-degree+1  
		frac=i/float(window)  
		gauss=1/(numpy.exp((4*(frac))**2))  
		weightGauss.append(gauss)  
	weight=numpy.array(weightGauss)*weight  
	smoothed=[0.0]*(len(list)-window)  
	for i in range(len(smoothed)):  
		smoothed[i]=sum(numpy.array(list[i:i+window])*weight)/sum(weight)  
	return smoothed


def __dosPlot(oneCalc,ID,yLim=[-10,10],LSDA=False,postfix=''):
	'''
	Function to take the data generated for the DOS and plot them with energy shifted relative to the Fermi Energy.                       

	Arguments:
              oneCalc (dict): Single calculation that is the value of the dictionary of dictionaries of calculations
              ID (str): ID of the calculation
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)

        Returns:
              None

	'''


	

	'''extracts HOMO from nscf calculation output file as input to the plotting'''
	print 'Plotting DOS'
	try:
		Efermi=AFLOWpi.retr._getEfermi(oneCalc,ID)
		if type(Efermi)!=type(0.5):
			LSDA=True
			
	except:
		Efermi=0.0
	subdir=oneCalc['_AFLOWPI_FOLDER_']

	"""get the path to the subdirectory of the calc that you are making plots for"""
	filedos = os.path.join(subdir,'dos_%s.out'%ID)


        if postfix!='':
            postfix='_'+postfix
	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
	fileplot = os.path.join(subdir,'DOS_%s_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,postfix))

	#to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=20)             #to set the font size

	"""get the path to the subdirectory of the calc that you are making plots for"""
        filedos = os.path.join(subdir,'%s_dos.dat'%ID)
        
        width = 20
	height = 14
	pylab.figure(figsize=(width, height)) #to adjust the figure size

     ####################################
     #to plot the DOS
	try:
		data = open(filedos,'r').readlines()
	except Exception:
		logging.warning("output from dos calculation not found. Are you sure you ran ppDOS and it completed properly?")
		print 'Are you sure that you ran ppDos and that it completed properly?'
		return
		
	en = []
	enup = []
	endown = []
	dos = []
	dosdw = []
	for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
		try:
			if LSDA==True:
				enup.append(float(data[i].split()[0])-Efermi[0])
				endown.append(float(data[i].split()[0])-Efermi[1])
				
			else:
				en.append(float(data[i].split()[0])-Efermi) #to shift all the y values with respect to the Fermi level
			
			dos.append(float(data[i].split()[1]))
			dosdw.append(-1*float(data[i].split()[2]))
		
		except Exception, e:
			pass
	if LSDA==True:
		enup  = map(float,enup)
		endown= map(float,endown)
	else:
		endos=map(float,en)  #to convert the list x from float to numbers
	floatdos=map(float,dos)
	floatdosDOWN=map(float,dosdw)
	enshift = numpy.array(endos) #to treat the list b as an array?
  

	ax2=pylab.subplot(111)
	inDict=AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
	if "degauss" in inDict["&system"].keys():
		floatdos,floatdosDOWN,enshift=floatdos,floatdosDOWN,enshift
	else:
		floatdos,floatdosDOWN,enshift=__smoothGauss(floatdos),__smoothGauss(floatdosDOWN),__smoothGauss(enshift)

#        plot(f,enshift,'k') #to plot the original data
       	pylab.plot(enshift,floatdos,'k') #to plot the smoothed dat
	dosMAX = 1.0*max([dos[k] for k in range(len(dos)) if en[k] < yLim[1] and en[k] > yLim[0]])
	if not LSDA:
		pylab.plot(enshift,floatdos,'k-') #to plot the smoothed data
		pylab.ylim(0,dosMAX) # scales DOS to larges value of DOS in the given energy range
		
	else:
		dosMIN = 1.0*min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if en[k] < yLim[1] and en[k] > yLim[0]])
		enshiftup  = numpy.array(enup)
		enshiftdown= numpy.array(endown)
		floatdosup=map(float,dos)
		pylab.plot(enshiftdown,floatdosDOWN,'k-')
		pylab.plot(enshiftup,floatdosUP,'k-')

		pylab.ylim(dosMIN,dosMAX) # scales DOS to larges value of DOS in the given energy range
		pylab.axhline(0.0, color = 'k', linewidth = 1.3) #line separating up and down spin
#	else:
#		pylab.ylim(0,dosMAX) # scales DOS to larges value of DOS in the given energy range

	pylab.xlim(yLim[0],yLim[1])
	pylab.axvline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line
	pylab.xticks(numpy.arange(yLim[0],yLim[1]+1,2))
     ##############################
     #to increase the linewidth of the axis
	ax2.spines['bottom'].set_linewidth(1.5)
	ax2.spines['left'].set_linewidth(1.5)
	ax2.spines['right'].set_linewidth(1.5)
	ax2.spines['top'].set_linewidth(1.5)


     #other useful options for the frame! :D
	ax2.yaxis.set_ticks_position('left')
	pylab.ylabel('Density of States (States/eV)')
	ax2.axes.yaxis.set_label_position('right')	
	pylab.xlabel('E(eV)')

	figtitle = ''
	compoundName = AFLOWpi.retr._getStoicName(oneCalc,strip=True)
        compoundNameLatex = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)

     #to set a centerd title, with a bigger font
	figtitle = r'Density of States: %s' % (compoundNameLatex,)
	t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=14,horizontalalignment='center') #[x,y]

	matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')
	try:
			AFLOWpi.retr._moveToSavedir(fileplot)
	except Exception,e:
		pass

	AFLOWpi.retr._moveToSavedir(fileplot)
        pyplot.cla()
	pyplot.clf()
	pyplot.close()


def __combinePDOS(dosfiles):
	mat=[]  # matrix with total sum of ldos
	for i in range(len(dosfiles)):
		mati=[] # temporal matrix for each DOS file "i"
		k=0
		kresolved = False

		if not kresolved:
			with open(dosfiles[i],'r') as pDOSFile:
				pDOSString = pDOSFile.readlines()
				for line in pDOSString:
					'''projwfc.x outputs ####### when the value of the kpoint is less than -100 so we skip those'''

					try:
                                            if len(line.strip())!=0:
                                                mati.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])	  
			
					except:
                                            try:
                                                if len(line.strip())!=0:
                                                    mati.append([float(line.split()[0]),float(line.split()[1])])
                                            except:
                                                pass


						
				if mat == []: # if it is the first dos file, copy total matrix (mat) = the first dos files's data
					mat=mati[:]
				else:
					for j in range(len(mati)): # if it is not the first file, sum values
						try:
							mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2]]  
						except Exception,e:
                                                    try:
                                                        mat[j]=[mat[j][0],mat[j][1]+mati[j][1]]
                                                    except Exception,e:
                                                        AFLOWpi.run._fancy_error_log(e)
							

##############################################################################################################################################

		if kresolved:
			logging.warning('k resolved not supported')
			print 'k resolved not supported'
			raise SystemExit
			with open(dosfiles[i],'r') as pDOSFile:
				pDOSString = pDOSFile.readlines()
				for line in pDOSString:

					'''projwfc.x outputs ####### when the value of the kpoint is less than -100 so we skip those'''
					if len(line) > 10 and len(re.findall('[*]+',line)) == 0 and line.split()[0] != "#":

						ik = line.split()[0]


						ik = int(ik)


						if ik > k:  #if it is a different k block

							k=int(line.split()[0])
							oldmat=[] # temporal matrix for each k-point


						if ik == 1:
							mati.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]) # append: energy, ldosup, ldosdw
						elif ik == k and k > 1:
							oldmat.append([float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])
						elif len(line) < 5 and k > 1:  #if blank line, sum k-frame to the total
							for j in range(len(oldmat)):  
								mati[j]=[mati[j][0],mati[j][1]+oldmat[j][1],mati[j][2]+oldmat[j][2]]

				if mat == []: # if it is the first dos file, copy total matrix (mat) = the first dos files's data
					mat=mati[:]
				else:
					for j in range(len(mati)): # if it is not the first file, sum values
						mat[j]=[mat[j][0],mat[j][1]+mati[j][1],mat[j][2]+mati[j][2]]  


	return mat





def __sumpdos(oneCalc,ID,TB=False):
	'''
	Takes the output files from projwfx.x that is called in the ppDOS function and sums 
	the projected orbitals across the files of each orbital for each atomic species
	and outputs the summed data in files named <species>_<orbital>.sumpdos.


	Arguments:
	      oneCalc (dict): dictionary of one calculation
              ID (str): ID of the calculation

        Keyword Arguments:
              None

        Returns:
              None

	'''

	subdir=oneCalc['_AFLOWPI_FOLDER_']
	

	'''get a list of all the pdos files to be combined then plotted'''
	atomList = []
	pDOSDict = {}
	for k,v in oneCalc.iteritems():
		if re.match(r'_AFLOWPI_[A-Z][0-9]*_',k):
			atomList.append(v)

	
	'''
	just the possible names of orbitals that projwfc.x will output
	skips over orbitals for atomic species that are not there
	'''
	orbitalList = ['s','p','d','f','All']

	'''
	globs for the file names that match the expression to sort out files for
	the different atomic species and each of their orbitals

	it then gives a list of the files for each orbital for each species in the 
	calculation and saves them as a dictionary with:
	the atom species as key
	a list of dictionaries of each atomic orbital for the species as the value
	those dictionaries have the orbital as key 
	and the value is a list of the files obtained by glob
	'''

        if TB==True:
            glob_ID =  AFLOWpi.prep._return_ID(oneCalc,ID,step_type='PAO-TB',last=True,straight=False)
            glob_ID +='_TB'
        else:
            glob_ID =  AFLOWpi.prep._return_ID(oneCalc,ID,step_type='dos',last=True,straight=False)

	byAtomDict={}
	for atom in atomList:
		byAtom=[]
		for orbital in orbitalList:
			pDOSFiles= glob.glob(os.path.join(subdir,'%s.pdos_atm*(%s)*(%s)' % (glob_ID,atom,orbital)))

			if len(pDOSFiles):
				byAtom.append({'%s' % orbital:pDOSFiles})
		byAtomDict[atom] = byAtom

	for atom,orbital in byAtomDict.iteritems():			
			for orbitalDict in range(len(orbital)):
				for orbitalName,fileList in orbital[orbitalDict].iteritems():
					data = __combinePDOS(fileList)
					with open(os.path.join(subdir,'%s_%s.sumpdos' % (atom,orbitalName)),'wb') as outputFile:
						cPickle.dump(data,outputFile)


	byAtom={}
	for atom in atomList:	
		pDOSFiles= glob.glob(os.path.join(subdir,'%s.pdos_atm*(%s)_wfc*' % (glob_ID,atom)))
		if len(pDOSFiles):
			pDOSDict['%s_All'% (atom)] =  pDOSFiles


	'''Cycle through all of the atoms in a single calculation and sum over the pdos files
	in the calculation directory and saves  '''

	for atom,files in pDOSDict.iteritems():		
		data = __combinePDOS(files)
		with open(os.path.join(subdir,'%s.sumpdos' % (atom)),'wb') as outputFile:
					cPickle.dump(data,outputFile)

###################################################################################################################






def __plotByAtom(maxNum,speciesNum,fig,atom,oneCalc,ID,yLim=[-10,10],LSDA=False,ax=None,TB=False):
	"""
	Function to take the data files generated by the sumpdos function and
        plots them with energy shifted relative to the Fermi Energy.
	
	Arguments:
	      atom (str): Species of atom you are wanting to plot ('All' will
                          plot the combined pdos for all atoms in the system)
              oneCalc (dict): single calculation that is the value of the dictionary
                              of dictionaries of calculations
	
	Keyword Arguments:
    	      yLim (list): the limits of the ordinate on the plot (default: [-10,10]

        Returns:
              ax2 (matplotlib.pyplot.axis): returns an axes object with the plotted proj. DOS
              
	"""
	try:
            if TB==True:
		    fermi_ID = '%s_WanT_dos'%calcID
		    print fermi_ID
		    Efermi=AFLOWpi.retr._getEfermi(oneCalc,fermi_ID,directID=True)

            else:
		    Efermi=AFLOWpi.retr._getEfermi(oneCalc,ID)
		    if type(Efermi)!=type(0.5):
			    Efermi=Efermi[0]

	except Exception,e:

		Efermi=0.0
        try:
            ax2=ax[speciesNum]
        except:
            ax2=ax

	'''extracts HOMO from nscf calculation output file as input to the plotting'''

	"""get the path to the subdirectory of the calc that you are making plots for"""	
	subdir = oneCalc['_AFLOWPI_FOLDER_']
	
	try:
		filePlotName = 'PDOS_%s%s_%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atom)
	except IOError:
                filePlotName = 'PDOS_%s%s_%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atom)

	fileplot = os.path.join(subdir,filePlotName)

	#to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=9)             #to set the font size

        width = 20
	height = 14.0

     #to set a centerd title with a bigger font
     ####################################
     #to plot the DOS


	def getPlotData(sumpdosFile):
		with open(sumpdosFile,'rb') as dataFile:
			data = cPickle.load(dataFile)

		en = []
		pdos = []
		ldos = []
		ldosDOWN = []

		for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
			try:
                                one_en = float(data[i][0])-Efermi
                                if one_en > yLim[0]*1.05 and one_en < yLim[1]*1.05:
                                    en.append(float(data[i][0])-Efermi) #to shift all the y values with respect to the Fermi level
                                    ldos.append(float(data[i][1]))
                                    pdos.append(float(data[i][2]))
                                    ldosDOWN.append(-1*float(data[i][2]))
			except Exception, e:
				pass


		endos=map(float,en)  #to convert the list x from float to numbers
		floatdos=map(float,ldos)
		floatdosDOWN=map(float,ldosDOWN)
		enshift = numpy.array(endos) #to treat the list b as an array?
		inDict=AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])

		if "degauss" in inDict["&system"].keys():
			return enshift,floatdos,floatdosDOWN
		else:
			return __smoothGauss(enshift),__smoothGauss(floatdos),__smoothGauss(floatdosDOWN)



	maxDOS=0
	minDOS=0
	orbitalList = ['s','p','d','f',"All"]
	pDOSFiles= glob.glob(os.path.join(subdir,'%s_*.sumpdos' % (atom)))

	pDOSNoPath = []
	
	for item in reversed(orbitalList):
		if os.path.exists(os.path.join(subdir,'%s_%s.sumpdos' % (atom,item))):
			pDOSNoPath.append(os.path.join(subdir,'%s_%s.sumpdos' % (atom,item)))

	color = 'k'
	for filepath in pDOSNoPath:
		filename =  filepath.split('/')[-1]
		orbitalName = filename.split('_')[-1].split('.')[0]
		species =  filename.split('_')[0].split('.')[0]
			
		if orbitalName != 'All':
			#print 'Plotting %s orbital of %s for %s' % (orbitalName,atom,__getStoicName(oneCalc))
			logging.info('Plotting %s orbital of %s for %s' % (orbitalName,atom,AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
		else:
			#print 'Plotting %s orbital of  %s' % (orbitalName,__getStoicName(oneCalc))
			logging.info('Plotting %s orbital of %s' % (orbitalName,AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
		if orbitalName == 's':
			color = 'g'
		elif orbitalName == 'p':
			color = 'b'
		elif orbitalName == 'd':
			color = 'r'
		elif orbitalName == 'f':
			color = 'c'
		elif orbitalName == 'All':
			color = 'k'
		'''gets the energy and the DOS for the orbital for the atom'''
		enshift, floatdos,floatdosDOWN = getPlotData(filepath)

		''' scales DOS to larges value of DOS in the given energy 
		range and finds the largest DOS between the different orbitals'''
		if max([floatdos[k] for k in range(len(floatdos)) if (enshift[k] < yLim[1] and enshift[k] > yLim[0])]) > maxDOS:
                    maxDOS = max([floatdos[k] for k in range(len(floatdos)) if (enshift[k] < yLim[1] and enshift[k] > yLim[0] )])

		try:
			if min([floatdosDOWN[k] for k in range(len(floatdos)) if (enshift[k]<yLim[1] and enshift[k]>yLim[0])])<minDOS:
                            minDOS = min([floatdosDOWN[k] for k in range(len(floatdos)) if (enshift[k] < yLim[1] and enshift[k] > yLim[0] )])
		except:
			minDOS=0
		
		if not LSDA:
			minDOS=0
			
#		if orbitalName != 'All':

		ax2.plot(enshift,floatdos,color+'-',label=orbitalName)
		if LSDA:
			ax2.plot(enshift,floatdosDOWN,color+'-',label=orbitalName)
	
		       	
	handles, labels = ax2.get_legend_handles_labels()

	if LSDA:
		ax2.legend(handles[::-2], labels[::-2])
		pylab.ylim(minDOS,maxDOS) # scales DOS to larges value of DOS in the given energy range
		pylab.axhline(0.0, color = 'k', linewidth = 1.3) #line to separate up and down spin
	else:
		ax2.legend(handles[::-1], labels[::-1])
		pylab.ylim(0,maxDOS) # scales DOS to larges value of DOS in the given energy range




	try:
		max_x =  max(enshift)
	except:
		pass
	try:
		min_x =  min(enshift)
	except:
		pass

	pylab.xlim(yLim[0],yLim[1])

        #plot the ticks only on the bottom plot of the figure
        ax2.set_xticks(numpy.arange(yLim[0],yLim[1]+1,2))


	ax2.axvline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line

     ##############################
     #to increase the linewidth of the axis
	ax2.spines['bottom'].set_linewidth(1.5)
	ax2.spines['left'].set_linewidth(1.5)
	ax2.spines['right'].set_linewidth(1.5)
	ax2.spines['top'].set_linewidth(1.5)


     #other useful options for the frame! :D
	ax2.yaxis.set_ticks_position('left')

	for i in ax2.get_yticklabels():
		i.set_visible(False)
	try:
		speciesYPos = yLim[0]+(abs(yLim[0]-yLim[1])/50.0)
	except:
		speciesYPos = yLim[0]+(abs(yLim[0]-yLim[1])/50.0)

	ax2.axes.yaxis.set_label_position('right')	

	figtitle=''
	compoundName = AFLOWpi.retr._getStoicName(oneCalc,strip=True)
        compoundNameLatex = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)
        #to set a centerd title, with a bigger font
	figtitle = r'Orbital Projected DOS: %s' % (compoundNameLatex)

	try:
		return ax2
	except:
		return ''



def dos(calcs,yLim=[-10,10],runlocal=False,postfix=''):
	'''
	Generates DOS plots for the calculations in the dictionary of dictionaries of calculations

	Arguments:
              calcs (dict): dictionary of dictionaries representing the set of calculations
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)
              runlocal (bool): Do the plotting right now or if False do it when the calculations are running
              postfix (str): Postfix to the filename of the plot
              tight_banding (bool): Whether to treat the input data as from Quantum Espresso or WanT bands.x

        Returns:
              None
              
	'''

	if runlocal:
		for ID,oneCalc in calcs.iteritems():
			__dosPlot(oneCalc,ID,yLim,postfix=postfix)
	else:
		for ID,oneCalc in calcs.iteritems():
			AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__dosPlot(oneCalc,ID,[%s,%s],postfix='%s')" % (yLim[0],yLim[1],postfix))




def __plotByAtom(maxNum,speciesNum,fig,atom,oneCalc,ID,yLim=[-10,10],LSDA=False,ax=None,TB=False):
	"""
	Function to take the data files generated by the sumpdos function and
        plots them with energy shifted relative to the Fermi Energy.
	
	Arguments:
	      atom (str): Species of atom you are wanting to plot ('All' will
                          plot the combined pdos for all atoms in the system)
              oneCalc (dict): single calculation that is the value of the dictionary
                              of dictionaries of calculations
	
	Keyword Arguments:
    	      yLim (list): the limits of the ordinate on the plot (default: [-10,10]

        Returns:
              ax2 (matplotlib.pyplot.axis): returns an axes object with the plotted proj. DOS
              
	"""
	try:
            if TB==True:
		    fermi_ID = '%s_WanT_dos'%calcID
		    print fermi_ID
		    Efermi=AFLOWpi.retr._getEfermi(oneCalc,fermi_ID,directID=True)
            else:
		    Efermi=AFLOWpi.retr._getEfermi(oneCalc,ID)
		    if type(Efermi)!=type(0.5):
			    Efermi=Efermi[0]
	except:
		Efermi=0.0

	print 'Efermi/HOMO-LUMO = %s' % Efermi

        try:
            ax2=ax[speciesNum]
        except:
            ax2=ax

	'''extracts HOMO from nscf calculation output file as input to the plotting'''

	"""get the path to the subdirectory of the calc that you are making plots for"""	
	subdir = oneCalc['_AFLOWPI_FOLDER_']
	
	try:
		filePlotName = 'PDOS_%s%s_%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atom)
	except IOError:
                filePlotName = 'PDOS_%s%s_%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atom)

	fileplot = os.path.join(subdir,filePlotName)

	#to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=9)             #to set the font size

        width = 20
	height = 14.0

     #to set a centerd title with a bigger font
     ####################################
     #to plot the DOS


	def getPlotData(sumpdosFile):
		with open(sumpdosFile,'rb') as dataFile:
			data = cPickle.load(dataFile)

		en = []
		pdos = []
		ldos = []
		ldosDOWN = []

		for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
			try:
                                one_en = float(data[i][0])-Efermi
                                if one_en > yLim[0]*1.05 and one_en < yLim[1]*1.05:
                                    en.append(float(data[i][0])-Efermi) #to shift all the y values with respect to the Fermi level
                                    ldos.append(float(data[i][1]))
                                    pdos.append(float(data[i][2]))
                                    ldosDOWN.append(-1*float(data[i][2]))
			except Exception, e:
				pass


		endos=map(float,en)  #to convert the list x from float to numbers
		floatdos=map(float,ldos)
		floatdosDOWN=map(float,ldosDOWN)
		enshift = numpy.array(endos) #to treat the list b as an array?


		inDict=AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
		if "degauss" in inDict["&system"].keys() or TB==True:

			return enshift,floatdos,floatdosDOWN
		else:
			return __smoothGauss(enshift),__smoothGauss(floatdos),__smoothGauss(floatdosDOWN)



	maxDOS=0
	minDOS=0
	orbitalList = ['s','p','d','f',"All"]
	pDOSFiles= glob.glob(os.path.join(subdir,'%s_*.sumpdos' % (atom)))

	pDOSNoPath = []
	
	for item in reversed(orbitalList):
		if os.path.exists(os.path.join(subdir,'%s_%s.sumpdos' % (atom,item))):
			pDOSNoPath.append(os.path.join(subdir,'%s_%s.sumpdos' % (atom,item)))

	color = 'k'
	for filepath in pDOSNoPath:
		filename =  filepath.split('/')[-1]
		orbitalName = filename.split('_')[-1].split('.')[0]
		species =  filename.split('_')[0].split('.')[0]
			
		if orbitalName != 'All':
			#print 'Plotting %s orbital of %s for %s' % (orbitalName,atom,__getStoicName(oneCalc))
			logging.info('Plotting %s orbital of %s for %s' % (orbitalName,atom,AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
		else:
			#print 'Plotting %s orbital of  %s' % (orbitalName,__getStoicName(oneCalc))
			logging.info('Plotting %s orbital of %s' % (orbitalName,AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
		if orbitalName == 's':
			color = 'g'
		elif orbitalName == 'p':
			color = 'b'
		elif orbitalName == 'd':
			color = 'r'
		elif orbitalName == 'f':
			color = 'c'
		elif orbitalName == 'All':
			color = 'k'
		'''gets the energy and the DOS for the orbital for the atom'''
		enshift, floatdos,floatdosDOWN = getPlotData(filepath)

		''' scales DOS to larges value of DOS in the given energy range and finds the largest DOS between the different orbitals'''
		if max([floatdos[k] for k in range(len(floatdos)) if (enshift[k] < yLim[1] and enshift[k] > yLim[0])]) > maxDOS:
                    maxDOS = max([floatdos[k] for k in range(len(floatdos)) if (enshift[k] < yLim[1] and enshift[k] > yLim[0] )])

		try:
			if min([floatdosDOWN[k] for k in range(len(floatdos)) if (enshift[k]<yLim[1] and enshift[k]>yLim[0])])<minDOS:
                            minDOS = min([floatdosDOWN[k] for k in range(len(floatdos)) if (enshift[k] < yLim[1] and enshift[k] > yLim[0] )])
		except:
			minDOS=0
		
		if not LSDA:
			minDOS=0
			
#		if orbitalName != 'All':
		ax2.plot(enshift,floatdos,color+'-',label=orbitalName)

		if LSDA:
			ax2.plot(enshift,floatdosDOWN,color+'-',label=orbitalName)
	
		       	
	handles, labels = ax2.get_legend_handles_labels()

	if LSDA:
		ax2.legend(handles[::-2], labels[::-2])
		pylab.ylim(minDOS,maxDOS) # scales DOS to larges value of DOS in the given energy range
		pylab.axhline(0.0, color = 'k', linewidth = 1.3) #line to separate up and down spin
	else:
		ax2.legend(handles[::-1], labels[::-1])
		pylab.ylim(0,maxDOS) # scales DOS to larges value of DOS in the given energy range

	try:
		max_x =  max(enshift)
	except:
		pass
	try:
		min_x =  min(enshift)
	except:
		pass
	pylab.xlim(yLim[0],yLim[1])

        #plot the ticks only on the bottom plot of the figure

        ax2.set_xticks(numpy.arange(yLim[0],yLim[1]+1,2))


	ax2.axvline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line

     ##############################
     #to increase the linewidth of the axis
	ax2.spines['bottom'].set_linewidth(1.5)
	ax2.spines['left'].set_linewidth(1.5)
	ax2.spines['right'].set_linewidth(1.5)
	ax2.spines['top'].set_linewidth(1.5)


     #other useful options for the frame! :D
	ax2.yaxis.set_ticks_position('left')

	for i in ax2.get_yticklabels():
		i.set_visible(False)
	try:
		speciesYPos = yLim[0]+(abs(yLim[0]-yLim[1])/50.0)
	except:
		speciesYPos = yLim[0]+(abs(yLim[0]-yLim[1])/50.0)

	ax2.axes.yaxis.set_label_position('right')	

	figtitle=''
	compoundName = AFLOWpi.retr._getStoicName(oneCalc,strip=True)
        compoundNameLatex = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)
        #to set a centerd title, with a bigger font
	figtitle = r'Orbital Projected DOS: %s' % (compoundNameLatex)

	try:
		return ax2
	except:
		return ''


def opdos(calcs,yLim=[-10,10],runlocal=False,postfix='',scale=False,tight_binding=False):
	'''
	Generates electronic band structure plots for the calculations in the dictionary of dictionaries
	of calculations with the option to have a DOS or PDOS plot to accompany it.

	Arguments:
              calcs (dict): dictionary of dictionaries representing the set of calculations
              
	Keyword Arguments:
              yLim (list): List or tuple of two integers for max and min range in horizontal axis of DOS plot
              LSDA (bool): To plot DOS as spin polarized or not (calculation must have been done as spin polarized)
              runlocal (bool): Do the plotting right now or if False do it when the calculations are running
              postfix (str): Postfix to the filename of the plot
              tight_banding (bool): Whether to treat the input data as from Quantum Espresso or WanT bands.x

        Returns:
              None

	'''

	for ID,oneCalc in calcs.iteritems():	
		if runlocal:
			AFLOWpi.plot.__opdos(oneCalc,ID,yLim,postfix=postfix,scale=scale,tight_binding=tight_binding)
		else:
			AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__opdos(oneCalc,ID,[%s,%s],postfix='%s',scale=%s,tight_binding=%s)" % (yLim[0],yLim[1],postfix,scale,tight_binding))



def __opdos(oneCalc,ID,yLim,postfix='',scale=False,tight_binding=False):

	logging.info('Plotting Orbital Projected DOS')
	logging.info('summing pdos for %s' % oneCalc['_AFLOWPI_FOLDER_'].split('/')[-1])

        
        LSDA=False
        nspin = int(AFLOWpi.scfuj.chkSpinCalc(oneCalc,ID))
        if nspin!=1:
            LSDA=True



	__sumpdos(oneCalc,ID,TB=tight_binding)

	if postfix!='':
		postfix='_'+postfix

        atomList=[]
	atomPlotArray = collections.OrderedDict()
	perSpecies = collections.OrderedDict()
        for k,v in oneCalc.iteritems():
                if re.match(r'_AFLOWPI_[A-Z][0-9]*_',k):
                        atomList.append(v)

	atomList = sorted(list(set(atomList)))


	#to set figure size and default fonts
	matplotlib.rc("font", family="serif")      #to set the font type
	matplotlib.rc("font", size=9)             #to set the font size


        width = 20
	height = 14

        ##sometimes smearing will give negative pdos vals
        ##and we dont want that to show if not lsda
        figtitle='Orbital Proj. DOS: %s'%(AFLOWpi.retr._getStoicName(oneCalc,strip=True))
        if not LSDA:
            pyplot.gca().set_ylim(bottom=0)

	fig, ax_list = pylab.subplots(len(atomList),figsize=(14.0,math.floor(3.5*len(atomList))))

	atomAxisArray = collections.OrderedDict()
	frame1 = pyplot.gca()

	old_high=0.0
	old_low =0.0
	t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=20,horizontalalignment='center',) #[x,y]

	for species in range(len(atomList)):
		ax = __plotByAtom(len(atomList),species,fig,atomList[species],oneCalc,ID,yLim,LSDA=LSDA,ax=ax_list,TB=tight_binding)

		ax.axes.set_ylabel('Density of States (States/eV)')
		ax.axes.set_title(atomList[species],weight='bold',loc='left',fontsize=18)

		ax.axes.set_autoscalex_on(False)
		
		ylim=ax.axes.get_ylim()
		new_low=ylim[0]
		if new_low<old_low:
			old_low=new_low
		new_high=ylim[1]
		if new_high>old_high:
			old_high=new_high

                if not LSDA:
                    old_low=0.0


                if species!=len(atomList)-1:
                    ax.axes.set_xticklabels([])

                try:
                    ax_list[species]=ax
                except:
                    ax_list=[ax]

		atomAxisArray.update({species:ax})

        for species,ax in atomAxisArray.iteritems():
            for ax in fig.get_axes():

                ax.axes.set_ylim([1.05*old_low,1.05*old_high])
		ax.axes.set_xlim(yLim[0],yLim[1])

	subdir=oneCalc['_AFLOWPI_FOLDER_']

	"""get the path to the subdirectory of the calc that you are making plots for"""

	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
	
	fileplot = os.path.join(subdir,'PDOS_%s_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,postfix))

	fig.savefig(fileplot,bbox_inches='tight')
	try:
		AFLOWpi.retr._moveToSavedir(fileplot)
	except:
		pass

        pyplot.cla()
	pyplot.clf()
	pyplot.close()

