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
                                            'yl':r'$\sigma/\tau$ $(10^{21}$ $m/s/\Omega)$ ',
                                            'fp':'CONDUCTION',

                                            }
        trans_plot_dict['seebeck']       = {'pf':'seebeck_',
                                            'ft':'$Seebeck$:',
                                            'lc':'Seebeck ',
                                            'xl':'$\mu$ (eV)',
                                            'yl':'S $(10^{-3} $ $V/K)$',
                                            'fp':'SEEBECK',
                                            }

        # trans_plot_dict['sig_seebeck']   = {'pf':'sigma_seebeck_',
        #                                     'ft':'$\sigma S$:',
        #                                     'lc':'\sigma S ',
        #                                     'xl':'$\mu$ (eV)',
        #                                     'yl':'$\sigma S$ ($Vm/\Omega/K)$ ',
        #                                     'fp':'SIGMA_SEEBECK',

        #                                     }
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
		type_list=['kappa','seebeck','cond']

        for Type in type_list:
            try:
                max_val=0.0
                min_val=0.0


                file_name=[]
                file_name_up=[]
                file_name_down=[]
                for ID_ent in ID_list:

                    search=oneCalc['_AFLOWPI_FOLDER_']+'/%s_PAOpy_%s*.dat'%(ID_ent,trans_plot_dict[Type]['pf']) 

                    file_name.extend(glob.glob(search))
		    search = oneCalc['_AFLOWPI_FOLDER_']+'/%s_PAOpy_%s*up*.dat'%(ID_ent,trans_plot_dict[Type]['pf']) 
                    file_name_up.extend(glob.glob(search))
		    search = oneCalc['_AFLOWPI_FOLDER_']+'/%s_PAOpy_%s*down*.dat'%(ID_ent,trans_plot_dict[Type]['pf']) 

                    file_name_down.extend(glob.glob(search))

		

                spin_polarized=False
                if len(file_name_down)!=0:
                        spin_polarized=True
		
		if len(file_name_down) == 0 and len(file_name) == 0:
			continue

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
				    y_vals=numpy.asarray(y_vals)

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
				    y_vals=numpy.asarray(y_vals)


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

#			    x_vals,y_vals = AFLOWpi.plot.__smoothGauss(x_vals,degree=5),AFLOWpi.plot.__smoothGauss(y_vals,degree=5)

			    print y_vals
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

#			    if epsilon==True:
#				    if min_x<=0.2:
#					    min_x=0.2

#			    if Type=='kappa':
#				    y_vals=numpy.asarray(y_vals)/1.0e17
#			    elif Type=='sigma_seebeck':
#				    y_vals=numpy.asarray(y_vals)/1.0e-3
#			    elif Type=='cond':
#				    y_vals=numpy.asarray(y_vals)/1.0e20



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
                    postfix_append='_'+postfix
		else:
                    postfix_append=''

                fileplot = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_%s_%s%s.pdf'%(trans_plot_dict[Type]['fp'],compoundName,ID,postfix_append))
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
	print ep_data_file
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
			val_array.append(sum([float_array[1],float_array[2],float_array[3]])/3.0)
		except:
			pass


		en_array=[mult_x*x for x in en_array]
		val_array=[mult_y*y for y in val_array]

		
	return en_array,val_array	




###############################################################################################################################
