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
def transport_plots(calcs,runlocal=False,postfix=''):


	if runlocal:
		for ID,oneCalc in calcs.iteritems():
                    AFLOWpi.plot.__transport_plot(oneCalc,ID,postfix=postfix)
	else:
		for ID,oneCalc in calcs.iteritems():
			AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__transport_plot(oneCalc,ID,postfix='%s')" %postfix)



def __transport_plot(oneCalc,ID,nm=False,postfix=''):
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
                                            'lc':'Cond. ',
                                            'xl':'$\mu$ (eV)',
                                            'yl':'$\sigma$ $(10^{20}$ $m/\Omega)$ ',
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
        trans_plot_dict['epsilon_i']     = {'pf':'epsilon_imag_',
                                            'ft':'$\epsilon_{imaginary}$:',
                                            'lc':'\epsilon_{i} ',
                                            'xl':'$\hbar\omega$ (eV)',
                                            'yl':'$\epsilon_{imag}$ (au)',
                                            'fp':'EPSILON_IMAG',

                                            }
        trans_plot_dict['epsilon_r']     = {'pf':'epsilon_real_',
                                            'ft':'$\epsilon_{real}$:',
                                            'lc':'\epsilon_{r} ',
                                            'xl':'$\hbar\omega$ (eV)',
                                            'yl':'$\epsilon_{real}$ (au)',
                                            'fp':'EPSILON_REAL',

                                            }
        


####################################################################################################################
        ##needed to make sure the real and imaginary parts of epsilon have the same scaling
        min_val=0.0
        max_val=0.0




        file_name = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_epsilon*.dat'%(ID, ))




        try:                            
            for i in range(len(file_name)):
                en,val = read_transport_datafile(file_name[i])

                set_max = max(val)
                set_min = min(val)
                
                if set_max>max_val:
                    max_val=set_max
                if set_min<min_val:
                    min_val=set_min

        except Exception,e: 
            print e
        
        max_val_ep =  max_val
        min_val_ep =  min_val

        ID_list =  AFLOWpi.prep._return_ID(oneCalc,ID,step_type='transport',last=True,straight=True)

####################################################################################################################

        type_list=['kappa','seebeck','sig_seebeck','cond','ZT','epsilon_i','epsilon_r']
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

                    file_name_up.extend(glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_%sup_*.dat'%(ID_ent,trans_plot_dict[Type]['pf']) ))
                    file_name_down.extend(glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/%s_WanT_%sdown_*.dat'%(ID_ent,trans_plot_dict[Type]['pf']) ))



                spin_polarized=False
                if len(file_name_down)!=0:
                        spin_polarized=True


                if spin_polarized==False:
                        ############################################################################################
                        ##Non Spin Polarized Case
                        ############################################################################################
                        data_by_temp={}
                        for i in range(len(file_name)):
                            temperature = file_name[i][:-4].split('_')[-1][:-1]
                            data_by_temp[temperature] = read_transport_datafile(file_name[i])                 



                if spin_polarized==True:
                        #############################################################################################
                        ##Spin Polarized Case
                        #############################################################################################
                        data_by_temp_dn={}
                        data_by_temp_up={}
                        for i in range(len(file_name_down)):
                            temperature = file_name_down[i][:-4].split('_')[-1][:-1]
                            data_by_temp_dn[temperature] = read_transport_datafile(file_name_down[i]) 
                            data_by_temp_up[temperature] = read_transport_datafile(file_name_up[i]) 

                #############################################################################################

                #############################################################################################
                width = 20
                height = 14
                if spin_polarized==True:
                    height = 28
                pylab.figure(figsize=(width, height)) #to adjust the figure size

                markers=["o","s","8","D","H","1"]
                colors =['b','g','c','r','m','0.75','y','orange']
                lines  =['-','--',':',]                                   

                if spin_polarized==False:
                        ax2=pylab.subplot(111)

			ax2.tick_params(axis='both', which='major', labelsize=21)

                        sorted_all =  sorted([[float(x[0]),x[1]] for x in data_by_temp.items()])

                        for temp_index in range(len(sorted_all)):                    
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

                            max_x=max(x_vals)
                            min_x=min(x_vals)
                            for i in [y_vals]:
                                set_max=max(i)
                                if set_max>max_val:
                                    max_val=set_max
                            for i in [y_vals]:
                                set_min=min(i)
                                if set_min<min_val:
                                    min_val=set_min
			    
                            ax2.plot(x_vals,y_vals,label=label_text,color=color_choice,linestyle=ls_choice,linewidth=2)


                if spin_polarized==True:
                    ax1=pylab.subplot(211)
                    ax2=pylab.subplot(212)

		    ax1.tick_params(axis='both', which='minor', labelsize=21)
		    ax2.tick_params(axis='both', which='major', labelsize=21)

                    sorted_dn =  sorted([[float(x[0]),x[1]] for x in data_by_temp_dn.items()])
                    sorted_up =  sorted([[float(x[0]),x[1]] for x in data_by_temp_up.items()])

                    for temp_index in range(len(sorted_dn)):                    
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

                            max_x=max(x_vals)
                            min_x=min(x_vals)

                            for i in [y_vals]:
                                set_max=max(i)
                                if set_max>max_val:
                                    max_val=set_max
                            for i in [y_vals]:
                                set_min=min(i)
                                if set_min<min_val:
                                    min_val=set_min

                            ax1.plot(x_vals,y_vals,label=label_text,color=color_choice,linestyle=ls_choice,linewidth=2)
                            pyplot.xlabel(trans_plot_dict[Type]['xl'],{'fontsize':22})
                            pyplot.ylabel(trans_plot_dict[Type]['yl'],{'fontsize':22})
                            

                    pylab.axhline(0.0, color = 'k',linestyle='dashed', linewidth = 1.3)
                    ax1.legend(loc=0,fontsize=22)

		    if Type in ['epsilon_i','epsilon_r','ZT']:
			    ax1.yaxis.set_ticks([0],)

                    pyplot.xlabel(trans_plot_dict[Type]['xl'],{'fontsize':22})
                    pyplot.ylabel(trans_plot_dict[Type]['yl'],{'fontsize':22})
 
                    for temp_index in range(len(sorted_up)):                    
                            label_text='$'+trans_plot_dict[Type]['lc']+'^{up}$ @ %sK'%int(sorted_up[temp_index][0])
                            color_choice  = colors[temp_index%len(colors)]
                            marker_choice = markers[temp_index%len(markers)]
                            ls_choice  = lines[temp_index%len(lines)]
                            x_vals=sorted_up[temp_index][1][0]
                            y_vals=sorted_up[temp_index][1][1]

			    if Type=='kappa':
				    y_vals=numpy.asarray(y_vals)/1.0e17
			    elif Type=='sigma_seebeck':
				    y_vals=numpy.asarray(y_vals)/1.0e-3
			    elif Type=='cond':
				    y_vals=numpy.asarray(y_vals)/1.0e20

                            max_x=max(x_vals)
                            min_x=min(x_vals)

                            for i in [y_vals]:
                                set_max=max(i)
                                if set_max>max_val:
                                    max_val=set_max
                            for i in [y_vals]:
                                set_min=min(i)
                                if set_min<min_val:
                                    min_val=set_min
    
                            ax2.plot(x_vals,y_vals,label=label_text,color=color_choice,linestyle=ls_choice,linewidth=2)

                ax2.axhline(0.0, color = 'k',linestyle='dashed', linewidth = 1.3) #line separating up and down spni

                #we set these earlier
                if Type in ['epsilon_i','epsilon_r']:
                    max_val=max_val_ep
                    min_val=min_val_ep

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

			if Type in ['epsilon_i','epsilon_r','ZT']:
				ax1.yaxis.set_ticks([0.0])

                        ax1.xaxis.set_ticks([])
                        ax1.set_xlim([min_x,max_x])
                    except:
                        pass

                    pyplot.xlabel(trans_plot_dict[Type]['xl'],{'fontsize':22})
                    pyplot.ylabel(trans_plot_dict[Type]['yl'],{'fontsize':22})

                except Exception,e:
                    print e
                    pass
                ax2.set_ylim([min_val*mult_min,max_val*mult_max])
                ax2.legend(loc=0,fontsize=22)

                ax2.set_xlim([min_x,max_x])
		if Type in ['epsilon_i','epsilon_r','ZT']:
			ax2.yaxis.set_ticks([0.0])

                pyplot.xlabel(trans_plot_dict[Type]['xl'],{'fontsize':22})
                pyplot.ylabel(trans_plot_dict[Type]['yl'],{'fontsize':22})

                compoundName = AFLOWpi.retr._getStoicName(oneCalc,strip=True)
                compoundNameLatex = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)

                ax2.legend(loc=0,fontsize=22)
		if Type in ['epsilon_i','epsilon_r','ZT']:
			ax2.yaxis.set_ticks([0])

#                ax2.axes.yaxis.set_ticks_position('right')

                pylab.axhline(0.0, color = 'k',linestyle='dashed', linewidth = 1.3) #line separating up and down spin

                figtitle = r'%s %s' % (trans_plot_dict[Type]['ft'],compoundNameLatex)

                t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=26,horizontalalignment='center') #[x,y] 

                if postfix!='':
                    postfix='_'+postfix

                fileplot = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_%s_%s%s.pdf'%(trans_plot_dict[Type]['fp'],compoundName,ID,postfix))
                matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')


            except Exception,e:
                print e
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

	lines=ep_data_string.split('\n')[1:]
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
