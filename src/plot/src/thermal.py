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
import os
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot
from matplotlib import pylab
import numpy as np
import io
import logging 
import scipy
import matplotlib.gridspec as gridspec


def plot_lattice_TC(calcs,postfix="",temp_range=None):
        for ID,oneCalc in list(calcs.items()):
                AFLOWpi.plot._plot_lattice_TC(oneCalc,ID,temp_range=temp_range)

def plot_gruneisen(calcs,postfix="",w_range=None,grun_range=None):
        for ID,oneCalc in list(calcs.items()):
                AFLOWpi.plot.__gruneisen_of_omega_ap(oneCalc,ID,w_range=w_range,grun_range=grun_range)

def __plot_gruneisen(oneCalc,ID,postfix="",THz=False,w_range=None):

        
        calcCopy = oneCalc
        calcID = ID

        calcID = AFLOWpi.prep._return_ID(oneCalc,calcID,step_type='phonon',last=True)        


        if postfix!='':
                postfix='_'+postfix

        subdir=oneCalc['_AFLOWPI_FOLDER_']
        fileplot = os.path.join(subdir,'GRUN_PATH_%s_%s%s.png' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,postfix))

        """get the path to the subdirectory of the calc that you are making plots for"""

        filebands = os.path.join(subdir,'%s.phBAND.gp'%calcID)




        '''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''   


       #to set figure size and default fonts
        matplotlib.rc("font", family="serif")      #to set the font type
        matplotlib.rc("font", size=20)             #to set the font size

        
        width = 20
        height = 8
        fig = pylab.figure(figsize=(width, height))#to adjust the figure size
     

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
                    except Exception as e:
                        pass

                k_y = np.asarray(y).T.tolist()

                for entry in k_y:
                    k_x.append(x)


        except Exception as e:
                AFLOWpi.run._fancy_error_log(e)        
                logging.warning("output from bands calculation not found. Are you sure you ran ppBands and it completed properly?")
                print("Are you sure you ran ppBands and it completed properly?")
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

#       for j in range(len(k_y[0])):
#               if k_y[0][j]==0.0 and k_y[1][j]==0.0 and k_y[1][j]==0.0:
#                       for i in range(len(k_y)):
#                               k_y[i][j]=0.0
#       gamma_index=
#       for i in range(len(k_y[0])):
        
####################################################################################
#       gamma_index=[]
#       for sym in range(len(SymPrint)):
#               if SymPrint[sym]=='$\\Gamma$':
#                       gamma_index.append(symIndex[sym])


####################################################################################

        #set frequency of optical branches at gamma to either one before or one after 
        #so that the LOTO splitting at gamma doesn't result in an outlier point

        gamma_index=[]
        by_k = np.asarray(k_y).T
        for i in range(len(by_k)):
                num_zero =  len(by_k[i])-np.count_nonzero(by_k[i])
                if num_zero>=3:
                        gamma_index.append(i)

        for j in gamma_index:
                for i in range(len(k_y)):
                        k_y[i][j]=0.0

                        
        gdat = np.loadtxt(os.path.join(subdir,"%s.phSCATTER.band"%calcID)).T
        print((gdat.shape))
#       gdat = np.sqrt(gdat)
        gdat[np.where(gdat>15)]=15.0
        gdat[np.where(gdat<-15)]=-15.0

        cmap_neg = pylab.cm.get_cmap('YlGnBu_r')

        gdat =gdat

#       gdat =gdat/np.amax(gdat)*255


        colors=gdat




        k_x = np.asarray(k_x)
        k_y = np.asarray(k_y)
        k_x = np.ravel(k_x)
        k_y = np.ravel(k_y)

        colors = np.ravel(colors)



#       colors[np.where(k_y<1.0)]=1.e-16
        csort = np.argsort(colors)
#       print csort.shape

        k_y = k_y[csort]#[::-1]
        k_x = k_x[csort]#[::-1]

        colors=colors[csort]#[::-1]

        print((colors.shape))
#       print colors

        normalize_pos = matplotlib.colors.Normalize(vmin=np.amin(colors), 
                                                    vmax=np.amax(colors))

#       sizes[np.where(sizes<0.1)]=0.3



        pylab.scatter(k_x,k_y,s=10.0,marker="o",c=colors,cmap="cool",norm=normalize_pos)






        min_val=np.amin(k_y)
        max_val=np.amax(k_y)



            

            


        if THz==True:
            pylab.ylabel('Frequency (THz)')
        else:
            pylab.ylabel('Frequency (cm$^{-1}$)')
        pylab.xlim(np.amin(k_x),np.amax(k_x)) 
        
        if w_range is not None:
                min_val=w_range[0]/1.1
                max_val=w_range[1]/1.1

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
        buf = io.StringIO(bandSym)

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
                        except Exception as e:
                                print(e)
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
                except Exception as e:
                        pass
     #Print path labels to band structure x-axis
        try:
                pylab.xticks([a[index] for index in symIndex],SymPrint)
        except Exception as e:
                pass
                return
        pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Femi level line
        locs, labels = pylab.xticks()



##########################################################################################################

##########################################################################################################


        if w_range is not None:
                min_val=w_range[0]/1.1
                max_val=w_range[1]/1.1

        ax1.set_ylim(1.1*min_val,max_val*1.1) 


        pylab.axhline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line



        ax1.set_position([0.07,0.1,0.895,0.95]) #[left,bottom,width,height]             
#       ax1.set_position([0.07,0.1,0.925,0.9]) #[left,bottom,width,height]              

        neg_cax = fig.add_axes([0.97,0.25,0.02,0.6])


                
        cbar = matplotlib.colorbar.ColorbarBase(neg_cax, cmap="cool",
                                                norm=normalize_pos,
                                                extend="neither",
                                                spacing='uniform',
                                                orientation='vertical')

        
#       cbar.ax.tick_params('both',width=2,length=29)
#       cbarneg.ax.tick_params('both',width=2,length=29)
#       cbar.ax.tick_params(direction='in')

        
        
        
        cbar.ax.set_title('\n  $\gamma$',size=32 ,va="bottom")
                        

        labsti = [t.get_text() for t in cbar.ax.get_yticklabels() ]
        labsti[-1]+="+"
        print(labsti)
        cbar.update_ticks()
        tick_pos = cbar.get_ticks()

        cbar.set_ticklabels(labsti,update_ticks=False)
        cbar.update_ticks()



#       figtitle = 'Phonon Dispersion: %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)) 
#       t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=24,horizontalalignment='center') #[x,y]

        matplotlib.pyplot.savefig(fileplot,bbox_inches='tight',dpi=500)
        AFLOWpi.plot._copy_to_fig_dir(oneCalc,fileplot)

        pyplot.cla()
        pyplot.clf()
        pyplot.close()



def __gruneisen_of_omega(oneCalc,ID,projected=True):
    matplotlib.rc("font", size=24)             #to set the font size

    width = 10.0
    height = 24.0
    pylab.figure(figsize=(width, height))#to adjust the figure size
    gs = gridspec.GridSpec(3,1)
    therm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_up')
    for branch in range(3):
            therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER_%s.gp'%(therm_ID,branch))

            therm_data=[]

            with open(therm_file_name,"r") as fo:
                data = fo.read()

            data=data.split('\n')
            therm_data = np.asarray([list(map(float,x.split())) for x in data if len(x.strip())!=0])

            therm_data=np.asarray(therm_data)
            therm_data[:,0]=np.around(therm_data[:,0],decimals=2)
#           therm_data[:,1]=np.around(therm_data[:,1],decimals=2)

        #    therm_data = np.unique(b).view(therm_data.dtype).reshape(-1, therm_data.shape[1])

            therm_data  = np.vstack({tuple(row) for row in therm_data})
        #   therm_data[:,0] = np.ma.masked_where(therm_data[:,0] <= 1.0, therm_data, copy=False)
            for i in range(len(therm_data)):
                    if therm_data[i][0]<1.1 or therm_data[i][1]>10000.0:
                            therm_data[i][1]=0.0


            x = np.sort(np.unique(therm_data[:,0]))
            y = np.zeros_like(x)
            for i in range(x.shape[0]):
                    inds = np.where(therm_data[:,0]==x[i])[0]
                    y[i] = np.mean(therm_data[inds,1])
            from scipy.signal import savgol_filter
            y1 = savgol_filter(y, 199, 2)

            therm_data[:,1] = np.ma.masked_equal(therm_data[:,1],0.0)
            therm_data = np.ma.masked_equal(therm_data,0.0)
            therm_data[:,1] = np.ma.masked_greater(therm_data[:,1],1000.0)

            print((therm_data.shape))

            ax1 = pylab.subplot(gs[branch])
            pylab.ylabel('$\gamma^{2}$')
            pylab.xlabel('$\omega$ $(cm^{-1})$')
            pylab.plot(therm_data[:,0],therm_data[:,1]**2,'k',linestyle=' ',marker='o',fillstyle='none')

            
            pylab.axhline(np.mean(therm_data[:,1]),color='b',linewidth=2.0)
            pylab.plot(x,y1,'r',marker='',linewidth=2.0)
            figtitle = '$Gr\ddotuneisen$ $Parameter:$ %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)) 
            t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=24,horizontalalignment='center') #[x,y]
            pylab.xlim([0.0,np.amax(therm_data[:,0])])

    ext=AFLOWpi.plot._get_plot_ext_type()

    fileplot = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'SCATTER_%s_%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,ext))

    matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')
    AFLOWpi.plot._copy_to_fig_dir(oneCalc,fileplot)

    pyplot.cla()
    pyplot.clf()
    pyplot.close()

    

# def __plot_thermal_conductivity(oneCalc,ID):

#     therm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal')
#     therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_thermal_cond.dat'%therm_ID)

#     with open(therm_file_name,'r') as tcfo:
#         therm_data_str = tcfo.read()


import  matplotlib.pyplot

def _plot_lattice_TC(oneCalc,ID,temp_range=[80.0,800.0]):

        fname = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],)

        therm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_up')
        therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_thermal_cond.dat'%therm_ID)

        data=np.loadtxt(therm_file_name, dtype=np.float, comments='#',skiprows=1)
        matplotlib.pyplot.plot(data[:,0],data[:,1])

        dc = scipy.where(np.logical_and(data[:,0]>temp_range[0],data[:,0]<temp_range[1]))
        data_cut=data[dc]
        max_y = np.amax(data_cut[:,1])*1.1
        matplotlib.pyplot.ylim([0.0,max_y])

        matplotlib.pyplot.xlim(temp_range)
        matplotlib.pyplot.xlabel('T (K)')

        figtitle = 'Lattice Thermal Conductivity: %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)) 
        t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=14,horizontalalignment='center') #[x,y]

        matplotlib.pyplot.ylabel(r'$\kappa_{lat}$ $(\frac{W}{m\cdot K})$')

        ext=AFLOWpi.plot._get_plot_ext_type()

        fig_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"LATTICE_TC_%s_%s.%s"%(AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,ext))
        matplotlib.pyplot.savefig(fig_file_name,bbox_inches='tight')
        AFLOWpi.plot._copy_to_fig_dir(oneCalc,fig_file_name)




def __gruneisen_of_omega_ap(oneCalc,ID,w_range=None,grun_range=None,label_map={}):
    matplotlib.rc("font", size=20)             #to set the font size

    therm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_up')

    therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER.ap'%therm_ID)

    therm_data=[]

    with open(therm_file_name,"r") as fo:
        data = fo.read()

    data=data.split('\n')
    labs = data[0].split()[2:]


    therm_data = np.asarray([list(map(float,x.split())) for x in data[1:] if len(x.strip())!=0])
    
    therm_data=np.asarray(therm_data)
#    therm_data[:,0]=np.around(therm_data[:,0],decimals=1)
#    therm_data[:,1:]=np.around(therm_data[:,1:],decimals=1)


    
    therm_data  = np.vstack({tuple(row) for row in therm_data})

    therm_data[np.where(therm_data[:,0]<1.1)[0]]=0.0
            



#    therm_data = np.ma.masked_equal(therm_data,0.0)
    therm_data[np.where(np.abs(therm_data[:,1:])>100.0)[0],1:]=0
    

    width  = 18.0
    height = (len(labs)-1)*4.0
    

    fig = pylab.figure(figsize=(width, height))#to adjust the figure size
    gs = gridspec.GridSpec(len(labs),1)
    gs.update(hspace=0.10,top=0.915)

    
    dw =np.linspace(0.0,np.amax(therm_data[:,0]),250,endpoint=True)

    therm_data[:,1] = np.nan_to_num(therm_data[:,1])

    avg=np.zeros((dw.shape[0]-1,2))
    for i in range(dw.shape[0]-1):
            inds=np.where(np.logical_and(therm_data[:,0]<dw[i+1],therm_data[:,0]>=dw[i]))[0]
            if len(inds)==0:
                    try:
                            avg[i,1] = avg[i-1,1]
                    except: avg[i,1] =0.0
            else:
                    avg[i,1] =  np.nan_to_num(np.mean(therm_data[inds,1]**2)**(0.5))
                    
            avg[i,0] = dw[i]
            print((avg[i]))
            



    color_cycle=['k','r', 'g', 'b', 'orange','c', 'm', 'y','k']
                    
                    

    if w_range is None:
            w_range=[0.0,np.amax(therm_data[:,0])]
    if grun_range is None:
            if np.amin(therm_data[:,1])<0.0:
#                   grun_range=[np.amin(therm_data[:,1])*1.05,np.amax(therm_data[:,1])*1.05]
                    grun_range=[0.0,np.amax(therm_data[:,2:])*1.05]
            else:
                    grun_range=[np.amin(therm_data[:,1])*0.95,np.amax(therm_data[:,1])*1.05]
#           grun_range=[0.0,np.amax(therm_data[:,1])*1.05]


    for i in range(1,len(labs)):
            if labs[i] in list(label_map.keys()):
                    new_lab = label_map[labs[i]]
            else:
                    new_lab = labs[i]
            ax = pylab.subplot(gs[i])
            print((labs[i],(np.mean(therm_data[:,i+1]))**(0.5)))
            ax.plot(therm_data[:,0],np.abs(therm_data[:,i+1]),linestyle=' ',
                    marker='o',fillstyle='none',label=new_lab,
                    color=color_cycle[i%len(color_cycle)])
#           if labs[i]=="Total":
#                   ax.plot(avg[:,0],avg[:,1],color="r")
            if i<(len(labs)-1):
                    ax.set_xticklabels([])

            ax.set_ylabel('$\|\gamma\|_{i \\alpha}$',fontsize=24)

            ax.set_ylim(grun_range)
            ax.set_xlim(w_range)
            ax.legend(fontsize=24,loc=1)
            ax.axhline(0.0,color="k",ls="--")

    ax.set_xlabel('$\omega$ $(cm^{-1})$',fontsize=24)
    figtitle = '$Gr\ddotuneisen$ $Parameter:$ %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)) 
    t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=28,horizontalalignment='center') #[x,y]

    fileplot = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'SCATTER_ap_%s_%s.png' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,))
        
    matplotlib.pyplot.savefig(fileplot,bbox_inches='tight',layout='tight',dpi=300)
    AFLOWpi.plot._copy_to_fig_dir(oneCalc,fileplot)
    matplotlib.pyplot.close()

    fileplot = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'SCATTER_%s_%s.png' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,))

    fig = pylab.figure(figsize=(12, 10))#to adjust the figure size
    ax = pylab.subplot(111)

    therm_data[:,1] = np.ma.masked_equal(therm_data[:,1],0.0)
    ax.plot(therm_data[:,0],therm_data[:,1]**2,linestyle=' ',
            marker='o',fillstyle='none',color="k")


    ax.set_ylabel('$\gamma^{2}$',fontsize=24)

    sl = np.where(np.logical_and(therm_data[:,0]>=w_range[0],
                                 therm_data[:,0]<=w_range[1]))

    llim=np.amin(therm_data[sl,1]**2)*1.05
    if llim>0:
            llim=0
    ulim=np.amax(therm_data[sl,1]**2)*1.05
    if ulim<0:
            ulim=0

    ax.set_ylim([llim,ulim])
#    ax.set_ylim(grun_range)
    ax.set_xlim(w_range)
#    ax.set_xlim([0.0,180.0])
    
    ax.axhline(0.0,color="k",ls="--")

    ax.set_xlabel('$\omega$ $(cm^{-1})$',fontsize=24)
    figtitle = '$Gr\ddotuneisen$ $Parameter:$ %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)) 
    t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=28,horizontalalignment='center') #[x,y]
    matplotlib.pyplot.savefig(fileplot,bbox_inches='tight',layout='tight',dpi=300)
    AFLOWpi.plot._copy_to_fig_dir(oneCalc,fileplot)

    pyplot.cla()
    pyplot.clf()
    pyplot.close()

    fileplot = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'SCATTER2_%s_%s.png' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,))

    fig = pylab.figure(figsize=(12, 10))#to adjust the figure size
    ax = pylab.subplot(111)


#    avg = np.ma.masked_equal(avg,0.0)
#    ax.plot(therm_data[:,0],therm_data[:,1]**2,linestyle=' ',
#           marker='o',fillstyle='none',color="k")
    ax.plot(avg[:,0],avg[:,1]**2,linewidth=4,label="average",
            marker='',color="r")
#           marker='o',fillstyle='none',color="k")
    ax.legend()

    ax.set_ylabel('$\gamma^{2}$',fontsize=24)




    w_range=[0,np.amax(avg[:,0])]
    ax.set_ylim(grun_range)
    ax.set_ylim([0.0,np.amax(avg[:,1]**2)*1.05])
    ax.set_xlim(w_range)
    


    ax.set_xlabel('$\omega$ $(cm^{-1})$',fontsize=24)
    figtitle = '$Gr\ddotuneisen$ $Parameter:$ %s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)) 
    t = pylab.gcf().text(0.5,0.92, figtitle,fontsize=28,horizontalalignment='center') #[x,y]
    matplotlib.pyplot.savefig(fileplot,bbox_inches='tight',layout='tight',dpi=300) 
    AFLOWpi.plot._copy_to_fig_dir(oneCalc,fileplot)           

    pyplot.cla()
    pyplot.clf()
    pyplot.close()

