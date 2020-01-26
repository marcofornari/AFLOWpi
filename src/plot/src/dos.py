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
import pickle
import collections
import math

def __smoothGauss(list,strippedXs=False,degree=2):  

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
        print('Plotting DOS')
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

        exten=AFLOWpi.plot._get_plot_ext_type()
        fileplot = os.path.join(subdir,'DOS_%s_%s%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,postfix,exten))

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
                print('Are you sure that you ran ppDos and that it completed properly?')
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
                
                except Exception as e:
                        pass
        if LSDA==True:
                enup  = list(map(float,enup))
                endown= list(map(float,endown))
        else:
                endos=list(map(float,en))  #to convert the list x from float to numbers
        floatdos=list(map(float,dos))
        floatdosDOWN=list(map(float,dosdw))
        enshift = numpy.array(endos) #to treat the list b as an array?
  

        ax2=pylab.subplot(111)
        inDict=AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
#       if "degauss" in inDict["&system"].keys():
        floatdos,floatdosDOWN,enshift=floatdos,floatdosDOWN,enshift
#       else:
#               floatdos,floatdosDOWN,enshift=__smoothGauss(floatdos),__smoothGauss(floatdosDOWN),__smoothGauss(enshift)

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
                floatdosup=list(map(float,dos))
                pylab.plot(enshiftdown,floatdosDOWN,'k-')
                pylab.plot(enshiftup,floatdosUP,'k-')

                pylab.ylim(dosMIN,dosMAX) # scales DOS to larges value of DOS in the given energy range
                pylab.axhline(0.0, color = 'k', linewidth = 1.3) #line separating up and down spin
#       else:
#               pylab.ylim(0,dosMAX) # scales DOS to larges value of DOS in the given energy range
        
        pylab.xlim(yLim[0],yLim[1])
        pylab.axvline(0.0, color = 'k', linestyle='dashed', linewidth = 1.3) #Fermi level line
#       pylab.xticks(numpy.arange(yLim[0],yLim[1]+1,1))
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
        except Exception as e:
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
                                                except Exception as e:
                                                    try:
                                                        mat[j]=[mat[j][0],mat[j][1]+mati[j][1]]
                                                    except Exception as e:
                                                        AFLOWpi.run._fancy_error_log(e)
                                                        

##############################################################################################################################################

                if kresolved:
                        logging.warning('k resolved not supported')
                        print('k resolved not supported')
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
        the projected orbitals across the files of each orbital for each atomic speci5Bes
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
        atomList=list(AFLOWpi.retr._getAtomNum(oneCalc['_AFLOWPI_INPUT_'],strip=False).keys())


        '''
        just the possible names of orbitals that projwfc.x will output
        skips over orbitals for atomic species that are not there
        '''
        orbitalList = ['s','p','d','f','All']
        orbitalList = ['s','p','d','f']

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
                        pDOSFiles= glob.glob(os.path.join(subdir,'%s.pdos_atm*(%s)*wfc*%s*' % (glob_ID,atom,orbital)))

                        if len(pDOSFiles):
                                byAtom.append({'%s' % orbital:pDOSFiles})
                byAtomDict[atom] = byAtom

        for atom,orbital in list(byAtomDict.items()):                   
                        for orbitalDict in range(len(orbital)):
                                for orbitalName,fileList in list(orbital[orbitalDict].items()):
                                        data = __combinePDOS(fileList)
                                        with open(os.path.join(subdir,'%s_%s.sumpdos' % (atom,orbitalName)),'wb') as outputFile:
                                                pickle.dump(data,outputFile,protocol=0)
#                                       numpy.savetxt(os.path.join(subdir,'%s_%s.sumpdos.txt' % (atom,orbitalName)),data)

        byAtom={}
        for atom in atomList:   
                pDOSFiles= glob.glob(os.path.join(subdir,'%s.pdos_atm*(%s)_wfc*' % (glob_ID,atom)))
                if len(pDOSFiles):
                        pDOSDict['%s_All'% (atom)] =  pDOSFiles


        '''Cycle through all of the atoms in a single calculation and sum over the pdos files
        in the calculation directory and saves  '''

        for atom,files in list(pDOSDict.items()):               
                data = __combinePDOS(files)
                with open(os.path.join(subdir,'%s.sumpdos' % (atom)),'wb') as outputFile:
                                        pickle.dump(data,outputFile,protocol=0)
#               numpy.savetxt(os.path.join(subdir,'%s.sumpdos.txt' % (atom)),data)

                        

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

                    Efermi=AFLOWpi.retr._getEfermi(oneCalc,fermi_ID,directID=True)

            else:
                    Efermi=AFLOWpi.retr._getEfermi(oneCalc,ID)
                    if type(Efermi)!=type(0.5):
                            Efermi=Efermi[0]

        except Exception as e:

                Efermi=0.0
        try:
            ax2=ax[speciesNum]
        except:
            ax2=ax

        '''extracts HOMO from nscf calculation output file as input to the plotting'''

        """get the path to the subdirectory of the calc that you are making plots for"""        
        subdir = oneCalc['_AFLOWPI_FOLDER_']
        
        exten=AFLOWpi.plot._get_plot_ext_type()
        try:
                filePlotName = 'PDOS_%s%s_%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atom,exten)
        except IOError:
                filePlotName = 'PDOS_%s%s_%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atom,exten)

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
                        data = pickle.load(dataFile)

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
                        except Exception as e:
                                pass


                endos=list(map(float,en))  #to convert the list x from float to numbers
                floatdos=list(map(float,ldos))
                floatdosDOWN=list(map(float,ldosDOWN))
                enshift = numpy.array(endos) #to treat the list b as an array?
                inDict=AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])

#               if "degauss" in inDict["&system"].keys():
                return enshift,floatdos,floatdosDOWN
#               else:
#                       return __smoothGauss(enshift),__smoothGauss(floatdos),__smoothGauss(floatdosDOWN)



        maxDOS=0
        minDOS=0
        orbitalList = ['s','p','d','f',"All"]
        orbitalList = ['s','p','d','f',]
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
                        
#               if orbitalName != 'All':

                ax2.plot(enshift,floatdos,color+'-',label=orbitalName)
                if LSDA:
                        ax2.plot(enshift,floatdosDOWN,color+'-',label=orbitalName)
        
                        
        handles, labels = ax2.get_legend_handles_labels()

        if LSDA:
                ax2.legend(handles[::-2], labels[::-2],loc=1)
                pylab.ylim(minDOS,maxDOS) # scales DOS to larges value of DOS in the given energy range
                pylab.axhline(0.0, color = 'k', linewidth = 1.3) #line to separate up and down spin
        else:
                ax2.legend(handles[::-1], labels[::-1],loc=1)
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
#        ax2.set_xticks(numpy.arange(yLim[0],yLim[1]+1,1))


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
                for ID,oneCalc in list(calcs.items()):
                        __dosPlot(oneCalc,ID,yLim,postfix=postfix)
        else:
                for ID,oneCalc in list(calcs.items()):
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

                    Efermi=AFLOWpi.retr._getEfermi(oneCalc,fermi_ID,directID=True)
            else:
                    Efermi=AFLOWpi.retr._getEfermi(oneCalc,ID)
                    if type(Efermi)!=type(0.5):
                            Efermi=Efermi[0]
        except:
                Efermi=0.0

#       print 'Efermi/HOMO-LUMO = %s' % Efermi

        try:
            ax2=ax[speciesNum]
        except:
            ax2=ax

        '''extracts HOMO from nscf calculation output file as input to the plotting'''

        """get the path to the subdirectory of the calc that you are making plots for"""        
        subdir = oneCalc['_AFLOWPI_FOLDER_']
        
        exten=AFLOWpi.plot._get_plot_ext_type()
        try:
                filePlotName = 'PDOS_%s%s_%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atom,exten)
        except IOError:
                filePlotName = 'PDOS_%s%s_%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,atom,exten)

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
                        data = pickle.load(dataFile)

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
                        except Exception as e:
                                pass


                endos=list(map(float,en))  #to convert the list x from float to numbers
                floatdos=list(map(float,ldos))
                floatdosDOWN=list(map(float,ldosDOWN))
                enshift = numpy.array(endos) #to treat the list b as an array?


                inDict=AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])
#               if "degauss" in inDict["&system"].keys() or TB==True:
                return enshift,floatdos,floatdosDOWN
#               else:
#                       return __smoothGauss(enshift),__smoothGauss(floatdos),__smoothGauss(floatdosDOWN)



        maxDOS=0
        minDOS=0
        orbitalList = ['s','p','d','f',"All"]
        orbitalList = ['s','p','d','f',]
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
                        
#               if orbitalName != 'All':
                ax2.plot(enshift,floatdos,color+'-',label=orbitalName)

                if LSDA:
                        ax2.plot(enshift,floatdosDOWN,color+'-',label=orbitalName)
        
                        
        handles, labels = ax2.get_legend_handles_labels()

        if LSDA:
                ax2.legend(handles[::-2], labels[::-2],loc=1)
                pylab.ylim(minDOS,maxDOS) # scales DOS to larges value of DOS in the given energy range
                pylab.axhline(0.0, color = 'k', linewidth = 1.3) #line to separate up and down spin
        else:
                ax2.legend(handles[::-1], labels[::-1],loc=1)
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

#        ax2.set_xticks(numpy.arange(yLim[0],yLim[1]+1,2))


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

def apdos(calcs,yLim=[-10,10],runlocal=False,postfix='',scale=False,tight_binding=False):
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

        for ID,oneCalc in list(calcs.items()):  
                if runlocal:
                        AFLOWpi.plot.__apdos(oneCalc,ID,yLim,postfix=postfix,scale=scale,tight_binding=tight_binding)
                else:
                        AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__apdos(oneCalc,ID,[%s,%s],postfix='%s',scale=%s,tight_binding=%s)" % (yLim[0],yLim[1],postfix,scale,tight_binding))


def __apdos(oneCalc,ID,yLim,postfix='',scale=False,tight_binding=False,label_map={}):

        AFLOWpi.plot.__sumpdos(oneCalc,ID,TB=tight_binding)     
        subdir = oneCalc['_AFLOWPI_FOLDER_']


        print(('Plotting atomic projected DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True))))
        logging.info('Plotting atomic projected DOS of %s ' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True)))
        petype = AFLOWpi.plot._get_plot_ext_type()
        fileplot = os.path.join(subdir,'APDOS_%s_%s%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,postfix,petype)) 

        LSDA=False
        nspin = int(AFLOWpi.scfuj.chkSpinCalc(oneCalc,ID))
        if nspin!=1:
            LSDA=True

        try:
            if tight_binding==True:
                    fermi_ID = '%s_WanT_dos'%calcID

                    Efermi=AFLOWpi.retr._getEfermi(oneCalc,fermi_ID,directID=True)
            else:
                    Efermi=AFLOWpi.retr._getEfermi(oneCalc,ID)
                    if type(Efermi)!=type(0.5):
                            Efermi=Efermi[0]
        except Exception as e:
                Efermi=0.0


        dos_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='bands',last=True)

        width = 12
        height = 8
        pylab.figure(figsize=(width, height)) #to adjust the figure size

        ax2=pylab.subplot(111)



        def getPlotData(sumpdosFile):
                try:
                        with open(sumpdosFile,'rb') as dataFile:
                                data = pickle.load(dataFile)

                except Exception as e:

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

                        except Exception as e:
                                pass

                endos=list(map(float,en))  #to convert the list x from float to numbers
                floatdos=list(map(float,ldos))
                floatdosDOWN=list(map(float,ldosDOWN))
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




                ''' scales DOS to larges value of DOS in the given energy range and finds the largest DOS between the different orbitals'''
                try:
                        if max([floatdos[k] for k in range(len(floatdos)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ]) > maxDOS:
                                maxDOS = max([floatdos[k] for k in range(len(floatdos)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ])

                        if min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ]) < minDOS:
                                minDOS = min([floatdosDOWN[k] for k in range(len(floatdosDOWN)) if enshift[k] < yLim[1] and enshift[k] > yLim[0] ])
                except:
                    pass


                if species in list(label_map.keys()):
                        species_lab = label_map[species]
                else:
                        species_lab = species

                if species=='TOTAL':    
                        pylab.plot(enshift,floatdos,'-',label=species_lab,color='k',linewidth=2)       
                else:
                        pylab.plot(enshift,floatdos,'-',label=species_lab,linewidth=2,color=colors[fi%len(colors)])                        




                if LSDA:
                        if species=='TOTAL':
                                pylab.plot(enshift,floatdosDOWN,'-',label=species_lab,color='k',linewidth=2) 
                        else:
                                pylab.plot(enshift,floatdosDOWN,'-',label=species_lab,linewidth=2,color=colors[fi%len(colors)])               

        handles, labels = ax2.get_legend_handles_labels()

        pylab.axvline(0.0, color = 'k', linewidth = 2.0,linestyle='--') #line separating up and down spin

        if LSDA:
                ax2.legend(handles[::2], labels[::2],fontsize=14,loc=1)
                dosRange=max([numpy.abs(minDOS),numpy.abs(maxDOS)])
                pylab.ylim(1.1*minDOS,1.1*maxDOS) # scales DOS to larges value of DOS in the given energy range 

                pylab.axhline(0.0, color = 'k', linewidth = 2.0) #line separating up and down spin
        else:
                ax2.legend(handles, labels,fontsize=14,loc=1)
                pylab.ylim(0,1.1*maxDOS) # scales DOS to larges value of DOS in the given energy range

        pylab.xlim(yLim)



        ax2.spines['bottom'].set_linewidth(1.5)
        ax2.spines['left'].set_linewidth(1.5)
        ax2.spines['right'].set_linewidth(1.5)
        ax2.spines['top'].set_linewidth(1.5)






        ax2.yaxis.set_ticks_position('left')
        pylab.ylabel('Density of States (States/eV)',fontsize=20)
        pylab.xlabel('E (eV)',fontsize=20)

        ax2.axes.set_yticklabels([])
#        ax2.axes.xaxis.set_label_position('top')
        locs, labels = pylab.xticks()


        #save fig
        matplotlib.pyplot.savefig(fileplot,bbox_inches='tight')

        try:
                AFLOWpi.retr._moveToSavedir(fileplot)
        except Exception as e:
                pass

        pyplot.cla()
        pyplot.clf()
        pyplot.close()



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

        for ID,oneCalc in list(calcs.items()):  
                if runlocal:
                        AFLOWpi.plot.__opdos(oneCalc,ID,yLim,postfix=postfix,scale=scale,tight_binding=tight_binding)
                else:
                        AFLOWpi.prep._addToBlock(oneCalc,ID,'PLOT',"AFLOWpi.plot.__opdos(oneCalc,ID,[%s,%s],postfix='%s',scale=%s,tight_binding=%s)" % (yLim[0],yLim[1],postfix,scale,tight_binding))


def __opdos(oneCalc,ID,yLim,postfix='',scale=False,tight_binding=False,label_map={}):

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

        atomList=list(AFLOWpi.retr._getAtomNum(oneCalc['_AFLOWPI_INPUT_'],strip=False).keys())



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

                spec_lab=atomList[species]
                if spec_lab in list(label_map.keys()):
                        spec_lab=label_map[spec_lab]

                ax.axes.set_title(spec_lab,weight='bold',loc='left',fontsize=18)

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

        for species,ax in list(atomAxisArray.items()):
            for ax in fig.get_axes():

                ax.axes.set_ylim([1.05*old_low,1.05*old_high])
                ax.axes.set_xlim(yLim[0],yLim[1])
                if  LSDA:
                        ax.axes.axhline(0.0,linewidth=1.3,color="k")

        subdir=oneCalc['_AFLOWPI_FOLDER_']
                
        """get the path to the subdirectory of the calc that you are making plots for"""

        '''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
        exten=AFLOWpi.plot._get_plot_ext_type()

        fileplot = os.path.join(subdir,'PDOS_%s_%s%s.%s' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,postfix,exten))

        fig.savefig(fileplot,bbox_inches='tight')
        try:
                AFLOWpi.retr._moveToSavedir(fileplot)
        except:
                pass

        pyplot.cla()
        pyplot.clf()
        pyplot.close()
