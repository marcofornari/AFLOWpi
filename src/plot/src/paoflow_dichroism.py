
import AFLOWpi
import numpy as np
import os
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import glob

def __plot_dichroism(oneCalc,ID,spin=False,real=False,en_range=None):


    extension='png'
    subdir = oneCalc['_AFLOWPI_FOLDER_']

    if real:
        if spin:
            dem_files = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/SCDr_*')
        else:
            dem_files = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/MCDr_*')
    else:
        if spin:
            dem_files = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/SCDi_*')
        else:
            dem_files = glob.glob(oneCalc['_AFLOWPI_FOLDER_']+'/MCDi_*')


    for dat_file in dem_files:

        width = 9
        height = 6
        plt.figure(figsize=(width, height))#to adjust the figure size

        dat_name = dat_file[:-4]
        dat_split_name = dat_name.split('_')

        ipol=dat_split_name[-1][0]
        jpol=dat_split_name[-1][1]
        if spin:
            spol=dat_split_name[-2]
        else:
            spol=''

        dat = np.loadtxt(dat_file)

        if en_range is not None:
            dat=dat[np.where(np.logical_and(dat[:,0]>en_range[0],dat[:,0]<en_range[1]))]

        if spin:
            if real:
                lab = r"Re[$\omega\sigma^{%s}_{%s%s}$]"%(spol,ipol,jpol)
            else:
                lab = r"Im[$\omega\sigma^{%s}_{%s%s}$]"%(spol,ipol,jpol)
            plt.plot(dat[:,0],dat[:,1],label=lab,color="m")
        else:
            if real:
                lab = r"Re[$\omega\sigma_{%s%s}$]"%(ipol,jpol)
            else:
                lab = r"Im[$\omega\sigma_{%s%s}$]"%(ipol,jpol)
            plt.plot(dat[:,0],dat[:,1],label=lab,color="g")


        plt.ylim([np.amin(dat[:,1])*1.05,np.amax(dat[:,1])*1.05])
        plt.xlim([np.amin(dat[:,0]),np.amax(dat[:,0])])
        plt.legend(loc=1)
        plt.axhline(0.0,ls="-",color="k")

        plt.xlabel(r"$\hbar \omega$ $(eV)$")
        chem_name = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)

        if spin:
            if real:                
                plt.ylabel(r"Re[$\omega\sigma^{%s}_{%s%s}$] ($10^{29}sec^{-2}$)"%(spol,ipol,jpol))
                fileplot = os.path.join(subdir,'SPIN_DICHRO_REAL_%s_%s%s_%s_%s.%s'%(spol,ipol,jpol,AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,extension))
            else:
                plt.ylabel(r"Im[$\omega\sigma^{%s}_{%s%s}$] ($10^{29}sec^{-2}$)"%(spol,ipol,jpol))
                fileplot = os.path.join(subdir,'SPIN_DICHRO_IMAG_%s_%s%s_%s_%s.%s'%(spol,ipol,jpol,AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,extension))
            plt.title('Spin Circular Dichroism: %s'%chem_name)


        else:
            if real:
                plt.ylabel(r"Re[$\omega\sigma_{%s%s}$] ($10^{29}sec^{-2}$)"%(ipol,jpol))
                fileplot = os.path.join(subdir,'MAG_DICHRO_REAL_%s%s_%s_%s.%s'%(ipol,jpol,AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,extension))
            else:
                plt.ylabel(r"Im[$\omega\sigma_{%s%s}$] ($10^{29}sec^{-2}$)"%(ipol,jpol))
                fileplot = os.path.join(subdir,'MAG_DICHRO_IMAG_%s%s_%s_%s.%s'%(ipol,jpol,AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,extension))
            plt.title('Magnetic Circular Dichroism: %s'%chem_name)





        matplotlib.pyplot.savefig(fileplot,bbox_inches='tight',layout="tight")
        plt.close()
