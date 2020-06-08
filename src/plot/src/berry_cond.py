
import AFLOWpi
import numpy as np
import os
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import glob 

def __plot_berry_cond(oneCalc,ID,spin=False,en_range=None):



    extension = AFLOWpi.plot._get_plot_ext_type()

    subdir = oneCalc['_AFLOWPI_FOLDER_']

    if spin:
        dem_files = glob.glob(subdir+'/%s_shcEf_*'%ID)
    else:
        dem_files = glob.glob(subdir+'/%s_ahcEf_*'%ID)


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
            lab=r"$\sigma^{%s}_{%s%s}$"%(spol,ipol,jpol)
            plt.plot(dat[:,0],dat[:,1],label=lab,color="r")
        else:
            lab=r"$\sigma_{%s%s}$"%(ipol,jpol)
            plt.plot(dat[:,0],dat[:,1],label=lab,color="c")


        chem_name = AFLOWpi.retr._getStoicName(oneCalc,strip=True,latex=True)

        plt.ylim([np.amin(dat[:,1])*1.05,np.amax(dat[:,1])*1.05])
        plt.xlim([np.amin(dat[:,0]),np.max(dat[:,0])])
        plt.legend(loc=1)
        plt.axhline(0.0,ls="--",color="k")
        plt.axvline(0.0,ls="-",color="k")

        plt.xlabel(r"$\hbar \omega$ $(eV)$")

        if spin:
            plt.ylabel(r"$\sigma^{%s}_{%s%s}$ $(\Omega^{-1} cm^{-1})$"%(spol,ipol,jpol))

            if AFLOWpi.plot._get_title_option():
                plt.title('Spin Hall Conductivity: %s'%chem_name)

            fileplot = os.path.join(subdir,'SPIN_HALL_COND_%s_%s%s_%s_%s.%s' % (spol,ipol,jpol,AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,extension))
        else:
            plt.ylabel(r"$\sigma_{%s%s}$ $(\Omega^{-1} cm^{-1})$"%(ipol,jpol))
            if AFLOWpi.plot._get_title_option():
                plt.title('Anomalous Hall Conductivity: %s'%chem_name)

            fileplot = os.path.join(subdir,'ANOM_HALL_COND_%s%s_%s_%s.%s' % (ipol,jpol,AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID,extension))


        matplotlib.pyplot.savefig(fileplot,bbox_inches='tight',layout="tight")
        AFLOWpi.plot._copy_to_fig_dir(oneCalc,fileplot)
        plt.close()
