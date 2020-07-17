import AFLOWpi
import logging
import os

def _post_proc(oneCalc,ID,plot_num,**kwargs):
    if "data_postfix" in list(kwargs.keys()):
        data_postfix="_"+kwargs["data_postfix"]
        del kwargs["data_postfix"]
    else:
        data_postfix=""
    data_postfix=""

    pp_input='''&inputpp
prefix = '%s',
outdir = './'
filplot = 'pp_temp'
plot_num=%s
/
&PLOT
nfile = 1
filepp(1) = 'pp_temp'
weight(1) = 1.0
fileout = '%s_pp_%s%s.dat'
''' %(oneCalc["_AFLOWPI_PREFIX_"],plot_num,ID,plot_num,data_postfix)

    pp_input+=",\n".join(["%s=%s"%(k,v) for k,v in kwargs.items()])+"\n/\n"

    infn=os.path.join(oneCalc["_AFLOWPI_FOLDER_"],"%s_pp_%02d.in"%(ID,plot_num))
    with open(infn,"w") as ofo:
        ofo.write(pp_input)

    execPrefix = AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
    AFLOWpi.run._oneRun("",oneCalc,"%s_pp_%02d"%(ID,plot_num),engine='espresso',calcType='custom',execPath='./pp.x',execPrefix=execPrefix ) 
