import AFLOWpi
import os 
import re
import glob
import numpy as np

def _prep_berry(oneCalc,ID,gdir,kp_mult):

    inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    inputDict["&control"]["lberry"] = ".true."
    inputDict["&control"]["gdir"] = str(gdir)
    inputDict["&control"]["calculation"] = "'nscf'"


    grid = list(map(int,inputDict["K_POINTS"]["__content__"].split()))

    grid[gdir-1] = int(kp_mult*grid[gdir-1])
    try:
        if gdir>1:
            grid[gdir-2] = int(grid[gdir-2]/kp_mult)
    except: pass

    inputDict["K_POINTS"]["__content__"] = " ".join(map(str,grid))

    inputDict["&control"]["nppstr"] = "%s"%grid[gdir-1]

    inputString = AFLOWpi.retr._joinInput(inputDict)
    ID_gdir = ID+"_gdir%s"%gdir

    ifn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID_gdir)

    if not os.path.exists(ifn):
        with open(ifn,'w') as newIn:
            newIn.write(inputString)

    return oneCalc,ID_gdir

def setup_berry(calcs,kp_factor):
    add = "oneCalc,ID_gdir1 = AFLOWpi.prep._prep_berry(oneCalc,ID,1,%s)"%kp_factor
    AFLOWpi.prep.addToAll_(calcs,block='PREPROCESSING',addition=add)
    AFLOWpi.run._skeletonRun(calcs,calcType="gdir1",execPath='"./pw.x"',execPostfix="-northo 1")
    add = "        oneCalc,ID_gdir2 = AFLOWpi.prep._prep_berry(oneCalc,ID,2,%s)"%kp_factor
    AFLOWpi.prep.addToAll_(calcs,block='RUN',addition=add)
    AFLOWpi.run._skeletonRun(calcs,calcType="gdir2",execPath='"./pw.x"',execPostfix="-northo 1")
    add = "        oneCalc,ID_gdir3 = AFLOWpi.prep._prep_berry(oneCalc,ID,3,%s)"%kp_factor
    AFLOWpi.prep.addToAll_(calcs,block='RUN',addition=add)
    AFLOWpi.run._skeletonRun(calcs,calcType="gdir3",execPath='"./pw.x"',execPostfix="-northo 1")


def _pull_pol_berry(oneCalc,ID):

    afd=oneCalc["_AFLOWPI_FOLDER_"]

    direc=[1,2,3]

    tmp_pol=[]
    pdir=[]
    for i in range(len(direc)):
        fn=os.path.join(afd,"%s_gdir%s.out"%(ID,direc[i]))
        with open(fn) as ifo:
            ifs=ifo.read()

        pdir.append(list(map(float,re.findall(r"The polarization direction is:  \((.*)\)",ifs)[0].split(","))))
        tmp_pol.append(float(re.findall("P =\s*([-\d.]+)",ifs)[2]))

    tmp_pol = np.array(tmp_pol)
    pdir    = np.array(pdir)

    tmp_pol=pdir.T @ tmp_pol

    dst_nm="_".join(ID.split("_")[:2])+"_pol"
    savep=os.path.join(os.path.dirname(afd),dst_nm)

    with open(savep,"w") as ifo:
        ifs=ifo.write("% 16.8f % 16.8f % 16.8f\n"%(tmp_pol[0],tmp_pol[1],tmp_pol[2]))
    

def _pull_stress_piezo(oneCalc,ID):

    #pull stress from output
    stress=AFLOWpi.retr.getStress(oneCalc,ID)

    #convert stribg to 6 element array
    stress=np.array([list(map(float,i.split()[3:])) for i in stress.split("\n")[3:-1]])

    #convert from kbar to N/m^2
    stress*=1.e8


    afd=oneCalc["_AFLOWPI_FOLDER_"]


    fn=os.path.join(afd,"%s.out"%(ID))
    with open(fn) as ifo:
        ifs=ifo.read()


    dst_nm="_".join(ID.split("_")[:2])+"_stress"
    savep=os.path.join(os.path.dirname(afd),dst_nm)

    with open(savep,"w") as ifo:
        ifs=ifo.write("% 12.8e % 12.8e % 12.8e % 12.8e % 12.8e % 12.8e\n"%(stress[0,0],stress[1,1],stress[2,2],stress[1,2],stress[0,2],stress[0,1]))

def _read_piezo_dat(oneCalc,ID):

    afd=oneCalc["_AFLOWPI_FOLDER_"]

    fil = sorted(glob.glob("%s/%s_ELASTIC/Dst*pol"%(afd,ID)))

    dinfo=[]
    cvl=[]

    for i in range(len(fil)):
        fn=os.path.basename(fil[i])[3:].split("_")[:2]
        dinfo.append(list(map(int,fn)))

        dat=np.loadtxt(fil[i])
        cvl.append(dat)


    dinfo=np.array(dinfo)-1
    cvl=np.array(cvl)

    ndist=np.unique(dinfo[:,0]).size
    diter=np.unique(dinfo[:,1]).size


    cvl_sort=np.zeros((ndist,diter,3))

    for i in range(cvl.shape[0]):
        cvl_sort[dinfo[i,0],dinfo[i,1]]=cvl[i]

    return cvl_sort

def _read_piezo_stress_dat(oneCalc,ID):

    afd=oneCalc["_AFLOWPI_FOLDER_"]

    fil = sorted(glob.glob("%s/%s_ELASTIC/Dst*stress"%(afd,ID)))

    dinfo=[]
    cvl=[]

    for i in range(len(fil)):
        fn=os.path.basename(fil[i])[3:].split("_")[:2]
        dinfo.append(list(map(int,fn)))

        dat=np.loadtxt(fil[i])
        cvl.append(dat)


    dinfo=np.array(dinfo)-1
    cvl=np.array(cvl)

    ndist=np.unique(dinfo[:,0]).size
    diter=np.unique(dinfo[:,1]).size


    cvl_sort=np.zeros((ndist,diter,6))

    for i in range(cvl.shape[0]):
        cvl_sort[dinfo[i,0],dinfo[i,1]]=cvl[i]


    return cvl_sort

    

def _calc_piezo_tensor(oneCalc,ID):
    pols   = AFLOWpi.prep._read_piezo_dat(oneCalc,ID)
    stress = AFLOWpi.prep._read_piezo_stress_dat(oneCalc,ID)
    print(pols)
    print(stress)
    order=7
    if pols.shape[1]<6:
        order=1
    elif pols.shape[1]<8:
        order=3
    elif pols.shape[1]<10:
        order=5

    stress/=1.e8

    res=np.zeros((pols.shape[0],6,3))

    for s_dir in range(6):
        for p_dir in range(3):
            for dst in range(pols.shape[0]):
                res[dst,s_dir]=np.polyfit(stress[dst,:,s_dir],pols[dst,:,p_dir],order)[-2]


    print(res)
    
    
