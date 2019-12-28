import numpy as np
import os
import re 
import AFLOWpi
import scipy.linalg

def load_dyn(fn,nat,nq):

    dyn = np.zeros((nq,3,3,nat,nat),dtype=complex)

    with open(fn,"r") as fo:
        fs=fo.readlines()

    counter=-1.0
    tcounter=0
    for line in fs:
        try:
            try:
                dat = list(map(float,line.split()))
            except:
                counter+=0.5                
            if len(dat)==2:
                n=int(dat[0])-1
                m=int(dat[1])-1
                            
            ipol = int(tcounter)%3

            dyn[int(counter),:,ipol,n,m]=dat[0]+dat[1]*1.0j,dat[2]+dat[3]*1.0j,dat[4]+dat[5]*1.0j
            tcounter+=1
        
        except Exception as e: pass
        if int(counter)==nq:
            break

#    dyn = np.transpose(dyn,axes=(0,1,3,2,4))
    dyn = np.transpose(dyn,axes=(0,3,2,4,1))
    dyn=dyn.reshape((nq,nat*3,nat*3),order="C")

    
    return dyn


def _delta_dyn(oneCalc,ID,dat_type="DOS"):    

    norm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='phonon')
    expn_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_up')
    cont_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_dn')
    cent_diff=True
    if cont_ID is None:
        cont_ID = norm_ID
        cent_diff=False

    expn_vol_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_relax_up')
    cont_vol_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_relax_dn')
    if cont_vol_ID is None:
        cont_vol_ID = norm_ID
    
    norm_vol = AFLOWpi.retr.getCellVolume(oneCalc,norm_ID,string=False,conventional=False)
    cont_vol = AFLOWpi.retr.getCellVolume(oneCalc,cont_vol_ID,string=False,conventional=False)
    expn_vol = AFLOWpi.retr.getCellVolume(oneCalc,expn_vol_ID,string=False,conventional=False)

    dV = (expn_vol-cont_vol)/norm_vol
    nat = int(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['nat'])

    
    data_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.ph%s.gp'%(norm_ID,dat_type))
    tdata= np.loadtxt(data_file_name,usecols=(1,2))
    nq=tdata.shape[0]

    tdata = None





    fn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.ph%s.dyn"%(norm_ID,dat_type))
    print("loading dynmat for V")
    dyn = AFLOWpi.retr.load_dyn(fn,nat,nq)

    fn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.ph%s.dyn"%(expn_ID,dat_type))
    print("loading dynmat for V+dV")
    delDq = AFLOWpi.retr.load_dyn(fn,nat,nq)


    if cont_ID!=norm_ID:
        fn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.ph%s.dyn"%(cont_ID,dat_type))
        print("loading dynmat for V-dV")
        delDq -= AFLOWpi.retr.load_dyn(fn,nat,nq)
    else:
        delDq -= dyn


    vq   = np.zeros_like(dyn)
    eig  = np.zeros((dyn.shape[0],nat*3),dtype=float)
    grun = np.zeros_like(eig)

    n,n = np.diag_indices(nat*3)



    for i in range(dyn.shape[0]):
        eig[i],vq[i] = scipy.linalg.eigh(dyn[i])
        grun[i]      = -1.0*np.real(np.conj(vq[i].T).dot((delDq[i]/dV).dot(vq[i])))[n,n]

    
    


    PF = 1.0/(2.0*eig)
    PF_inf_ind = np.where(np.isinf(PF**2))
    PF[PF_inf_ind] = 0.0



    
    grun      =  np.nan_to_num(grun*PF)


    grun_i = np.mean(grun[:,:3]**2,axis=0)**(0.5)


    if dat_type=="DOS":
        return np.ravel(grun,order="C"),grun_i

    else:
        fn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.phSCATTER.band"%norm_ID)
        np.savetxt(fn,grun)
