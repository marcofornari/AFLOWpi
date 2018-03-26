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
import re
import numpy as np
import scipy.integrate
import collections
import os


def _increase_celldm1(oneCalc,ID,amount):
    if ID in oneCalc['prev']:
        return oneCalc,ID

    #change prefix in input file so it doesn't conflict with later calculations after thermal

    if amount > 1.0:
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','prefix','"%s_volup"'%oneCalc['_AFLOWPI_PREFIX_'])
        oneCalc["_AFLOWPI_PREFIX_"]=oneCalc["_AFLOWPI_PREFIX_"]+"_volup"
    else:
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','prefix','"%s_voldn"'%oneCalc['_AFLOWPI_PREFIX_'])
        oneCalc["_AFLOWPI_PREFIX_"]=oneCalc["_AFLOWPI_PREFIX_"]+"_voldn"


    inp_file = AFLOWpi.retr._getInputFileString(oneCalc,ID)
    celldm1 = float(AFLOWpi.retr._splitInput(inp_file)['&system']['celldm(1)'])

    new_celldm1=celldm1*amount


    oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','celldm(1)',new_celldm1)

    return oneCalc,ID


def _get_gruneisen_ap(oneCalc,ID):
    norm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='phonon')
    cont_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_dn')
    expn_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_up')

    if cont_ID == None:
        cont_ID = norm_ID

    cont_vol_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_relax_dn')
    if cont_vol_ID == None:
        cont_vol_ID = norm_ID
    
    expn_vol_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_relax_up')

    norm_vol = AFLOWpi.retr.getCellVolume(oneCalc,norm_ID,string=False,conventional=False)
    cont_vol = AFLOWpi.retr.getCellVolume(oneCalc,cont_vol_ID,string=False,conventional=False)
    expn_vol = AFLOWpi.retr.getCellVolume(oneCalc,expn_vol_ID,string=False,conventional=False)

    bohr2meter=5.29177e-11

    extension="eig.ap"

    norm_freq,q_point_old,labels =AFLOWpi.retr._get_ph_dos_data_ap(oneCalc,norm_ID)
    cont_freq,q_point_old,labels =AFLOWpi.retr._get_ph_dos_data_ap(oneCalc,cont_ID)
    expn_freq,q_point_old,labels =AFLOWpi.retr._get_ph_dos_data_ap(oneCalc,expn_ID)





    # de=0.5
    # dmax = np.amax(norm_freq)
    # dmin = np.amin(norm_freq)
    # nbin = np.ceil(dmax/de+1.0)
    # freqs= np.linspace(dmin,dmax,nbin)
    # dos_data = np.zeros((freqs.shape[0]))
    # data = np.ravel(norm_freq)



    # for w in xrange(freqs.shape[0]):
    #     dos_data[w] = np.sum(np.exp(-((freqs[w]-data)/de)**2))/float(norm_freq.shape[0])

    # print np.cumsum(dos_data)

    # raise SystemExit

    grun=[]
    q_point=[]
    omega=[]

    grun = np.zeros((norm_freq.shape[0],norm_freq.shape[1]+1),dtype=float)

    grun[:,0] = norm_freq[:,0]



    grun[:,1],grun_i = AFLOWpi.retr._delta_dyn(oneCalc,ID)
#    grun[:,1],_ = AFLOWpi.retr._get_gruneisen(oneCalc,ID,band=False)


    grun[:,2:] = grun[:,1][:,None]*norm_freq[:,1:]

    labels = np.asarray(labels[3:])

    ap_grun = np.zeros((grun.shape[0],np.unique(labels).shape[0]+1),dtype=float)


    ap_grun[:,0] = grun[:,0]
    ap_grun[:,1] = grun[:,1]

    unique_labels = np.unique(labels[1:])

    for lab in xrange(1,unique_labels.shape[0]+1):
        lab_val = np.unique(labels)[lab-1]
        lab_ind = np.where(labels==lab_val)[0]
        ap_grun[:,lab+1] = np.sum(norm_freq[:,lab_ind],axis=1)*ap_grun[:,1]
        

    unique_labels = unique_labels.tolist()

    temp = ['Freq','Total']
    temp.extend(unique_labels)
    unique_labels=temp


    grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER.ap'%ID)
    np.savetxt(grun_file_name,ap_grun,fmt='% 4.8f',header=' '.join(unique_labels))
    
    return grun_i

        
###################################################################################################
###################################################################################################    
###################################################################################################

def _get_gruneisen(oneCalc,ID,band=True):
    norm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='phonon')
    expn_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_up')

    cont_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_dn')
    if cont_ID == None:
        cont_ID = norm_ID

    expn_vol_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_relax_up')
    cont_vol_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_relax_dn')
    if cont_vol_ID == None:
        cont_vol_ID = norm_ID
    




    norm_vol = AFLOWpi.retr.getCellVolume(oneCalc,norm_ID,string=False,conventional=False)
    cont_vol = AFLOWpi.retr.getCellVolume(oneCalc,cont_vol_ID,string=False,conventional=False)
    expn_vol = AFLOWpi.retr.getCellVolume(oneCalc,expn_vol_ID,string=False,conventional=False)

    bohr2meter=5.29177e-11
 
    if band==True:
       print band
       raise SystemExit
       extension='phBAND.gp'
    else:
        extension="phDOS.gp"


    norm_freq,q_point_old = AFLOWpi.retr._get_ph_dos_data(oneCalc,norm_ID,extension=extension)
    cont_freq,q_point_old = AFLOWpi.retr._get_ph_dos_data(oneCalc,cont_ID,extension=extension)
    expn_freq,q_point_old = AFLOWpi.retr._get_ph_dos_data(oneCalc,expn_ID,extension=extension)

    grun=[]
    q_point=[]
    omega=[]


    PF = np.nan_to_num(1.0*norm_vol/norm_freq)
#    PF = np.nan_to_num((expn_vol+cont_vol)/(expn_freq+cont_freq))

    PF_inf_ind = np.where(np.isinf(PF**2))
    PF[PF_inf_ind] = 0.0
    
    grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'PF.dat')
    np.savetxt(grun_file_name,np.concatenate([np.ravel(norm_freq,order="C")[:,None],np.ravel(PF,order="C")[:,None]],axis=1),fmt='% 4.8f')




    grun = -(expn_freq-cont_freq)/(expn_vol-cont_vol) * PF
#    grun = -(expn_freq-cont_freq)/(expn_vol-cont_vol) * 1.0

    av = np.sqrt(np.mean(grun**2,axis=0))


    grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER_old.gp'%ID)
    np.savetxt(grun_file_name,np.concatenate([np.ravel(norm_freq,order="C")[:,None],np.ravel(grun,order="C")[:,None]],axis=1),fmt='% 4.8f')
    grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER_0.gp'%ID)
    np.savetxt(grun_file_name,np.concatenate([norm_freq[:,0:1],grun[:,0:1]**2],axis=1),fmt='% 4.8f')
    grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER_1.gp'%ID)
    np.savetxt(grun_file_name,np.concatenate([norm_freq[:,1:2],grun[:,1:2]**2],axis=1),fmt='% 4.8f')
    grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER_2.gp'%ID)
    np.savetxt(grun_file_name,np.concatenate([norm_freq[:,2:3],grun[:,2:3]**2],axis=1),fmt='% 4.8f')

    
    
    

    theta = np.amax(norm_freq[:,:3],axis=0)*1.2398e-4*1.16045221e4


    return np.ravel(grun,order="C"),theta

###################################################################################################
###################################################################################################    
###################################################################################################

def _get_ph_dos_data(oneCalc,ID,extension='phBAND.gp',postfix=''):

    data_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s%s.%s'%(ID,postfix,extension))


    #data =np.loadtxt(data_file_name,dtype=np.float64,)
    data = []

    with open(data_file_name,'r') as fo:
        fs=fo.read()
    fs=fs.split('\n')

    for line in fs:
        if len(line.strip())!=0:
            try:
                data.append(map(float,line.split()[1:]))
            except: pass
    data = np.asarray(data)


    return data,None


def _get_ph_dos_data_ap(oneCalc,ID,postfix=''):

    extension='eig.ap'

    data_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s%s.%s'%(ID,postfix,extension))


    #data =np.loadtxt(data_file_name,dtype=np.float64,)
    data = []

    with open(data_file_name,'r') as fo:
        fs=fo.read()
    fs=fs.split('\n')
    labels=fs[0].split()
    fs=fs[1:]
    for line in fs:
        if len(line.strip())!=0:
            dat_temp = map(float,line.split())
            temp_one = [dat_temp[3]]
            temp_one.extend(dat_temp[4:])
            data.append(temp_one)
    data = np.asarray(data)

    ret_dat=np.zeros(data.shape)

    ret_dat[:,0]=data[:,0]
    for i in range(1,ret_dat.shape[1]):
        ret_dat[:,i] = data[:,i]

    return ret_dat,ret_dat[1],labels


def _get_ph_band_data(oneCalc,ID,extension='phBAND.gp',postfix=''):

    data_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s%s.%s'%(ID,postfix,extension))


    #data =np.loadtxt(data_file_name,dtype=np.float64,)
    data = []

    with open(data_file_name,'r') as fo:
        fs=fo.read()
    fs=fs.split('\n')
    labels=fs[0]
    fs=fs[1:]
    for line in fs:
        if len(line.strip())!=0:
            data.append(map(float,line.split()))
    data = np.asarray(data)

    return data[:,1:],data[:,0]



def _get_debye_temp(v_i,cell_vol,nat):

    #convert to meters
    bohr2meter=5.29177e-11
    V=cell_vol*bohr2meter**3.0

    #get num atoms in cell


    #some constants
    h_bar=1.0545718*10**-34 
    k_b=1.38064852*10**-23

    #get v_debye for each branch
    v_s_TA       = v_i[0]
    v_s_TA_prime = v_i[1]
    v_s_LA       = v_i[2]


    #calculate debye temperature for each branch
    wo_v           = h_bar/k_b#*(6.0*np.pi**2.0*nat/V)**(1.0/3.0)
    debye_TA       = wo_v*v_s_TA
    debye_TA_prime = wo_v*v_s_TA_prime
    debye_LA       = wo_v*v_s_LA

    return debye_TA,debye_TA_prime,debye_LA


def _therm_pp(__submitNodeName__,oneCalc,ID,run_matdyn=True):


    grun_i=[0.0,0.0,0.0]

    #for plotting

    #for DC model
    _,theta_i = AFLOWpi.retr._get_gruneisen(oneCalc,ID,band=False)


    grun_i = AFLOWpi.retr._get_gruneisen_ap(oneCalc,ID)
#    raise SystemExit

    norm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='phonon')

    v_i = AFLOWpi.run.do_sound_velocity(__submitNodeName__,oneCalc,norm_ID,dk_theta=0.1,dk_phi=0.2,dk_r=0.0125,
                                        r_max=0.05,theta_range=[-np.pi/2.0,np.pi/2.0],phi_range=[0.0,2.0*np.pi],
                                        origin=[0.0,0.0,0.0],nspin=1,kpi=0,read_S=False,shift=0.0,run_matdyn=run_matdyn)


    #get volume of original cell
    cell_vol = AFLOWpi.retr.getCellVolume(oneCalc,norm_ID,string=False,conventional=False)
    #convert to meters
    bohr2meter=5.29177e-11
    V=cell_vol*bohr2meter**3.0

    #get num atoms in cell
    N = float(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['nat'])




    cell_mass = AFLOWpi.retr._get_cell_mass(oneCalc,ID)
    
    #convert amu to kg
    M=cell_mass*1.66054e-27
    Mass=M/float(N)
    Vol=V/float(N)

    TEMP = np.linspace(5.0,1205.0,1201.0,endpoint=True)

    therm_cond_data_str="T       Total        TA           TA'          LA"    

    for T in TEMP:
        total,TA_cont,TA_prime_cont,LA_cont = AFLOWpi.retr._do_therm(v_i,theta_i,grun_i,Mass,Vol,T)
        
        therm_cond_data_str+='\n%7.1f %12.3f %12.3f %12.3f %12.3f'%(T,total,TA_cont,TA_prime_cont,LA_cont)

    therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_thermal_cond.dat'%ID)
    with open(therm_file_name,'w') as tcfo:
        tcfo.write(therm_cond_data_str) 

######################################################################################################################                        

    theta_a = [i*N**(-1.0/3.0) for i in theta_i]
    
    therm_stat_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_therm_stats.dat'%ID)
    therm_tex_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_therm_stats.tex'%ID)
    comp_name=AFLOWpi.retr._getStoicName(oneCalc,strip=True)
    total,TA_cont,TA_prime_cont,LA_cont = AFLOWpi.retr._do_therm(v_i,theta_i,grun_i,Mass,Vol,300.0)
    try:
        with open(therm_tex_file,"w") as ofo:
            ofo.write("%5.5s % 12.3f& % 12.3f& % 12.3f& % 12.3f& % 12.3f& % 12.3f& % 12.3f& % 12.3f& % 12.3f& % 12.3f&% 12.3f& % 12.3f& % 12.3f& % 12.3f& % 12.3f & % 12.3f \\ \n"%\
                          (comp_name,grun_i[0],grun_i[1],grun_i[2],(grun_i[0]+grun_i[1]+grun_i[2])/3,
                           v_i[0],v_i[1],v_i[2],(v_i[0]+v_i[1]+v_i[2])/3,
                           theta_i[0],theta_i[1],theta_i[2],(theta_i[0]+theta_i[1]+theta_i[2])/3,
                           TA_cont,TA_prime_cont,LA_cont,total))
    except Exception,e:
        print e
    with open(therm_stat_file,"w") as ofo:
        ofo.write("--------------------------------------------------------------\n")
        ofo.write("%s\n"%comp_name)
        ofo.write("--------------------------------------------------------------\n")
        ofo.write("                  TA           TA'            LA           AVG\n")
        ofo.write("--------------------------------------------------------------\n")
        ofo.write("Grun       % 9.3f     % 9.3f     % 9.3f     % 9.3f\n"%(grun_i[0],grun_i[1],grun_i[2],(grun_i[0]+grun_i[1]+grun_i[2])/3))
        ofo.write("v_s        % 9.3f     % 9.3f     % 9.3f     % 9.3f\n"%(v_i[0],v_i[1],v_i[2],(v_i[0]+v_i[1]+v_i[2])/3))
        ofo.write("T_Debye    % 9.3f     % 9.3f     % 9.3f     % 9.3f\n"%(theta_i[0],theta_i[1],theta_i[2],(theta_i[0]+theta_i[1]+theta_i[2])/3))
#        ofo.write("T_Debye_a  % 9.3f     % 9.3f     % 9.3f     % 9.3f\n"%(theta_a[0],theta_a[1],theta_a[2],(theta_a[0]+theta_a[1]+theta_a[2])/3))
        ofo.write("--------------------------------------------------------------\n")
        ofo.write("                  TA           TA'            LA           TOT\n")
        ofo.write("C_v @ 300K % 9.3f     % 9.3f     % 9.3f     % 9.3f\n"%(TA_cont,TA_prime_cont,LA_cont,total))
        ofo.write("--------------------------------------------------------------\n")






######################################################################################################################                        
######################################################################################################################                        

def _do_therm(v_i,theta_i,grun_i,Mass,Vol,T):
    ######################################################################################################################
    def calc_tau_N_TA(x,grun,velocity,vol,mass,T):
        #Relax Temp Scattering TA'                                                   
        h_bar = np.float64(1.0545718*(10.0**(-34.0))) # hbar
        k_b   = np.float64(1.38064852*(10.0**(-23.0))) #boltzmann                                               

        # k_b/h_bar done this way to avoid overflow errors 
        one_over_tau_N_TA  = k_b*(k_b/h_bar)**4.0
        one_over_tau_N_TA *= grun**2.0
        one_over_tau_N_TA *= vol/mass
        one_over_tau_N_TA /= velocity**5.0

        one_over_tau_N_TA *= x
        one_over_tau_N_TA *= T**5.0

        tau_N_TA  = 1.0/one_over_tau_N_TA

        return tau_N_TA
    ######################################################################################################################
    def calc_tau_N_LA(x,grun,velocity,vol,mass,T):
        #Relax Temp Scattering LA                                                                                                         
        h_bar = np.float64(1.0545718*(10.0**(-34.0))) # hbar
        k_b   = np.float64(1.38064852*(10.0**(-23.0))) #boltzmann                                                      
 
        one_over_tau_N_LA  = k_b*(k_b/h_bar)**4.0
        one_over_tau_N_LA *= grun**2.0
        one_over_tau_N_LA *= vol/mass
        one_over_tau_N_LA /= velocity**5.0

        one_over_tau_N_LA *= x**2.0
        one_over_tau_N_LA *= T**5.0

        tau_N_LA  = 1.0/one_over_tau_N_LA

        return tau_N_LA
    ######################################################################################################################
    def calc_tau_U(x,grun,velocity,debye_temp,mass,T):
        #Relax Temp Umklamp LA                                                                                     
        h_bar = np.float64(1.0545718*(10.0**(-34.0))) # hbar
        k_b   = np.float64(1.38064852*(10.0**(-23.0))) #boltzmann                                                                

        one_over_tau_U = k_b*(k_b/h_bar)
        one_over_tau_U *= grun**2.0
        one_over_tau_U /= mass
        one_over_tau_U /= velocity**2.0
        one_over_tau_U /= debye_temp


        one_over_tau_U *= x**2.0
        one_over_tau_U *= T**3.0

        one_over_tau_U *= np.exp(-1.0*debye_temp/(3.0*T))

        tau_U = 1.0/one_over_tau_U

        return tau_U
    ######################################################################################################################
    def first_int(x,grun,velocity,debye_temp,vol,mass,T,trans):
        #do first integral in DB model
        tau_U=calc_tau_U(x,grun,velocity,debye_temp,mass,T)
        if trans:
            tau_N=calc_tau_N_TA(x,grun,velocity,vol,mass,T)
        else:
            tau_N=calc_tau_N_LA(x,grun,velocity,vol,mass,T)

        Tc = 1.0/(1.0/tau_U + 1.0/tau_N)

        sol  = Tc
        sol *= x**4.0
        sol *= np.exp(x)/(np.exp(x)-1.0)**2.0


        return sol
    ######################################################################################################################      
    def second_int(x,grun,velocity,debye_temp,vol,mass,T,trans):
        #do second integral in DB model
        tau_U=calc_tau_U(x,grun,velocity,debye_temp,mass,T)
        if trans:
            tau_N=calc_tau_N_TA(x,grun,velocity,vol,mass,T)
        else:
            tau_N=calc_tau_N_LA(x,grun,velocity,vol,mass,T)

        Tc = 1.0/(1.0/tau_U + 1.0/tau_N)

        #Tc/tau_N
        Tc_mod = tau_U/(tau_U + tau_N)


        sol  = Tc_mod
        sol *= x**4.0
        sol *= np.exp(x)/(np.exp(x)-1.0)**2.0


        return sol
    ######################################################################################################################
    def third_int(x,grun,velocity,debye_temp,vol,mass,T,trans):
        #do third integral in DB model
        tau_U=calc_tau_U(x,grun,velocity,debye_temp,mass,T)
        if trans:
            tau_N=calc_tau_N_TA(x,grun,velocity,vol,mass,T)
        else:
            tau_N=calc_tau_N_LA(x,grun,velocity,vol,mass,T)


        #Tc/(tau_U*tau_N)
        Tc_mod = 1.0/(tau_U + tau_N)

        sol  = Tc_mod
        sol *= x**4.0
        sol *= np.exp(x)/(np.exp(x)-1.0)**2.0

        return sol



    ######################################################################################################################
    #define constants
    h_bar = np.float64(1.0545718*(10.0**(-34.0))) # hbar
    k_b   = np.float64(1.38064852*(10.0**(-23.0))) #boltzmann  

    #define a constant from other constants
    C_PHON_CONST  = 1.0/3.0
    C_PHON_CONST *= k_b*(k_b/h_bar)**3.0
    C_PHON_CONST *= T**3.0
    C_PHON_CONST /= 2.0*np.pi**2.0

    #speed of sound for each accoustic branch
    vel_TA    = v_i[0]
    vel_TA1   = v_i[1]
    vel_LA    = v_i[2]
    #scaling constants for each accoustic branch
    C_TA      = C_PHON_CONST/vel_TA
    C_T1      = C_PHON_CONST/vel_TA1
    C_LA      = C_PHON_CONST/vel_LA
    #debye temp for each accoustic branch
    theta_TA  = theta_i[0]
    theta_TA1 = theta_i[1]
    theta_LA  = theta_i[2]
    #gruneisen parameter for each accoustic branch
    grun_TA   = grun_i[0]
    grun_TA1  = grun_i[1]
    grun_LA   = grun_i[2]

    #upper limit of integration for a given T for the integrals for each accoustic branch
    max_TA    = theta_TA/T
    max_TA1   = theta_TA1/T
    max_LA    = theta_LA/T
    
    
    #calculate the three integrals for TA phonon and find lattice k for TA
    TA_1 =  scipy.integrate.quad(first_int, 0.00,max_TA,(grun_TA,vel_TA,theta_TA,Vol,Mass,T,True),limit=100)[0]
    TA_2 =  scipy.integrate.quad(second_int,0.00,max_TA,(grun_TA,vel_TA,theta_TA,Vol,Mass,T,True),limit=100)[0]
    TA_3 =  scipy.integrate.quad(third_int, 0.00,max_TA,(grun_TA,vel_TA,theta_TA,Vol,Mass,T,True),limit=100)[0]
    klattice_TA  = C_TA * (TA_1 + TA_2*(TA_2/TA_3))           

    #calculate the three integrals for TA' phonon and find lattice k for TA'
    T1_1 = scipy.integrate.quad(first_int, 0.00,max_TA1,(grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True),limit=100)[0]
    T1_2 = scipy.integrate.quad(second_int,0.00,max_TA1,(grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True),limit=100)[0]
    T1_3 = scipy.integrate.quad(third_int, 0.00,max_TA1,(grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True),limit=100)[0]
    klattice_T1  = C_T1 * (T1_1 + T1_2*(T1_2/T1_3))

    #calculate the three integrals for LA phonon and find lattice k for LA
    LA_1 =  scipy.integrate.quad(first_int, 0.00,max_LA,(grun_LA,vel_LA,theta_LA,Vol,Mass,T,False),limit=100)[0]
    LA_2 =  scipy.integrate.quad(second_int,0.00,max_LA,(grun_LA,vel_LA,theta_LA,Vol,Mass,T,False),limit=100)[0]
    LA_3 =  scipy.integrate.quad(third_int, 0.00,max_LA,(grun_LA,vel_LA,theta_LA,Vol,Mass,T,False),limit=100)[0]
    klattice_LA  =  C_LA * (LA_1 + LA_2*(LA_2/LA_3))

    # num_samples=1000
    # samples = np.linspace(1.e-8,max_TA,num_samples)
    # ##################################################################################                 
    # #calculate the three integrals for TA phonon and find lattice k for TA
    # TA_1 =  np.sum( first_int(samples,grun_TA,vel_TA,theta_TA,Vol,Mass,T,True))*max_TA/num_samples
    # TA_2 =  np.sum(second_int(samples,grun_TA,vel_TA,theta_TA,Vol,Mass,T,True))*max_TA/num_samples
    # TA_3 =  np.sum( third_int(samples,grun_TA,vel_TA,theta_TA,Vol,Mass,T,True))*max_TA/num_samples
    # klattice_TA  = C_TA * (TA_1 + TA_2*(TA_2/TA_3))           
    # ##################################################################################                
    # #calculate the three integrals for TA1 phonon and find lattice k for TA1
    # TA1_1 =  np.sum( first_int(samples,grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True))*max_TA1/num_samples
    # TA1_2 =  np.sum(second_int(samples,grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True))*max_TA1/num_samples
    # TA1_3 =  np.sum( third_int(samples,grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True))*max_TA1/num_samples
    # klattice_T1  = C_T1 * (TA1_1 + TA1_2*(TA1_2/TA1_3))
    # ##################################################################################                
    # #calculate the three integrals for LA phonon and find lattice k for LA
    # LA_1 =  np.sum( first_int(samples,grun_LA,vel_LA,theta_LA,Vol,Mass,T,True))*max_LA/num_samples
    # LA_2 =  np.sum(second_int(samples,grun_LA,vel_LA,theta_LA,Vol,Mass,T,True))*max_LA/num_samples
    # LA_3 =  np.sum( third_int(samples,grun_LA,vel_LA,theta_LA,Vol,Mass,T,True))*max_LA/num_samples
    # klattice_LA  =  C_LA * (LA_1 + LA_2*(LA_2/LA_3))
    # ##################################################################################                
    
    #sum over all contibutions for lattice k at Temp T 
    klattice_total  = 0.0
    klattice_total += klattice_TA
    klattice_total += klattice_T1
    klattice_total += klattice_LA

    return klattice_total,klattice_TA,klattice_T1,klattice_LA
