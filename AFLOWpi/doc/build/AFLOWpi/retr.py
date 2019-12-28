import AFLOWpi
import re
import numpy
import scipy.integrate
import collections
import os
from functools import reduce


def _increase_celldm1(oneCalc,ID,amount):
    if ID in oneCalc['prev']:
        return oneCalc,ID

    #change prefix in input file so it doesn't conflict with later calculations after thermal

    oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&control','prefix','"%s_vol"'%oneCalc['_AFLOWPI_PREFIX_'])
    oneCalc["_AFLOWPI_PREFIX_"]=oneCalc["_AFLOWPI_PREFIX_"]+"_vol"

    inp_file = AFLOWpi.retr._getInputFileString(oneCalc,ID)
    celldm1 = float(AFLOWpi.retr._splitInput(inp_file)['&system']['celldm(1)'])

    new_celldm1=celldm1*amount


    oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','celldm(1)',new_celldm1)

    return oneCalc,ID


def _get_gruneisen(oneCalc,ID,band=True):
    norm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='phonon')
    expn_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal')

    expn_vol_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='thermal_relax')

    norm_vol = AFLOWpi.retr.getCellVolume(oneCalc,norm_ID,string=False,conventional=False)
    expn_vol = AFLOWpi.retr.getCellVolume(oneCalc,expn_vol_ID,string=False,conventional=False)

    bohr2meter=5.29177e-11
    norm_vol*=bohr2meter**3.0
    expn_vol*=bohr2meter**3.0


#    print expn_vol/norm_vol
#    print 
 
    if band==True:
       print(band)
       raise SystemExit
       extension='phBAND.gp'
    else:
#        extension=''
        extension='eig.ap'

    norm_freq,q_point_old = AFLOWpi.retr._get_ph_dos_data(oneCalc,norm_ID,extension=extension)
    expn_freq,q_point_old = AFLOWpi.retr._get_ph_dos_data(oneCalc,expn_ID,extension=extension)



    grun=[]
    q_point=[]
    omega=[]

    for i in range(len(norm_freq)):
        if band:
            q_point.append(q_point_old[i])
        for j in range(len(norm_freq[i])):
            try:
                    deriv  = (expn_freq[i][j]-norm_freq[i][j])/(expn_vol-norm_vol)
                    deriv *= -1.0*norm_vol/norm_freq[i][j]
                    if not numpy.isnan(deriv) and not numpy.isinf(deriv):
                        try:
                            grun[i].append(deriv)
                        except Exception as e:
                            grun.append([])
                            grun[i].append(deriv)
                        try:
                            omega[i].append(norm_freq[i][j])
                        except Exception as e:
                            omega.append([])
                            omega[i].append(norm_freq[i][j])
                    else:

                        try:
                            grun[i].append(0.0)
                        except Exception as e:
                            grun.append([])
                            grun[i].append(0.0)
                        try:
                            omega[i].append(norm_freq[i][j])
                        except Exception as e:
                            omega.append([])
                            omega[i].append(norm_freq[i][j])
                        continue
            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)

#                print e


#    print len(grun)
#    print len(q_point_old)
#    print len(q_point)


#   raise SystemExit
    if band:
        grun_data_str='%s          %s                  %s                 %s'%('q','TA',"TA'",'LA')
        for i in range(len(grun)):
            try:
                if len(grun[i]) == len(grun[12]):
                    grun_data_str+='\n%10.8f '%q_point[i]

                    for j in range(len(grun[i])):
                        grun_data_str+='%16.16f '%grun[i][j]
            except:
                continue
    else:
        grun_data_str=""
        for i in range(len(grun)):
            try:
                for j in range(len(grun[i])):
#                    print grun[i][j]
#                    print grun[i][j]**2.0
                    grun_data_str+='%16.16f %16.16f\n'%(omega[i][j],grun[i][j]**2.0)
            except:
                continue
 #   if len(grun)==3:
 #       pass
    # else:
    #     for i in range(len(grun[0][3:-1])):
    #         grun_data_str+='Optical             '
    #     for i in range(len(q_point)):
    #         grun_data_str+='\n%10.8f %16.16f %16.16f  %16.16f'%(q_point[i],grun[0][i],grun[1][i],grun[2][i])
    #         for j in range(len(grun[i][3:-1])):
    #             grun_data_str+='%16.16f'%(grun[i][j])
        

    if band:
        grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phGRUN.gp'%ID)
    else:
        grun_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phSCATTER.gp'%ID)

    with open(grun_file_name,'w') as gpfo:
        gpfo.write(grun_data_str)
    
    av_TA       = sum(grun[0])/len(grun[0])
    av_TA_prime = sum(grun[1])/len(grun[1])
    av_LA      = sum(grun[2])/len(grun[2])
        
    return [av_TA,av_TA_prime,av_LA]
        
    
            

#    frequency = 
#    delta_vol_frequency= 
def _get_ph_dos_data(oneCalc,ID,extension='phBAND.gp',postfix=''):

    data_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s%s.%s'%(ID,postfix,extension))


    #data =numpy.loadtxt(data_file_name,dtype=numpy.float64,)
    data = []

    with open(data_file_name,'r') as fo:
        fs=fo.read()
    fs=fs.split('\n')
    labels=fs[0]
    fs=fs[1:]
    for line in fs:
        if len(line.strip())!=0:
            dat_temp = list(map(float,line.split()))
            temp_one = [dat_temp[3]]
            temp_one.extend(dat_temp[4:])
            data.append(temp_one)
    data = numpy.asarray(data)
#    print data

    ret_dat=numpy.zeros(data.shape)
    print(ret_dat)
    ret_dat[:,0]=data[:,0]
    for i in range(1,ret_dat.shape[1]):
        ret_dat[:,i] = data[:,0]*data[:,i]
    print(ret_dat)
    return ret_dat,ret_dat[1]


def _get_ph_band_data(oneCalc,ID,extension='phBAND.gp',postfix=''):

    data_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s%s.%s'%(ID,postfix,extension))


    #data =numpy.loadtxt(data_file_name,dtype=numpy.float64,)
    data = []

    with open(data_file_name,'r') as fo:
        fs=fo.read()
    fs=fs.split('\n')
    labels=fs[0]
    fs=fs[1:]
    for line in fs:
        if len(line.strip())!=0:
            data.append(list(map(float,line.split())))
    data = numpy.asarray(data)
    print(data)
    return data[:,1:],data[:,0]



    

def _get_gamma_velocity(freq,q):
    return (freq[0]-freq[10])/(q[0]-q[10])

def _get_debye_freq(oneCalc,ID):


    path_str = AFLOWpi.run._phonon_band_path(oneCalc,ID)

    path_pts_list = [int(i.split()[1]) for i in path_str.split('\n')[1:] if len(i.strip())!=0]
#    print path_pts_list
    freq,q_vals = AFLOWpi.retr._get_ph_dos_data(oneCalc,ID)
    TA       = freq[:,0]
    TA_prime = freq[:,1]
    LA       = freq[:,2]
#    q_vals   = freq[-1]

    total=0
#    for i in range(len(path_pts_list)):
    path_index= []
    for i in range(len(path_pts_list)):
        if path_pts_list[i]==0:
            continue
        if path_pts_list[i+1]==0:
            path_index.append([sum(path_pts_list[:i]),sum(path_pts_list[:i+1])+1])
        else:
            path_index.append([sum(path_pts_list[:i]),sum(path_pts_list[:i+1])])
    


    
#    index = for i in path_pts_list
#    print path_index
#    print TA
    average_freq_TA       = numpy.average([max(TA[i[0]:i[1]]) for i in path_index])
    average_freq_TA_prime = numpy.average([max(TA_prime[i[0]:i[1]]) for i in path_index])
    average_freq_LA       = numpy.average([max(LA[i[0]:i[1]]) for i in path_index])
    return average_freq_TA,average_freq_TA_prime,average_freq_LA

def _get_debye_temp(oneCalc,ID):
    #get the frequencies for TA, TA', and LA
    frequencies,q_vals = AFLOWpi.retr._get_ph_dos_data(oneCalc,ID)

    #q vals so we can find v_debye near gamma
#    q_vals = frequencies[-1]

    
#    frequencies[0] = [i*0.0299792458 for i in frequencies[0]]
#    frequencies[1] = [i*0.0299792458 for i in frequencies[1]]
#    frequencies[2] = [i*0.0299792458 for i in frequencies[2]]
    #get volume of original cell
    norm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='phonon')
    cell_vol = AFLOWpi.retr.getCellVolume(oneCalc,norm_ID,string=False,conventional=False)
    #convert to meters
    bohr2meter=5.29177e-11
    V=cell_vol*bohr2meter**3.0

    #get num atoms in cell
    N = float(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['nat'])

    #some constants
    h_bar=1.0545718*10**-34 
    k_b=1.38064852*10**-23

    #get v_debye for each branch
    v_i = AFLOWpi.retr._get_debye_freq(oneCalc,ID)
    v_s_TA       = v_i[0]
    v_s_TA_prime = v_i[1]
    v_s_LA       = v_i[2]


    #calculate debye temperature for each branch
    wo_v           = 2.0*numpy.pi*h_bar/k_b#*(6.0*numpy.pi**2.0*N/V)**(1.0/3.0)
    debye_TA       = wo_v*v_s_TA
    debye_TA_prime = wo_v*v_s_TA_prime
    debye_LA       = wo_v*v_s_LA


#    print debye_TA
#    print debye_TA_prime
#    print debye_LA

    return debye_TA,debye_TA_prime,debye_LA



def _therm_pp(oneCalc,ID):
#    grun_i  = AFLOWpi.retr._get_gruneisen(oneCalc,ID)
    grun_i=[0.0,0.0,0.0]
    AFLOWpi.retr._get_gruneisen(oneCalc,ID,band=False)

    theta_i = AFLOWpi.retr._get_debye_temp(oneCalc,ID)
#    print theta_i
    #get volume of original cell
    norm_ID  = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='phonon')
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

    #get the frequencies for TA, TA', and LA
    frequencies,q_vals = AFLOWpi.retr._get_ph_dos_data(oneCalc,ID)

    v_i=AFLOWpi.retr._get_debye_freq(oneCalc,ID)

    print(Vol)
    print(Mass)
    print(grun_i)
    print(v_i)
    print(theta_i)


    therm_cond_data_str="T       Total        TA           TA'          La"
    TEMP = numpy.linspace(0.0,2000.0,400.0)

    Vol   = 375.60264000000006 # volume of one cell in angstrom^3                                                      
    Vol  *= 10.0**(-30.0) # convert angstrom to m^3                                                                    
    Vol  /= 8.0 # do volume per atom                                                                                   

    Mass   = 63.5463*3.0  # 3 Cu                                                                                       
    Mass += 78.9718*4.0  # 4 Se                                                                                        
    Mass += 121.760*1.0  # 1 Antimony                                                                                  
    Mass *= 1.66054e-27 # convert amu to kg                                                                            
    Mass /= 8.0        # average mass per atom      

    v_i     = [1485.0, 1699.0, 3643.0] #TA, TA', LA                                                                         
    theta_i = [60.0,   65.0,   78.0  ] #TA, TA', LA                                                                      
    grun_i  = [1.27,   1.14,   1.26, ] #TA, TA', LA     


    for T in TEMP:
        total,TA_cont,TA_prime_cont,LA_cont = AFLOWpi.retr._do_therm(v_i,theta_i,grun_i,Mass,Vol,T)
        
        therm_cond_data_str+='\n%7.1f %12.3f %12.3f %12.3f %12.3f'%(T,total,TA_cont,TA_prime_cont,LA_cont)

    therm_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_thermal_cond.dat'%ID)
    with open(therm_file_name,'w') as tcfo:
        tcfo.write(therm_cond_data_str) 


def _test_therm_pp():

    therm_cond_data_str="T       Total        TA           TA'          La"
#   TEMP = numpy.linspace(1.0,2001.0,2000.0)
    TEMP = [80.0,
            87.3333333333333,
            94.6666666666667 ,
            102,
            109.333333333333,
            116.666666666667,
            124,
            131.333333333333,
            138.666666666667,
            146,
            153.333333333333,
            160.666666666667,
            168,
            175.333333333333,
            182.666666666667,
            190,
            197.333333333333,
            204.666666666667,
            212,
            219.333333333333,
            226.666666666667,
            234,
            241.333333333333,
            248.666666666667,
            256,
            263.333333333333,
            270.666666666667,
            278,
            285.333333333333,
            292.666666666667,
            300,
            ]

    Vol   = 375.60264000000006 # volume of one cell in angstrom^3                                                      
    Vol  *= 10.0**(-30.0) # convert angstrom to m^3                                                                    
    Vol  /= 8.0 # do volume per atom                                                                                   

    Mass   = 63.5463*3.0  # 3 Cu                                                                                       
    Mass += 78.9718*4.0  # 4 Se                                                                                        
    Mass += 121.760*1.0  # 1 Antimony                                                                                  
    Mass *= 1.66054e-27 # convert amu to kg                                                                            
    Mass /= 8.0        # average mass per atom      

    v_i     = [1485.0, 1699.0, 3643.0] #TA, TA', LA                                                                         
    theta_i = [60.0,   65.0,   78.0  ] #TA, TA', LA                                                                      
    grun_i  = [1.27,   1.14,   1.26, ] #TA, TA', LA     


    for T in TEMP:
        total,TA_cont,TA_prime_cont,LA_cont = AFLOWpi.retr._do_therm(v_i,theta_i,grun_i,Mass,Vol,T)
        
        therm_cond_data_str+='\n%7.1f %12.3f %12.3f %12.3f %12.3f'%(T,total,TA_cont,TA_prime_cont,LA_cont)

    therm_file_name = os.path.join('./','test_thermal_cond.dat')
    with open(therm_file_name,'w') as tcfo:
        tcfo.write(therm_cond_data_str) 


######################################################################################################################                                                                                               
######################################################################################################################                                                                                               

def _do_therm(v_i,theta_i,grun_i,Mass,Vol,T):
    ######################################################################################################################
    def calc_tau_N_TA(x,grun,velocity,vol,mass,T):
        #Relax Temp Scattering TA'                                                                                                                                          
        h_bar = numpy.float64(1.0545718*(10.0**(-34.0))) # hbar
        k_b   = numpy.float64(1.38064852*(10.0**(-23.0))) #boltzmann                                               

        one_over_tau_N_TA  = k_b**4.0
        one_over_tau_N_TA *= grun**2.0
        one_over_tau_N_TA *= vol
        one_over_tau_N_TA /= mass
        one_over_tau_N_TA /= h_bar**3.0
        one_over_tau_N_TA /= velocity**5.0

        one_over_tau_N_TA *= k_b/h_bar
        one_over_tau_N_TA *= x
        one_over_tau_N_TA *= T**5.0

        tau_N_TA  = 1.0/one_over_tau_N_TA

        return tau_N_TA
    ######################################################################################################################
    def calc_tau_N_LA(x,grun,velocity,vol,mass,T):
        #Relax Temp Scattering LA                                                                                                         
        h_bar = numpy.float64(1.0545718*(10.0**(-34.0))) # hbar
        k_b   = numpy.float64(1.38064852*(10.0**(-23.0))) #boltzmann                                                                                                                                                                                                      
        one_over_tau_N_LA  = k_b**3.0
        one_over_tau_N_LA *= grun**2.0
        one_over_tau_N_LA *= vol
        one_over_tau_N_LA /= mass
        one_over_tau_N_LA /= h_bar**2.0
        one_over_tau_N_LA /= velocity**5.0

        one_over_tau_N_LA *= (k_b/h_bar)**2.0
        one_over_tau_N_LA *= x**2.0
        one_over_tau_N_LA *= T**5.0

        tau_N_LA  = 1.0/one_over_tau_N_LA

        return tau_N_LA
    ######################################################################################################################                                                                                   
    def calc_tau_U(x,grun,velocity,debye_temp,mass,T):
        #Relax Temp Umklamp LA                                                                                     
        h_bar = numpy.float64(1.0545718*(10.0**(-34.0))) # hbar
        k_b   = numpy.float64(1.38064852*(10.0**(-23.0))) #boltzmann                                                                                            

        one_over_tau_U  = h_bar
        one_over_tau_U *= grun**2.0
        one_over_tau_U /= mass
        one_over_tau_U /= velocity**2.0
        one_over_tau_U /= debye_temp

        one_over_tau_U *= (k_b/h_bar)**2.0
        one_over_tau_U *= x**2.0
        one_over_tau_U *= T**3.0

        one_over_tau_U *= numpy.exp(-1.0*debye_temp/(3.0*T))

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
        sol *= numpy.exp(x)/(numpy.exp(x)-1.0)**2.0


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
        sol *= numpy.exp(x)/(numpy.exp(x)-1.0)**2.0


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
        sol *= numpy.exp(x)/(numpy.exp(x)-1.0)**2.0

        return sol



    ######################################################################################################################
    #define constants
    h_bar = numpy.float64(1.0545718*(10.0**(-34.0))) # hbar
    k_b   = numpy.float64(1.38064852*(10.0**(-23.0))) #boltzmann  

    #define a constant from other constants
    C_PHON_CONST  = 1.0/3.0
    C_PHON_CONST *= k_b**4.0
    C_PHON_CONST *= T**3.0
    C_PHON_CONST /= 2.0*numpy.pi**2.0
    C_PHON_CONST /= h_bar**3.0

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
    print(max_TA,max_TA1,max_LA)
    print() 
    print() 
    
    ##################################################################################
    #calculate the three integrals for TA phonon and find lattice k for TA
    TA_1 =  scipy.integrate.quad(first_int, 0.00,max_TA,(grun_TA,vel_TA,theta_TA,Vol,Mass,T,True),
                                 epsabs=0, epsrel=1.49e-6)[0]
    TA_2 =  scipy.integrate.quad(second_int,0.00,max_TA,(grun_TA,vel_TA,theta_TA,Vol,Mass,T,True),
                                 epsabs=0, epsrel=1.49e-6)[0]
    TA_3 =  scipy.integrate.quad(third_int, 0.00,max_TA,(grun_TA,vel_TA,theta_TA,Vol,Mass,T,True),
                                 epsabs=0, epsrel=1.49e-6)[0]

    klattice_TA  = C_TA * (TA_1 + TA_2**2.0/TA_3)           
    #calculate the three integrals for TA' phonon and find lattice k for TA'
    T1_1 = scipy.integrate.quad(first_int, 0.00,max_TA1,(grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True),
                                epsabs=0, epsrel=1.49e-6)[0]
    T1_2 = scipy.integrate.quad(second_int,0.00,max_TA1,(grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True),
                                epsabs=0, epsrel=1.49e-6)[0]
    T1_3 = scipy.integrate.quad(third_int, 0.00,max_TA1,(grun_TA1,vel_TA1,theta_TA1,Vol,Mass,T,True),
                                epsabs=0, epsrel=1.49e-6)[0]

    klattice_T1  = C_T1 * (T1_1 + T1_2**2.0/T1_3)           
    #calculate the three integrals for LA phonon and find lattice k for LA
    LA_1 =  scipy.integrate.quad(first_int, 0.00,max_LA,(grun_LA,vel_LA,theta_LA,Vol,Mass,T,False),
                                 epsabs=0, epsrel=1.49e-6)[0]
    LA_2 =  scipy.integrate.quad(second_int,0.00,max_LA,(grun_LA,vel_LA,theta_LA,Vol,Mass,T,False),
                                 epsabs=0, epsrel=1.49e-6)[0]
    LA_3 =  scipy.integrate.quad(third_int, 0.00,max_LA,(grun_LA,vel_LA,theta_LA,Vol,Mass,T,False),
                                 epsabs=0, epsrel=1.49e-6)[0]

#   LA_3 =  scipy.integrate.quad(third_int, 0.001,max_LA,(grun_LA,vel_LA,theta_LA,Vol,Mass,T,False),epsabs=1.49e-20, epsrel=1.49e-20)[
#   print LA_3.message
#   LA_3 = LA_3[0]


    klattice_LA  =  C_LA * (LA_1 + LA_2**2.0/LA_3)
    ##################################################################################                                           
    #sum over all contibutions for lattice k at Temp T 
    klattice_total  = 0.0
    klattice_total += klattice_TA
    klattice_total += klattice_T1
    klattice_total += klattice_LA

    return klattice_total,klattice_TA,klattice_T1,klattice_LA

#def  _get_long_phon(oneCalc,ID,optical=True):

#     band_file_name = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phBAND.gp'%ID)

#     with open(band_file_name,'r') as bdfo:
#         bdfs = bdfo.read()



#     gg=re.findall('^\s*0.00.*\n',bdfs,re.M)
#     split_at_gamma = gg[0].split()[1:]

#     long_phon_ind = [i+1 for i in range(len(split_at_gamma)) if abs(float(split_at_gamma[i]))==0.00]
# #    print long_phon_ind 
# #    print long_phon_ind 
# #    print long_phon_ind 

#     phon=[list(),list(),list()]
#     q_vals=[]
#     opt=[[]]
# #     for i in bdfs.split('\n'):
# # #        print i

# #         per_q = [float(j) for j in i.split() ]
# #         try:
# # #            phon[0].append(per_q[long_phon_ind[0]])
# # #            phon[1].append(per_q[long_phon_ind[1]])
# # #            phon[2].append(per_q[long_phon_ind[2]]) 
# #             q_vals.append(per_q[0])
# #             phon[0].append(per_q[1])
# #             phon[1].append(per_q[2])
# #             phon[2].append(per_q[3]) 


# #             if optical==True:
# #                 for k in range(len(per_q))[4:]:
# #                     try:
# # #                        print k-4
# #                         opt[k-4].append(per_q[k])
# #                     except Exception as e:
# #                         opt[k-4]=[]
# #                         opt[k-4]=[per_q[k]]

# #         except:
# #             pass


# #     v0 = _get_gamma_velocity(phon[0],q_vals)
# #     v1 = _get_gamma_velocity(phon[1],q_vals)
# #     v2 = _get_gamma_velocity(phon[2],q_vals)
# # #    grad_first_phon  = numpy.gradient(phon[0])
# # #    grad_second_phon  =numpy.gradient(phon[1])
# # #    grad_third_phon  = numpy.gradient(phon[2])


# #   #  sums= {0:sum(grad_first_phon[:30]),1:sum(grad_second_phon[:30]),2:sum(grad_third_phon[:30]),}
# #     sums = {0:v0,1:v1,2:v2}
# # #    print grad_first_phon[:50]
# # #    print grad_second_phon[:50]
# # #    print grad_third_phon[:50]


# #     sorted_ind = [i[0] for i in  sorted(sums.items())]

# #     #convert cm^-1 to rad/s
#      sf = 0.0299792458*10.0**12
# # #    print sf
# # #    print sf
# # #    print sf
# #     #scale frequencies to meters
#      phon[sorted_ind[0]]=[i*sf for i in phon[sorted_ind[0]]]
#      phon[sorted_ind[1]]=[i*sf for i in phon[sorted_ind[1]]]
#      phon[sorted_ind[2]]=[i*sf for i in phon[sorted_ind[2]]]
# # #    print q_vals
# # #    return [phon[0],phon[1],phon[2],q_vals]
# #     if optical==True:
# #         ret_list = [phon[sorted_ind[0]],phon[sorted_ind[1]],phon[sorted_ind[2]],]
# #         ret_list.extend(opt)
# #         ret_list.append(q_vals)
# #         return ret_list
# #     else:
# #         return [phon[sorted_ind[0]],phon[sorted_ind[1]],phon[sorted_ind[2]],q_vals]


import AFLOWpi
import os
import numpy

def gap_size(oneCalc,ID):
    try:
        bandgap_type=gap_type(oneCalc,ID)

        if bandgap_type not in ['p-type','n-type','insulator']:
            print('conductor..no gap')
            return 0.0
        elif bandgap_type=='p-type':
            en,dos=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=[0.1,40.0])
            start=0.1
            end=0.1
            for i in range(len(dos))[:-2]:
#            for i in range(len(dos)):
                if dos[i]>0.00001:
                    end=en[i]
                    break
            return end-start


        elif bandgap_type=='n-type':
            en,dos=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=[-40.0,-0.1])
            start=-0.1
            end=0.1
            for i in reversed(list(range(len(dos)))[:-2]):
#            for i in reversed(range(len(dos))):
                if dos[i]>0.00001:
                    end=en[i]
                    break
            return numpy.abs(end-start)


        elif bandgap_type=='insulator':
            en,dos_down=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=[-40.0,-0.1])
            en,dos_up=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=[0.1,40.0])
            start=-0.1
            end=0.1
            for i in reversed(list(range(len(dos_down)))[:-2]):
                if dos_down[i]>0.00001:
                    end=en[i]
                    break
            for i in range(len(dos_up)[-2]):
                if dos_down[i]>0.00001:
                    start=en[i]
                    break

            return numpy.abs(end-start)



                       
    except Exception as e:
        print(e)
        return 0.0

def gap_type(oneCalc,ID):
    try:
        dos_range = 0.1
        range_up=[0.05,dos_range]
        range_down=[-1.0*dos_range,-0.05]
        en_above,dos_above=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=range_up)
        en_below,dos_below=AFLOWpi.retr._get_dos(oneCalc,ID,dos_range=range_down)
        dos_above_found=False
        dos_below_found=False

        for i in reversed(list(range(len(dos_above)))[2:]):

            if dos_above[i]>0.0001:
                dos_above_found=True
        
        for i in reversed(list(range(len(dos_below)))[2:]):

            if dos_below[i]>0.0001:
                dos_below_found=True

        if dos_above_found==True and dos_below_found==True:
            return 'conductor'
        elif dos_above_found==False and dos_below_found==True:
            return 'p-type'
        elif dos_above_found==True and dos_below_found==False:
            return 'n-type'
        elif dos_above_found==False and dos_below_found==False:
            return 'insulator'
    except:
        return 'None'


def _get_dos(oneCalc,ID,LSDA=False,dos_range=[-0.1,0,1],normalize=True):
    try:
        dos_ID = AFLOWpi.prep._return_ID(oneCalc,ID,step_type='dos',last=True)

	'''extracts HOMO from nscf calculation output file as input to the plotting'''

	try:
		Efermi=AFLOWpi.retr._getEfermi(oneCalc,ID)
		if type(Efermi)!=type(0.5):
			LSDA=True
			
	except:
		Efermi=0.0
	subdir=oneCalc['_AFLOWPI_FOLDER_']

	"""get the path to the subdirectory of the calc that you are making plots for"""


	'''name of file of the DOS plots is dosBandPlot_<_AFLOWPI_PREFIX_>'''
	fileplot = os.path.join(subdir,'DOS_%s%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),oneCalc['_AFLOWPI_PREFIX_']))

	#to set figure size and default fonts


	"""get the path to the subdirectory of the calc that you are making plots for"""
        filedos = os.path.join(subdir,'%s_dos.dat'%dos_ID)
        
	try:
		data = open(filedos,'r').readlines()
	except Exception:
            pass

	en = []
	enup = []
	endown = []
	dos = []
	dosdw = []
        scaling_dosup=[]
        scaling_downdw=[]
        scaling_dos=[]
	for i in range(1, len(data)):      #append DOS data to en,dos,dosw lists
		try:
			if LSDA==True:
                            val_up   = float(data[i].split()[0])-Efermi[0]
                            val_down = float(data[i].split()[0])-Efermi[1]
 #                           if val_up < dos_range[1] and valZ_up > dos_range[0]:
                            enup.append(val_up)
#                            if val_down <dos_range[1] and val_down > dos_range[0]:
                            endown.append(val_down)
				
			else:
                            val=float(data[i].split()[0])-Efermi

                        try:
                            scaling_dosdw.append(-1*float(data[i].split()[2]))
                            scaling_dosup.append(float(data[i].split()[1]))

                            if val_down <dos_range[1] and val_down > dos_range[0]:
                                dosdw.append(-1*float(data[i].split()[2]))
                                en_down.append(val_down) #to shift all the y values with respect to the Fermi level
                            if val_up <dos_range[1] and val_up > dos_range[0]:
                                dos.append(float(data[i].split()[1]))
                                en_up.append(val_up) #to shift all the y values with respect to the Fermi level
                        except:
                            scaling_dos.append(float(data[i].split()[1]))
                            if val > dos_range[0] and val <dos_range[1]:
                                dos.append(float(data[i].split()[1]))
                                en.append(val) #to shift all the y values with respect to the Fermi level
		
		except Exception as e:
			pass
	if LSDA==True:
		enup  = list(map(float,enup))
		endown= list(map(float,endown))

                floatdosDOWN=list(map(float,dosdw))
                floatdos=list(map(float,dos))

#                renormalize_dw=1.0/sum(scaling_dosdw)
#                renormalize_up=1.0/sum(scaling_dosup)
                renormalize_dw=1.0/sum(dosdw)
                renormalize_up=1.0/sum(dosup)

                array_dosup = numpy.asarray(floatdos)*renormalize_up
                floatdos=array_dos.tolist()

                array_dosdw = numpy.asarray(floatdosDOWN)*renormalize_dw
                floatdosDOWN=array_dosdw.tolist()

	else:
		endos=list(map(float,en))  #to convert the list x from float to numbers
#                renormalize=1.0/sum(scaling_dos)
                renormalize=1.0/sum(dos)                
                floatdos=list(map(float,dos))
                array_dos = numpy.asarray(floatdos)*renormalize
                floatdos=array_dos.tolist()

	enshift = numpy.array(endos) #to treat the list b as an array?

        
        if LSDA==True:
            return enshift,floatdos,floatdosDOWN
        else:
            return enshift,floatdos

    except:
        return None

import numpy
import AFLOWpi
import copy
import decimal


def supercell(inputString,numX=1,numY=1,numZ=1):
    inputString = AFLOWpi.prep._transformInput(inputString)    
    return AFLOWpi.retr._constructSupercell(inputString,numX=numX,numY=numY,numZ=numZ)




def _expandBoundaries(labels,symMatrix,numX,numY,numZ,beginX=0,beginY=0,beginZ=0):
    superList=[]
    try:
        symMatrix=symMatrix.getA()
    except:
        pass

    for entry in range(len(symMatrix)):
        for x in range(0,numX):
            first  = (symMatrix[entry][0]+x)/float(abs(numX))
            second = symMatrix[entry][1]
            third  = symMatrix[entry][2]
            superList.append([labels[entry],first,second,third])

    newPos=[x[1:] for x in superList]
    labels=[x[0] for x in superList]
    superList=[]
    #################################################################
    #################################################################
    for entry in range(len(newPos)):
        for y in range(0,numY):
            first  = newPos[entry][0]
            second = (newPos[entry][1]+y)/float(abs(numY))
            third  = newPos[entry][2]
            superList.append([labels[entry],first,second,third])

    labels=[x[0] for x in superList]
    newPos=[x[1:] for x in superList]
    superList=[]
    #################################################################
    #################################################################
    for entry in range(len(newPos)):
        for z in range(0,numZ):
            first  = newPos[entry][0]
            second = newPos[entry][1]
            third  = (newPos[entry][2]+z)/float(abs(numZ))
            superList.append([labels[entry],first,second,third])

    labels=[x[0] for x in superList]
    newPos=[x[1:] for x in superList]

    orig_list=[]
    orig_atom_ss_index = [(x-1)*numX*numY*numZ for x in range(1,len(symMatrix)+1)]


    for i in range(len(orig_atom_ss_index)):
        
        popped=superList.pop(orig_atom_ss_index[i])
        superList.insert(i,popped)

    symMatrix= numpy.array([x[1:] for x in superList])
    labels = numpy.array([x[0] for x in superList])

    return labels,symMatrix

def _constructSupercell(inputString,numX=1,numY=1,numZ=1,stringOrMatrix='String',newVectors=True):
    splitInput =  AFLOWpi.retr._splitInput(inputString)
    
    if '{crystal}' != splitInput['ATOMIC_POSITIONS']['__modifier__']:
        logging.error('unit in AFLOWpi.retr._constructSupercell not for ATOMIC_POSITIONS MUST BE {crystal}')
        return inputString
    cellParamMatrix = AFLOWpi.retr.getCellMatrixFromInput(inputString)

    labels =  AFLOWpi.retr._getPosLabels(inputString)
    symMatrix = AFLOWpi.retr._getPositions(inputString)

    symMatrix_orig=copy.deepcopy(symMatrix)
    labels_orig=copy.deepcopy(labels)

    coordold,flags = AFLOWpi.retr.detachPosFlags(AFLOWpi.qe.regex.atomic_positions(inputString))

    outputString=''
    superList=[]

    if newVectors==True:
        labels,symMatrix=AFLOWpi.retr._expandBoundaries(labels,symMatrix,numX,numY,numZ)
    else:
        labels,symMatrix=AFLOWpi.retr._expandBoundariesNoScale(labels,symMatrix,numX,numY,numZ)

    for entry in range(len(symMatrix)):
        posLineStr = ' '.join(['%20.14f' % (decimal.Decimal(str(numpy.around(i,9)))) for i in symMatrix[entry]])+'\n'
        outputString+='%4s %8s' % (labels[entry],posLineStr)

    
    splitInput['&system']['nat']=str(len(labels))

    splitInput['&system']['celldm(1)']=str(float(splitInput['&system']['celldm(1)'])*numX)
    if 'celldm(2)' in list(splitInput['&system'].keys()):
        scaleY = float(numY)/float(numX)
        splitInput['&system']['celldm(2)']=str(float(splitInput['&system']['celldm(2)'])*scaleY)
    if 'celldm(3)' in list(splitInput['&system'].keys()):
        scaleZ = float(numZ)/float(numX)
        splitInput['&system']['celldm(3)']=str(float(splitInput['&system']['celldm(3)'])*scaleZ)


    splitInput['ATOMIC_POSITIONS']['__content__']=outputString

    returnString=AFLOWpi.retr._joinInput(splitInput)

    return returnString



import AFLOWpi
import os
import shutil
import subprocess

def atomicDistances(calcs,runlocal=False,inpt=False,outp=True):
    engineDir  = AFLOWpi.prep._ConfigSectionMap("prep",'engine_dir')	
    if os.path.isabs(engineDir) == False:
        configFileLocation = AFLOWpi.prep._getConfigFile()
        engineDir =  os.path.join(configFileLocation, enginedir)
    distXPath = os.path.join(engineDir,'dist.x')
    try:
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_exec').lower()!='false':
            AFLOWpi.prep.totree(distXPath,calcs)
        for ID,oneCalc in calcs.items():
            if runlocal:
                if outp==True:
                    AFLOWpi.retr._getDist(oneCalc,ID,outp=True)
                if inpt==True:
                    AFLOWpi.retr._getDist(oneCalc,ID,outp=False)
            else:
                if outp==True:
                    AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','AFLOWpi.retr._getDist(oneCalc,ID,outp=True)\n')
                if inpt==True:
                    AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','AFLOWpi.retr._getDist(oneCalc,ID,outp=False)\n')

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)


def _getDist(oneCalc,ID,outp=True):
    try:
        if outp==True:
            AFLOWpi.retr._writeInputFromOutput(oneCalc,ID)
            os.rename(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_new.in'),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_dist.in')) 
        else:
            shutil.copyfile('%s.in'%ID,'%s_dist.in'%ID)

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    #fix to clear restart_mode in input file to avoid error with dist.x
    dist_ID='%s_dist'%ID        
    infil = AFLOWpi.retr._getInputFileString(oneCalc,dist_ID)
    inDict=AFLOWpi.retr._splitInput(infil)
    inDict['&control']['restart_mode']='"from_scratch"'
    infil = AFLOWpi.retr._joinInput(inDict)
    dest_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_dist.in')

    with open(dest_file,'w') as fo:
        fo.write(infil)

    if AFLOWpi.prep._ConfigSectionMap('prep','copy_exec').lower()=='false':
        engineDir = AFLOWpi.prep._ConfigSectionMap('prep','engine_dir')
        if os.path.isabs(enginedir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            enginedir =  os.path.join(configFileLocation, enginedir)

        execPath=os.path.join(engineDir,'dist.x')
    else:
        execPath='./dist.x'
    try:
        subprocess.Popen('%s < %s_dist.in > /dev/null' % (execPath,ID),cwd=oneCalc['_AFLOWPI_FOLDER_'],shell=True).wait()
        if outp==True:
            os.rename(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'dist.out'),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_output_dist.out'))
        else:
            os.rename(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'dist.out'),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_input_dist.out'))

    except Exception as e:
        logging.error('dist.x did not run properly')
        print('ERROR: dist.x did not run properly')


import AFLOWpi.run
import AFLOWpi.prep
import os
import datetime
import pickle
import logging 
import re
import numpy
import copy
import AFLOWpi.retr
import traceback
import sys 
import AFLOWpi.qe
import decimal
import contextlib
import io

@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = io.StringIO()
    yield
    sys.stdout = save_stdout


##############################################################################################################################
##############################################################################################################################
## auto k point
##############################################################################################################################
##############################################################################################################################
def _find_numkpoints(outputFile):
    '''
    DEFUNCT <TAGGED FOR REMOVAL>

    Arguments:


    Keyword Arguments:


    Returns:


    '''

    kpointNumRegex = re.compile(r'\s*number\sof\sk\spoints\s*=\s*(\d*).*\n')
    try:
        return int(kpointNumRegex.findall(outputFile)[-1])
    except Exception as e:
        return 1

def _get_pool_num(oneCalc,ID):
    '''
    Gets number of pools requested for this particular execution for engine.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
           None

    Returns:
          npool (int): number of pools requested for this particular execution for engine

    '''

    hard_limit=4
    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out'%ID),'r') as outputFile:
            outputFileString=outputFile.read()
        n_reduced_k = AFLOWpi.retr._find_numkpoints(outputFileString)
        execPrefix=AFLOWpi.prep._ConfigSectionMap('run','exec_prefix')
        if n_reduced_k!=1:
            num_mpi_procs = int(re.findall(r'np\s*(\d*)',execPrefix)[-1])
            npool=1
            for i in reversed(list(range(1,n_reduced_k+1))):
                if num_mpi_procs%i == 0 and i<hard_limit+1:
                    npool=i
                    break
        return int(npool)

    except Exception as e:
        return 1




def detachPosFlags(positionString):
    '''
    Detach the position flags from the string of the atomic positions

    Arguments:
          positionString (str): string of the atomic positions

    Keyword Arguments:
          None

    Returns:
          positions (list): atomic positions as a list
          flags (list): a list of the flags for by each position

    '''

    posArray = positionString.split('\n')
    posArrayStripped=[]
    flagArray=[]
    
    for pos in range(len(posArray)):
        if len(posArray[pos])!=0:
            if len(posArray[pos].split())<5:
                posArray[pos]+=' 1 1 1'
            try:

                posArrayStripped.append(' '.join(posArray[pos].split()[:-3]))
                flagArray.append(' '.join(posArray[pos].split()[-3:]))
            except:
                pass
    flags='\n'.join(flagArray)
    positions='\n'.join(posArrayStripped)

    return positions,flags




def attachPosFlags(positionString,flagString):
    '''
    Reattaches the flags to the end of the atomic position

    Arguments:
          positionString (str): string of the atomic positions
          flagString (str): string of the flags

    Keyword Arguments:
          None

    Returns:
          positionString (str): Positions with flags attached

    '''

    positionStringSplit = positionString.split('\n')
    flagStringSplit = flagString.split('\n')
    combined = ['%s %s' % (positionStringSplit[x],flagStringSplit[x]) for x in range(len(flagStringSplit))]
    positionString = '\n'.join(combined)

    return positionString

def checkStatus(PROJECT,SET='',config='',step=0,status={},negate_status=False):
    '''
    function that loads the calclogs of each of the calculations in your run and 
    displays status information about them.

    Arguments:
          PROJECT        (str): Project name
    Keyword Arguments:
          SET            (str): Set name
          config         (str): Config file used
          step           (int): Which number step to seee (default all of them)
          status        (dict): dictionary containing the status and the value 
                                you want to see (ex. {'Error':True})
          negate_status (bool): Whether to take the calculations with that status or without it

    Returns:
          calcsList (list): list of each step of the calculations you
                           selected that satisfay the status criteria

    '''

    calcsList=[]    
    if step==0:    
        while True:
            step+=1
            try:
                one_step = AFLOWpi.prep.loadlogs(PROJECT, SET,'step_%02d' % step,config=config,suppress_warning=True)
                calcsList.append(one_step)

            except:
                break
    else:
        calcsList=[AFLOWpi.prep.loadlogs(PROJECT, SET,'step_%02d' % step,config=config)]

    outString=''

    

    for step in range(len(calcsList)):

        origLength=len(calcsList[step])
        string_prev=''
        header = ['Folder'.ljust(8),'ID'.ljust(25)]
        for ID,oneCalc in calcsList[step].items():
            string = ['%-8s' % x for x in list(calcsList[step][ID]['__status__'].keys())]
            if len(string_prev)>len(string):
                string=string_prev
                string_prev=string
                header.extend(string)
            
        calcCopy=copy.deepcopy(calcsList[step])
        if len(calcsList[step])!=0:
            for ID,oneCalc in calcsList[step].items():
                try:
                    for k,v in status.items():
                        if negate_status:
                            if oneCalc['__status__'][k]==v:
                                del calcCopy[ID]                                
                        else:
                            if oneCalc['__status__'][k]!=v:
                                del calcCopy[ID]                                
                except:
                    pass

            headerString=' | '.join(header)
            outString+=  '-'*(len(headerString))+'\n'
            outString+=  'STEP: %s' % (step+1)+'  | %s/%s'%(len(calcCopy),origLength)+'\n'
            outString+=  '-'*(len(headerString))+'\n'
            outString+=  headerString+'\n'
            outString+=  '-'*(len(headerString))+'\n'
            for ID,oneCalc in calcCopy.items():
                stringStatusList = ['%-8s' % x for x in list(oneCalc['__status__'].values())]
                string=[os.path.basename(oneCalc['_AFLOWPI_FOLDER_'].split('_')[-1]).ljust(8),ID.ljust(25)]
                string.extend(stringStatusList)
                outString+= ' | '.join(string)+'\n'
            #copy back for return of subset
            calcsList[step]=copy.deepcopy(calcCopy)


    return calcsList
            
def _getOutputString(oneCalc,ID):
    '''
    Gets string of output for that particular step in the workflow for a single calcualtion

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
           None

    Returns:
          outFileString (str): output of that particular step in the workflow

    '''

    outFilePath =  os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out'  % ID)
    try:
        with open(outFilePath,'r') as outFile:
            outFileString = outFile.read()
        return outFileString
    except:
        logging.warning('could not get output file: %s' % outFilePath)
        print('could not get output file: %s' % outFilePath)

def getCellVolume(oneCalc,ID,conventional=True,string=True):
    '''
    Gets the cell volume from output if avaiable and if not gets it from the input

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step


    Keyword Arguments:
          conventional (bool): return the volume of the primitive or conventional lattice 
          string (bool): to return as a string or as a float

    Returns:
          vol (str): Volume of the cell 

    '''


    try:
        print(ID)
        outFileString = AFLOWpi.retr._getOutputString(oneCalc,ID)
        vol = float(re.findall(r'unit-cell volume\s*=\s*([0-9.-]*)',outFileString)[-1])




    except:
        try:
            input_str= AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],ID)['_AFLOWPI_INPUT_']
            cell_vec = AFLOWpi.retr.getCellMatrixFromInput(input_str)
            vol =  AFLOWpi.retr.getCellVolumeFromVectors(cell_vec)

        except Exception as e:
            print(e)
            raise SystemExit
            logging.warning('could not get volume from output')
            print('could not get volume from output')

    if conventional==True:
        ibrav=int(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['ibrav'])
        if ibrav in [3,7,9,11,13]:
            vol*=2.0
        if ibrav in [2,10]:
            vol*=4.0

    if string==True:
        return str(vol)
    else:
        return vol


def getCellVolumeFromVectors(cellInput):
    '''
    Calculates cell volume from basis vectors

    Arguments:
          cellInput (numpy.matrix): basis set defining the cell

    Keyword Arguments:
          None

    Returns:
          vol (float): volume of the cell (may be not scaled. it depends on your input matrix)

    '''
    vol = numpy.cross(cellInput[0],cellInput[1]).dot(cellInput[2].T).getA()[0][0]
    return vol

def _getInputFileString(oneCalc,ID):
    '''
    Gets the string of the input that step of the calculation

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          None

    Returns:
          inFileString (str): input string to that calculation step

    '''

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'r') as inFile:
        inFileString = inFile.read()
    
    return inFileString 


def _getInitialInputString(oneCalc):
    '''
    Gets initial input to the workflow 

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation

    Keyword Arguments:
          None

    Returns:
          inFileString (str): initial input to the workflow 

    '''

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],oneCalc['_AFLOWPI_PREFIX_'][1:]+'.in'),'r') as inFile:
        inFileString = inFile.read()
        
    return inFileString 


def _getInitialOutputString(oneCalc):
    '''
    Gets initial output to the workflow 

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation

    Keyword Arguments:
          None

    Returns:
          inFileString (str): initial output to the workflow 

    '''


    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],oneCalc['_AFLOWPI_PREFIX_'][1:]+'.out'),'r') as outFile:
        outFileString = outFile.read()
        
    return outFileString 



def _getSymList(ID,oneCalc):
    '''
    Gets symmetry operations from output if available

    Arguments:
          ID (str): ID string for the particular calculation and step
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
       
    Keyword Arguments:
          None

    Returns:
          joinedList (str): a string of the sym ops pulled from the engine output

    '''

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],oneCalc['_AFLOWPI_PREFIX_'][1:]+'.out'),'r') as outFile:
        outFileString = outFile.read()
    joinedList = ''.join(set(re.findall(r'isym.*\n',outFileString)))
    joinedList = "Symmetry Operations:\n"+joinedList
    return joinedList

def _moveToSavedir(filePath):
    '''
    Move a file to the savedir specified in the config file

    Arguments:
          filePath (str): path of the file to be copied

    Keyword Arguments:
          None

    Returns:
          None

    '''

    try:
        if AFLOWpi.prep._ConfigSectionMap('prep','save_dir') != '':
            try:
                savedir = AFLOWpi.prep._ConfigSectionMap('prep','save_dir')
                workdir = AFLOWpi.prep._ConfigSectionMap('prep','work_dir')

                try:
                    workdir=os.path.realpath(workdir)
                except:
                    pass

                zeroBaseDir = os.path.dirname(os.path.dirname(filePath))
                zeroBaseDirName = os.path.basename(os.path.dirname(filePath))
                firstBaseDir = os.path.dirname(zeroBaseDir)
                try:
                    secondBaseDir = os.path.dirname(firstBaseDir)
                    secondBaseDir = os.path.realpath(secondBaseDir)
                    firstBaseDir =  os.path.realpath(firstBaseDir)

                except:
                    secondBaseDir = ''


                if os.path.isabs(savedir) == False:
                    configFileLocation = AFLOWpi.prep._getConfigFile()
                    savedir =  os.path.join(configFileLocation, savedir)


                if not os.path.exists(savedir):
                    try:
                        os.mkdir(savedir)
                    except:
                        logging.warning('Could not transfer %s to %s. %s does not exist and could not be created. ' % (filePath,savedir,savedir))

                if os.path.realpath(firstBaseDir)==os.path.realpath(workdir):
                    firstLevel = os.path.basename(os.path.dirname(zeroBaseDir))
                    firstLevelNew = os.path.join(savedir,firstLevel)
                    secondLevelNew=''
                    savedir = os.path.join(firstLevelNew,zeroBaseDirName)
                elif os.path.realpath(secondBaseDir)==os.path.realpath(workdir):

                    firstLevel = os.path.basename(os.path.dirname(zeroBaseDir))
                    secondLevel = os.path.basename(os.path.dirname(firstBaseDir))
                    firstLevelNew = os.path.join(savedir,secondLevel)
                    secondLevelNew = os.path.join(firstLevelNew,firstLevel)

                    savedir = os.path.join(secondLevelNew,zeroBaseDirName,)
                else:
                    firstLevelNew=savedir
                    secondLevelNew=''



                if not os.path.exists(firstLevelNew):
                    try:
                        os.mkdir(firstLevelNew)
                    except:
                        logging.warning('Could not transfer %s to %s. %s does not exist and could not be created. ' % (filePath,firstLevelNew,firstLevelNew))


                if secondLevelNew!='':
                    if not os.path.exists(secondLevelNew):
                        try:
                            os.mkdir(secondLevelNew)
                        except:
                            logging.warning('Could not transfer %s to %s. %s does not exist and could not be created. ' % (filePath,secondLevelNew,secondLevelNew))

                if os.path.isfile(savedir):
                    savedir+='_SAVEDIR'
                if not os.path.exists(savedir):
                    try:
                        os.mkdir(savedir)
                    except:
                        logging.warning('Could not transfer %s to %s. %s does not exist and could not be created. ' % (filePath,savedir,savedir))

                os.system('cp  %s %s/' % (filePath,os.path.abspath(savedir)))
            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

def grabEnergy(oneCalc,ID):
    '''
    Grabs energy for oneCalc and adds the keyword 'Energy' with the value of the energy grabbed from the output

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          None

    Returns:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
                          with the 'Energy' keyword added
     
    '''

    newDict=grabEnergyOut({ID:oneCalc})
    return newDict[ID]['Energy']

def grabEnergyOut(calcs):
    '''
    Goes in every subdirectory of the calculation and searches
    for the final energy of the calculation and returns a new 
    copy of the input dictionary that includes the final energy.

    Arguments:
	  calcs (dict): dictionary of dictionaries of calculations

    Keyword Arguments:
          None

    Returns:
	  calcs (dict): dictionary of dictionaries of calculations with energy added

    '''

    calcs1 = copy.deepcopy(calcs)
    energyRegex = re.compile(r'(?:(?:(?:(?:\!\s+)total)|(?:Final)) en\w+\s*=\s+(.+?)Ry)',re.MULTILINE)
    
    for ID,oneCalc in calcs1.items():
        try:
            if os.path.exists(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID)):
                with file(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID),'r') as outFile:
                    outFileString=outFile.read()
                    #searches for either total energy or final energy in the scf output and 
                    #adds it to the output dictionary
                    finalEnergy=energyRegex.findall(outFileString)

                    if len(finalEnergy):
                        energyFloat  = float(finalEnergy[-1])

                        calcs[ID]['Energy']= energyFloat
                    else: #if the energy can not be found the test entry is deleted from the output dictionary
                        outCalcPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % oneCalc['_AFLOWPI_PREFIX_'][1:])
                        logging.warning('could not get energy. check output file: %s' % outCalcPath)
                        print('could not get energy. check output file: %s' % outCalcPath)
                        calcs[ID]['Energy']=0.0
            else:
                
                outCalcPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % oneCalc['_AFLOWPI_PREFIX_'][1:])
                logging.warning('could not get energy. check output file: %s' % outCalcPath)
                print('could not get energy. check output file: %s' % outCalcPath)
                calcs[ID]['Energy']=0.0
        except:
            outCalcPath = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % oneCalc['_AFLOWPI_PREFIX_'][1:])
            logging.warning('could not get energy. check output file: %s' % outCalcPath)
            print('could not get energy. check output file: %s' % outCalcPath)
    return calcs



def getForce(oneCalc,ID,string=True):
    '''
    Gets last entry of the total force in the calculation from the engine output.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
           None

    Returns:
            force_string (str): Total force from engine output

    '''

    output = AFLOWpi.retr._getOutputString(oneCalc,ID)
    forceRegex = re.compile(r'(Total force =\s+[0-9.]+)')
    try:
        force_string = forceRegex.findall(output)[-1]

        if string==False:
            return float(force_string.strip().split()[-1])
        else:
            return force_string
    except:
        return ''


def getStress(oneCalc,ID):
    stressRegex = re.compile(r'(\s+total\s+stress.+(?:\(kbar\).+\n)(?:.+\n)*)')
    output = AFLOWpi.retr._getOutputString(oneCalc,ID)
    try:
        stress_string = stressRegex.findall(output)[-1]
        return stress_string
    except:
        return ''


def getCellOutput(oneCalc,ID):
    '''
    retreives information about the structure and its chemistry and prints
    it into a file in the AFLOWpi folder inside the project/set directories

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          None

    Returns:
          outputStr (str): Output of the report for the calculation

    '''

    outputStr='********************************************************\n'
    outputStr+=ID+'\n'
    folder=oneCalc['_AFLOWPI_FOLDER_']

    try:
        scfOutput = '%s.out' % ID
        with open(os.path.join(folder,scfOutput),'r') as outFile:

            lines = outFile.read()
    except Exception as e:
        print("No caclulation output available for %s. Are you sure the test ran properly?" % scfOutput)

        return

    try:
        scfInput = '%s.in' % ID
        with open(os.path.join(folder,scfInput),'r') as inFile:

            inLines = inFile.read()
    except Exception as e:


        print("No caclulation output available for %s. Are you sure the test ran properly?" % scfOutput)
        return
    retrDict = {}
    try:
        outputStr+=os.path.abspath(os.path.join(folder,scfInput))+'\n'
        outputStr+=__getStoicName(oneCalc)+'\n'
        outputStr+='********************************************************\n\n'        
        coordRegex = re.compile(r'(Begin final coordinates(?:.*\n)+End final coordinates)',re.MULTILINE)
        energyRegex = re.compile(r'(?:(?:(?:(?:\!\s+)total)|(?:Final)) en\w+\s*=\s+(.+?)Ry)',re.MULTILINE)

        stressRegex = re.compile(r'(\s+total\s+stress.+(?:\(kbar\).+\n)(?:.+\n)*)')
        forceRegex = re.compile(r'(Total force =\s+[0-9.]+)')



        coord   = coordRegex.findall(lines)
        energyArr = energyRegex.findall(lines)
        stressArr = stressRegex.findall(lines)
        forceArr  = forceRegex.findall(lines)
    except Exception as e:
        pass


    try:
        symList = AFLOWpi.retr._getSymList(ID,oneCalc)
        outputStr+=symList+'\n\n'
    except:
        pass
    try:
        outputStr+='Total Energy '+energyArr[-1]+' Ry\n\n'
        retrDict['energy']= energyArr[-1]
    except Exception as e:
        retrDict['energy']= ''

        pass
    try:
        outputStr+=forceArr[-1]+'\n\n'
        retrDict['force'] = forceArr[-1]
    except Exception as e:
        retrDict['force']= ''
        pass
    try:
        outputStr+=stressArr[-1]+'\n\n'
        retrDict['stress']= stressArr[-1]
    except Exception as e:
        retrDict['stress']= ''
        pass

    ###################
    ####CELL PARAMS####
    ###################
    atomicCELLPARAMRegex = re.compile(r"Begin final coordinates\n.+\n+.+\n(.+)\n(.+)\n(.+)\n",re.MULTILINE)
    initialParamsRegex = re.compile(r"crystal axes: \(cart. coord. in units of alat\)\n(.+)\n(.+)\n(.+)\n",re.MULTILINE)

    inputStr=''
    atomCELLPARAMSortList = []
    try:
        alatStr =  re.findall(r'alat\s*=\s*[0-9.-]*',lines)[-1]
    except:
        alatStr=''
    try:
        try:
            atomicCELLPARAMArr = atomicCELLPARAMRegex.findall(lines)
            atomicCELLSTRING = '\n'.join(atomicCELLPARAMArr[0])
        except:
            atomicCELLSTRING='Failed to get final coordinated. May have had convergence issue'
        try:
            atomicCELLPARAMArr = atomicCELLPARAMArr[-1]
        except:
            atomicCELLPARAMArr = initialParamsRegex.findall(lines)
            try:
                atomicCELLPARAMArr = atomicCELLPARAMArr[-1]
            except:
                pass
        if not len(atomicCELLPARAMArr) > 1:
            try: 
                splitAtomicCELLPARAM = atomicCELLPARAMArr.split('\n')[1:]
            except:
                pass
        else:
            splitAtomicCELLPARAM = atomicCELLPARAMArr
            testDictString=''
        try:
            outputStr+= alatStr+'\nCELL_PARAMETERS\n'+atomicCELLSTRING+'\n'
            inputStr+=alatStr+'\nCELL_PARAMETERS\n'+atomicCELLSTRING
            testDictString+=atomicCELLSTRING
        except:
            pass
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        testDictString = ''
    try:
        retrDict['CELL_PARAMETERS']=testDictString                
    except:
        retrDict['CELL_PARAMETERS'] = ''

    ######################
    ###ATOMIC POSITIONS###
    ######################
    atomicNLRegex      = re.compile(r'(ATOMIC_POSITIONS\s*.+\n)')
    atomicPosRegex     = AFLOWpi.qe.regex.atomic_positions(lines,'content','regex')
    try:
        try:
                inputStr+=atomicNLRegex.findall(lines)[-1].strip()+'\n'
        except ValueError:
                pass

        atomPosSortList = []
        try:
                atomicPosArr = atomicPosRegex.findall(lines)[-1]
                outputStr+="ATOMIC_POSITIONS\n"+atomicPosArr+'\n\n'
                splitAtomicPos = atomicPosArr.split('\n')
        except ValueError:
                pass

        for item in splitAtomicPos:
                if item!='':
                        atomSplit = item.split()
                        formattedAtomPos = atomSplit[0]+' '+' '.join(['%g' % (decimal.Decimal(str(atomSplit[i]))) for i in range(1,len(atomSplit))])+'\n'
                        atomPosSortList.append(formattedAtomPos)
        atomPosSortList = sorted(atomPosSortList)
        testDictString=''
        for formattedAtomPos in atomPosSortList:
                testDictString+=formattedAtomPos
                inputStr+=formattedAtomPos

        retrDict['ATOMIC_POSITIONS']=testDictString

    except:
        retrDict['ATOMIC_POSITIONS']=''
    try:
        baseDir='/'.join(folder.split('/')[:-1])+'/AFLOWpi/summary.log'
    except:
        baseDir='./'

    with open(baseDir,'a+') as outFile:
        outFile.write(outputStr)

    return outputStr

def _getCellParams(oneCalc,ID):
    '''
    Reads the output from the SCF or relax calculation to get the primitive cell parameters
    produced by the calculation. If it fails it defaults to the input lattice vectors

    Arguments:
          oneCalc (dict): a dictionary of a single calculation
          ID (str): ID string for the particular calculation and step          

    Keyword Arguments:
          None

    Returns:
          alat (float): the scale for the lattice vectors
          cell_matrix (numpy.array): the unscaled primitive lattice vectors

    '''

    try:
        scfOutput = '%s.out' % oneCalc['_AFLOWPI_PREFIX_'][1:]
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],scfOutput),'r') as outFile:
                lines = outFile.read()




        cellParams = AFLOWpi.qe.regex.cell_parameters(lines,'content')
        splitInput = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
        alat=float(splitInput['&system']['celldm(1)'])
        paramMatrix=_cellStringToMatrix(cellParams)

        if len(cellParams):

                try:
                    return float(alat[0]),paramMatrix
                except:
                    return float(alat),paramMatrix
        elif len(initialParams):
                paramMatrixList = []
                for params in initialParams[0]:
                        paramArray = params.split()

                        paramMatrixList.append([paramArray[3],paramArray[4],paramArray[5]])
                paramMatrix =  numpy.array(paramMatrixList,dtype='|S10')
                paramMatrix = paramMatrix.astype(numpy.float)

                return float(alat[0]),paramMatrix

        else:
                print('No card!')
                return AFLOWpi.retr.getCellMatrixFromInput(oneCalc['_AFLOWPI_INPUT_'])

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
     
        paramMatrix = AFLOWpi.retr.getCellMatrixFromInput(oneCalc['_AFLOWPI_INPUT_'])

        alat = AFLOWpi.retr._getAlatFromInput(oneCalc['_AFLOWPI_INPUT_']) 
        paramMatrix/=alat

        return alat,paramMatrix 

def get_parameters(oneCalc,ID,conventional=True):
    if conventional:
        cell = AFLOWpi.retr._getConventionalCellFromInput(oneCalc,ID)[2]

    else:
        cell = AFLOWpi.retr.getCellMatrixFromInput(oneCalc['_AFLOWPI_INPUT_'])

    return AFLOWpi.retr.free2abc(cell,cosine=False,bohr=False,string=False)
        
def _getAlatFromInput(inputString):
    '''
    Grabs the value of alat from the QE pwscf input

    Arguments:
          inputString

    Keyword Arguments:
          None

    Returns:
          alat (float): the alat scale for the primitive lattice vectors

    '''

    splitInput = AFLOWpi.retr._splitInput(inputString)
    try:
        return float(splitInput['&system']['celldm(1)'])
    except:
        pass
    try:
        return float(splitInput['&system']['A'])
    except:
        pass
    try:
        return float(splitInput['&system']['a'])
    except:
        return 0

def getRecipParams(oneCalc):
    '''
    reads the output from the SCF or relax calculation to get the 
    reciprocal cell parameters produced by the calculation.

    Arguments:
          oneCalc (dict): a dictionary of a single calculation

    Keyword Arguments:
          None

    Returns:
          alat (float): alat multiplier
          paramMatrix (numpy.matrix): matrix of reciprocal lattice vectors
    '''
    scfOutput = '%s.out' % oneCalc['_AFLOWPI_PREFIX_'][1:]
    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],scfOutput),'r') as outFile:
            lines = outFile.read()
            re0 = re.compile(r"^.+reciprocal axes: \(\cart. coord. in units 2 .+\)\s*\n^.+\((.+)\)\s*\n^.+\((.+)\)\s*\n^.+\((.+)\)\s*\n",re.MULTILINE)
            re2 = re.compile(r"Begin final coordinates\n.+\n+.+alat= ([\d.]+).+\n",re.MULTILINE)

            alat = re2.findall(lines)
            cellParams = re0.findall(lines)

            if len(cellParams):
                    paramMatrixList = []
                    for params in cellParams[-1]:
                            paramArray = params.split()

                            paramMatrixList.append(paramArray)

                    paramMatrixList = [paramMatrixList[0],paramMatrixList[1],paramMatrixList[2]]
                    paramMatrix =  numpy.matrix(paramMatrixList,dtype='|S10')
                    paramMatrix = paramMatrix.astype(numpy.float)

                    return alat,paramMatrix

            else:
                    print('No card!')
                    return alat,[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
	

def _getAtomNum(inputString,strip=False):
    '''
    Gets the number of atoms from a QE pwscf input file

    Arguments:
          inputString (str): String of a QE pwscf input file

    Keyword Arguments:
          strip (bool): If True then Cr2,Cr45,Cr would be all considered the same species
    
    Returns:
          numOfEach (collections.OrderedDict): dictionary with keys for the species labels and values the number

    '''

    nullstr=''
    ##for the string in the dictionary of dictionaries
    inputSplit = inputString.split('\n')

    '''find all matches to get each of the atoms in the ATOMIC_POSITIONS LIST'''
    atomList =AFLOWpi.qe.regex.atomic_positions(inputString).split('\n')
    atomList = [x.strip(' \r').split()[0] for x in atomList if len(x.strip(' \r'))!=0]

    if strip==True:
        atomList = [x.strip('1234567890') for x in atomList]
    '''put those matches into a dictionary of the species and its count'''
    numOfEach = OrderedDict((i,atomList.count(i)) for i in atomList)

    return numOfEach

from collections import OrderedDict
def _getStoicName(oneCalc,strip=False,latex=False,order=True):
    '''
    Determines the name of the compound by looking at the input and
    finds the stoichiometric number of species in the compound

    Arguments:
          oneCalc (dict): one calculation that is a dictionary of values associated with that calculation

    Keyword Arguments:
          strip (bool): strip species number when it equals 1
          latex (bool): output string in latex format

    Returns:
          name (str): chemical namex

    '''

    inputString = oneCalc['_AFLOWPI_INPUT_']
    numOfEach = AFLOWpi.retr._getAtomNum(inputString,strip=strip) 
    '''get GCD of it'''
    def GCD(a, b):

            if b == 0:
                    return a
            else:
                    return GCD(b, a % b)
    '''cycle through the number of each species for each calculation
    to get the GCD of all the number of species in the calculation'''
    stoicGCD = reduce(GCD, list(numOfEach.values()))

    numOfEachCopy = OrderedDict(numOfEach)
    '''go through and divide the number of each species by the GCD
    and update the dictionary with {species:number of them in the cell}'''

    for species,num in numOfEach.items():
            numOfEachCopy[species] = numOfEach[species] / stoicGCD

    '''builds name for printing in order of elements are listed in the 
    ATOMIC POSTIONS from top of that list is first atom down to bottom 
    of the list is the last species in the name'''
    name = ''

    if latex==True:
        name+='$'

    if order==True:
        iterator = sorted(numOfEachCopy.items())
    else:
        iterator = list(numOfEachCopy.items())

    for key,value in iterator:
        number=value
        if strip==True and value==1:
            number=''

        if latex==True and value!=1:
            name+="%s_{%s}" % (key,number)

        else:
            name+="%s%s" % (key,number)
    '''returns the name of the compound ready to print with LaTeX rendering'''

    if latex==True:
        name+='$'

    return name


def _getPathFromFile(oneCalc):
    '''
    Reads the input file for the band structure calculation and retrieves the 
    k point path from the file

    Arguments:
          oneCalc (dict): one calculation that is a dictionary of values associated with that calculation

    Keyword Arguments:
          None

    Returns:
          path (str): k point path for bands in the bands pwscf input file

    '''

    calcID  = oneCalc[0]
    oneCalc = oneCalc[1]
    bandsIn = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in' % calcID)
    try:
            with open(bandsIn,'r') as bandsInFile:
                    bandsInString = bandsInFile.read()

            path = AFLOWpi.qe.regex.k_points(bandsInString)		



            return path
    except Exception as e:
            print(e)
            print('Did you run ppBands and did it complete properly?')
            return 

###############################################################################
###############################################################################
def _getHighSymPoints(oneCalc,ID=None):
    '''
    Searching for the ibrav number in the input file for the calculation
    to determine the path for the band structure calculation

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
                
    Keyword Arguments:
          ID (str): ID string for the particular calculation and step          

    Returns:
          special_points (list): list of the HSP names
          band_path (str): path in string form

    '''



    ibrav = 0
    ibravRegex = re.compile('ibrav[\s]*=[\s]*([\d]+)\s*[,\n]*')

    ibrav=int(ibravRegex.findall(oneCalc['_AFLOWPI_INPUT_'])[-1])

    if ibrav==0:
            return
    if ibrav < 0:
        print('Lattice type %s is not implemented' % ibrav)
        logging.error('The ibrav value from expresso has not yet been implemented to the framework')
        raise Exception


    alat,cellOld = AFLOWpi.retr._getCellParams(oneCalc,ID)        
    a,b,c,alpha,beta,gamma =  free2abc(cellOld,cosine=False,bohr=False,string=False)
###############################################################################
###############################################################################
    if   ibrav==1:  ibrav_var =  'CUB'
    elif ibrav==2:  ibrav_var =  'FCC'
    elif ibrav==3:  ibrav_var =  'BCC'
    elif ibrav==4:  ibrav_var =  'HEX'
    elif ibrav==6:  ibrav_var =  'TET'
    elif ibrav==8:  ibrav_var =  'ORC'
    elif ibrav==9:  ibrav_var =  'ORCC'
    elif ibrav==11: ibrav_var =  'ORCI'
    elif ibrav==5:
        if alpha < numpy.pi/2.0:   ibrav_var =  'RHL1'
        elif alpha > numpy.pi/2.0: ibrav_var =  'RHL2'
    elif ibrav==7:
        if(c < a):   ibrav_var =  'BCT1'
        elif(c > a): ibrav_var =  'BCT2'
        else:        ibrav_var =  'BCC'
    elif ibrav==10:
        if (1.0/a**2 >1.0/b**2+1.0/c**2):  ibrav_var =  'ORCF1'
        elif (1.0/a**2<1.0/b**2+1.0/c**2): ibrav_var =  'ORCF2'
        else:                              ibrav_var =  'ORCF2'
    elif(int(ibrav)==14):
        minAngle = numpy.amin([alpha,beta,gamma])
        maxAngle = numpy.amax([alpha,beta,gamma])
        if alpha==90.0 or beta==90.0 or gamma==90.0:
            if alpha>=90.0 or beta>=90.0 or gamma>=90.0: ibrav_var =  'TRI2A'
            if alpha<=90.0 or beta<=90.0 or gamma<=90.0: ibrav_var =  'TRI2B'
        elif minAngle>90.0:                              ibrav_var =  'TRI1A'
        elif maxAngle<90:                                ibrav_var =  'TRI1B'
###############################################################################
###############################################################################
    if ibrav_var=='CUB':
        band_path = 'gG-X-M-gG-R-X|M-R'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                           'M'   : (0.5, 0.5, 0.0),
                           'R'   : (0.5, 0.5, 0.5),
                           'X'   : (0.0, 0.5, 0.0)}
                           
    if ibrav_var=='FCC':
        default_band_path = 'gG-X-W-K-gG-L-U-W-L-K|U-X'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                          'K'    : (0.375, 0.375, 0.750),
                          'L'    : (0.5, 0.5, 0.5),
                          'U'    : (0.625, 0.250, 0.625),
                          'W'    : (0.5, 0.25, 0.75),
                          'X'    : (0.5, 0.0, 0.5)}
                          
    if ibrav_var=='BCC':
        band_path = 'gG-H-N-gG-P-H|P-N'
        special_points = {'gG'   : (0, 0, 0),
                          'H'    : (0.5, -0.5, 0.5),
                          'P'    : (0.25, 0.25, 0.25,), 
                          'N'    : (0.0, 0.0, 0.5)}
            
    if ibrav_var=='HEX':
        band_path = 'gG-M-K-gG-A-L-H-A|L-M|K-H'
        special_points = {'gG'   : (0, 0, 0),
                          'A'    : (0.0, 0.0, 0.5),
                          'H'    : (1.0/3.0, 1.0/3.0, 0.5),
                          'K'    : (1.0/3.0, 1.0/3.0, 0.0),
                          'L'    : (0.5, 0.0, 0.5),
                          'M'    : (0.5, 0.0, 0.0)}
        
    if ibrav_var=='RHL1':
        eta1 = 1.0 + 4.0*numpy.cos(alpha)
        eta2 = eta1+1.0
        eta=eta1/eta2
        nu =0.75-eta/2.0
        band_path = 'gG-L-B1|B-Z-gG-X|Q-F-P1-Z|L-P'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                          'B'    : (eta, 0.5, 1.0-eta),
                          'B1'   : (0.5, 1.0-eta, eta-1.0),
                          'F'    : (0.5, 0.5, 0.0),
                          'L'    : (0.5, 0.0, 0.0),
                          'L1'   : (0.0, 0.0, -0.5),
                          'P'    : (eta, nu, nu),
                          'P1'   : (1.0-nu, 1.0-nu, 1.0-eta),
                          'P2'   : (nu, nu, eta-1.0),
                          'Q'    : (1.0-nu, nu, 0.0),
                          'X'    : (nu, 0.0, -nu),
                          'Z'    : (0.5, 0.5, 0.5)}
         
    if ibrav_var=='RHL2':
        eta=1.0/(2*numpy.tan(alpha/2.0)**2)
        nu =0.75-eta/2.0
        band_path = 'gG-P-Z-Q-gG-F-P1-Q1-L-Z'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                          'F'    : (0.5, -0.5, 0.0),
                          'L'    : (0.5, 0.0, 0.0),
                          'P'    : (1.0-nu, -nu, 1.0-nu),
                          'P1'   : (nu, nu-1.0, nu-1.0),
                          'Q'    : (eta, eta, eta),
                          'Q1'   : (1.0-eta, -eta, -eta),
                          'Z'    : (0.5, -0.5, 0.5)} 

    if ibrav_var=='TET':
        band_path = 'gG-X-M-gG-Z-R-A-Z|X-R|M-A'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                          'A'    : (0.5, 0.5, 0.5),
                          'M'    : (0.5, 0.5, 0.0),
                          'R'    : (0.0, 0.5, 0.5),
                          'X'    : (0.0, 0.5, 0.0),
                          'Z'    : (0.0, 0.0, 0.5)}

    if ibrav_var=='BCT1':
       eta = (1+c**2/a**2)/4
       band_path = 'gG-X-M-gG-Z-P-N-Z1-M|X-P'
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'M'     : (-0.5, 0.5, 0.5),
                         'N'     : (0.0, 0.5, 0.0),
                         'P'     : (0.25, 0.25, 0.25),
                         'X'     : (0.0, 0.0, 0.5),
                         'Z'     : (eta, eta, -eta),
                         'Z1'    : (-eta, 1.0-eta, eta)}
         
    if ibrav_var=='BCT2':
       band_path = 'gG-X-Y-gS-gG-Z-gS1-N-P-Y1-Z|X-P'
       eta = (1 + a**2/c**2)/4.0
       zeta = a**2/(2.0*c**2)
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'N'     : (0.0, 0.5, 0.0),
                         'P'     : (0.25, 0.25, 0.25),
                         'gS'    : (-eta, eta, eta),
                         'gS1'   : (eta, 1-eta, -eta),
                         'X'     : (0.0, 0.0, 0.5),
                         'Y'     : (-zeta, zeta, 0.5),
                         'Y1'    : (0.5, 0.5, -zeta),
                         'Z'     : (0.5, 0.5, -0.5)}
         
    if ibrav_var=='ORC':
         band_path = 'gG-X-S-Y-gG-Z-U-R-T-Z|Y-T|U-X|S-R'
         special_points = {'gG'  : (0.0, 0.0, 0.0),
                           'R'   : (0.5, 0.5, 0.5),
                           'S'   : (0.5, 0.5, 0.0),
                           'T'   : (0.0, 0.5, 0.5),
                           'U'   : (0.5, 0.0, 0.5),
                           'X'   : (0.5, 0.0, 0.0),
                           'Y'   : (0.0, 0.5, 0.0),
                           'Z'   : (0.0, 0.0, 0.5)}

    if ibrav_var=='ORCF1':
       band_path = 'gG-Y-T-Z-gG-X-A1-Y|T-X1|X-A-Z|L-gG'
       eta = (1+a**2/b**2+a**2/c**2)/4
       zeta =(1+a**2/b**2-a**2/c**2)/4
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'A'     : (0.5, 0.5 + zeta, zeta),
                         'A1'    : (0.5, 0.5-zeta, 1.0-zeta),
                         'L'     : (0.5, 0.5, 0.5),
                         'T'     : (1.0, 0.5, 0.5),
                         'X'     : (0.0, eta, eta),
                         'X1'    : (1.0, 1.0-eta, 1.0-eta),
                         'Y'     : (0.5, 0.0, 0.5),
                         'Z'     : (0.5, 0.5, 0.0)}

    if ibrav_var=='ORCF2':
       band_path = 'gG-Y-C-D-X-gG-Z-D1-H-C|C1-Z|X-H1|H-Y|L-gG'
       eta = (1+a**2/b**2-a**2/c**2)/4
       phi = (1+c**2/b**2-c**2/a**2)/4
       delta=(1+b**2/a**2-b**2/c**2)/4
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'C'     : (0.5, 0.5-eta, 1.0-eta),
                         'C1'    : (0.5, 0.5+eta, eta),
                         'D'     : (0.5-delta, 0.5, 1.0-delta),
                         'D1'    : (0.5+delta, 0.5, delta),
                         'L'     : (0.5, 0.5, 0.5),
                         'H'     : (1.0-phi, 0.5-phi, 0.5),
                         'H1'    : (phi, 0.5+phi, 0.5),
                         'X'     : (0.0, 0.5, 0.5),
                         'Y'     : (0.5, 0.0, 0.5),
                         'Z'     : (0.5, 0.5, 0.0),}

    if ibrav_var=='ORCF3':
       band_path = 'gG-Y-T-Z-gG-X-A1-Y|X-A-Z|L-R'
       eta =(1+a**2/b**2+a**2/c**2)/4
       zeta=(1+a**2/b**2-a**2/c**2)/4
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'A'     : (0.5, 0.5 + zeta, zeta),
                         'A1'    : (0.5, 0.5-zeta, 1.0-zeta),
                         'L'     : (0.5, 0.5, 0.5),
                         'T'     : (1.0, 0.5, 0.5),
                         'X'     : (0.0, eta, eta),
                         'X1'    : (1.0, 1.0-eta, 1.0-eta),
                         'Y'     : (0.5, 0.0, 0.5),
                         'Z'     : (0.5, 0.5, 0.0)}

    if ibrav_var=='ORCC':
       band_path = 'gG-X-S-R-A-Z-gG-Y-X1-A1-T-Y|Z-T'
       zeta=(1+a**2/b**2)/4.0
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'A'     : (zeta, zeta, 0.5),
                         'A1'    : (-zeta, 1.0-zeta, 0.5),
                         'R'     : (0.0, 0.5, 0.5),
                         'S'     : (0.0, 0.5, 0.0),
                         'T'     : (-0.5, 0.5, 0.5),
                         'X'     : (zeta, zeta, 0.0),
                         'X1'    : (-zeta, 1.0-zeta, 0.0),
                         'Y'     : (-0.5, 0.5, 0.0),
                         'Z'     : (0.0, 0.0, 0.5)}
         
    if ibrav_var=='ORCI':
         band_path = 'gG-X-L-T-W-R-X1-Z-gG-Y-S-W|L1-Y|Y1-Z'
         chi = (1.0 + (a/c)**2)/4.0
         eta = (1.0 + (b/c)**2)/4.0
         delta = (b*b - a*a)/(4*c*c)
         mu = (b*b + a*a)/(4*c*c)
         special_points = {'gG'   : (0, 0, 0),
                           'L'    : (-mu, mu, 0.5-delta),
                           'L1'   : (mu, -mu, 0.5+delta),
                           'L2'   : (0.5-delta, 0.5+delta, -mu),
                           'R'    : (0.0, 0.5, 0.0),
                           'S'    : (0.5, 0.0, 0.0),
                           'T'    : (0.0, 0.0, 0.5),
                           'W'    : (0.25,0.25,0.25),
                           'X'    : (-chi, chi, chi),
                           'X1'   : (chi, 1.0-chi, -chi),
                           'Y'    : (eta, -eta, eta),
                           'Y1'   : (1.0-eta, eta, -eta),
                           'Z'    : (0.5, 0.5, -0.5)}
   
    if ibrav_var=='TRI1A':        
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG' 
        special_points = {'gG'    : (0.0,0.0,0.0),
                          'L'     : (0.5,0.5,0.0),
                          'M'     : (0.0,0.5,0.5),
                          'N'     : (0.5,0.0,0.5),
                          'R'     : (0.5,0.5,0.5),
                          'X'     : (0.5,0.0,0.0),
                          'Y'     : (0.0,0.5,0.0),
                          'Z'     : (0.0,0.0,0.5),}
        
    if ibrav_var=='TRI2A':        
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG'
        special_points = {'gG'    : (0.0,0.0,0.0),
                          'L'     : (0.5,0.5,0.0),
                          'M'     : (0.0,0.5,0.5),
                          'N'     : (0.5,0.0,0.5),
                          'R'     : (0.5,0.5,0.5),
                          'X'     : (0.5,0.0,0.0),
                          'Y'     : (0.0,0.5,0.0),
                          'Z'     : (0.0,0.0,0.5),}
 
    if ibrav_var=='TRI1B':        
        band_path = "X-gG-Y|L-gG-Z|N-gG-M|R-gG"
        special_points = {'gG'    : ( 0.0, 0.0,0.0),
                          'L'     : ( 0.5,-0.5,0.0),
                          'M'     : ( 0.0, 0.0,0.5),
                          'N'     : (-0.5,-0.5,0.5),
                          'R'     : ( 0.0,-0.5,0.5),
                          'X'     : ( 0.0,-0.5,0.0),
                          'Y'     : ( 0.5, 0.0,0.0),
                          'Z'     : (-0.5, 0.0,0.5),}

    if ibrav_var=='TRI2B':        
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG'
        special_points = {'gG'    : ( 0.0, 0.0,0.0),
                          'L'     : ( 0.5,-0.5,0.0),
                          'M'     : ( 0.0, 0.0,0.5),
                          'N'     : (-0.5,-0.5,0.5),
                          'R'     : ( 0.0,-0.5,0.5),
                          'X'     : ( 0.0,-0.5,0.0),
                          'Y'     : ( 0.5, 0.0,0.0),
                          'Z'     : (-0.5, 0.0,0.5),}
            
    if ibrav_var=='MCLC':
       sin_gamma = numpy.sin(numpy.arccos(cos_gamma)) 
       mu        = (1+(b/a)**2.0)/4.0
       delta     = b*c*cos_gamma/(2.0*a**2.0)
       xi        = mu -0.25*+(1.0 - b*cos_gamma/c)*(4.0*(sin_gamma**2.0))
       eta       = 0.5 + 2.0*xi*c*cos_gamma/b
       phi       = 1.0 + xi - 2.0*mu
       psi       = eta - 2.0*delta

    if ibrav_var=='MCLC2':        
        pass

    if ibrav_var=='MCLC5':        
        pass

    '''
    def MCL(cellOld):
       a1 = cellOld[0][0]
       b1 = (cellOld[1][0]**2 + cellOld[1][1]**2)**(0.5)
       c1 = cellOld[2][2]
       gamma = numpy.arctan(cellOld[1][1]/cellOld[1][0])
       myList = [a1, b1, c1]
       c = max(myList)
       a = min(myList)
       myList.remove(c)
       myList.remove(a)
       b = myList[0]
       alpha = gamma

       eta = (1 - b*numpy.cos(alpha)/c)/(2*numpy.sin(alpha)**2)
       nu = 0.5 - eta*c*numpy.cos(alpha)/b
       special_points = {
           'gG'    : (0.0, 0.0, 0.0),
           'A'    : (0.5, 0.5, 0.0),
           'C'    : (0.0, 0.5, 0.5),
           'D'    : (0.5, 0.0, 0.5),
           'D1'    : (0.5, 0.0, -0.5),
           'E'    : (0.5, 0.5, 0.5),
           'H'    : (0.0, eta, 1.0-nu),
           'H1'    : (0.0, 1.0-eta, nu),
           'H2'    : (0.0, eta, -nu),
           'M'    : (0.5, eta, 1.0-nu),
           'M1'    : (0.5, 1.0-eta, nu),
           'M2'    : (0.5, eta, -nu),
           'X'    : (0.0, 0.5, 0.0),
           'Y'    : (0.0, 0.0, 0.5),
           'Y1'    : (0.0, 0.0, -0.5),
           'Z'    : (0.5, 0.0, 0.0)
         }

       default_band_path = 'gG-Y-H-C-E-M1-A-X-H1|M-D-Z|Y-D'
       band_path = default_band_path
       return special_points, band_path
    '''
###############################################################################
###############################################################################           
    aflow_conv = numpy.identity(3)
    qe_conv    = numpy.identity(3)

    if ibrav==2:
        aflow_conv = numpy.asarray([[ 0.0, 1.0, 1.0],[ 1.0, 0.0, 1.0],[ 1.0, 1.0, 0.0]])/2.0                       
        qe_conv    = numpy.asarray([[-1.0, 0.0, 1.0],[ 0.0, 1.0, 1.0],[-1.0, 1.0, 0.0]])/2.0
    if ibrav==3:
        aflow_conv = numpy.asarray([[-1.0, 1.0, 1.0],[ 1.0,-1.0, 1.0],[ 1.0, 1.0,-1.0]])/2.0                    
        qe_conv    = numpy.asarray([[ 1.0, 1.0, 1.0],[-1.0, 1.0, 1.0],[-1.0,-1.0, 1.0]])/2.0                    
    if ibrav==7:
        aflow_conv = numpy.asarray([[-1.0, 1.0, 1.0],[ 1.0,-1.0, 1.0],[ 1.0, 1.0,-1.0]])/2.0
        qe_conv    = numpy.asarray([[ 1.0,-1.0, 1.0],[ 1.0, 1.0, 1.0],[-1.0,-1.0, 1.0]])/2.0
    if ibrav==9:
        aflow_conv = numpy.asarray([[ 1.0,-1.0, 0.0],[ 1.0, 1.0, 0.0],[ 0.0, 0.0, 2.0]])/2.0
        qe_conv    = numpy.asarray([[ 1.0, 1.0, 0.0],[-1.0, 1.0, 0.0],[ 0.0, 0.0, 2.0]])/2.0
    if ibrav==10:
        aflow_conv = numpy.asarray([[ 0.0, 1.0, 1.0],[ 1.0, 0.0, 1.0],[ 1.0, 1.0, 0.0]])/2.0
        qe_conv    = numpy.asarray([[ 1.0, 0.0, 1.0],[ 1.0, 1.0, 0.0],[ 0.0, 1.0, 1.0]])/2.0  
    if ibrav==11:
        aflow_conv = numpy.asarray([[-1.0, 1.0, 1.0],[ 1.0,-1.0, 1.0],[ 1.0, 1.0,-1.0]])/2.0
        qe_conv    = numpy.asarray([[ 1.0, 1.0, 1.0],[-1.0, 1.0, 1.0],[-1.0,-1.0, 1.0]])/2.0

                                   
    for k,v in special_points.items():
        second = (aflow_conv*numpy.linalg.inv(qe_conv))*numpy.matrix(v).T
        special_points[k]=tuple(second.flatten().tolist()[0])


    return special_points, band_path	
#############################################################################################################
#############################################################################################################
#############################################################################################################










def writeInputFromOutput(calcs,replace=False,runlocal=False):
    '''


    Arguments:
	  calcs (dict): dictionary of dictionaries of calculations

    Keyword Arguments:
          replace (bool): If true replace the input with updated atomic positions and cell parameters
          runlocal (bool): run local or write to <ID>.py

    Returns:
          None

    '''

    for ID,oneCalc in calcs.items():
        if runlocal==False:
            AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','AFLOWpi.retr._writeInputFromOutput(oneCalc,ID,replace=%s)\n' % replace)
        else:
            AFLOWpi.retr._writeInputFromOutput(oneCalc,ID,replace=replace)
        
def _writeInputFromOutputString(oneCalc,ID):
    '''
    Generates an input of that step updated with the positions and lattice vectors of its output

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step          

    Keyword Arguments:
          None

    Returns:
          newInput (str): string of that step's input updated with the positions and lattice vectors of its output

    '''

    alatRegex = re.compile(r'(?:CELL_PARAMETERS)\s*\(\s*alat\s*=\s*([0-9.]*)\s*\)',re.MULTILINE)

    try:
        if os.path.exists(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.out')):
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.out'),'r') as outputFile:
                engineOutput=outputFile.read()
        else:
            
            return
    except:
        return
    try:
        coordRegex= AFLOWpi.qe.regex.atomic_positions(engineOutput,'content','regex')
        atomCoord = coordRegex.findall(outputFile)

    except:
        pass
    try:
        alat = alatRegex.findall(engineOutput)[-1]
    except IndexError as e:
        pass
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        print(e)
    
    newInput = ''

    inFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in')
    outFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_new.in')
    with open(inFileName,'r') as inFile:
        inFileString = inFile.read()

    try:
        inFileStringOrig = inFileString
        newInput = inFileString
        nameListRegex = re.compile(r'(&[A-Za-z]+)|(?:\s*(\S*?)\s*=\s*(\S+?)(?:(?:\s*\,\s*)|(?:\s*\n)))')
        namelist = nameListRegex.findall(inFileString)
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        
    try:
        cellParamRegex=AFLOWpi.qe.regex.cell_parameters(engineOutput,'content','regex')
        CELL_PARAMETERS = cellParamRegex.findall(engineOutput)[-1]        
    except:
        pass
        
    try:
        
        tokenizedInput =  AFLOWpi.retr._splitInput(newInput)
        tokenCopy = copy.deepcopy(tokenizedInput)

        for nameList,nameListDict in tokenizedInput.items():
            if type(nameListDict)==type(OrderedDict({'someDict':42})):
                for parameter,value in nameListDict.items():
                    if parameter.upper() =='A' or parameter.upper() =='B' or parameter.upper() =='C' or parameter.upper() =='COSBC' or parameter.upper() =='COSAC' or parameter.upper() =='COSAB':
                        del tokenCopy[nameList][parameter]
                    elif parameter=='ibrav':
                        tokenCopy[nameList][parameter]='0'
            if nameList=='&system':
                tokenCopy[nameList]['celldm(1)']=str(alat)


        try:
            if len(coordRegex.findall(newInput))!=0 and len(coordRegex.findall(engineOutput))!=0:

                coords,flags = detachPosFlags(AFLOWpi.qe.regex.atomic_positions(oneCalc['_AFLOWPI_INPUT_']))
                tokenCopy['ATOMIC_POSITIONS']['__content__']=attachPosFlags(atomCoord,flags)
            else:
                pass
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)

        newInput = AFLOWpi.retr._joinInput(tokenCopy)        
        
        if len(cellParamRegex.findall(engineOutput))!=0:
            newInput = re.sub('ibrav\s*=\s*[0-9]*','ibrav = 0',newInput)
            if len(re.findall(r'(?:celldm\(\s*[2-6]\s*\))\s*',newInput)) !=0:
                newInput = re.sub(r'(?:celldm\(\s*[2-6]\s*\))\s*=\s*[0-9.]*[,]*','' ,newInput)
                newInput = re.sub(r'(?:celldm\(\s*1\s*\))\s*=\s*[0-9.]*[,]*','celldm(1)=%s' % alat ,newInput)

            if len(re.findall('CELL_PARAMETERS',newInput))!=0:
                cellParamRegex = re.compile(r'(?:CELL_PARAMETERS)\s+.+\n((?:[-0-9|.]+)\s+(?:[-0-9|.]+)\s+(?:[-0-9|.]+)\s*\n*)+(?=(?:[A-Z|_|\s]+\n)|)',re.MULTILINE)
                newInput = cellParamRegex.sub(CELL_PARAMETERS,newInput)
            else:
                newInput+='\nCELL_PARAMETERS\n%s' % CELL_PARAMETERS
        else:
            pass
    except IndexError:
        pass
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    newInput = AFLOWpi.prep.remove_blank_lines(newInput)
    return newInput


def _prefixFromInput(inputString):
    '''
    DEFUNCT <CONSIDER FOR REMOVAL>

    Arguments:


    Keyword Arguments:


    Returns:


    '''

    prefix = AFLOWpi.retr._splitInput(inputString)['&control']['prefix']
    prefix=prefix.replace("'","")
    prefix=prefix.replace('"','')

    return prefix
    

def _writeInputFromOutput(oneCalc,ID,replace=False):
    '''
    Writes an input file of that step updated with the positions and lattice vectors of its output

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step          

    Keyword Arguments:
          replace (bool): if True then replace the input with the updated one

    Returns:
          None

    '''

    outFileName = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_new.in')
    newInput = AFLOWpi.retr._writeInputFromOutputString(oneCalc,ID)
    if newInput!=None:
        with open(outFileName,'w') as outFile:
            outFile.write(newInput)


        if replace == True:
            os.rename(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_old.in'))
            os.rename(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_new.in'),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'))

        else:
            pass
        



def _writeEfermi(oneCalc,ID):
    '''
    Grabs the fermi enery or HOMO energy from the output files of the dos calculations
    and converts into rydberg then writes it to file <ID>.efermi

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step          

    Keyword Arguments:
          None

    Returns:
          None

    '''

    Efermi = 0.0

    #if it's not a MP grid don't pull the fermi level
    outfiles = [os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out' % ID)]
    try:

            for outfile in outfiles:
                    with open(outfile,'r') as outFileObj:
                            outFileString = outFileObj.read()	

                            HOMOList = re.findall(r'highest occupied, lowest unoccupied level \(ev\):\s*(.+?)\s',outFileString)	      
                            EfermiList = re.findall(r'\s*the Fermi energy is\s*([.0-9]+)\s*ev',outFileString)	       
                            justHOMOList = re.findall(r'highest occupied level\s*\(ev\)\:\s*([.0-9]+)\s*',outFileString)

                            spin_polarized=re.compile('the spin up/dw Fermi energies are\s*([-]*[0-9.])\s*([-]*[0-9.])')
                            if len(HOMOList):
                                    Efermi=float(HOMOList[-1])
                                    logging.info('highest occupied, lowest unoccupied level: %s eV' % Efermi)
                            elif len(EfermiList):
                                    Efermi=float(EfermiList[-1])
                                    logging.info('the Fermi energy is: %s eV' % Efermi)
                            elif len(justHOMOList):
                                    Efermi=float(justHOMOList[-1])
                                    logging.info('highest occupied level: %s eV' % Efermi)

                            elif len(spin_polarized.findall(outFileString))!=0:
                                Efermi=spin_polarized.findall(outFileString)[-1]



    except Exception as e:
            print(e)
            logging.info('could not find eFermi/HOMO-LUMO in %s checking to see if EFERMI/HOMOLUMO is in calc dictionary from previous calculation.')

    'look to see if efermi is there from dos calc and if not add efermi as part of the dos calc'
    try:
        Efermi=oneCalc['EFERMI']
    except:
        if Efermi!=0.0:
            oneCalc['EFERMI']=Efermi
    '''save value of efermi'''
    fermi_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.efermi'%ID)

    with open(fermi_file,'w') as outFileObj:
            outFileObj.write(str(Efermi))
        

def _getEfermi(oneCalc,ID,directID=False):
    '''
    Grabs fermi level from <ID>.efermi file. if there is none returns 0.0

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step          

    Keyword Arguments:
          None

    Returns:
          efermi (float): fermi level 
   
    '''
    if directID==True:
        glob_ID=ID
    else:
        glob_ID=  AFLOWpi.prep._return_ID(oneCalc,ID,step_type='dos',last=True,straight=False)

    fermi_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.efermi'%glob_ID)
    try:
        with open(fermi_file,'r') as outFileObj:
            efermi = float(outFileObj.read())
            return efermi
    except Exception as e:
        print(e)
        return 0.0


def _joinInput(inputDict):
    '''
    Joins input tokenized by AFLOWpi.retr._splitInput and returns a string of a QE pwscf input file

    Arguments:
          inputDict (dict): tokenized input file

    Keyword Arguments:
          None

    Returns:
          newInputString (str): a string of a QE pwscf input file

    '''

    newInputString = ''
    for namelist,parameters in inputDict.items():
        if '&' in namelist:        
            newInputString+=namelist+'\n'
            for parameter,value in inputDict[namelist].items():
                newInputString+= '   %s = %s,\n' % (parameter.lower(),value)

            newInputString+='/\n'
        else:
            try:
                newInputString+=namelist
                if namelist=='ATOMIC_POSITIONS' and inputDict[namelist]['__modifier__'].strip()=='':
                    inputDict[namelist]['__modifier__']='{crystal}'
                    
                newInputString+=' '+inputDict[namelist]['__modifier__'].lower()+'\n'                
                newInputString+=inputDict[namelist]['__content__']+'\n'
            except:
                pass

    newInputString = AFLOWpi.prep.remove_blank_lines(newInputString)

    return newInputString


def _getPath(dk, oneCalc,ID=None,points=False):
    '''
    Get path between HSP

    Arguments:
          dk (float): distance between points 
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
       

    Keyword Arguments:
          ID (str): ID string for the particular calculation and step          
          points (bool): Give the path as points or in aflow convention

    Returns:
          numPointStr (str): path between HSP

    '''

    def kdistance(hs, p1, p2):
        g = numpy.dot(hs.T, hs)
        p1, p2 = numpy.array(p1), numpy.array(p2)
        d = p1 - p2
        dist2 = numpy.dot(d.T, numpy.dot(g, d).T)
        return numpy.sqrt(dist2)

    def getSegments(path):
        segments = path.split('|')
        return segments

    def getPoints(pathSegment):
        pointsList = pathSegment.split('-')
        return pointsList
    def getNumPoints(path):
        list1 = getSegments(path)
        numPts = 0
        for index in (list1):
            numPts += len(getPoints(index))
        return numPts

    ibrav = 0
    ibravRegex = re.compile('ibrav[\s]*=[\s]*([\d]+)\s*[,\n]*')
    ibrav = int(ibravRegex.findall(oneCalc['_AFLOWPI_INPUT_'])[0])	

    if ibrav==0:
            return
    if ibrav<0:
        print('Lattice type %s is not implemented' % ibrav)
        logging.error('The ibrav value from expresso has not yet been implemented to the framework')
        raise Exception

    alat,cell = AFLOWpi.retr._getCellParams(oneCalc,ID)



    totalK=0
    special_points, band_path = AFLOWpi.retr._getHighSymPoints(oneCalc,ID=ID)
    hs = 2*numpy.pi*numpy.linalg.inv(cell)  # reciprocal lattice
    segs = getSegments(band_path)

    numPointsStr = str(getNumPoints(band_path)) + '\n' #starts the string we need to make
    if points==True:
        numPointsStr=''

    path_points=''
    for index in segs:

        a = getPoints(index) #gets the points in each segment of path spereated by |            
        point1 = None
        point2 = None

        for index2 in range(len(a)-1):
            try:
                point1 = a[index2]
                point2 = a[index2+1]
                p1 = special_points[point1] 
                p2 = special_points[point2]

                newDK = (2*numpy.pi/alat)*dk
                numK = int(numpy.ceil((kdistance(hs, p1, p2)/newDK)))
                totalK+=numK

                numK = str(numK)

                if points==True:
                    a0 = numpy.linspace(p1[0],p2[0],numK).astype(numpy.float16)
                    a1 = numpy.linspace(p1[1],p2[1],numK).astype(numpy.float16)
                    a2 = numpy.linspace(p1[2],p2[2],numK).astype(numpy.float16)



                    for x in range(len(a0)):
                        lbl=point1
                        if x==(int(numK)-1):
                            lbl=point2
                        numPointsStr += '%8.7f %8.7f %8.7f !! %s\n'%(a0[x],a1[x],a2[x],lbl)



                else:
                    sp1 = special_points[point1]
                    sp2 = special_points[point2]

                    numPointsStr += '%s %s %s %s ! %s\n' % (sp1[0],sp1[1],sp1[2],numK,point1.rjust(2))
                    if(index2 == len(a)-2):
                        numPointsStr += '%s %s %s %s ! %s\n' % (sp2[0],sp2[1],sp2[2],str(0),point2.rjust(2))

            except Exception as e:
                print(e)
    if points==True:

        numPointsStr=numPointsStr.replace('!!','%5.4f !!' % (2.0/float(totalK)))
        numPointsStr='%s\n'%totalK+numPointsStr
    return numPointsStr


def _joinMatrixLabels(labels,matrix):
    '''
    Joins a list of the atomic position species labels with a list or array of their atomic positions

    Arguments:
          labels (list): list or array of the atomic positions species labels
          matrix (numpy.matrix): the positions in matrix, array, or list form

    Keyword Arguments:
          None

    Returns:
          outString (str): a string of the atomic positions joined with their species labels

    '''
    
    matrixString = AFLOWpi.retr._cellMatrixToString(matrix)

    outString = '\n'.join( [' '.join(x) for x  in zip(labels,matrixString.split('\n'))])
    return outString

def _orderSplitInput(inputCalc):
    '''
    Order tokenized input as QE needs it

    Arguments:
          inputCalc (calc): tokenized QE pwscf input

    Keyword Arguments:
          None

    Returns:
          newOrderedDict (dict): the ordered tokenized QE pwscf input

    '''

    newOrderedDict=OrderedDict()
    inputOrder = ['&control','&system','&electrons','&ions','&cell','ATOMIC_SPECIES','ATOMIC_POSITIONS','K_POINTS']
    for item in inputOrder:
        if item in list(inputCalc.keys()):
            newOrderedDict.update({item:inputCalc[item]})
        else:
            newOrderedDict.update({item:OrderedDict()})
    for key in list(inputCalc.keys()):
        if key not in inputOrder:
            newOrderedDict.update({key:inputCalc[key]})

    return newOrderedDict

def _splitInput(inFileString):
    '''
    Tokenizes the QE pwscf input file

    Arguments:
          inFileString (str): string of a QE pwscf input file

    Keyword Arguments:
          None

    Returns:
          inputDict (dict): the tokenized input

    '''

    inputDict=OrderedDict()
    nameListRegex = re.compile(r'\s*(&[A-Za-z]+)|(?:\s*(\S*?)\s*=\s*(\S+?)(?:(?:\s*\,\s*)|(?:\s*\n)))')
    namelist = nameListRegex.findall(inFileString)

    for item in range(len(namelist)):
        if len(namelist[item][0])!=0:
            lastOne = namelist[item][0].lower()
            inputDict[lastOne]=OrderedDict()
        else:
            try:                    
                if lastOne not in list(inputDict.keys()):
                    inputDict[lastOne.lower()]=OrderedDict()
                inputDict[lastOne.lower()][namelist[item][1].lower()]=namelist[item][2]
            except Exception as e:
                AFLOWpi.run._fancy_error_log(e)


    atomicSpecRegex = re.compile(r'ATOMIC_SPECIES\s*\n((?:(?:\s*[A-Za-z0-9]+)\s+(?:[0-9.]+)\s+(?:.+)\n*)+)(?=(?:[A-Z_\s]+\n)|)',(re.MULTILINE))
    try: 
        inputDict['ATOMIC_SPECIES']=OrderedDict()

        inputDict['ATOMIC_SPECIES']['__content__']=atomicSpecRegex.findall(inFileString)[-1]
        inputDict['ATOMIC_SPECIES']['__modifier__']=''

    except:
        inputDict['ATOMIC_SPECIES']['__content__']='\n\n'
        inputDict['ATOMIC_SPECIES']['__modifier__']=''



    modifier_regex = AFLOWpi.qe.regex.atomic_positions('','modifier','regex')
    if len(modifier_regex.findall(inFileString)):

        paramUnits=modifier_regex.findall(inFileString)[-1]

        if len(paramUnits.strip()):
            paramUnits = paramUnits.strip()
        else:
            paramUnits=' '
        paramUnits = paramUnits.lower()
        if len(paramUnits.strip()):
            paramUnits = '{%s}' % paramUnits
        else:
            paramUnits=''
    else:
        paramUnits=''

    coordRegex = AFLOWpi.qe.regex.atomic_positions(inFileString,'content','regex')

    try:
        coords = coordRegex.findall(inFileString)[-1]

        inputDict['ATOMIC_POSITIONS']=OrderedDict()
        inputDict['ATOMIC_POSITIONS']['__content__']=coords
        inputDict['ATOMIC_POSITIONS']['__modifier__']=paramUnits

    except:
        pass

    cellParamRegex = AFLOWpi.qe.regex.cell_parameters('','content','regex')

    try:
        unitsRegex = AFLOWpi.qe.regex.cell_parameters('','modifier','regex')

        if len(cellParamRegex.findall(inFileString)):
            paramUnits=unitsRegex.findall(inFileString)[-1].strip()
            paramUnits = paramUnits.lower()
            if paramUnits in ['cubic','crystal','angstrom','bohr']:
                if len(paramUnits.strip()):
                    paramUnits = '{%s}' % paramUnits
                else:
                    paramUnits=''
            else:
                paramUnits=''


        cellParams = cellParamRegex.findall(inFileString)
        if len(cellParams)!=0:
            cellParams=cellParams[-1]
            inputDict['CELL_PARAMETERS']=OrderedDict()
            inputDict['CELL_PARAMETERS']['__content__']=cellParams
            inputDict['CELL_PARAMETERS']['__modifier__']=paramUnits



    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    try:
        kpoints = AFLOWpi.qe.regex.k_points(inFileString)

        try:            
            modifier=AFLOWpi.qe.regex.k_points(inFileString,'modifier')
            if len(modifier.strip()):
                modifier = '{%s}' % modifier
            else:
                modifier=''

            kpointInput = modifier+'\n'+kpoints
            inputDict['K_POINTS']=OrderedDict()
            inputDict['K_POINTS']['__modifier__']=modifier
            inputDict['K_POINTS']['__content__']=kpoints
        except Exception as e:
            modifier=''
            print(e)


    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)


    return inputDict


def getBravaisLatticeName(bravaisLatticeNumber):
    '''
    Returns the name of the crystal system from the bravias lattice number

    Arguments:
          braviasLatticeNumber (int): the number of the bravais lattice in QE convention

    Keyword Arguments:
          None

    Returns:
          A string of the name of the crystal system

    '''        

    if ibrav in (1,2,3):
        return 'cubic'
    elif ibrav in (4):
        return 'rhombohedral'
    elif ibrav in (5,-5):
        return 'hexagonal'
    elif ibrav in (6,7):
        return 'tetragonal'
    elif ibrav in (8,9,-9,10,11):
        return 'orthorhombic'
    elif ibrav in (12,-12,13,14):
        return 'triclinic'
    else:
        logging.error('Not a valid bravais lattice number')
        print('Not a valid bravais lattice number')
        return ''


def pw2cif(calcs,inpt=True,outp=True,runlocal=False,outputFolder=None,filePrefix=''):
    '''
    Writes a simple CIF file for engine input and outputs for viewing the structure

    Arguments:
	  calcs (dict): dictionary of dictionaries of calculations

    Keyword Arguments:
          inpt (bool): Do CIF for input
          outp (bool): Do CIF for output
          outputFolder (str): Output directory for the CIF
          filePrefix (str): Optional prefix to the CIF filenames

    Returns:
          None

    '''        
    if inpt==True and outp==True:
        inpt_or_outp="'both'"
    elif inpt==True and output==False:
        inpt_or_outp="'False'"
    elif inpt==False and output==True:
        inpt_or_outp="'output'"
    else:
        inpt_or_outp='""'

    try:
        for ID,oneCalc in calcs.items():
            if runlocal:
                inOrOut=eval(inpt_or_outp)
                try:
                    AFLOWpi.retr._pw2cif(oneCalc,ID,inOrOut=inOrOut,outputFolder=outputFolder,filePrefix=filePrefix)
                except Exception as e:
                    AFLOWpi.run._fancy_error_log(e)
                    continue
            else:
                inOrOut=inpt_or_outp
                if outputFolder==None:

                    AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','AFLOWpi.retr._pw2cif(oneCalc,ID,inOrOut=%s,outputFolder=None,filePrefix="%s")\n' % (inOrOut,filePrefix))
                else:
                    AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','AFLOWpi.retr._pw2cif(oneCalc,ID,inOrOut=%s,outputFolder="%s",filePrefix="%s")\n' % (inOrOut,outputFolder,filePrefix))

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

def _pw2cif(oneCalc,ID,inOrOut='input',outputFolder=None,filePrefix=''):
    '''
    Generates the CIF from input or output of single calculation at that step

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step         

    Keyword Arguments:
          inOrOut (str): which structs to do CIF for.. 'input' 'output' or 'both'
          outputFolder (str): Output directory for the CIF
          filePrefix (str): Optional prefix to the CIF filenames

    Returns:
          None

    '''        

    def atomPosPretty(atomicPositionString,ibrav):
            
            a,b,c,alpha,beta,gamma = (0.0,0.0,0.0,0.0,0.0,0.0)
            atomPosSortList = []
            try:
                splitAtomicPos = atomicPositionString.split('\n')                
            except ValueError:
                pass


            for item in splitAtomicPos:
                
                    if len(item.strip())!=0:
                            atomSplit = item.split()
                            formattedAtomPos = atomSplit[0]+' '+' '.join(['%14.12f' % float(atomSplit[i]) for i in range(1,len(atomSplit))])+'\n'
                            atomPosSortList.append(formattedAtomPos)


            atomPosSortCifList = []
            atomPosSortCifListSecond = []
            countDict = {}
            for item in atomPosSortList:
                splitItem = item.split()
                firstItem = splitItem[0]
                if firstItem not in list(countDict.keys()):
                    countDict[firstItem] = 1
                else:
                    count = countDict[firstItem]+1
                    countDict[firstItem]=count

            countDictOrig = copy.deepcopy(countDict)
            origDictSecond = copy.deepcopy(countDict)

            for item1 in atomPosSortList:
                try:
                    splitItem = item1.split()
                    firstItem = splitItem[0]

                    
                    secondItem = '%s%s' % (firstItem,countDictOrig[firstItem]-countDict[firstItem]+1)

                    

                    try:
                        entry = '%s 1.0 %04f %04f %04f Biso 1.0 %s' % (secondItem,float(splitItem[1]),float(splitItem[2]),float(splitItem[3]),firstItem)
                        atomPosSortCifList.append(entry)
                    except:
                        pass
                except Exception as e:
                    AFLOWpi.run._fancy_error_log(e)

                    
                """
                decriment the counter for the atom in question
                """
                countDict[firstItem] = countDict[firstItem]-1                

########################################################################################
########################################################################################
########################################################################################
            loopblock = """loop_
      _atom_site_label
      _atom_site_occupancy
      _atom_site_fract_x
      _atom_site_fract_y
      _atom_site_fract_z
      _atom_site_thermal_displace_type
      _atom_site_B_iso_or_equiv
      _atom_site_type_symbol
"""
            
            outputString = loopblock+'\n'.join(atomPosSortCifList)

            return outputString

    if inOrOut != 'input' and inOrOut != 'output' and inOrOut != 'both':
        print("options for variable 'inOrOut' must be 'in','out',or 'both'")
        logging.warning("options for variable 'inOrOut' must be 'in','out',or 'both'")
        return
    if filePrefix!='':
        filePrefix=filePrefix+'_'
    input_suffix=''
    cif_suffix=''
    cif_prefix=''

    cif_prefix=__getStoicName(oneCalc)+'_'
    outFileString=''
    inputFileString  = oneCalc['_AFLOWPI_INPUT_']
    
    ibravRegex = re.compile('ibrav[\s]*=[\s]*([\d]+)\s*[,\n]*')
    a,b,c,alpha,beta,gamma = (0.0,0.0,0.0,0.0,0.0,0.0)

    if inOrOut.lower()=='input' or inOrOut.lower()=='both':

        labels,atomicPositions,cellParamVec = AFLOWpi.retr._getConventionalCellFromInput(oneCalc,ID)

        input_suffix='.in'
        cif_suffix='.in.cif'
        

        inputDict = AFLOWpi.retr._splitInput(inputFileString)
        ibravRegex= re.compile(r'ibrav\s*=\s*([0-9]{1,2})')

        try:
            ibravOrig= ibravRegex.findall(inputFileString)[-1]
            ibrav=int(ibravOrig)
            if ibrav==5:

                labels,atomicPositions,cellParamVec = AFLOWpi.retr._rho2hex(cellParamVec,atomicPositions,labels)
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)


        a,b,c,alpha,beta,gamma =  free2abc(cellParamVec,cosine=False,degrees=True,string=False,)
        
        outFileString+="""data_global
_chemical_name '%s'

_cell_length_a %s
_cell_length_b %s
_cell_length_c %s
_cell_angle_alpha %s
_cell_angle_beta %s
_cell_angle_gamma %s
""" % (cif_prefix[:-1],a,b,c,alpha,beta,gamma)

        atomicPositions = AFLOWpi.retr._cellMatrixToString(atomicPositions)
        atomicPositions = '\n'.join( [' '.join(x) for x  in zip(labels,atomicPositions.split('\n'))])        
        outFileString += atomPosPretty(atomicPositions,ibravOrig)
        outFileString += '\n'

        if outputFolder == None:
            outputFolder = oneCalc['_AFLOWPI_FOLDER_']
        with open(os.path.join(outputFolder,filePrefix+'STRUCT_'+cif_prefix+ID+cif_suffix),'w') as outFile:
            outFile.write(outFileString)

    if inOrOut.lower()=='output' or inOrOut.lower()=='both':
        input_suffix='.out'
        cif_suffix='.out.cif'
        inputDict = {}
        outFileString=''




        try:
#            labels,positions,cell = AFLOWpi.retr._getConventionalCellFromOutput(oneCalc,ID)
            labels=AFLOWpi.retr._getPosLabels(inputString)
            alat,cellParamMatrix = AFLOWpi.retr._getCellParams(oneCalc,ID)
            cell=cellParamMatrix*float(alat)
            positions = AFLOWpi.retr.getPositionsFromOutput(oneCalc,ID)
            

            ibrav = AFLOWpi.retr.getIbravFromVectors(cell)        
            ibrav=int(ibrav)
            positions= positions.getA()

        except Exception as e:

            AFLOWpi.run._fancy_error_log(e)

            try:
                with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'r') as inFileObj:
                     exceptionInputFileString  = inFileObj.read()



            except:
                print("Error finding ATOMIC_POSITIONS from output")
                logging.warning("Error finding ATOMIC_POSITIONS from output")
                return
        try:
            cellParamVec = AFLOWpi.retr._getCellOutDim(oneCalc,ID)
        except:
            return


        a = cellParamVec['a']
        b  = cellParamVec['b']
        c  = cellParamVec['c']
        alpha  = cellParamVec['alpha']
        beta = cellParamVec['beta']

        gamma = cellParamVec['gamma']
        a,b,c,alpha,beta,gamma =  free2abc(cell,cosine=False,degrees=True,string=False,)        



        outFileString+="""data_global
_chemical_name '%s'
_cell_length_a %s
_cell_length_b %s
_cell_length_c %s
_cell_angle_alpha %s
_cell_angle_beta %s
_cell_angle_gamma %s
""" % (cif_prefix[:-1],a,b,c,alpha,beta,gamma)
        
        positions= AFLOWpi.retr._joinMatrixLabels(labels,positions)


        outFileString += atomPosPretty(positions,ibrav)
        outFileString += '\n'

        if outputFolder == None:
            outputFolder = oneCalc['_AFLOWPI_FOLDER_']
        joiner = filePrefix+'STRUCT_'+cif_prefix+ID+cif_suffix

        with open(os.path.join(outputFolder,joiner),'w') as outFile:
            outFile.write(outFileString)



def celldm2params(a,b,c,alpha,beta,gamma,k,v):
    '''
    DEFUNCT <Marked for removal>

    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    '''not sure what this is. probably is not needed anymore much and should be removed'''
    if re.match(r'\s*celldm.*1.*',k):
        a=float(v)
    if re.match(r'\s*celldm.*2.*',k):
        b=float(v)*a
    if re.match(r'\s*celldm.*3.*',k):
        c=float(v)*a
    if re.match(r'\s*celldm.*4.*',k):
        alpha=numpy.arccos(float(v))*180.0/numpy.pi
    if re.match(r'\s*celldm.*5.*',k):
        beta=numpy.arccos(float(v))*180.0/numpy.pi
    if re.match(r'\s*celldm.*6.*',k):
        gamma=numpy.arccos(float(v))*180.0/numpy.pi



    a=numpy.around(a,decimals=5)
    b=numpy.around(b,decimals=5)
    c=numpy.around(c,decimals=5)
    alpha=numpy.around(alpha,decimals=5)
    beta=numpy.around(beta,decimals=5)
    gamma=numpy.around(gamma,decimals=5)


    return a,b,c,alpha,beta,gamma


def inputDict2params(inputDict):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    ibrav=0
    for k,v in sorted(inputDict.items()):

        if re.match(r'ibrav',k):
            ibrav=int(v)

    a,b,c,alpha,beta,gamma,k = (0.0,0.0,0.0,0.0,0.0,0.0,k)
    for k,v in sorted(inputDict.items()):
        a,b,c,alpha,beta,gamma = celldm2params(a,b,c,alpha,beta,gamma,k,v)

    if ibrav == 1:
        b,c = (a,a)
        alpha,beta,gamma = (90.0,90.0,90.0)
    if ibrav == 2: 
        b,c = (a,a)
        alpha,beta,gamma = (90.0,90.0,90.0)
    if ibrav == 3: 
        b,c = (a,a)
        alpha,beta,gamma = (90.0,90.0,90.0)
    if ibrav == 4:
        b,c = (a,a)
        beta = alpha
        gamma = 120.0
    if ibrav == 5 or ibrav == -5:
        b = a
#        beta = alpha
#        gamma = 120.0
    if ibrav == 6:
        b = a
        alpha,beta,gamma = (90.0,90.0,90.0)
    if ibrav == 7:
        b = a
        alpha,beta,gamma = (90.0,90.0,90.0)
    if ibrav == 8:
        alpha,beta,gamma = (90.0,90.0,90.0)
    if ibrav == 9 or ibrav == -9:
        alpha,beta,gamma = (90.0,90.0,90.0)
    if ibrav == 10:
        alpha,beta,gamma = (90.0,90.0,90.0)                    
    if ibrav == 11:
        alpha,beta,gamma = (90.0,90.0,90.0)


    a=numpy.around(a,decimals=5)
    b=numpy.around(b,decimals=5)
    c=numpy.around(c,decimals=5)
    alpha=numpy.around(alpha,decimals=5)
    beta=numpy.around(beta,decimals=5)
    gamma=numpy.around(gamma,decimals=5)

    return a,b,c,alpha,beta,gamma

def _getCellOutDim(oneCalc,ID,cosine=True,degrees=True):
    '''
    Returns the cell parameters a,b,c,alpha,beta,gamma of the primitive lattice from the output.
    If the primitive lattice does not change during the calculation it takes it from the input

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step         

    Keyword Arguments:
          None

    Returns:
         out_params (dict): dictionary of a,b,c,alpha,beta,gamma of the primitive lattice from the output

    '''        

    ibravRegex= re.compile(r'ibrav\s*=\s*([0-9]{1,2})')
    inputFileString  = oneCalc['_AFLOWPI_INPUT_']
    try:
        ibravOrig= ibravRegex.findall(inputFileString)[-1]
    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.out'),'r') as inFileObj:
        inFileString  = inFileObj.read()

    cellParamRegex = re.compile(r'(?:CELL_PARAMETERS)\s*\(alat\s*=\s*[0-9.]*\)\n([.0-9\s-]*)',re.MULTILINE)
    try:

        splitParams =  cellParamRegex.findall(inFileString)[-1].split('\n')
    except Exception as e:
        try:
            cellParamRegex = re.compile(r'crystal axes: \(cart. coord. in units of alat\)\n\s*a\(1\)\s*=\s*\(\s*([0-9.\s]*)\)\s*\n\s*a\(2\)\s*=\s*\(\s*([0-9.\s]*)\)\s*\n\s*a\(3\)\s*=\s*\(\s*([0-9.\s]*)\)\s*\n')
            splitParams =  cellParamRegex.findall(inFileString)[-1]
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)

    alat,cellOld = AFLOWpi.retr._getCellParams(oneCalc,ID)

    paramList = []
    alat,cellParaMatrix=__getCellParams(oneCalc,ID)
    cellParaMatrix*=alat

    alatRegex = re.compile(r'(?:CELL_PARAMETERS)\s*\(alat\s*=\s*([0-9.]*)\)',re.MULTILINE)

    splitInput  = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

    string= free2abc(cellParaMatrix)
    cellParams=[float(x.split('=')[1]) for x in string.split(',')]
    for x in range(len(cellParams)):
        if cellParams[x]<0.0000000001:
            cellParams[x]=0.0



        

    ibravOrig = int(splitInput['&system']['ibrav'])


    celldmDict = {}
    celldmDict['ibrav']=ibravOrig

    out_params = {'a':cellParams[0],'b':cellParams[1],'c':cellParams[2],'alpha':cellParams[3],'beta':cellParams[4],'gamma':cellParams[5]}

    if cosine==False:
        out_params['alpha']=numpy.arccos(out_params['alpha'])
        out_params['beta']=numpy.arccos(out_params['beta'])
        out_params['gamma']=numpy.arccos(out_params['gamma'])
        
        if degrees==True:
            out_params['alpha']*=180.0/numpy.pi
            out_params['beta']*=180.0/numpy.pi
            out_params['gamma']*=180.0/numpy.pi


    for k,v in out_params.items():
        out_params[k]=float('%10.5e'%numpy.around(v,decimals=5))

    return out_params

def _getInputParams(oneCalc,ID):
    '''
    Returns the cell parameters a,b,c,alpha,beta,gamma of the primitive lattice from the input.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step         

    Keyword Arguments:
          None

    Returns:
         out_params (dict): dictionary of a,b,c,alpha,beta,gamma of the primitive lattice from the input

    '''        
        

    inputFileString  = oneCalc['_AFLOWPI_INPUT_']

    inputDict = AFLOWpi.retr._splitInput(inputFileString)['&system']
    
    for key in list(inputDict.keys()):
        if re.match(r'celldm',key): 
            pass
        elif re.match(r'ibrav',key): 
            pass
        else:
            del inputDict[key]

    a,b,c,alpha,beta,gamma = inputDict2params(inputDict)

    BohrToAngstrom = 0.529177249
    a*= BohrToAngstrom
    b*= BohrToAngstrom
    c*= BohrToAngstrom

    return {'a':float(a),'b':float(b),'c':float(c),'alpha':float(alpha),'beta':float(beta),'gamma':float(gamma)}

def celldm2free(ibrav=None,celldm1=None,celldm2=None,celldm3=None,celldm4=None,celldm5=None,celldm6=None,returnString=True):
    '''
    Converts  QE's celldm format into QE convention primitive lattice vectors 

    Arguments:
          None

    Keyword Arguments:
           ibrav (int): bravais lattice type by QE convention
           celldm1 (float): a
           celldm2 (float): b/a
           celldm3 (float): c/a
           celldm4 (float): Cosine of the angle between axis b and c
           celldm5 (float): Cosine of the angle between axis a and c
           celldm6 (float): Cosine of the angle between axis a and b
           returnString (bool): return vectors as a string or as a numpy matrix

    Returns:
         Either a numpy.matrix or a string of the primitive lattice vectors generated from the celldm input

    '''        


    if int(ibrav)==0 or abs(int(ibrav))>14:
        logging.error('ibrav = %s not a valid option' % ibrav)
        print('ibrav = %s not a valid option' % ibrav)

    ibrav=int(ibrav)




    if ibrav==1:
        matrix = numpy.matrix(((celldm1,0      ,0),
                               (0      ,celldm1,0),
                               (0     ,0       ,celldm1)))
    

    if ibrav==2:
        matrix = numpy.matrix(((-celldm1  ,0      ,celldm1),
                               ( 0        ,celldm1,celldm1),
                               (-celldm1  ,celldm1,0      )))/2.0
        
    if ibrav==3:
        matrix = numpy.matrix((( celldm1, celldm1, celldm1),
                               (-celldm1, celldm1, celldm1),
                               (-celldm1,-celldm1, celldm1)))/2.0
        
    if ibrav==4:
        matrix = numpy.matrix((( 1.        ,0.               , 0.      ),
                               (-0.5       ,numpy.sqrt(3)/2  , 0.      ),
                               ( 0.        ,0.               , celldm3)))*celldm1

        # if numpy.abs(celldm4-celldm5)<0.01 and numpy.abs(celldm4-celldm6)<0.01 and numpy.abs(celldm6-celldm5)<0.01:
        #     c=celldm4
        #     tx=numpy.sqrt((1-c)/2)
        #     ty=numpy.sqrt((1-c)/6)
        #     tz=numpy.sqrt((1+2*c)/3)

        #     matrix = celldm1*numpy.matrix(((tx  ,-ty  , tz),
        #                                    (0   , 2*ty, tz),
        #                                    (-tx ,-ty  , tz)))

    if ibrav==5:
        print(celldm1,celldm2,celldm3,celldm4,celldm5,celldm6, end=' ')
        c=celldm4
        tx=numpy.sqrt((1-c)/2)
        ty=numpy.sqrt((1-c)/6)
        tz=numpy.sqrt((1+2*c)/3)

        matrix = celldm1*numpy.matrix(((tx  ,-ty  , tz),
                                       (0   , 2*ty, tz),
                                       (-tx ,-ty  , tz)))




        # if numpy.abs(celldm4+0.5)<0.01 or numpy.abs(celldm5+0.5)<0.01 or numpy.abs(celldm6+0.5)<0.01:
        #     matrix = numpy.matrix((( 1.        ,0.               , 0.      ),
        #                            (-0.5       ,numpy.sqrt(3)/2  , 0.      ),
        #                            ( 0.        ,0.               , celldm3)))*celldm1



                    
    if ibrav==-5:
        c=celldm4
        tx=numpy.sqrt((1-c)/2)
        ty=numpy.sqrt((1-c)/6)
        tz=numpy.sqrt((1+2*c)/3)
        u = tz - 2*numpy.sqrt(2)*ty
        v = tz + numpy.sqrt(2)*ty
        matrix = celldm1/numpy.sqrt(3)*numpy.matrix(((u,v,v),
                                 (v,u,v),
                                                     (v,v,u)))

    if ibrav==6:
        matrix=celldm1*numpy.matrix(((1,0,0),
                                     (0,1,0),
                                     (0,0,celldm3)))

    if ibrav==7:
        matrix = numpy.matrix((
                ( 1.0,-1.0,celldm3),
                ( 1.0, 1.0,celldm3),
                (-1.0,-1.0,celldm3),
                ))*celldm1
        matrix/=2.0

    if ibrav==8:
        matrix = celldm1*numpy.matrix(((1.0, 0.0     , 0.0    ,),
                                       (0.0, celldm2 , 0.0    ,),
                                       (0.0, 0.0     , celldm3,)))
    try:
        a=celldm1
    except:
        pass
    try:
        b=celldm2*celldm1
    except:
        pass
    try:
        c=celldm3*celldm1
    except:
        pass

    if ibrav==9:
        matrix=numpy.matrix((( a/2, b/2,0),
                             (-a/2, b/2,0),
                             ( 0,   0,  c)))

    if ibrav==-9:
        matrix=numpy.matrix(((a/2,-b/2,0),
                             (a/2, b/2,0),
                             (0,   0,  c)))

    if ibrav==10:
        matrix = numpy.matrix(((a,   0.,  c   ),
                               (a,   b,   0.  ),  
                               (0.,  b,   c   )))/2.0




    if ibrav==11:
        matrix=numpy.matrix((( a/2, b/2, c/2), 
                             (-a/2, b/2, c/2), 
                             (-a/2,-b/2, c/2)))

    if ibrav==12:
        gamma=numpy.arccos(celldm4)
        matrix = numpy.matrix(((a,           0,           0),
                               (b*numpy.cos(gamma),b*numpy.sin(gamma),0),
                               (0,           0,           c)))

    if ibrav==13:
        print(celldm1,celldm2,celldm3,celldm4,celldm5,celldm6, end=' ')
        gamma=numpy.arccos(celldm4)
        matrix=numpy.matrix(((a/2,                0,                -c/2),
                             (b*numpy.cos(gamma), b*numpy.sin(gamma),0),
                             (a/2,                0,                 c/2)))

    if ibrav==14:
        alpha=numpy.arccos(celldm4)
        beta=numpy.arccos(celldm5)
        gamma=numpy.arccos(celldm6)

        twoone = c*(numpy.cos(alpha)-numpy.cos(beta)*numpy.cos(gamma))/numpy.sin(gamma)
        insert =  1 + 2*celldm4*celldm5*celldm6-celldm4**2-celldm5**2-celldm6**2
        twotwo = c*numpy.sqrt(insert)/numpy.sin(gamma)

        matrix=numpy.matrix(((a,                  0,                  0),
                             (b*numpy.cos(gamma), b*numpy.sin(gamma), 0),
                             (c*numpy.cos(beta),  twoone,             twotwo)))


    if returnString==True:
        vecString = '''{:^6.10f} {:^6.10f} {:^6.10f}
{:6^.10f} {:6^.10f} {:6^.10f}
{:6^.10f} {:6^.10f} {:6^.10f}
'''.format(float(matrix.item(0,0)),float(matrix.item(0,1)),float(matrix.item(0,2)),float(matrix.item(1,0)),float(matrix.item(1,1)),float(matrix.item(1,2)),float(matrix.item(2,0)),float(matrix.item(2,1)),float(matrix.item(2,2)))

        return vecString
    else:
        return matrix

def _pw2aflowConvention(cellParamMatrix,symMatrix):
    '''
    DEFUNCT <CONSIDER FOR REMOVAL>

    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    '''transform the positions first with the pre-transformed vectors'''
    newSymMatrix       = AFLOWpi.retr._pw2aflowPositions(cellParamMatrix,symMatrix)
    '''now get the transformed cell vectors'''
    newCellParamMatrix = AFLOWpi.retr._pw2aflowConventionVec(cellParamMatrix)
    
    return newCellParamMatrix,newSymMatrix

def _pw2aflowConventionVec(cellParamMatrix):
    '''
    DEFUNCT <CONSIDER FOR REMOVAL>

    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    returnMatrix = copy.deepcopy(cellParamMatrix)
    returnMatrix*= 0.529177249
    temp=[]
    temp.append(numpy.copy(returnMatrix[0]))
    temp.append(numpy.copy(returnMatrix[1]))
    temp.append(numpy.copy(returnMatrix[2]))

    lengths = ((numpy.sqrt(temp[0].dot(temp[0].T)),0),(numpy.sqrt(temp[1].dot(temp[1].T)),1),(numpy.sqrt(temp[2].dot(temp[2].T)),2))
    lengths = sorted(lengths)

    tempIndex= [x[1] for x in lengths]

    returnMatrix[0] = numpy.copy(cellParamMatrix[tempIndex[0]])
    returnMatrix[1] = numpy.copy(cellParamMatrix[tempIndex[1]])
    returnMatrix[2] = numpy.copy(cellParamMatrix[tempIndex[2]])

    return returnMatrix



def _pw2aflowPositions(cellParamMatrix,symMatrix):
    '''
    DEFUNCT <CONSIDER FOR REMOVAL>    

    Arguments:


    Keyword Arguments:


    Returns:


    '''        


    cellMatrix = copy.deepcopy(cellParamMatrix)
    cellMatrix*= 0.529177249

    temp=[]
    tempPos=[]
    temp.append(numpy.copy(cellMatrix[0]))
    temp.append(numpy.copy(cellMatrix[1]))
    temp.append(numpy.copy(cellMatrix[2]))

    lengths = ((numpy.sqrt(temp[0].dot(temp[0].T)),0),(numpy.sqrt(temp[1].dot(temp[1].T)),1),(numpy.sqrt(temp[2].dot(temp[2].T)),2))
    lengths = sorted(lengths)

    tempIndex= [x[1] for x in lengths]
    '''transpose so we can switch the vertical vectors easier'''
    transposedPositions=symMatrix.T
    returnMatrix = copy.deepcopy(transposedPositions)

    returnMatrix[0] = numpy.copy(transposedPositions[tempIndex[0]])
    returnMatrix[1] = numpy.copy(transposedPositions[tempIndex[1]])
    returnMatrix[2] = numpy.copy(transposedPositions[tempIndex[2]])
        
    '''transpose the positions again so they're back to normal'''
    return returnMatrix.T


def _free2celldm(cellparamatrix,ibrav=0,primitive=True):
    '''
    Convert lattice vectors to celldm

    Arguments:
          cellparamatrix (numpy.matrix): matrix of cell vectors

    Keyword Arguments:
          ibrav (int): Overrides the bravais lattice automatically detected 
                       (must be in QE convention if primitive)
          primitive (bool): If True it will treat cellparamatrix as primitive lattice vectors

    Returns:
         paramDict (dict): dictionary of the celldm generated by the input matrix

    '''        

    splitString=  [x.split('=') for x in  AFLOWpi.retr.free2ibrav(cellparamatrix,ibrav=ibrav,primitive=primitive).split()]
    paramDict=OrderedDict()
    for item in splitString:
        paramDict[item[0]]=item[1].strip(' ,')
        
    newParamDict=OrderedDict()
    for k,v in paramDict.items():
        key=k.replace('(','').replace(')','')
        if k=='ibrav':
            newParamDict[key]=int(v)
        else:
            newParamDict[key]=numpy.around(numpy.float(v),decimals=5)

    return paramDict

def free2ibrav(cellparamatrix,ibrav=0,primitive=True):
    '''
    Convert lattice vectors to celldm

    Arguments:
          cellparamatrix (numpy.matrix): matrix of cell vectors

    Keyword Arguments:
          ibrav (int): Overrides the bravais lattice automatically detected 
                       (must be in QE convention if primitive)
          primitive (bool): If True it will treat cellparamatrix as primitive lattice vectors

    Returns:
         ibravStr (str): string of the celldm used for QE pwscf input

    '''        
    
    if primitive==True:
        cellParamMatrixConv = AFLOWpi.retr._prim2ConvVec(cellparamatrix,ibrav=ibrav)        

    else:
        cellParamMatrixConv = cellparamatrix



    if ibrav==0:
        ibrav = int(AFLOWpi.retr.getIbravFromVectors(cellparamatrix))

    conv=cellParamMatrixConv.T.getA()*0.529177249
    
    import numpy as np
    try:

        

        a = np.sqrt(cellParamMatrixConv[0].dot(cellParamMatrixConv[0].T))
        b = np.sqrt(cellParamMatrixConv[1].dot(cellParamMatrixConv[1].T))
        c = np.sqrt(cellParamMatrixConv[2].dot(cellParamMatrixConv[2].T))

    except:
        cellParamMatrixConv = np.array(cellParamMatrixConv)
        a = np.sqrt(cellParamMatrixConv[0].dot(cellParamMatrixConv[0]))
        b = np.sqrt(cellParamMatrixConv[1].dot(cellParamMatrixConv[1]))
        c = np.sqrt(cellParamMatrixConv[2].dot(cellParamMatrixConv[2]))

    celldm_1 = a
    celldm_2 = b/a
    celldm_3 = c/a

    celldm_6 = cellParamMatrixConv[1].dot(cellParamMatrixConv[2].T)/(b*c)
    celldm_5 = cellParamMatrixConv[0].dot(cellParamMatrixConv[2].T)/(a*c)
    celldm_4 = cellParamMatrixConv[0].dot(cellParamMatrixConv[1].T)/(a*b)


    ibravStr = ""
    celldm_1=numpy.around(celldm_1,decimals=5)
    celldm_2=numpy.around(celldm_2,decimals=5)
    celldm_3=numpy.around(celldm_3,decimals=5)
    celldm_4=numpy.around(celldm_4,decimals=5)
    celldm_5=numpy.around(celldm_5,decimals=5)
    celldm_6=numpy.around(celldm_6,decimals=5)

    

    try:
        if ibrav==1 or ibrav==2 or ibrav==3:
                ibravStr = "ibrav=%s, celldm(1)=%f\n"%(ibrav,celldm_1)

        elif ibrav==4 or ibrav==6 or ibrav==7:
            ibravStr = "ibrav=%s, celldm(1)=%f, celldm(3)=%f\n"%(ibrav,celldm_1, celldm_3) 
        elif ibrav==5:
                ibravStr = "ibrav=%s, celldm(1)=%f, celldm(4)=%f\n"%(ibrav,celldm_1, celldm_4) 


        elif ibrav==8 or ibrav==9 or ibrav==-9 or ibrav ==10 or ibrav==11:
            ibravStr = "ibrav=%s, celldm(1)=%f, celldm(2)=%f, celldm(3)=%f\n"%(ibrav,celldm_1, celldm_2, celldm_3)
        elif ibrav==12:
            ibravStr = "ibrav=%s, celldm(1)=%f, celldm(2)=%f, celldm(3)=%f, celldm(4)=%f\n"%(ibrav,celldm_1, celldm_2, celldm_3, celldm_4)
        elif ibrav==-12:
                ibravStr = "ibrav=%s, celldm(1)=%f, celldm(2)=%f, celldm(3)=%f, celldm(5)=%f\n"%(ibrav,celldm_1, celldm_2, celldm_3, celldm_5)
        elif ibrav==13:
                ibravStr = "ibrav=%s, celldm(1)=%f, celldm(2)=%f, celldm(3)=%f, celldm(4)=%f\n"%(ibrav,celldm_1, celldm_2, celldm_3, celldm_4)
        elif ibrav==14:
                ibravStr = "ibrav=%s, celldm(1)=%f, celldm(2)=%f, celldm(3)=%f, celldm(4)=%f celldm(5)=%f celldm(6)=%f\n"%(ibrav,celldm_1, celldm_2, celldm_3, celldm_4, celldm_5, celldm_6)
        else:
            logging.error('ibrav %s is not a valid option.' % ibrav)
            print('ibrav %s is not a valid option.' % ibrav)
            return 

    except Exception as e:

        AFLOWpi.run._fancy_error_log(e)

    return ibravStr




def free2abc(cellparamatrix,cosine=True,degrees=True,string=True,bohr=False):
    '''
    Convert lattice vectors to a,b,c,alpha,beta,gamma of the primitive lattice

    Arguments:
          cellparamatrix (numpy.matrix): matrix of cell vectors

    Keyword Arguments:
          cosine (bool): If True alpha,beta,gamma are cos(alpha),cos(beta),cos(gamma),
          degrees (bool): If True return alpha,beta,gamma in degrees; radians if False
          string (bool): If True return a,b,c,alpha,beta,gamma as a string; if False return as a list
          bohr (bool): If True return a,b,c in bohr radii; if False return in angstrom

    Returns:
         paramArray (list): a list of the parameters a,b,c,alpha,beta,gamma generated from the input matrix

    '''        

    ibrav = getIbravFromVectors(cellparamatrix)
    try:
        cellparamatrix=cellparamatrix.getA()
    except Exception as e:
        pass
#        print e
    try:
        a = numpy.sqrt(cellparamatrix[0].dot(cellparamatrix[0].T))
        b = numpy.sqrt(cellparamatrix[1].dot(cellparamatrix[1].T))
        c = numpy.sqrt(cellparamatrix[2].dot(cellparamatrix[2].T))
    except:
        cellparamatrix = numpy.array(cellparamatrix)
        a = numpy.sqrt(cellparamatrix[0].dot(cellparamatrix[0]))
        b = numpy.sqrt(cellparamatrix[1].dot(cellparamatrix[1]))
        c = numpy.sqrt(cellparamatrix[2].dot(cellparamatrix[2]))

    degree2radian = numpy.pi/180
    alpha,beta,gamma=(0.0,0.0,0.0)


    alpha = numpy.arccos(cellparamatrix[1].dot(cellparamatrix[2].T)/(b*c))
    beta  = numpy.arccos(cellparamatrix[0].dot(cellparamatrix[2].T)/(a*c))
    gamma = numpy.arccos(cellparamatrix[0].dot(cellparamatrix[1].T)/(a*b))

    if numpy.abs(alpha)<0.000001:
        alpha=0.0
    if numpy.abs(beta)<0.000001:
        beta=0.0
    if numpy.abs(gamma)<0.000001:
        gamma=0.0


    AngstromToBohr = 1.88971616463207
    BohrToAngstrom = 1/AngstromToBohr
    if bohr==False:
        a*=BohrToAngstrom
        b*=BohrToAngstrom
        c*=BohrToAngstrom

        a=float('%10.5e'%numpy.around(a,decimals=5))
        b=float('%10.5e'%numpy.around(b,decimals=5))
        c=float('%10.5e'%numpy.around(c,decimals=5))

    if cosine==True:
        cosBC=numpy.cos(alpha)
        cosAC=numpy.cos(beta)
        cosAB=numpy.cos(gamma)
        paramArray = [a,b,c,cosBC,cosAC,cosAB]

        param_list=[]
        for v in range(len(paramArray)):
            param_list.append(float('%10.5e'%numpy.around(paramArray[v],decimals=5)))
        paramArray=tuple(param_list)

        returnString = 'a=%s,b=%s,c=%s,cos(alpha)=%s,cos(beta)=%s,cos(gamma)=%s' % tuple(paramArray)

    if degrees==True:
        alpha/=degree2radian
        beta/= degree2radian
        gamma/=degree2radian

    if cosine!=True:
        paramArray = (a,b,c,alpha,beta,gamma)

        param_list=[]
        for v in range(len(paramArray)):
            param_list.append(float('%10.5e'%numpy.around(paramArray[v],decimals=5)))
        paramArray=tuple(param_list)
            
        returnString = 'A=%s,B=%s,C=%s,alpha=%s,beta=%s,gamma=%s' % tuple(paramArray)

    if string==True:
        return returnString
    else:

        return paramArray
            
 

def celldm2abc(ibrav=None,celldm1=None,celldm2=None,celldm3=None,celldm4=None,celldm5=None,celldm6=None,cosine=True,degrees=False):
    '''
    Convert celldm for espresso into A,B,C,and angles alpha,beta,gamma

    Arguments:
          cellparamatrix (numpy.matrix): matrix of cell vectors

    Keyword Arguments:
          cosine (bool): If True alpha,beta,gamma are cos(alpha),cos(beta),cos(gamma),
          degrees (bool): If True return alpha,beta,gamma in degrees; radians if False

    Returns:
         paramArray (list): a list of the parameters a,b,c,alpha,beta,gamma generated from the input matrix

    '''        

    if ibrav==1 or ibrav==2 or ibrav==3:
        celldm2=1.0
        celldm3=1.0
        celldm4=0.0
        celldm5=0.0
        celldm6=0.0


    elif ibrav==8 or ibrav==9 or ibrav==10 or ibrav==11:
        celldm4=0.0
        celldm5=0.0
        celldm6=0.0
    elif ibrav==12 or ibrav==13:
        celldm5=0.0
        celldm6=0.0


        
    elif ibrav==6 or ibrav==7:
        celldm2=1.0
        celldm4= 0.0
        celldm5= 0.0
        celldm6= 0.0

    elif ibrav==4:
        celldm2= 1.0
        celldm4= 0.0
        celldm5= 0.0
        celldm6=-0.5

    elif ibrav==5:
        celldm2=1.0
        celldm3=1.0
        celldm5=celldm4
        celldm6=celldm4
        
    else:

        print('ibrav=%s not supported' % ibrav)
        logging.warning('ibrav=%s not supported' % ibrav)
        return

    a=celldm1
    b=celldm2*celldm1
    c=celldm3*celldm1
    if cosine==False:
        alpha= numpy.arccos(celldm4)
        beta=  numpy.arccos(celldm5)
        gamma= numpy.arccos(celldm6)
        if numpy.abs(alpha)<0.00000001:
            alpha=0.0
        if numpy.abs(beta)<0.00000001:
            beta=0.0
        if numpy.abs(gamma)<0.00000001:
            gamma=0.0
    else:
        alpha= celldm4
        beta = celldm5
        gamma= celldm6
    degree2radian = numpy.pi/180
    if degrees==True:
        alpha/= degree2radian
        beta /= degree2radian
        gamma/= degree2radian


    a=numpy.around(a,decimals=5)
    b=numpy.around(b,decimals=5)
    c=numpy.around(c,decimals=5)
    alpha=numpy.around(alpha,decimals=5)
    beta=numpy.around(beta,decimals=5)
    gamma=numpy.around(gamma,decimals=5)

    return a,b,c,alpha,beta,gamma

def abc2celldm(a=None,b=None,c=None,alpha=None,beta=None,gamma=None,ibrav=None):
    '''
    Convert a,b,c,alpha,beta,gamma into celldm for QE

    Arguments:
          None

    Keyword Arguments:
          a (float): length of a
          b (float): length of b
          c (float): length of c
          alpha (float): Angle between axis b and c
          beta  (float): Angle between axis a and c
          gamma (float): Angle between axis a and b
          ibrav (int): ibrav to be used to convert for defaults

    Returns:
          celldm_array (array): an array of (ibrav,celldm(1),celldm(2),celldm(3),celldm(4),celldm(5),celldm(6))

    '''        


    if ibrav==1 or ibrav==2 or ibrav==3 or ibrav==6 or ibrav==7 or ibrav == 8 or ibrav==9 or ibrav==-9 or ibrav==10 or ibrav==11:
        alpha,beta,gamma = (90.0,90.0,90.0)
    if ibrav == 1 or ibrav==2 or ibrav==3:
        b,c = (a,a)
    if ibrav == 4:
        b = a
        alpha = 90.0
        beta  = 90.0
        gamma = 120.0
    if ibrav == 5 or ibrav == -5:
        b = a
        c = b
#        beta = alpha
#        gamma = alpha
    if ibrav == 6 or ibrav==7:
        b = a
    if ibrav==12 or ibrav==13:
        alpha=90.0
        beta=90.0

    celldm1 = a
    celldm2 = b/a
    celldm3 = c/a
    celldm4 = numpy.cos(gamma*numpy.pi/180.0)
    celldm5 = numpy.cos(beta*numpy.pi/180.0)
    celldm6 = numpy.cos(alpha*numpy.pi/180.0)


    celldm1=numpy.around(celldm1,decimals=5)
    celldm2=numpy.around(celldm2,decimals=5)
    celldm3=numpy.around(celldm3,decimals=5)
    celldm4=numpy.around(celldm4,decimals=5)
    celldm5=numpy.around(celldm5,decimals=5)
    celldm6=numpy.around(celldm6,decimals=5)


    if numpy.abs(celldm4)<0.000000001:
       celldm4=0.0
    if numpy.abs(celldm5)<0.000000001:
       celldm5=0.0
    if numpy.abs(celldm6)<0.000000001:
       celldm6=0.0


    return ibrav,celldm1,celldm2,celldm3,celldm4,celldm5,celldm6

def abcVol(a=None,b=None,c=None,alpha=None,beta=None,gamma=None,ibrav=None):
    '''
    Get volume from a,b,c,alpha,beta,gamma

    Arguments:
          None

    Keyword Arguments:
          a (float): length of a
          b (float): length of b
          c (float): length of c
          alpha (float): Angle between axis b and c
          beta  (float): Angle between axis a and c
          gamma (float): Angle between axis a and b
          ibrav (int): ibrav to be used to convert for defaults

    Returns:
          cell_vol (float): volume of the cell

    '''        

    cellParamMatrix = abc2free(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,ibrav=ibrav,returnString=False)

    return AFLOWpi.retr.getCellVolumeFromVectors(cellParamMatrix)

def abc2free(a=None,b=None,c=None,alpha=None,beta=None,gamma=None,ibrav=None,returnString=True):
    '''
    Converts a,b,c,alpha,beta,gamma into QE convention primitive lattice vectors 

    Arguments:
          None

    Keyword Arguments:
          a (float): length of a
          b (float): length of b
          c (float): length of c
          alpha (float): Angle between axis b and c
          beta  (float): Angle between axis a and c
          gamma (float): Angle between axis a and b
          ibrav (int): bravais lattice type by QE convention
          returnString (bool): return vectors as a string or as a numpy matrix

    Returns:
          Either a numpy.matrix or a string of the primitive lattice vectors generated from the celldm input

    '''

    ibrav,celldm1,celldm2,celldm3,celldm4,celldm5,celldm6 = abc2celldm(a,b,c,alpha,beta,gamma,ibrav=ibrav)
    cell_vectors = AFLOWpi.retr.celldm2free(ibrav,celldm1,celldm2,celldm3,celldm4,celldm5,celldm6,returnString=returnString)
    return cell_vectors

# def ibrav2String(ibrav=None,celldm1=None,celldm2=None,celldm3=None,celldm4=None,celldm5=None,celldm6=None):
#     '''
#     DEFUNCT <CONSIDER FOR REMOVAL>

#     Arguments:


#     Keyword Arguments:


#     Returns:


#     '''        

#     ibravStr = ""

#     try:
#         if ibrav == 1 or ibrav == 2 or ibrav ==3:
#                 ibravStr = "ibrav=1, celldm(1)=%f\n"%celldm1
#         elif ibrav == 5:
#                 ibravStr = "ibrav=5, celldm(1)=%f, celldm(4)=%f\n"%(celldm1, celldm2) 
#         elif ibrav == 4 or ibrav == 6 or ibrav ==7:
#                 ibravStr = "ibrav=2, celldm(1)=%f, celldm(3)=%f\n"%(celldm1, celldm3)
#         elif ibrav == 8 or ibrav == 9 or ibrav == 10 or ibrav == 11:
#                 ibravStr = "ibrav=2, celldm(1)=%f, celldm(2)=%f, celldm(3)=%f\n"%(celldm1, celldm2, celldm3)
#         elif ibrav == 12 or ibrav == 13:
#                 ibravStr = "ibrav=2, celldm(1)=%f, celldm(2)=%f, celldm(3)=%f, celldm(4)=%f\n"%(celldm1, celldm2, celldm3, celldm4)
#         elif ibrav == 14:
#                 ibravStr = "ibrav=2, celldm(1)=%f, celldm(2)=%f, celldm(3)=%f, celldm(4)=%f celldm(5)=%f celldm(6)=%f\n"%(celldm1, celldm2, celldm3, celldm4, celldm5, celldm6)
#         else:
#             logging.error('ibrav %s is not a valid option.' % ibrav)
#             print 'ibrav %s is not a valid option.' % ibrav
#             return 

#     except Exception as e:
#         AFLOWpi.run._fancy_error_log(e)

#     return ibravStr



def getPositionsFromOutput(oneCalc,ID):
    '''
    Get the atomic positions from the output. If atomic positons can not be read from
    output then the the positions from the input are returned.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          None

    Returns:
          pos (numpy.matrix): atomic positions

    '''        

    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.out'),'r') as outFileObj:
            outFileString  = outFileObj.read()

        pos= AFLOWpi.retr._getPositions(outFileString)

        return pos

        if len(pos)==0:
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'r') as outFileObj:
                outFileString  = outFileObj.read()
            pos = AFLOWpi.retr._getPositions(outFileString)

            return pos 
    except:
        try:
            with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'r') as outFileObj:
                outFileString  = outFileObj.read()

            pos = AFLOWpi.retr._getPositions(outFileString)

            return pos 
        except Exception as e:
            AFLOWpi.run._fancy_error_log(e)



def getCellMatrixFromInput(inputString,string=False,scale=True):
    '''
    Get the primitive cell vectors from the input.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          string (bool): Return a string or matrix
          scale (bool): scale the vectors by alat

    Returns:
          cellParamMatrix (numpy.matrix): primitive lattice params

    '''        

    splitInput = AFLOWpi.retr._splitInput(inputString)
    if 'CELL_PARAMETERS' in list(splitInput.keys()):
        '''sanitize the modifier'''
        modifier=splitInput['CELL_PARAMETERS']['__modifier__'].strip('(){}')
        cellTimes=1.0
        if scale==True:
            if modifier=='angstrom':
                cellTimes=1.0/0.529177249
        cellParamString = splitInput['CELL_PARAMETERS']['__content__']
        cellParamMatrix = AFLOWpi.retr._cellStringToMatrix(cellParamString)*cellTimes

        return cellParamMatrix
    celldm2freeDict={'ibrav':splitInput['&system']['ibrav'],'returnString':string}

    for items in list(splitInput['&system'].keys()):

        if len(re.findall('celldm',items)):
            inputVar = items.replace(')','')
            inputVar = inputVar.replace('(','')

            celldm2freeDict[inputVar]=float(splitInput['&system'][items])

    pos=AFLOWpi.retr.celldm2free(**celldm2freeDict)

    return pos


def _getCelldm2freeDict(free2celldm_output_dict):
    '''
    Cleans fortran array format into a format for python ex. celldm(1) to celldm1

    Arguments:
          free2celldm_output_dict (dict): output from AFLOWpi.retr.free2celldm

    Keyword Arguments:
          None

    Returns:
          output_dict (dict) dictionary with cleaned keywords

    '''        

    outputDict=OrderedDict()
    for k,v in free2celldm_output_dict.items():
        if len(re.findall('celldm',k)):
            inputVar = k.replace(')','')
            inputVar = inputVar.replace('(','')
            outputDict[inputVar]=float(v)
        if k=='ibrav':
            outputDict[k]=int(v)

    return outputDict

def getPositionsFromInput(oneCalc,ID):
    '''


    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:


    Returns:


    '''        

    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'r') as inFileObj:
            inFileString  = inFileObj.read()
    except:
        return
    return AFLOWpi.retr._getPositions(inFileString)


def _getPositions(inputString,matrix=True):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    try:
        atomicPos = AFLOWpi.qe.regex.atomic_positions(inputString)
    except:
        return
    atomicPos,flags = detachPosFlags(atomicPos)

    splitAtomicPos = atomicPos.split('\n')
    posArray=[]
    for atom in splitAtomicPos:
        splitAtom = atom.split()
        if len(splitAtom):
            posArray.append([float(x) for x in splitAtom[1:]])

    if matrix==False:
        arrayStr=''
        for entry in posArray:
            addition = ' '.join([str(x) for x in entry])+'\n'                

            arrayStr+=addition

        return arrayStr

    else:
        outputMatrix = numpy.matrix(posArray)
        return outputMatrix


def _getPosLabels(inputString):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    try:
        atomicPos     = AFLOWpi.qe.regex.atomic_positions(inputString)
    except:
        return

    splitAtomicPos = atomicPos.split('\n')
    labels=[]
    for atom in splitAtomicPos:
        splitAtom = atom.split()
        if len(splitAtom):
            labels.append(splitAtom[0])

    return labels


#################################################################################################################
#################################################################################################################
## prim to conv cell transform
#################################################################################################################
#################################################################################################################

def _rho2hex(cellParamMatrix,symMatrix,labels):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    conv2prim=numpy.matrix([[ 2.0  , 1.0 , 1.0 ],
                            [-1.0  , 1.0 , 1.0 ],
                            [-1.0  ,-2.0 , 1.0 ]])/3.0


    prim2conv=numpy.linalg.inv(conv2prim)
    hex_vec=prim2conv.dot(cellParamMatrix)
    first=numpy.matrix([[ 2.0  , 0.0 , 0.0 ],
                        [ 0.0  , 1.0 , 0.0 ],
                        [ 0.0  , 0.0 , 1.0 ]])/3.0
    second=numpy.matrix([[ 1.0  ,0.0 , 0.0 ],
                        [ 0.0  , 2.0 , 0.0 ],
                        [ 0.0  , 0.0 , 2.0 ]])/3.0

    first_copy=copy.deepcopy(symMatrix)
    second_copy=copy.deepcopy(symMatrix)
    for i in range(len(symMatrix)):
        first_copy[i]=first.dot(symMatrix[i].T).T
    for i in range(len(symMatrix)):
        second_copy[i]=first.dot(symMatrix[i].T).T


    sym_list=symMatrix.astype(numpy.float16).tolist()


    first_list = first_copy.astype(numpy.float16).tolist()

    second_list= second_copy.astype(numpy.float16).tolist()
    sym_list.extend(first_list)
    sym_list.extend(second_list)
    symMatrix=numpy.matrix(sym_list)

    for i in range(len(symMatrix)):
        symMatrix[i]=prim2conv.dot(symMatrix[i].T).T


    return labels,symMatrix,hex_vec

def _prim2ConvMatrix(cellParamMatrix,ibrav=0):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    try:
        if ibrav==0:
            ibrav = int(getIbravFromVectors(cellParamMatrix))

    except:
        ibrav=0
        cellParamMatrixCopy=copy.deepcopy(cellParamMatrix)
        return numpy.matrix([[ 1.0  , 0.0 , 0.0 ],
                             [ 0.0  , 1.0 , 0.0 ],
                             [ 0.0  , 0.0 , 1.0 ]])

#    '''conventional cells are primative cells so don't do anything'''
    if ibrav==1  or  ibrav==4  or  ibrav==6 or ibrav==8 or ibrav==12 or ibrav==14 or ibrav==None or ibrav==0:
        conv2prim=numpy.matrix([[ 1.0  , 0.0 , 0.0 ],
                                [ 0.0  , 1.0 , 0.0 ],
                                [ 0.0  , 0.0 , 1.0 ]])

    if ibrav==5:
        conv2prim=numpy.matrix([[ 1.0  , 0.0 , 0.0 ],
                                [ 0.0  , 1.0 , 0.0 ],
                                [ 0.0  , 0.0 , 1.0 ]])        


#    '''Face Centered Cubic and Orthorhombic'''
    if ibrav==10:
        conv2prim=numpy.matrix([
                [ 1.0  , 0.0 , 1.0 ],
                [ 1.0  , 1.0 , 0.0 ],
                [ 0.0  , 1.0 , 1.0 ],
                ])/2.0

    if ibrav==2:
        conv2prim=numpy.matrix([
                [-1.0  , 0.0 , 1.0 ],
                [ 0.0  , 1.0 , 1.0 ],
                [-1.0  , 1.0 , 0.0 ],
                ])/2.0


    if ibrav in [7]:
        conv2prim=numpy.matrix([
                [ 1.0  ,-1.0 , 1.0 ],
                [ 1.0  , 1.0 , 1.0 ],
                [-1.0  ,-1.0 , 1.0 ], 
                ])/2.0

    if ibrav in [3,11]:
        conv2prim=numpy.matrix([
                [ 1.0  , 1.0 , 1.0 ],
                [-1.0  , 1.0 , 1.0 ],
                [-1.0  ,-1.0 , 1.0 ], 
                ])/2.0
        

#    '''Base Centered Orthorhombic and Monoclinic C centered'''
    if ibrav==9:
        conv2prim=numpy.matrix([[ 0.5  , 0.5 , 0.0 ],
                                [-0.5  , 0.5 , 0.0 ],
                                [ 0.0  , 0.0 , 1.0 ]])

    if ibrav==13:
        conv2prim=numpy.matrix([[ 0.5  , 0.0 ,-0.5 ],
                                [ 0.0  , 1.0 , 0.0 ],
                                [ 0.5  , 0.0 , 0.5 ]])

    
    prim2conv=numpy.linalg.inv(conv2prim)

    return prim2conv




def _prim2ConvVec(cellParamMatrix,ibrav=0):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    prim2conv = AFLOWpi.retr._prim2ConvMatrix(cellParamMatrix,ibrav=ibrav)
    cellParamMatrixCopy=copy.deepcopy(cellParamMatrix)
#    conv=prim2conv.dot(cellParamMatrixCopy).astype(numpy.float32)
    return prim2conv

def _conv2PrimVec(cellParamMatrix,ibrav=0):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    prim2conv = numpy.linalg.inv(_prim2ConvMatrix(cellParamMatrix,ibrav=ibrav))
    cellParamMatrixCopy=copy.deepcopy(cellParamMatrix)
#    prim=prim2conv.dot(cellParamMatrixCopy).astype(numpy.float32)
    return prim2conv

def _getConventionalCellFromInput(oneCalc,ID):
    '''


    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:


    Returns:


    '''        

    inputString = oneCalc['_AFLOWPI_INPUT_']
    cellParamMatrix = AFLOWpi.retr.getCellMatrixFromInput(inputString)

    convCell= AFLOWpi.retr._prim2ConvVec(cellParamMatrix).getA()
    '''for aflow integration'''
        
    symMatrix = AFLOWpi.retr.getPositionsFromInput(oneCalc,ID)
    labels=AFLOWpi.retr._getPosLabels(inputString)
    labels,transformedPos = AFLOWpi.retr._prim2convPositions(labels,symMatrix,cellParamMatrix)
    convCell = numpy.matrix(convCell)

    return labels,transformedPos,convCell


def _getConventionalCellFromOutput(oneCalc,ID):
    '''


    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:


    Returns:


    '''        

    outputString = AFLOWpi.retr._getOutputString(oneCalc,ID)
    alat,cellParamMatrix = AFLOWpi.retr._getCellParams(oneCalc,ID)

    convCell= AFLOWpi.retr._prim2ConvVec(cellParamMatrix*alat).getA()

    '''for aflow integration'''


    symMatrix = AFLOWpi.retr.getPositionsFromOutput(oneCalc,ID)
    inputString=oneCalc['_AFLOWPI_INPUT_']
    labels=AFLOWpi.retr._getPosLabels(inputString)
    labels,transformedPos = AFLOWpi.retr._prim2convPositions(labels,symMatrix,cellParamMatrix)



    convCell = numpy.matrix(convCell)

    return labels,transformedPos,convCell

def transform_input_conv(oneCalc,ID):
    '''


    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:


    Returns:


    '''        

    labels,transformedPos,convCell = AFLOWpi.retr._getConventionalCellFromInput(oneCalc,ID)
    
    atomic_pos = AFLOWpi.retr._joinMatrixLabels(labels,transformedPos)
    celldmDict  = AFLOWpi.retr._free2celldm(convCell)
    
    splitInput = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

    for k,v in splitInput['&system'].items():
        try:
            splitInput['&system'][k]=celldmDict[k]
        except Exception as e:
            pass
    splitInput['ATOMIC_POSITIONS']['__content__']=atomic_pos
    splitInput['&system']['nat']=len(labels)
    newInput  = AFLOWpi.retr._joinInput(splitInput)

    return newInput



def _prim2convPositions(labels,symMatrix,cellParamMatrix):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    try:
        symMatrix=symMatrix.getA()
    except:
        pass

    symMatrix = AFLOWpi.retr._shiftAfter(symMatrix)
    returnMatrix=copy.deepcopy(symMatrix)
    ibrav = getIbravFromVectors(cellParamMatrix)

    #'''primative=conv'''





    if ibrav==1 or ibrav==6 or ibrav==8 or ibrav==12 or ibrav==14 or ibrav==5:
        return labels,symMatrix
    

    #'''base centered'''
    if ibrav==7:
        pass


        
    if ibrav==9 or ibrav==13:  

        pass

    #'''face centered'''        

    #'''body centered'''    

    #'''hexagonal'''    


    prim2Conv = AFLOWpi.retr._prim2ConvMatrix(cellParamMatrix)




    if ibrav==3 or ibrav==9:

        largerReturn=numpy.zeros(shape=(symMatrix.shape[0]*2,3))
        #double length of labels
        labels.extend(labels)

        lat1=copy.deepcopy(returnMatrix)
        lat2=copy.deepcopy(returnMatrix)


        returnMatrix = largerReturn


        returnMatrix = AFLOWpi.retr._shiftAfter(returnMatrix)

        for entry in range(len(symMatrix)):
            pos = symMatrix[entry]
            returnMatrix[entry] = prim2Conv.dot(pos.T).T
            returnMatrix[entry] = pos

#        for entry in range(len(lat1)):
#            largerReturn[entry]=AFLOWpi.retr._shiftAfter(lat1[entry])

        for entry in range(len(lat1)):
            lat1[entry] = AFLOWpi.retr._shiftX(returnMatrix[entry],0.5)
            lat1[entry] = AFLOWpi.retr._shiftY(lat1[entry],0.5)
            lat1[entry] = AFLOWpi.retr._shiftZ(lat1[entry],0.5)

#            largerReturn[entry+1*len(lat1)]= AFLOWpi.retr._shiftAfter(lat1[entry])
#            returnMatrix[entry+1*len(lat1)]= lat1[entry]







    if ibrav==2  or ibrav==10:
        largerReturn=numpy.zeros(shape=(symMatrix.shape[0]*4,3))
        #double the length of the labels twice
        labels.extend(labels)
        labels.extend(labels)

        lat1=copy.deepcopy(returnMatrix)
        lat2=copy.deepcopy(returnMatrix)
        lat3=copy.deepcopy(returnMatrix)
        lat4=copy.deepcopy(returnMatrix)

        for entry in range(len(lat1)):
            largerReturn[entry]=AFLOWpi.retr._shiftAfter(lat1[entry])
            
        for entry in range(len(lat2)):
            lat2[entry] = AFLOWpi.retr._shiftY(lat2[entry],0.5)
            lat2[entry] = AFLOWpi.retr._shiftZ(lat2[entry],0.5)
            largerReturn[entry+1*len(lat1)]= AFLOWpi.retr._shiftAfter(lat2[entry])

        for entry in range(len(lat3)):
            lat3[entry] = AFLOWpi.retr._shiftX(lat3[entry],0.5)
            lat3[entry] = AFLOWpi.retr._shiftZ(lat3[entry],0.5)
            largerReturn[entry+2*len(lat1)]= AFLOWpi.retr._shiftAfter(lat3[entry])

        for entry in range(len(lat4)):
            lat4[entry] = AFLOWpi.retr._shiftY(lat4[entry],0.5)
            lat4[entry] = AFLOWpi.retr._shiftZ(lat4[entry],0.5)
            largerReturn[entry+3*len(lat1)]= AFLOWpi.retr._shiftAfter(lat4[entry])

        
        largerReturn = AFLOWpi.retr._shiftAfter(numpy.matrix(largerReturn))
        returnMatrix = largerReturn
        

        for entry in range(len(largerReturn)):
            pos = largerReturn[entry]
            returnMatrix[entry] = prim2Conv.dot(pos.T).T


    returnMatrix=numpy.matrix(returnMatrix)

    returnMatrixList=returnMatrix
    returnMatrix = numpy.matrix(returnMatrix)


    return labels,returnMatrix


def _conv2primPositions(symMatrix,cellParamMatrix):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''        

    returnMatrix=copy.deepcopy(symMatrix)
    prim2Conv = AFLOWpi.retr._prim2ConvMatrix(cellParamMatrix)
    symMatrix=symMatrix.getA()

    for entry in range(len(symMatrix)):
        returnMatrix[entry] =  prim2Conv.dot(symMatrix[entry])


    return returnMatrix


def convertFCC(oneCalc,ID):
    '''


    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:


    Returns:


    '''        

    try:
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'r') as inFileObj:
            inputString  = inFileObj.read()
    except:
        return

    labels = AFLOWpi.retr._getPosLabels(inputString)
    cellVectors =  AFLOWpi.retr._getCellInput(oneCalc,ID,scaled=False)
    positions = getPositionsFromInput(oneCalc,ID)
    posNew=[]

    for atom in range(len(positions)):

        posTemp = AFLOWpi.retr._duplicateEdgeAtoms(positions[atom])
        posTemp =  AFLOWpi.retr._reduceDuplicates(posTemp)
        posNew.extend([labels[atom],x] for x in posTemp.tolist())

    labels = [x[0] for x in posNew]

    conv2prim=numpy.matrix([[-1.0  , 0.0 , 1.0 ],
                            [ 0.0  , 1.0 , 1.0 ],
                            [-1.0  , 1.0 , 0.0 ]])/2

    prim2conv =     numpy.linalg.inv(conv2prim)

    a = numpy.sqrt(cellVectors[0].dot(cellVectors[0].T)).getA()[0][0]
    b = numpy.sqrt(cellVectors[1].dot(cellVectors[1].T)).getA()[0][0]
    c = numpy.sqrt(cellVectors[2].dot(cellVectors[2].T)).getA()[0][0]

    conv =  prim2conv.dot(cellVectors)
    a =  a.getA()[0][0]
    b =  b.getA()[0][0]
    c =  c.getA()[0][0]
    mat = numpy.matrix([[a,0,0],
                        [0,b,0],
                        [0,0,c]],)
    for entry in range(len(positions)):
        positions[entry] = prim2conv.dot(positions[entry].T).T




##############################################################################################################
##
## NOT DONE
##
##############################################################################################################

# def get_aflow_transform(ibrav):

#     #SC ST SO RHO TRI SMCL                                                                                                                                                                                                                                                        
#     if ibrav in [1,4,5,6,8,12,14]:
#         cellParamMatrix = numpy.matrix([[ 1.0  , 0.0 , 0.0 ],
#                                         [ 0.0  , 1.0 , 0.0 ],
#                                         [ 0.0  , 0.0 , 1.0 ]])
#     #FC                                                                                                                                                                                                                                                                           
#     if ibrav in [2,10]:
#         cellParamMatrix = numpy.matrix([[ 0.0  , 1.0 , 1.0, ],
#                                         [ 1.0  , 0.0 , 1.0, ],
#                                         [ 1.0  , 1.0 , 0.0, ],])/2.0
#     #BC                                                                                                                                                                                                                                                                           
#     if ibrav in [3,7,11]:
#         cellParamMatrix = numpy.matrix([[-1.0  , 1.0 , 1.0 ],
#                                         [ 1.0  ,-1.0 , 1.0 ],
#                                         [ 1.0  , 1.0 ,-1.0 ]])/2.0
#     #ORCC                                                                                                                                                                                                                                                                         
#     if ibrav in [9]:
#         cellParamMatrix = numpy.matrix([[ 1.0  ,-1.0 , 0.0 ],
#                                         [ 1.0  , 1.0 , 0.0 ],
#                                         [ 0.0  , 0.0 , 2.0 ]])/2.0
#     #MCLC                                                                                                             

#     if ibrav in [13]:
#         cellParamMatrix = numpy.matrix([[ 1.0  , 1.0 , 0.0 ],
#                                         [-1.0  , 1.0 , 0.0 ],
#                                         [ 0.0  , 0.0 , 2.0 ]])/2.0
#     return cellParamMatrix

# def get_conv_aflow(cellParamMatrix,ibrav=0):

#     transform_matrix = get_aflow_transform(ibrav)
#     transform = numpy.linalg.inv(transform_matrix)
#     return transform.dot(cellParamMatrix)


# def QE_AFLOW(inputString):

#     cellParamMatrix = AFLOWpi.retr.getCellMatrixFromInput(inputString,scale=False)
#     to_prim=numpy.linalg.inv(AFLOWpi.retr._prim2ConvMatrix(cellParamMatrix,ibrav+1))

#     symMatrix = AFLOWpi.retr._getPositions(inputString).astype(numpy.float16)

#     input_dict = AFLOWpi.retr._splitInput(inputString)
#     cellParamMatrix = to_prim.dot(get_conv_aflow(cellParamMatrix,ibrav=ibrav+1))
#     cellParamMatrix = AFLOWpi.retr._cellMatrixToString(cellParamMatrix)

#     input_dict['CELL_PARAMETERS']['__content__']=cellParamMatrix

#     pos_conv_aflow = get_conv_aflow(symMatrix.T,ibrav=ibrav+1)
#     symMatrix = to_prim.dot(pos_conv_aflow).T
#     symMatrix = AFLOWpi.retr._shiftAfter(symMatrix)
#     labels = AFLOWpi.retr._getPosLabels(inputString)
#     symMatrix = AFLOWpi.retr._joinMatrixLabels(labels,symMatrix)

#     input_dict['ATOMIC_POSITIONS']['__content__']=symMatrix
#     input_dict['ATOMIC_SPECIES']['__content__']='\n'.join(list(set(labels)))+'\n'
#     inputString = AFLOWpi.retr._joinInput(input_dict)

#     return inputString

##############################################################################################################
##
## NOT DONE
##
##############################################################################################################

def getIbravFromVectors(cellVectors):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    a = numpy.sqrt(cellVectors[0].dot(cellVectors[0].T))
    b = numpy.sqrt(cellVectors[1].dot(cellVectors[1].T))
    c = numpy.sqrt(cellVectors[2].dot(cellVectors[2].T))


    alpha = cellVectors[1].dot(cellVectors[2].T)/(b*c)
    beta  = cellVectors[0].dot(cellVectors[2].T)/(a*c)
    gamma = cellVectors[0].dot(cellVectors[1].T)/(a*b)
    
    try:
        alphaDegrees = numpy.arccos(alpha).getA()[0][0]*180/numpy.pi
    except:
        alphaDegrees = numpy.arccos(alpha)*180/numpy.pi
    try:
        betaDegrees  = numpy.arccos(beta).getA()[0][0]*180/numpy.pi
    except:
        betaDegrees  = numpy.arccos(beta)*180/numpy.pi
    try:
        gammaDegrees = numpy.arccos(gamma).getA()[0][0]*180/numpy.pi
    except:
        gammaDegrees = numpy.arccos(gamma)*180/numpy.pi



    if numpy.abs(alphaDegrees)<0.0001:
        alphaDegrees=0.0
    if numpy.abs(betaDegrees)<0.0001:
        betaDegrees=0.0
    if numpy.abs(gammaDegrees)<0.0001:
        gammaDegrees=0.0

    

    # print a
    # print b
    # print c
    # print alphaDegrees
    # print betaDegrees
    # print gammaDegrees

    if numpy.abs(a-c)<0.01 and numpy.abs(a-b)<0.01 and numpy.abs(c-b)<0.01 and numpy.abs(alphaDegrees-betaDegrees)<0.01 and numpy.abs(gammaDegrees-betaDegrees)<0.01 and numpy.abs(gammaDegrees-alphaDegrees)<0.01:

        sortedVec =  AFLOWpi.retr._sortMatrix(cellVectors)

        reversedMat = numpy.matrix([x for x in reversed(sortedVec.getA())]).getA()
        '''in the case of FCC'''
        if  numpy.abs(alphaDegrees-60.0)<0.01 and numpy.abs(gammaDegrees-60.0)<0.01 and numpy.abs(betaDegrees-60.0)<0.01:
            return 2
            
            '''in the case of simple cubic'''
        elif reversedMat[0][0]==a and reversedMat[1][1]==b and reversedMat[2][2]==c:
            return 1
            '''in the case of BCC'''


        elif  numpy.abs(a-b)<0.01  and numpy.abs(b-c)<0.01 and numpy.abs(c-a)<0.01 and numpy.abs(alphaDegrees-betaDegrees)<0.01 and numpy.abs(alphaDegrees-gammaDegrees)<0.01 and numpy.abs(gammaDegrees-betaDegrees)<0.01 and not numpy.abs(gammaDegrees-90.0)<0.01:

            if numpy.abs(gammaDegrees-109.471225753)<0.01:
                return 3
            else:
                return 5            
            
        '''in the case of hexagonal'''
    elif numpy.abs(gammaDegrees-120.0)<0.01:
        return 4
        '''in the case of rhombohedral'''


    elif ((numpy.abs(a-b)<0.01  and numpy.abs(b-c)<0.01 and numpy.abs(c-a)<0.01) and (numpy.abs(alphaDegrees-109.471220)<0.01 and numpy.abs(betaDegrees-gammaDegrees)<0.01) or (numpy.abs(gammaDegrees-109.471220)<0.01 and numpy.abs(betaDegrees-alphaDegrees)<0.01) or (numpy.abs(betaDegrees-109.471220)<0.01 and numpy.abs(gammaDegrees-alphaDegrees)<0.01)):
        return 3

    
    elif (numpy.abs(a-b)<0.01 and not numpy.abs(c-b)<0.01) or (numpy.abs(c-b)<0.01 and not numpy.abs(a-b)<0.01) or (numpy.abs(c-a)<0.01 and not numpy.abs(c-b)<0.01):
        '''in the case of simple Tet'''
        if numpy.abs(alphaDegrees-90.0)<0.01 and numpy.abs(betaDegrees-90.0)<0.01 and numpy.abs(gammaDegrees-90.0)<0.01:
            return 6
            '''in the case of ORCC base cenetered ortho'''
        elif numpy.abs(alphaDegrees-90)<0.01 and numpy.abs(betaDegrees-90)<0.01 and not numpy.abs(gammaDegrees-90)<0.01:
            return 9
        elif  numpy.abs(gammaDegrees-90)<0.01 and numpy.abs(alphaDegrees-90)<0.01 and not numpy.abs(betaDegrees-90)<0.01:
            return 9
        elif  numpy.abs(betaDegrees-90)<0.01 and numpy.abs(gammaDegrees-90)<0.01 and not numpy.abs(alphaDegrees-90)<0.01:
            return 9


            'in the case of MCLC base centered monoclinic'
        elif numpy.abs(alphaDegrees-betaDegrees)<0.01 and numpy.abs(gammaDegrees-betaDegrees)>0.01 and not numpy.abs(gammaDegrees-90)<0.01:
            return 13
        elif  numpy.abs(gammaDegrees-alphaDegrees)<0.01 and numpy.abs(betaDegrees-alphaDegrees)>0.01 and not numpy.abs(betaDegrees-90)<0.01:
            return 13
        elif  numpy.abs(gammaDegrees-betaDegrees)<0.01 and numpy.abs(betaDegrees-alphaDegrees)>0.01 and not numpy.abs(alphaDegrees-90)<0.01:
            return 13



        '''in the case of BCO'''
    elif numpy.abs(a-b)<0.01  and numpy.abs(b-c)<0.01 and numpy.abs(c-a)<0.01:
        if (numpy.abs(alphaDegrees-gammaDegrees)<0.01 and not  numpy.abs(alphaDegrees-gammaDegrees)<0.01) or (numpy.abs(alphaDegrees-betaDegrees)<0.01 and not  numpy.abs(betaDegrees-gammaDegrees)<0.01) or (numpy.abs(betaDegrees-gammaDegrees)<0.01 and not  numpy.abs(alphaDegrees-gammaDegrees)<0.01):
            return 7

            '''in the case of ORCI'''
        elif numpy.abs(alphaDegrees-betaDegrees)>0.01 and numpy.abs(gammaDegrees-betaDegrees)>0.01 and numpy.abs(gammaDegrees-alphaDegrees)>0.01:
            return 11


        '''in the case of simple ortho'''
    elif  numpy.abs(a-b)>0.01  and numpy.abs(b-c)>0.01 and numpy.abs(c-a)>0.01:

        sortedVec =  AFLOWpi.retr._sortMatrix(cellVectors)
        reversedMat = numpy.matrix([x for x in reversed(sortedVec.getA())]).getA()

        if  numpy.abs(alphaDegrees-gammaDegrees)<0.01 and numpy.abs(gammaDegrees-betaDegrees)<0.01 and numpy.abs(gammaDegrees-90.0)<0.01:
            return 8



            'in the case of MCL simple monoclinic'
        elif numpy.abs(alphaDegrees-90)<0.01 and numpy.abs(betaDegrees-90)<0.01 and not numpy.abs(gammaDegrees-90)<0.01:
            return 12
        elif  numpy.abs(gammaDegrees-90)<0.01 and numpy.abs(alphaDegrees-90)<0.01 and not numpy.abs(betaDegrees-90)<0.01:
            return 12
        elif  numpy.abs(betaDegrees-90)<0.01 and numpy.abs(gammaDegrees-90)<0.01 and not numpy.abs(alphaDegrees-90)<0.01:
            return 12



        elif  numpy.abs(gammaDegrees-alphaDegrees)>0.01 and numpy.abs(betaDegrees-alphaDegrees)>0.01:
            ORCFtest = AFLOWpi.retr._convertCellBC(cellVectors)

            aTest = numpy.sqrt(ORCFtest[0].dot(ORCFtest[0].T))
            bTest = numpy.sqrt(ORCFtest[1].dot(ORCFtest[1].T))
            cTest = numpy.sqrt(ORCFtest[2].dot(ORCFtest[2].T))


            alphaTest = ORCFtest[1].dot(ORCFtest[2].T)/(bTest*cTest)
            betaTest  = ORCFtest[0].dot(ORCFtest[2].T)/(aTest*cTest)
            gammaTest = ORCFtest[0].dot(ORCFtest[1].T)/(aTest*bTest)

            alphaDegreesTest = numpy.arccos(alphaTest)*180/numpy.pi
            betaDegreesTest  = numpy.arccos(betaTest)*180/numpy.pi
            gammaDegreesTest = numpy.arccos(gammaTest)*180/numpy.pi
        

            'in the case of ORCF try to convert cell to ORCF and see if alpha,beta,gamma=90 degrees'
            if numpy.abs(alphaDegreesTest-90.0) < 0.01 and numpy.abs(betaDegreesTest-90.0) < 0.01 and numpy.abs(gammaDegreesTest-90.0) < 0.01:
                return 10
                'if not then this is triclinic '
            else:
                return 14
        else:
            raise TypeError
        

        'in the case of ORCI body centered ortho'

def _convertCellBC(cellVectors,toPrimOrConv='conv'):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    returnMatrix = copy.deepcopy(cellVectors)
    convertMatrix = numpy.matrix([[ 1, 1,-1 ],
                                  [-1, 1, 1 ],
                                  [ 1,-1, 1]])


    div2matrix = numpy.matrix([[0.5,0.0,0.0],
                               [0.0,0.5,0.0],
                               [0.0,0.0,0.5]])

    convertMatrix = convertMatrix.dot(div2matrix)
    
    if toPrimOrConv=='prim':
        convertMatrix=numpy.linalg.inv(convertMatrix)
    
    returnMatrix = convertMatrix.dot(cellVectors)
    return returnMatrix
    '''needs to be removed sometime soon'''
    '''needs to be removed sometime soon'''

def _cellStringToMatrix(cellParamString):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''
    
    cellParamSplit = cellParamString.split('\n')                                                                  
                                                                                                                   
    splitList = []                                                                                                    

    splitList = [split.split() for split in cellParamSplit if len(split)]
    symMatrix = numpy.matrix(splitList).astype(numpy.float32)

    return symMatrix

def _cellMatrixToString(cellMatrix,indent=True):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    try:
        cellMatrix=cellMatrix.getA()
    except:
        pass
    outputString = ''
    if indent:
        for entry in range(len(cellMatrix)):
            outputString+=' '.join(['%16.10f' % x for x in cellMatrix[entry]])+'\n'
    else:
        for entry in range(len(cellMatrix)):
            outputString+=' '.join(['%10.10f' % x for x in cellMatrix[entry]])+'\n'

    return outputString

def _getCellInput(oneCalc,ID,scaled=True):
    '''
    Gets the primitive cell vectors from the input file

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:
          scaled (bool): scale vectors by alat

    Returns:
          symMatrix (numpy.matrix): matrix representing primive lattice vectors

    '''

    try:
        cellParamMatrixDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
        ibrav=cellParamMatrixDict['&system']['ibrav']
        celldmDict={'ibrav':int(ibrav)}
        for y in  [(x[0],float(x[1])) for x in list(cellParamMatrixDict['&system'].items()) if len(re.findall('celldm',x[0])) !=0 ]:
            fixed = re.sub('[)(\s]', '', y[0])
            celldmDict[fixed]=y[1]
                       
        
        cellParams = celldm2free(**celldmDict)

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)
        print('no CELL_PARAMETERS in input')
        return

    symMatrix = AFLOWpi.retr._cellStringToMatrix(cellParams)

    a = numpy.sqrt(symMatrix[0].dot(symMatrix[0].T))

    if scaled:
        symMatrix=symMatrix/a

    return symMatrix


#####################################################################################################################
#####################################################################################################################
## cell and position transform
#####################################################################################################################
#####################################################################################################################

def _getSymmInput(oneCalc,ID):
    '''


    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step

    Keyword Arguments:


    Returns:


    '''        

    symMatrix = getPositionsFromInput(oneCalc,ID)
    orig = copy.deepcopy(symMatrix)
    return orig


def _duplicateEdgeAtoms(symMatrix):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    for i in range(3):
        arrayList=[]
        for arrayMat in symMatrix:
            arrayVec=arrayMat.getA()
            '''if atom is on side of cell make duplicate on the other side'''
            if arrayVec[0][0]==1:
                add=arrayVec[0]-[1,0,0]
                arrayList.append(add.tolist())

            if arrayVec[0][1]==1:
                add=arrayVec[0]-[0,1,0]
                arrayList.append(add.tolist())

            if arrayVec[0][2]==1:
                add=arrayVec[0]-[0,0,1]
                arrayList.append(add.tolist())

            arrayList.append(copy.deepcopy(arrayVec[0].tolist()))            

        formMatrix = numpy.matrix(arrayList).astype(float)
        symMatrix = numpy.matrix(numpy.round(formMatrix,decimals=10))

    for i in range(3):
        arrayList=[]
        for arrayMat in symMatrix:
            arrayVec=arrayMat.getA()
            '''if atom is on side of cell make duplicate on the other side'''
            if arrayVec[0][0]==0:
                add=arrayVec[0]+[1,0,0]
                arrayList.append(add.tolist())

            if arrayVec[0][1]==0:
                add=arrayVec[0]+[0,1,0]
                arrayList.append(add.tolist())

            if arrayVec[0][2]==0:
                add=arrayVec[0]+[0,0,1]
                arrayList.append(add.tolist())

            arrayList.append(copy.deepcopy(arrayVec[0].tolist()))            

        formMatrix = numpy.matrix(arrayList).astype(float)
        symMatrix = numpy.matrix(numpy.round(formMatrix,decimals=10))

    return symMatrix


def _getCartConvMatrix(symMatrix,cellMatrix):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    '''get the conventional matrix from the primitive including the positions'''
    a = numpy.sqrt(cellMatrix[0].dot(cellMatrix[0].T))
    b = numpy.sqrt(cellMatrix[1].dot(cellMatrix[1].T))
    c = numpy.sqrt(cellMatrix[2].dot(cellMatrix[2].T))


    cosA = cellMatrix[1].dot(cellMatrix[2].T)/(b*c)
    cosB  = cellMatrix[0].dot(cellMatrix[2].T)/(a*c)
    cosG = cellMatrix[0].dot(cellMatrix[1].T)/(a*b)


    alphaRad = numpy.arccos(cosA)
    betaRad  = numpy.arccos(cosB)
    gammaRad = numpy.arccos(cosG)
    if numpy.abs(alphaRad)<0.00000001:
        alphaRad=0.0
    if numpy.abs(betaRad)<0.00000001:
        betaRad=0.0
    if numpy.abs(gammaRad)<0.00000001:
        gammaRad=0.0
    sinA = numpy.sin(alphaRad)
    sinB = numpy.sin(betaRad)
    sinG = numpy.sin(gammaRad)

    
    cartMatrix = numpy.matrix([[ a, b*cosG, c*cosB                     ],
                               [ 0, b*sinG,-c*sinB*cosA                ],
                               [ 0, 0,      c*sinB*sinA                ],])


    return cartMatrix

def _convertCartesian(symMatrix,cellMatrix,scaleFactor=1.0):
    '''
    Converts atomic positions from crystal to cartesian coordinates

    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form

    Keyword Arguments:
          scaleFactor (int): scaling factor for output matrix

    Returns:
          in_cart (postions in cartiesian coordinates)

    '''

    try:
        cellMatrix = cellMatrix.getA()
    except:
        pass 

    try:
        symMatrix = symMatrix.getA()
    except:
        pass 

    cellMatrix=cellMatrix.astype(numpy.float)
    symMatrix=symMatrix.astype(numpy.float)
    returnMatrix=copy.deepcopy(symMatrix)

    for entry in range(len(returnMatrix)):
        returnMatrix[entry]=cellMatrix.dot(returnMatrix[entry].T).T

    in_cart = numpy.matrix(numpy.round(numpy.matrix(returnMatrix).astype(numpy.float),decimals=8))

    return in_cart

def _convertFractional(symMatrix,cellMatrix,scaleFactor=1):
    '''
    Converts atomic positions from cartesian to crystal coordinates

    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form

    Keyword Arguments:
          scaleFactor (int): scaling factor for output matrix

    Returns:
          in_cart (numpy.matrix): postions in crystal coordinates

    '''

    try:
        symMatrix=symMatrix.getA()
    except:
        pass
    try:
        cellMatrix=cellMatrix.getA()
    except:
        pass


    cartMatrix = AFLOWpi.retr._getCartConvMatrix(symMatrix,cellMatrix)

    cartMatrixInv=numpy.linalg.inv(cartMatrix).getA()
    cartMatrixInv=numpy.linalg.inv(cellMatrix).T

    for entry in range(len(symMatrix)):
        symMatrix[entry]=cartMatrixInv.dot(symMatrix[entry])

    returnMat = numpy.matrix(numpy.round(numpy.matrix(symMatrix).astype(numpy.float),decimals=10))

    return returnMat




def _expandBoundariesNoScale(labels,symMatrix,numX,numY,numZ,inList=False,expand=True):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form

    Keyword Arguments:


    Returns:


    '''

    symMatrix=symMatrix.getA()
    symMatrixList=symMatrix.tolist()
    superList=[]

    if expand==True:
        for item in range(len(symMatrixList)):
            superList.append([labels[item],symMatrixList[item][0],symMatrixList[item][1],symMatrixList[item][2],])

    for entry in range(len(symMatrix)):
        superList.append([labels[entry],symMatrix[entry][0]+numX,symMatrix[entry][1]+numY,symMatrix[entry][2]+numZ])


                    
    labels = [x[0] for x in superList]
    if inList==False:
        symMatrix= numpy.matrix([x[1:] for x in superList])
    else:
        symMatrix=[x[1:] for x in superList]

    return labels,symMatrix


def _shiftX(coordVec,shift):
    '''


    Arguments:
          coordVec (numpy.matrix): atomic positions in matrix form

    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)

    try:
        arrayVec = coordVec.getA()
    except:
        arrayVec=coordVec

    try:
        arrayVec[0][0]
        arrayVec=arrayVec[0]
    except:
        pass

    arrayVec[0]+=shift


    return numpy.matrix(arrayVec)

def _shiftY(coordVec,shift):
    '''


    Arguments:
          coordVec (numpy.matrix): atomic positions in matrix form

    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)

    try:
        arrayVec = coordVec.getA()
    except:
        arrayVec=coordVec

    try:
        arrayVec[0][0]
        arrayVec=arrayVec[0]
    except:
        pass



    arrayVec[1]+=shift

    return numpy.matrix(arrayVec)

def _shiftZ(coordVec,shift):
    '''


    Arguments:
          coordVec (numpy.matrix): atomic positions in matrix form

    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)

    try:
        arrayVec = coordVec.getA()
    except:
        arrayVec=coordVec

    try:
        arrayVec[0][0]
        arrayVec=arrayVec[0]
    except:
        pass

    arrayVec[2]+=shift


    return numpy.matrix(arrayVec)

def _checkEqualPoint(oldMatrix,newMatrix):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    return oldMatrix==newMatrix


def _shiftAfterRotation(coordVec):
    '''


    Arguments:
          coordVec (numpy.matrix): atomic positions in matrix form

    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)
    arrayVec = coordVec.getA()
    '''shift back to centering around the middle of cell'''
    arrayVec[0][0]+=0.5
    arrayVec[0][1]+=0.5
    arrayVec[0][2]+=0.5


    return numpy.matrix(arrayVec)




def getPointGroup(oneCalc,ID,source='input'):
    '''
    NOT COMPLETED.

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
          ID (str): ID string for the particular calculation and step


    Keyword Arguments:
          source (string): either 'input' or 'output' to get point group from input
                           or output atomic positions

    Returns:
          None

    '''

    if source!='input' and source!='output':
        print('source must equal "input" or "output"')
        logging.error('source must equal "input" or "output"')
        return 
    if source=='input':
        cellInput = AFLOWpi.retr._getCellInput(oneCalc,ID)
        symInput = AFLOWpi.retr._getSymmInput(oneCalc,ID)
    if source=='output':
        symInput = getPositionsFromOutput(oneCalc,ID)
        alat,symInput = AFLOWpi.retr._getCellParams(oneCalc,ID)


        
    ibrav = AFLOWpi.retr.getIbravFromVectors(cellInput)

        

    

    
def _sortMatrix(matrix):
    '''


    Arguments:
          matrix (numpy.matrix): sorts the matrix

    Keyword Arguments:


    Returns:


    '''

    sortedMatrix = matrix.tolist()
    uniqueSet=set()
    """if atom is at side of cell pull it back"""
    for arrayVec in sortedMatrix:
        '''add it to a set so we don't have repeats'''
        uniqueSet.add(tuple(arrayVec))
    sortedMatrix=list(uniqueSet)
    '''sort it so we can have a proper comparison'''
    sortedMatrix = sorted(sortedMatrix)
    return numpy.matrix(sortedMatrix)

def _shiftAfter(matrix):
    '''


    Arguments:
          matrix (numpy.matrix): sorts the matrix

    Keyword Arguments:


    Returns:


    '''

    arrayList=[]
    try:
        matrix=matrix.getA()
    except:
        matrix=matrix



    for arrayVec in matrix:
        arrayVec=arrayVec.tolist()

        '''shift positions in cell to account for sides of cell'''
        if arrayVec[0]>=1:
            arrayVec[0]-=numpy.floor(arrayVec[0])/1.0
        if arrayVec[1]>=1:
            arrayVec[1]-=numpy.floor(arrayVec[1])/1.0
        if arrayVec[2]>=1:
            arrayVec[2]-=numpy.floor(arrayVec[2])/1.0

        if arrayVec[0]<0:
            arrayVec[0]+=numpy.abs(numpy.floor(arrayVec[0]))/1.0
        if arrayVec[1]<0:
            arrayVec[1]+=numpy.abs(numpy.floor(arrayVec[1]))/1.0
        if arrayVec[2]<0:
            arrayVec[2]+=numpy.abs(numpy.floor(arrayVec[2]))/1.0
        arrayList.append(arrayVec)
    

    MATRIX=numpy.matrix(numpy.round(numpy.matrix(arrayList).astype(numpy.float32),decimals=8))

    return numpy.matrix(numpy.round(numpy.matrix(arrayList).astype(numpy.float32),decimals=8))


def compMatrices(matrix1,matrix2):
    '''
    Compare two numpy matrices to see if they are equal or not

    Arguments:
          matrix1 (numpy.matrix): first matrix
          matrix2 (numpy.matrix): second matrix

    Keyword Arguments:
          None

    Returns:
          equalBool (bool): True if they're equal False if they're not

    '''

    matrix1 =  AFLOWpi.retr._shiftCell(matrix1)
    matrix2 =  AFLOWpi.retr._shiftCell(matrix2)

    matrix1 = AFLOWpi.retr._reduceDuplicates(matrix1)
    matrix2 = AFLOWpi.retr._reduceDuplicates(matrix2)

    sortedMatrix1 = AFLOWpi.retr._sortMatrix(matrix1)
    sortedMatrix2 =  AFLOWpi.retr._sortMatrix(matrix2)

    try:
        equalBool = numpy.matrix(sortedMatrix1 == sortedMatrix2).all()
        if equalBool==False:
            matrix1 = AFLOWpi.retr._duplicateEdgeAtoms(matrix1)
            matrix2 = AFLOWpi.retr._duplicateEdgeAtoms(matrix2)
            matrix1 = AFLOWpi.retr._reduceDuplicates(matrix1)
            matrix2 = AFLOWpi.retr._reduceDuplicates(matrix2)

            sortedMatrix1 = AFLOWpi.retr._sortMatrix(matrix1)
            sortedMatrix2 =  AFLOWpi.retr._sortMatrix(matrix2)
            equalBool = numpy.matrix(sortedMatrix1 == sortedMatrix2).all()

        return equalBool
    except:

        return False

        
        


def invertXYZ(symMatrix,cellMatrix):    
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    symMatrix = AFLOWpi.retr._convertCartesian(symMatrix,cellMatrix)
    outputMatrix=[]
    for coordMatrix in symMatrix:
        outputMatrix.append(__invertXYZ(coordMatrix).tolist()[0])


    outputMatrix = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(numpy.float32),decimals=10))
    symMatrix = AFLOWpi.retr._convertFractional(outputMatrix,cellMatrix)

    return symMatrix


        

def invertX(symMatrix,cellMatrix):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    symMatrixCopy=copy.deepcopy(symMatrix)
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)
    outputMatrix=[]

    for coordMatrix in symMatrixCopy:
        outputMatrix.append(AFLOWpi.retr._invertX(coordMatrix).tolist()[0])

    outputMatrix = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))

    outputMatrix = AFLOWpi.retr._convertFractional(outputMatrix,cellMatrix)

    for coordMatrix in range(len(outputMatrix)):
        outputMatrix[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(outputMatrix[coordMatrix])


    return outputMatrix
def invertY(symMatrix,cellMatrix):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    symMatrixCopy=copy.deepcopy(symMatrix)
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)
    outputMatrix=[]

    for coordMatrix in symMatrixCopy:
        outputMatrix.append(AFLOWpi.retr._invertY(coordMatrix).tolist()[0])

    outputMatrix = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))

    outputMatrix = AFLOWpi.retr._convertFractional(outputMatrix,cellMatrix)

    for coordMatrix in range(len(outputMatrix)):
        outputMatrix[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(outputMatrix[coordMatrix])


    return outputMatrix

def invertZ(symMatrix,cellMatrix):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    symMatrixCopy=copy.deepcopy(symMatrix)
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)
    outputMatrix=[]

    for coordMatrix in symMatrixCopy:
        outputMatrix.append(AFLOWpi.retr._invertZ(coordMatrix).tolist()[0])

    outputMatrix = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))

    outputMatrix = AFLOWpi.retr._convertFractional(outputMatrix,cellMatrix)

    for coordMatrix in range(len(outputMatrix)):
        outputMatrix[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(outputMatrix[coordMatrix])


    return outputMatrix


        
def shiftX(symMatrix,cellMatrix,shift):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    symMatrixCopy=copy.deepcopy(symMatrix)
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)
    outputMatrix=[]

    a = numpy.sqrt(cellMatrix[0].dot(cellMatrix[0].T))

    for coordMatrix in symMatrixCopy:
        outputMatrix.append(AFLOWpi.retr._shiftX(coordMatrix,shift*a).tolist()[0])

    outputMatrix = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))
    symMatrixCopy = AFLOWpi.retr._convertFractional(outputMatrix,cellMatrix)

    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(symMatrixCopy[coordMatrix])


    return symMatrixCopy


def shiftY(symMatrix,cellMatrix,shift):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    symMatrixCopy=copy.deepcopy(symMatrix)
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)
    outputMatrix=[]

    b = numpy.sqrt(cellMatrix[1].dot(cellMatrix[1].T))

    for coordMatrix in symMatrixCopy:
        outputMatrix.append(AFLOWpi.retr._shiftY(coordMatrix,shift*b).tolist()[0])

    outputMatrix = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))
    symMatrixCopy = AFLOWpi.retr._convertFractional(outputMatrix,cellMatrix)

    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(symMatrixCopy[coordMatrix])


    return symMatrixCopy




def shiftZ(symMatrix,cellMatrix,shift):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    symMatrixCopy=copy.deepcopy(symMatrix)
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)
    outputMatrix=[]

    c = numpy.sqrt(cellMatrix[2].dot(cellMatrix[2].T))

    for coordMatrix in symMatrixCopy:
        outputMatrix.append(AFLOWpi.retr._shiftZ(coordMatrix,shift*c).tolist()[0])

    outputMatrix = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))
    symMatrixCopy = AFLOWpi.retr._convertFractional(outputMatrix,cellMatrix)

    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(symMatrixCopy[coordMatrix])


    return symMatrixCopy


        
def glideXshiftY(symMatrix,cellMatrix,shift):
    return shiftY(invertX(symMatrix,cellMatrix),cellMatrix,shift)
def glideXshiftZ(symMatrix,cellMatrix,shift):
    return shiftZ(invertX(symMatrix,cellMatrix),cellMatrix,shift)
def glideYshiftZ(symMatrix,cellMatrix,shift):
    return shiftZ(invertY(symMatrix,cellMatrix),cellMatrix,shift)
def glideYshiftX(symMatrix,cellMatrix,shift):
    return shiftX(invertY(symMatrix,cellMatrix),cellMatrix,shift)
def glideZshiftX(symMatrix,cellMatrix,shift):
    return shiftX(invertZ(symMatrix,cellMatrix),cellMatrix,shift)
def glideZshiftY(symMatrix,cellMatrix,shift):
    return shiftY(invertZ(symMatrix,cellMatrix),cellMatrix,shift)




def rotateAlpha(symMatrix,cellMatrix,angle):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    outputMatrix=[]
    symMatrixCopy=copy.deepcopy(symMatrix)
    'shift fractional coords to center around 0,0,0'
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    'convert to cartesian'
    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)

    'rotate in cartesian around z axis'
    for coordMatrix in symMatrixCopy:
        coordMatrix = AFLOWpi.retr._rotateAlpha(coordMatrix,angle)
        outputMatrix.append(coordMatrix.tolist()[0])

    symMatrixCopy = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))

    'transform back to fractional coords'
    symMatrixCopy = AFLOWpi.retr._convertFractional(symMatrixCopy,cellMatrix)



    'shift back to center around 0.5,0.5,0.5'
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(symMatrixCopy[coordMatrix])

    return symMatrixCopy

def rotateBeta(symMatrix,cellMatrix,angle):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form

    Keyword Arguments:


    Returns:


    '''

    outputMatrix=[]
    symMatrixCopy=copy.deepcopy(symMatrix)
    'shift fractional coords to center around 0,0,0'
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    'convert to cartesian'
    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)

    'rotate in cartesian around z axis'
    for coordMatrix in symMatrixCopy:
        coordMatrix = AFLOWpi.retr._rotateBeta(coordMatrix,angle)
        outputMatrix.append(coordMatrix.tolist()[0])

    symMatrixCopy = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))

    'transform back to fractional coords'
    symMatrixCopy = AFLOWpi.retr._convertFractional(symMatrixCopy,cellMatrix)



    'shift back to center around 0.5,0.5,0.5'
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(symMatrixCopy[coordMatrix])

    return symMatrixCopy

def rotateGamma(symMatrix,cellMatrix,angle):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form
          cellMatrix (numpy.matrix): primitive lattice vectors in matrix form


    Keyword Arguments:


    Returns:


    '''

    outputMatrix=[]
    symMatrixCopy=copy.deepcopy(symMatrix)
    'shift fractional coords to center around 0,0,0'
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftBeforeRotation(symMatrixCopy[coordMatrix])

    'convert to cartesian'
    symMatrixCopy = AFLOWpi.retr._convertCartesian(symMatrixCopy,cellMatrix)

    'rotate in cartesian around z axis'
    for coordMatrix in symMatrixCopy:
        coordMatrix = AFLOWpi.retr._rotateGamma(coordMatrix,angle)
        outputMatrix.append(coordMatrix.tolist()[0])

    symMatrixCopy = numpy.matrix(numpy.round(numpy.matrix(outputMatrix).astype(float),decimals=10))

    'transform back to fractional coords'
    symMatrixCopy = AFLOWpi.retr._convertFractional(symMatrixCopy,cellMatrix)



    'shift back to center around 0.5,0.5,0.5'
    for coordMatrix in range(len(symMatrixCopy)):
        symMatrixCopy[coordMatrix] = AFLOWpi.retr._shiftAfterRotation(symMatrixCopy[coordMatrix])

    return symMatrixCopy
        

def _invertXYZ(coordVec,xScale,yScale,zScale):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)
    arrayVec = coordVec.getA()
    arrayVec[0][0]=-arrayVec[0][0]
    arrayVec[0][1]=-arrayVec[0][1]
    arrayVec[0][2]=-arrayVec[0][2]

    return numpy.matrix(arrayVec)

def _invertX(coordVec,xScale=1):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)
    arrayVec = coordVec.getA()
    arrayVec[0][0]=-(arrayVec[0][0])
    return numpy.matrix(arrayVec)

def _invertY(coordVec,yScale=1):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)
    arrayVec = coordVec.getA()
    arrayVec[0][1]=-(arrayVec[0][1])
    return numpy.matrix(arrayVec)

def _invertZ(coordVec,zScale=1):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)
    arrayVec = coordVec.getA()
    arrayVec[0][2]=-(arrayVec[0][2])
    return numpy.matrix(arrayVec)


                      
def _rotateAlpha(coordVec,angle):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    angle = numpy.radians(angle,numpy.ndarray([angle]))[0]
    coordVec=copy.deepcopy(coordVec)

    rotationMatrix =     numpy.matrix([[       1         ,       0         ,    0               ],
                                       [       0          ,numpy.cos(angle),-numpy.sin(angle)   ],
                                       [       0          ,numpy.sin(angle), numpy.cos(angle)   ], ])



    coordVec = rotationMatrix.dot(coordVec.T).T


    return coordVec



def _rotateBeta(coordVec,angle):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    angle = numpy.radians(angle,numpy.ndarray([angle]))[0]
    coordVec=copy.deepcopy(coordVec)



    rotationMatrix =     numpy.matrix([[numpy.cos(angle) ,     0           ,numpy.sin(angle)    ],
                                       [     0           ,     1           ,    0               ],
                                       [-numpy.sin(angle),     0           ,numpy.cos(angle)    ], ])

    


    coordVec = rotationMatrix.dot(coordVec.T).T


    return coordVec


    


def _rotateGamma(coordVec,angle):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    angle = numpy.radians(angle,numpy.ndarray([angle]))[0]
    coordVec=copy.deepcopy(coordVec)


    rotationMatrix =     numpy.matrix([[numpy.cos(angle) ,-numpy.sin(angle),      0             ],
                                       [numpy.sin(angle) , numpy.cos(angle),      0             ],
                                       [        0         ,     0           ,     1             ], ])

    coordVec = rotationMatrix.dot(coordVec.T).T
    return coordVec

def _shiftBeforeRotation(coordVec):
    '''


    Arguments:


    Keyword Arguments:


    Returns:


    '''

    coordVec=copy.deepcopy(coordVec)
    arrayVec = coordVec.getA()
    '''shift to center all atoms around cartesian origin'''
    arrayVec[0][0]-=0.5
    arrayVec[0][1]-=0.5
    arrayVec[0][2]-=0.5
    return numpy.matrix(arrayVec)

def _shiftCell(symMatrix):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form

    Keyword Arguments:


    Returns:


    '''

    arrayList=[]
    for arrayMat in symMatrix:
        arrayVec = arrayMat.getA()
        '''shift positions in cell to account for sides of cell'''
        if arrayVec[0][0]>1:
            shift = numpy.floor(arrayVec[0][0])
            arrayMat-= numpy.array([1,0,0])*shift
        if arrayVec[0][1]>1:
            shift = numpy.floor(arrayVec[0][1])
            arrayMat-= numpy.array([0,1,0])*shift
        if arrayVec[0][2]>1:
            shift = numpy.floor(arrayVec[0][2])
            arrayMat-= numpy.array([0,0,1])*shift
    for arrayMat in symMatrix:
        arrayVec = arrayMat.getA()
        if arrayVec[0][0]<1:
            shift = numpy.abs(numpy.floor(arrayVec[0][0]))
            arrayMat+= numpy.array([1,0,0])*shift
        if arrayVec[0][1]<1:
            shift = numpy.abs(numpy.floor(arrayVec[0][1]))
            arrayMat+= numpy.array([0,1,0])*shift
        if arrayVec[0][2]<1:
            shift = numpy.abs(numpy.floor(arrayVec[0][2]))
            arrayMat+= numpy.array([0,0,1])*shift



        arrayList.append(arrayVec[0])


    outMatrix = numpy.matrix(numpy.round(numpy.matrix(arrayList).astype(float),decimals=10))


    return outMatrix



def shiftCell(symMatrix):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form

    Keyword Arguments:


    Returns:


    '''

    outMatrix = AFLOWpi.retr._shiftAfter(symMatrix)

    outMatrix = AFLOWpi.retr._reduceDuplicates(outMatrix)

    outMatrix = AFLOWpi.retr._duplicateEdgeAtoms(outMatrix)
    outMatrix = AFLOWpi.retr._reduceDuplicates(outMatrix)
    outMatrix = AFLOWpi.retr._shiftAfter(outMatrix)
    outMatrix = AFLOWpi.retr._reduceDuplicates(outMatrix)

    return outMatrix

def _reduceDuplicates(symMatrix):
    '''


    Arguments:
          symMatrix (numpy.matrix): atomic positions in matrix form

    Keyword Arguments:


    Returns:


    '''
    
    arrayListCopy=set()

    '''to get rid of duplicates'''
    for arrayMat in symMatrix:
        arrayVec=arrayMat.tolist()
        arrayListCopy.add(tuple(arrayVec[0]))
    arrayListCopy = list(arrayListCopy)
    outMatrix = numpy.matrix(numpy.round(numpy.matrix(arrayListCopy).astype(float),decimals=10))
    return outMatrix


def _get_cell_mass(oneCalc,ID):
    num_of_each = AFLOWpi.retr._getAtomNum(oneCalc['_AFLOWPI_INPUT_'],strip=True)

    total_mass=0.0
    for atom,number in num_of_each.items():
        total_mass+=float(AFLOWpi.prep._getAMass(atom))*float(number)

    return total_mass
#####################################################################################################################
#####################################################################################################################
## misc
#####################################################################################################################
#####################################################################################################################
def chemAsKeys(calcs):
    '''
    For sets of calcs that only differ by chemistry you can replace the
    hash by the chemical name. May not always work. Especially if there
    are two compounds with swapped positions of atoms of different elements

    Arguments:
	  calcs (dict): dictionary of dictionaries of calculations

    Keyword Arguments:
          None

    Returns:
	  calcsCopy (dict): returns new dictionay with keys as the chemical stoiciometry

    '''

    calcsCopy = OrderedDict()
    stoicNameList=[]
    for ID,oneCalc in OrderedDict(calcs).items():
        newKey=AFLOWpi.retr._getStoicName(oneCalc)
        stoicNameList.append(newKey)
        calcsCopy[newKey]=oneCalc


    if len(set(stoicNameList))!=len(calcs):
        logging.warning('cannot uniquely identify each calc with a corresponding chem. returning original calcs')
        print('cannot uniquely identify each calc with a corresponding chem. returning original calcs')
    

    return calcsCopy

import functools

import numpy as np
from numpy import linalg as la


__all__ = ('aflow_prim2conv', 'aflow_conv2prim', 'qe_prim2conv', 'qe_conv2prim',
           'conv_aflow2qe', 'conv_qe2aflow', 'prim_aflow2qe', 'prim_qe2aflow',
           'InvalidIbravError')


class InvalidIBravError(ValueError):
    """Supplied Bravais lattice index (ibrav) is not valid.
    
    Attributes:
        value: The invalid Bravais lattice index that caused this error.
    """
    
    def __init__(self, value):
        self.value = value
        super(ValueError, self).__init__('%s is not a valid Bravais lattice' % value)


def validate_ibrav(func):
    """Decorates a basis transformation function to validate the Bravais lattice
    index in the range [0, 14].
    
    Arguments:
        func: The decorated function of two arguments, cell and ibrav.
    
    Returns:
        Decorated function that raises `InvalidIBravError` when bad values are
        encountered.
    """
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        
        ibrav = args[1]
        
        try:
            assert ibrav in range(15)
        
        except AssertionError:
            raise InvalidIBravError(ibrav)
        
        return func(*args, **kwargs)
    
    return decorator


def matmul(a, b, *others):
    """Calculates the accumulated matrix or dot product of the arguments."""
    args = (a, b) + others
    return functools.reduce(np.dot, args)


class IBrav(object):
    """ibrav namespace"""
    FREE = 0
    CUB = 1
    FCC = 2
    BCC = 3
    HEX = 4
    RHL = 5
    TET = 6
    BCT = 7
    ORC = 8
    ORCC = 9
    ORCF = 10
    ORCI = 11
    MCL = 12
    MCLC = 13
    TRI = 14


def aflow_primitive_to_aflow_conventional_transform(ibrav):
    """
    Builds a transformation matrix for converting primtive cell vectors to
    conventional cell vectors using AFLOW convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    # face-centered lattices
    if ibrav in (IBrav.FCC, IBrav.ORCF):
        transform = np.matrix([[-1.,  1.,  1.],
                               [ 1., -1.,  1.],
                               [ 1.,  1., -1.]])
        
    # body-centered lattices
    elif ibrav in (IBrav.BCC, IBrav.BCT, IBrav.ORCI):
        transform = np.matrix([[ 0.,  1.,  1.],
                               [ 1.,  0.,  1.],
                               [ 1.,  1.,  0.]])
        
    # orthorhombic base-centered lattice
    elif ibrav == IBrav.ORCC:
        transform = np.matrix([[ 1.,  1.,  0.],
                               [-1.,  1.,  0.],
                               [ 0.,  0.,  1.]])
        
    # monoclinic base-centered lattice
    elif ibrav == IBrav.MCLC:
        transform = np.matrix([[ 1., -1.,  0.],
                               [ 1.,  1.,  0.],
                               [ 0.,  0.,  1.]])
        
    # lattices without centering share the same cell vectors between
    # their primitive cells and conventional cells, so no change is necessary
    else:
        transform = np.matrix(np.identity(3))
    
    return transform
    
    
def aflow_conventional_to_aflow_primitive_transform(ibrav):
    """
    Builds a transformation matrix for converting conventional cell vectors to
    primitive cell vectors using AFLOW convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    # primitive to conventional is the inverse of the transform we want to do
    inverse_transform = aflow_primitive_to_aflow_conventional_transform(ibrav)
    
    # invert the inverse transform
    return la.inv(inverse_transform)


def qe_primitive_to_qe_conventional_transform(ibrav):
    """
    Builds a transformation matrix for converting primitive cell vectors to
    conventional cell vectors using Quantum Espresso convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    if ibrav == IBrav.FCC:
        transform = np.matrix([[-1.,  1., -1.],
                               [-1.,  1.,  1.],
                               [ 1.,  1., -1.]])
        
    # body-centered lattices
    elif ibrav in (IBrav.BCC, IBrav.ORCI):
        transform = np.matrix([[ 1., -1.,  0.],
                               [ 0.,  1., -1.],
                               [ 1.,  0.,  1.]])
        
    # can't stuff BCT in with the other body-centered lattices
    # like in AFLOW because QE uniquely defines its lattice vectors
    elif ibrav == IBrav.BCT:
        transform = np.matrix([[ 1.,  0., -1.],
                               [-1.,  1.,  0.],
                               [ 0.,  1.,  1.]])
    
    # orthorhombic base-centered lattice
    elif ibrav == IBrav.ORCC:
        transform = np.matrix([[ 1., -1.,  0.],
                               [ 1.,  1.,  0.],
                               [ 0.,  0.,  1.]])
        
    # orthorhombic face-centered lattice
    elif ibrav == IBrav.ORCF:
        transform = np.matrix([[ 1.,  1., -1.],
                               [-1.,  1.,  1.],
                               [ 1., -1.,  1.]])
    
    # monoclinic base-centered lattice
    elif ibrav == IBrav.MCLC:
        transform = np.matrix([[ 1.,  0.,  1.],
                               [ 0.,  1.,  0.],
                               [-1.,  0.,  1.]])
        
    # lattices without centering share the same cell vectors between
    # their primitive cells and conventional cells, so no change is necessary
    else:
        transform = np.matrix(np.identity(3))
    
    return transform


def qe_conventional_to_qe_primitive_transform(ibrav):
    """
    Builds a transformation matrix for converting conventional cell vectors to
    primitive cell vectors using Quantum Espresso convention.
    """
    # primitive to conventional is the inverse of the transform we want to do
    inverse_transform = qe_primitive_to_qe_conventional_transform(ibrav)
    
    # invert the inverse transform
    return la.inv(inverse_transform)


def aflow_conventional_to_qe_conventional_transform(ibrav, angle=None):
    """
    Builds a transformation matrix for converting conventional cell vectors
    using AFLOW convention to conventional cell vectors using Quantum Espresso
    convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        angle (float): Characteristic angle for the lattice structure, only
                       necessary for rhombohedral (RHL) lattices
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    if ibrav == IBrav.HEX:
        transform = np.matrix([[ 1.,  1.,  0.],
                               [-1.,  0.,  0.],
                               [ 0.,  0.,  1.]])
    
    elif ibrav == IBrav.RHL:
        raise NotImplementedError
        
        assert angle is not None, "characteristic angle must be supplied for RHL"
        
        cos = np.cos(angle)
        cos_half = np.cos(angle / 2)
        sin_half = np.sin(angle / 2)
        
        aflow = np.matrix([
            [      cos_half, -sin_half,                                0],
            [      cos_half,  sin_half,                                0],
            [cos / cos_half,         0, np.sqrt(1 - (cos / cos_half)**2)]
        ])
        
        tx = np.sqrt((1 - cos) / 2)
        ty = np.sqrt((1 - cos) / 6)
        tz = np.sqrt((1 + 2*cos) / 3)
        
        qe = np.matrix([
            [ tx,  -ty, tz],
            [  0, 2*ty, tz],
            [-tx,  -ty, tz]
        ])
        
        transform = matmul(qe, la.inv(aflow))
    
    elif ibrav in (IBrav.MCL, IBrav.MCLC):
        raise NotImplementedError
    
    else:
        transform = np.identity(3)
    
    return transform


def qe_conventional_to_aflow_conventional_transform(ibrav, angle=None):
    """
    Builds a transformation matrix for converting conventional cell vectors
    using Quantum Espresso convention to conventional cell vectors using AFLOW
    convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        angle (float): Characteristic angle for the lattice structure, only
                       necessary for rhombohedral (RHL) lattices
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    inverse_transform = aflow_conventional_to_qe_conventional_transform(ibrav, angle)
    return la.inv(inverse_transform)


def aflow_primitive_to_qe_primitive_transform(ibrav, angle=None):
    """
    Builds a transformation matrix for converting primitive cell vectors using
    AFLOW convention to primitive cell vectors using Quantum Espresso convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    # successive transforms AFLOW prim -> AFLOW conv -> QE conv -> QE prim
    transforms = (aflow_primitive_to_aflow_conventional_transform(ibrav),
                  aflow_conventional_to_qe_conventional_transform(ibrav, angle),
                  qe_conventional_to_qe_primitive_transform(ibrav))
    
    return matmul(*reversed(transforms))


def qe_primitive_to_aflow_primitive_transform(ibrav, angle=None):
    """
    Builds a transformation matrix for converting primitive cell vectors using
    Quantum Espresso convention to primitive cell vectors using AFLOW convention.
    
    Arguments:
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        
    Returns:
        3x3 `numpy.matrix` transformation matrix.
    """
    # successive transforms QE prim -> QE conv -> AFLOW conv -> AFLOW prim
    transforms = (qe_primitive_to_qe_conventional_transform(ibrav),
                  qe_conventional_to_aflow_conventional_transform(ibrav, angle),
                  aflow_conventional_to_aflow_primitive_transform(ibrav))
    
    return matmul(*reversed(transforms))


@validate_ibrav
def aflow_primitive_to_aflow_conventional(cell, ibrav):
    """
    Converts primitive cell vectors to conventional cell vectors using AFLOW
    convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW primitive cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in AFLOW conventional lattice vectors.
    """
    transform = aflow_primitive_to_aflow_conventional_transform(ibrav)
    return matmul(transform, cell)
    

@validate_ibrav
def aflow_conventional_to_aflow_primitive(cell, ibrav):
    """
    Converts conventional cell vectors to primitive cell vectors using AFLOW
    convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in AFLOW primitive lattice vectors.
    """
    transform = aflow_conventional_to_aflow_primitive_transform(ibrav)
    return matmul(transform, cell)


@validate_ibrav
def qe_primitive_to_qe_conventional(cell, ibrav):
    """
    Converts primitive cell vectors to conventional cell vectors using Quantum
    Espresso convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 Quantum Espresso conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso primitive lattice
        vectors.
    """
    transform = qe_primitive_to_qe_conventional_transform(ibrav)
    return matmul(transform, cell)


@validate_ibrav
def qe_conventional_to_qe_primitive(cell, ibrav):
    """
    Converts conventional cell vectors to primitive cell vectors using Quantum
    Espresso convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 Quantum Espresso conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso primitive lattice
        vectors.
    """
    transform = qe_conventional_to_qe_primitive_transform(ibrav)
    return matmul(transform, cell)


@validate_ibrav
def aflow_conventional_to_qe_conventional(cell, ibrav, angle=None):
    """
    Converts conventional cell vectors using AFLOW convention to conventional
    cell vectors using Quantum Espresso convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        angle (float): Characteristic angle for the lattice structure, only
                       necessary for rhombohedral (RHL) lattices
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso conventional lattice
        vectors.
    """
    transform = aflow_conventional_to_qe_conventional_transform(ibrav, angle)
    return matmul(transform, cell)


@validate_ibrav
def qe_conventional_to_aflow_conventional(cell, ibrav, angle=None):
    """
    Converts conventional cell vectors using Quantum Espresso convention to
    conventional cell vectors using AFLOW convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 Quantum Espresso conventional cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
        angle (float): Characteristic angle for the lattice structure, only
                       necessary for rhombohedral (RHL) lattices
    
    Returns:
        3x3 `numpy.matrix` of the cell in AFLOW conventional lattice vectors.
    """
    transform = qe_conventional_to_aflow_conventional_transform(ibrav, angle)
    return matmul(transform, cell)


@validate_ibrav
def aflow_primitive_to_qe_primitive(cell, ibrav, angle=None):
    """
    Converts primitive cell vectors using AFLOW convention to primitive cell
    vectors using Quantum Espresso convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW primitive cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso primitive lattice
        vectors.
    """
    transform = aflow_primitive_to_qe_primitive_transform(ibrav, angle)
    return matmul(transform, cell)


@validate_ibrav
def qe_primitive_to_aflow_primitive(cell, ibrav, angle=None):
    """
    Converts primitive cell vectors using Quantum Espresso convention to
    primitive cell vectors using AFLOW convention.
    
    Arguments:
        cell (numpy.matrix): 3x3 AFLOW primitive cell vectors
        ibrav (int): Quantum Espresso Bravais lattice index of the crystal
                     structure
    
    Returns:
        3x3 `numpy.matrix` of the cell in Quantum Espresso primitive lattice
        vectors.
    """
    transform = qe_primitive_to_aflow_primitive_transform(ibrav, angle)
    return matmul(transform, cell)


# aliases
aflow_prim2conv = aflow_primitive_to_aflow_conventional
aflow_conv2prim = aflow_conventional_to_aflow_primitive
conv_aflow2qe = aflow_conventional_to_qe_conventional
conv_qe2aflow = qe_conventional_to_aflow_conventional
qe_prim2conv = qe_primitive_to_qe_conventional
qe_conv2prim = qe_conventional_to_qe_primitive
prim_aflow2qe = aflow2qe = aflow_primitive_to_qe_primitive
prim_qe2aflow = qe2aflow = qe_primitive_to_aflow_primitive

