import re
import copy
import logging 
import os
import AFLOWpi
import numpy

def tb_prep(oneCalc,ID):
#############################################################################################
	espressodir = AFLOWpi.prep._ConfigSectionMap('prep','engine_dir')
        if os.path.isabs(espressodir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            espressodir =  os.path.join(configFileLocation, espressodir)
#############################################################################################
	pdosExec = 'projwfc.x'
	pdosExec = os.path.join(espressodir,pdosExec)
        if os.path.exists(pdosExec)==False:
            print 'ProjectWFC executable not found. Check your config file to make sure engine_dir path is correct and that projwfc.x is in that directory..Exiting'
            logging.error('ProjectWFC executable not found. Check your config file to make sure engine_dir path is correct and that projwfc.x is in that directory..Exiting')
            raise SystemExit
	dosExec = 'dos.x'
	dosExec = os.path.join(espressodir,dosExec)
        if os.path.exists(dosExec)==False:
            print 'DOS executable not found. Check your config file to make sure engine_dir path is correct and that dos.x is in that directory..Exiting'
            logging.error('DOS executable not found. Check your config file to make sure engine_dir path is correct and that dos.x is in that directory..Exiting')
            raise SystemExit
	scfExec = 'pw.x'
	scfExec = os.path.join(espressodir,scfExec)
        if os.path.exists(scfExec)==False:
            print 'SCF executable not found. Check your config file to make sure engine_dir path is correct and that pw.x is in that directory..Exiting'
            logging.error('SCF executable not found. Check your config file to make sure engine_dir path is correct and that pw.x is in that directory..Exiting')
            raise SystemExit


	if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
		
                AFLOWpi.prep.totree(pdosExec,{ID:oneCalc})
                AFLOWpi.prep.totree(dosExec,{ID:oneCalc})
                AFLOWpi.prep.totree(scfExec,{ID:oneCalc})
	else:
                AFLOWpi.prep.totree(pdosExec,{ID:oneCalc},symlink=True)
                AFLOWpi.prep.totree(dosExec,{ID:oneCalc},symlink=True)
                AFLOWpi.prep.totree(scfExec,{ID:oneCalc},symlink=True)


def transport_prep(oneCalc,ID):

	"""
        sets up the environment do do scf->nscf->projwfc to get overlap for transport calcs
        """
        AFLOWpi.scfuj.tb_prep(oneCalc,ID)
	AFLOWpi.scfuj.want_dos_prep(oneCalc,ID)
#	AFLOWpi.scfuj.want_bands_prep(oneCalc,ID)
#	AFLOWpi.scfuj.want_epsilon_prep(oneCalc,ID)

#############################################################################################
def want_bands_prep(oneCalc,ID):
	wantdir = AFLOWpi.prep._ConfigSectionMap('prep','want_dir')
        if os.path.isabs(wantdir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            wantdir =  os.path.join(configFileLocation, wantdir)



	#copy WanT bands.x
	try:
            wantBandsExec = os.path.join(wantdir,'bands.x')
	    if os.path.exists(wantBandsExec)==False:
		    print 'WanT bands executable not found. Check your config file to make sure want_dir path is correct and that bands.x is in that directory..Exiting'
		    logging.error('WanT bands executable not found. Check your config file to make sure want_dir path is correct and that bands.x is in that directory..Exiting')
		    raise SystemExit

            if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
                AFLOWpi.prep.totree(wantBandsExec,{ID:oneCalc},rename='want_bands.x')

            else:
                AFLOWpi.prep.totree(wantBandsExec,{ID:oneCalc},rename='want_bands.x',symlink=True)

	except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


#############################################################################################
def want_dos_prep(oneCalc,ID):
	wantdir = AFLOWpi.prep._ConfigSectionMap('prep','want_dir')
        if os.path.isabs(wantdir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            wantdir =  os.path.join(configFileLocation, wantdir)




        ###copy WanT dos.x
	try:
       	    wantDOSExec = os.path.join(wantdir,'dos.x')

            if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
                AFLOWpi.prep.totree(wantDOSExec,{ID:oneCalc},rename='want_dos.x')

            else:
                AFLOWpi.prep.totree(wantDOSExec,{ID:oneCalc},rename='want_dos.x',symlink=True)

	except Exception,e:
		AFLOWpi.run._fancy_error_log(e)


#############################################################################################
def want_epsilon_prep(oneCalc,ID,en_range=[0.5,10.0],ne=95):
	wantdir = AFLOWpi.prep._ConfigSectionMap('prep','want_dir')
        if os.path.isabs(wantdir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            wantdir =  os.path.join(configFileLocation, wantdir)



        ###copy epsilon.x
	try:
   	    wantEpsilonExec = os.path.join(wantdir,'epsilon.x')
            if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
                AFLOWpi.prep.totree(wantEpsilonExec,{ID:oneCalc},rename='want_epsilon.x')

            else:
                AFLOWpi.prep.totree(wantEpsilonExec,{ID:oneCalc},rename='want_epsilon.x',symlink=True)

	except Exception,e:
            AFLOWpi.run._fancy_error_log(e)


def want_eff_mass_prep(oneCalc,ID):
	wantdir = AFLOWpi.prep._ConfigSectionMap('prep','want_dir')
        if os.path.isabs(wantdir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            configFileLocation = os.path.dirname(configFileLocation)
            wantdir =  os.path.join(configFileLocation, wantdir)



        ###copy effmass.x
	try:
   	    want_eff_mass_exec = os.path.join(wantdir,'effmass.x')
            if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
                AFLOWpi.prep.totree(want_eff_mass_exec,{ID:oneCalc},rename='effmass.x')

            else:
                AFLOWpi.prep.totree(want_eff_mass_exec,{ID:oneCalc},rename='effmass.x',symlink=True)

	except Exception,e:
            AFLOWpi.run._fancy_error_log(e)





def WanT_dos(oneCalc,ID=None,eShift=5.0,temperature=None,energy_range=[-21.0,21.0],boltzmann=True,k_grid=None,pdos=False,num_e=4001,cond_bands=True,fermi_surface=False,compute_ham=False,proj_thr=0.95):
	'''
		Make input files for  WanT bands calculation

		Arguments:
        	 - calc_copy -- dictionary of dictionaries of calculations
	'''
	output_calc = {}
        temp_dir= AFLOWpi.prep._get_tempdir()
	calc_copy=copy.deepcopy(oneCalc)
	subdir = calc_copy['_AFLOWPI_FOLDER_']
	nspin = AFLOWpi.scfuj.chkSpinCalc(calc_copy,ID)

	compute_ham_str=""
	if compute_ham==False:
		compute_ham_str="!"

        try:
            prefix = AFLOWpi.retr._prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])
        except:
            prefix = oneCalc['_AFLOWPI_PREFIX_']
                
#        if temp_dir!='./':
#            subdir=temp_dir

        
	if cond_bands!=0 and type(cond_bands) == type(452):
		nbnds=int(AFLOWpi.prep._num_bands(oneCalc,mult=False))+cond_bands
	
	elif cond_bands==False:
		nbnds=int(AFLOWpi.prep._num_bands(oneCalc,mult=False))
	# if unoccupied_bands==True:
	# 	if nbnd==None:
	# 		nbnds=int(AFLOWpi.prep._num_bands(oneCalc))
	# 	else:
	# 		nbnds=nbnd
	if cond_bands==True: 
		nbnds=int(AFLOWpi.prep._num_bands(oneCalc,mult=False)*1.75)
		       

        ###############################
        ###############################
        #hardcode for now 
        #hardcode for now 


        try:
		if eShift==-1.0:
			nscf_ID=ID+'_nscf'
			Efermi = AFLOWpi.retr._getEfermi(oneCalc,nscf_ID,directID=True)
			eShift=5.5
#			eShift=float(Efermi)+15.0
		
        except:
            eShift+4.0

	if temperature==None:
		temp_temp=300.0
	else:
		temp_temp = temperature

	temp_ev=float(temp_temp)/11604.505

	if k_grid==None:
		inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
		kpt_str  = inputDict['K_POINTS']['__content__']    
		k_grid = [int(numpy.ceil(float(x)*5)) for x in kpt_str.split()[:3]]

        nk1   =  k_grid[0]
        nk2   =  k_grid[1]
        nk3   =  k_grid[2]
        emin  =  energy_range[0]
        emax  =  energy_range[1]
        nbnd  =  nbnds
	delta = (emax-emin)/(num_e)
        ne    =  num_e
	if boltzmann==True:
		boltz_flag='.TRUE.'
	else:
		boltz_flag='.FALSE.'
	if pdos==True:
		proj_flag='.TRUE.'
	else:
		proj_flag='.FALSE.'
	if fermi_surface==True:
		ferm_surf_str='do_fermisurf=.TRUE.'
	else:
		ferm_surf_str=''
        #hardcode for now 
        #hardcode for now 
        ###############################
        ###############################
	inputfile = """ &INPUT							  
prefix 		= '%s_TB' 		  	                        
postfix 	= \'_WanT\'		        
work_dir	= \'./\'		 
%sdatafile_dft	= \'./%s_TB.save/atomic_proj.dat\'	 
nk(1)           = %d
nk(2)           = %d
nk(3)           = %d
emin            = %f
emax            = %f
atmproj_nbnd    = %d
ne              = %d
delta           = %f
!smearing_type   = 'mv'
do_orthoovp	=.True.
atmproj_do_norm =.TRUE.
atmproj_thr     = %f 
atmproj_sh      = %f         
%s                
"""%(ID,compute_ham_str,ID,nk1,nk2,nk3,emin,emax,nbnd,ne,delta,proj_thr,eShift,ferm_surf_str)	
	
	if pdos==True:
		inputfile+='projdos         = .TRUE.\n'
	if temperature!=None:
		inputfile += 'temperature     = %f\n'%temp_ev
	else:
		inputfile += 'temperature     = 1.0D-12\n'

	if nspin == 1:
		if boltzmann!=True:
			inputfile += "fileout	      = '%s_WanT_dos.dat'\n"%(ID)
			inputfile +='/\n'
		else:
			inputfile += "fileout	      = '%s_WanT_dos_%sK.dat'\n"%(ID,temperature)

			inputfile += 'do_boltzmann_conductivity = %s\n'  %boltz_flag


			inputfile += "fileout2        = '%s_WanT_cond_%sK.dat'\n"%(ID,temperature)
			inputfile += "fileout3        = '%s_WanT_seebeck_%sK.dat'\n"%(ID,temperature)
			inputfile += "fileout4        = '%s_WanT_sigma_seebeck_%sK.dat'\n"%(ID,temperature)
			inputfile += "fileout5        = '%s_WanT_kappa_%sK.dat'\n"%(ID,temperature)		
			inputfile += "fileout6        = '%s_WanT_ZetaT_%sK.dat'\n/ \n"%(ID,temperature)
		# Insert k-path 
#		inputfile += kpathStr

		

		calc_label = ID + "_WanT_dos"
		
	        a = calc_label+'.in'
	        new_inputfile = open(os.path.join('./',a),'w')
	        new_inputfile.write(inputfile)
	        new_inputfile.close()

	        output_calc = calc_copy
	        output_calc['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]

		single_output_calc = {calc_label:output_calc}

		return single_output_calc

	else:

		inputfile1 = inputfile

                inputfile  = re.sub('_WanT','_WanT_up',inputfile)
                inputfile1 = re.sub('_WanT','_WanT_dn',inputfile1)

		#SPIN UP  SPIN UP  SPIN UP
		if boltzmann!=True:
			inputfile += "fileout	      = '%s_WanT_dos_up.dat'\n"%(ID)
			inputfile += '\nspin_component	= "up"  \n'
			inputfile +='/\n'
		else:
			inputfile += "fileout	      = '%s_WanT_dos_up_%sK.dat'\n"%(ID,temperature)
			inputfile +=  """fileout = '%s_WanT_dos_up_%sK.dat'   \n"""  %(ID,temperature)		      
			inputfile += 'spin_component	= "up"  \n'
			inputfile += 'do_boltzmann_conductivity = %s\n'  %boltz_flag
#			inputfile += 'temperature     = %f\n'%temp_ev
			#Create two input files for spin-up and spin-down components.
			inputfile += "fileout2  = '%s_WanT_cond_up_%sK.dat'\n"%(ID,temperature)
			inputfile += "fileout3  = '%s_WanT_seebeck_up_%sK.dat'\n"%(ID,temperature)
			inputfile += "fileout4  = '%s_WanT_sigma_seebeck_up_%sK.dat'\n"%(ID,temperature)
			inputfile += "fileout5  = '%s_WanT_kappa_up_%sK.dat'\n"%(ID,temperature)		
			inputfile += "fileout6  = '%s_WanT_ZetaT_up_%sK.dat'\n / \n"%(ID,temperature)


		#SPIN DOWN  SPIN DOWN  SPIN DOWN
		if boltzmann!=True:
			inputfile1 += "fileout	      = '%s_WanT_dos_down.dat'\n"%(ID)
			inputfile1 += '\nspin_component	= "down"  \n'
			inputfile1 +='/\n'
		else:
			inputfile1 += "fileout	      = '%s_WanT_dos_down_%sK.dat'\n"%(ID,temperature)
			inputfile1 += 'do_boltzmann_conductivity = %s\n'  %boltz_flag
#			inputfile1 += 'temperature     = %f\n'%temp_ev
			inputfile1 +=  """fileout = '%s_WanT_dos_down_%sK.dat'   \nspin_component	= "down"  \n"""  %(ID,temperature)
			inputfile1 += "fileout2  = '%s_WanT_cond_down_%sK.dat'\n"%(ID,temperature)
			inputfile1 += "fileout3  = '%s_WanT_seebeck_down_%sK.dat'\n"%(ID,temperature)
			inputfile1 += "fileout4  = '%s_WanT_sigma_seebeck_down_%sK.dat'\n"%(ID,temperature)
			inputfile1 += "fileout5  = '%s_WanT_kappa_down_%sK.dat'\n"%(ID,temperature)		
			inputfile1 += "fileout6  = '%s_WanT_ZetaT_down_%sK.dat'\n / \n"%(ID,temperature)



                calc_label_up = ID + "_WanT_dos_up" 
                a = calc_label_up+'.in'
                new_inputfile = open(os.path.join(subdir,a),'w')
                new_inputfile.write(inputfile)
                new_inputfile.close()
                output_calc_up= calc_copy
                output_calc_up['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]


                calc_label_down = ID + "_WanT_dos_down" 
                a = calc_label_down+'.in'
                new_inputfile = open(os.path.join(subdir,a),'w')
                new_inputfile.write(inputfile1)
	        new_inputfile.close()
                output_calc_down = calc_copy
                output_calc_down['_AFLOWPI_INPUT_'] = inputfile1
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]


		output_calc = {calc_label_up:output_calc_up, calc_label_down:output_calc_down}		
		return output_calc


###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################

def WanT_epsilon(oneCalc,ID=None,eShift=5.0,temperature=300.0,energy_range=[0.01,8.01],ne=160,k_grid=None,compute_ham=False,proj_thr=0.95,proj_nbnd=None):
	'''
		Make input files for  WanT bands calculation

		Arguments:
        	 - calc_copy -- dictionary of dictionaries of calculations
	'''

	if energy_range[0]<=0.0:
		energy_range[0]=0.01

	compute_ham_str=""
	if compute_ham==False:
		compute_ham_str="!"

	output_calc = {}
        temp_dir= AFLOWpi.prep._get_tempdir()
	calc_copy=copy.deepcopy(oneCalc)
	subdir = calc_copy['_AFLOWPI_FOLDER_']
	nspin = AFLOWpi.scfuj.chkSpinCalc(calc_copy,ID)


        try:
            prefix = AFLOWpi.retr._prefixFromInput(oneCalc['_AFLOWPI_INPUT_'])
        except:
            prefix = oneCalc['_AFLOWPI_PREFIX_']
                
#        if temp_dir!='./':
#            subdir=temp_dir
	if proj_nbnd==None:
		nbnds=int(AFLOWpi.prep._num_bands(oneCalc,mult=False)*1.75)
	else:
		nbnds=int(AFLOWpi.prep._num_bands(oneCalc,mult=False))+proj_nbnd

        try:
            nscf_ID=ID+'_nscf'
            Efermi = AFLOWpi.retr._getEfermi(oneCalc,nscf_ID,directID=True)
            eShift=float(Efermi)+10.0
        except:
            eShift=5.0

#        eShift=5.0

        ###############################
        ###############################
        #hardcode for now 
        #hardcode for now 
	temp_ev=float(temperature)/11604.505

	if k_grid==None:
		inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
		kpt_str  = inputDict['K_POINTS']['__content__']    
		k_grid = [int(numpy.ceil(float(x)*5)) for x in kpt_str.split()[:3]]

        nk1   =  k_grid[0]
        nk2   =  k_grid[1]
        nk3   =  k_grid[2]
        emin  =  energy_range[0]
        emax  =  energy_range[1]
        nbnd  =  nbnds
	delta = (emax-emin)/(ne*2.0)
        #hardcode for now 
        #hardcode for now 
        ###############################
        ###############################
	inputfile = """ &INPUT							  
prefix		= '%s_TB' 		  	                        
postfix 	= \'_WanT\'		        
work_dir	= \'./\'		 
%sdatafile_dft	= \'./%s_TB.save/atomic_proj.dat\'	 
nk(1)           = %d
nk(2)           = %d
nk(3)           = %d
emin            = %f
emax            = %f
temperature     = %f
atmproj_nbnd    = %d
ne              = %d
delta           = %f
!smearing_type   = 'mv'
do_orthoovp	= .TRUE.
atmproj_do_norm =.TRUE.
atmproj_sh      = %f                         
atmproj_thr     = %f
"""%(ID,compute_ham_str,ID,nk1,nk2,nk3,emin,emax,temp_ev,nbnd,ne,delta,eShift,proj_thr)

	if nspin == 1:
		inputfile += "fileout	      = '%s_WanT_epsilon_imag.dat'\n"%(ID)
		inputfile += "fileout2	      = '%s_WanT_epsilon_real.dat'\n"%(ID)
		inputfile += "fileout3	      = '%s_WanT_epsilon_eels.dat'\n/ \n"%(ID)
		# Insert k-path 
#		inputfile += kpathStr

		calc_label = ID + "_WanT_epsilon"
		
	        a = calc_label+'.in'
	        new_inputfile = open(os.path.join('./',a),'w')
	        new_inputfile.write(inputfile)
	        new_inputfile.close()

	        output_calc = calc_copy
	        output_calc['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]

		single_output_calc = {calc_label:output_calc}

		return single_output_calc

	else:
		inputfile1 = inputfile
                inputfile  = re.sub('_WanT','_WanT_up',inputfile)
                inputfile1 = re.sub('_WanT','_WanT_dn',inputfile1)
		#Create two input files for spin-up and spin-down components.

		inputfile += """fileout = '%s_WanT_epsilon_up_imag.dat' \nspin_component	= "up"\n"""  %(ID,)
		inputfile += "fileout2  = '%s_WanT_epsilon_up_real.dat'\n"%(ID)
		inputfile += "fileout3  = '%s_WanT_epsilon_up_eels.dat'\n/ \n"%(ID)


		inputfile1 += """fileout = '%s_WanT_epsilon_down_imag.dat' \nspin_component	= "down"\n"""  %(ID,)
		inputfile1 += "fileout2  = '%s_WanT_epsilon_down_real.dat'\n"%(ID)
		inputfile1 += "fileout3  = '%s_WanT_epsilon_down_eels.dat'\n/ \n"%(ID)



		#Insert k-path
#		inputfile += kpathStr
#		inputfile1 += kpathStr

                calc_label_up = ID + "_WanT_epsilon_up" 
                a = calc_label_up+'.in'
                new_inputfile = open(os.path.join(subdir,a),'w')
                new_inputfile.write(inputfile)
                new_inputfile.close()
                output_calc_up= calc_copy
                output_calc_up['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]


                calc_label_down = ID + "_WanT_epsilon_down" 
                a = calc_label_down+'.in'
                new_inputfile = open(os.path.join(subdir,a),'w')
                new_inputfile.write(inputfile1)
	        new_inputfile.close()
                output_calc_down = calc_copy
                output_calc_down['_AFLOWPI_INPUT_'] = inputfile
                try:
                    output_calc['prev'].append(ID)
                except:
                    output_calc['prev']=[ID]


		output_calc = {calc_label_up:output_calc_up, calc_label_down:output_calc_down}		
		return output_calc




def run_transport(__submitNodeName__,oneCalc,ID,run_scf=True,run_transport_prep=True,run_bands=False,epsilon=False,temperature=300,en_range=[0.05,10.0],ne=1000,compute_ham=False,proj_thr=0.95,proj_sh=5.5,proj_nbnd=True):

	execPrefix = ''
	execPostfix = ''
        oneCalcID = ID

	compute_ham_str=""
	if compute_ham==False:
		compute_ham_str="!"

        def abortIFRuntimeError(subdir, ID):
            outfile = file(os.path.join(subdir, "%s.out"%ID)).read()
            errorList = re.findall(r'from (.*) : error #.*\n',outfile)
            if len(errorList) > 0 and 'zmat_hdiag' not in errorList:        
                logging.error("Error in %s.out -- ABORTING ACBN0 LOOP"%ID)
                print "Error in %s.out -- ABORTING ACBN0 LOOP"%ID                    
                raise SystemExit






        if '__runList__' not in oneCalc.keys():
            oneCalc['__runList__']=[]
#            if run_scf==False:
#                oneCalc['__runList__']=['scf']
	if len(oneCalc['__runList__'])==0:
		if run_transport_prep==False:
			oneCalc['__runList__']=['scf','nscf','pdos']
		if run_bands==False:
			oneCalc['__runList__'].append('bands')
	    
        
        config=None
	if config!=None:
		AFLOWpi.prep._forceGlobalConfigFile(config)
		logging.debug('forced config %s' % config)
	else:
		try:
			config = AFLOWpi.prep._getConfigFile()
			AFLOWpi.prep._forceGlobalConfigFile(config)
		except Exception,e:
			AFLOWpi.run._fancy_error_log(e)


	if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
            execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")

	else:
            execPrefix=''


	if AFLOWpi.prep._ConfigSectionMap("run","exec_postfix") != '':
		execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
	else:
		execPostfix=''


	if AFLOWpi.prep._ConfigSectionMap('run','engine') == '':
		engine = AFLOWpi.prep._ConfigSectionMap('run','engine')
	else:
		engine = 'espresso'


        subdir = oneCalc['_AFLOWPI_FOLDER_']
	oneCalc['_AFLOWPI_CONFIG_']=config

        if 'scf' not in oneCalc['__runList__']:

            try:
                npool=AFLOWpi.retr._get_pool_num(oneCalc,ID)        

                if npool!=1:
                    if len(re.findall(r'npool[s]*\s*(?:\d*)',execPostfix))!=0:
                        execPostfixPrime=re.sub(r'npool[s]*\s*(?:\d*)','npool %s'%npool,execPostfix)
              #          logging.debug(execPostfixPrime)

            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)


##################################################################################################################
            AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)


            oneCalc['__runList__'].append('scf')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,1.5,unoccupied_states=True)	

        else:
            '''if we are restarting from a job killed from going walltime 
            try to load ID_nscf and if we can't then just make a new one'''
            try:
                nscf_ID='%s_nscf' % ID
                nscf_calc = AFLOWpi.prep._loadOneCalc(oneCalc['_AFLOWPI_FOLDER_'],nscf_ID)                
                '''we have to make sure nscf step has the correct walltime and start time if it's a restart'''
                nscf_calc['__walltime_dict__']=oneCalc['__walltime_dict__']
            except Exception,e:
                try:
                    nscf_calc,nscf_ID= AFLOWpi.scfuj.nscf_nosym_noinv(oneCalc,ID,1.5,unoccupied_states=True)	

                except Exception,e:
                    AFLOWpi.run._fancy_error_log(e)


##################################################################################################################
        if 'nscf' not in oneCalc['__runList__']:

            try:
                npool=AFLOWpi.retr._get_pool_num(nscf_calc,nscf_ID)        

                if npool!=1:
                    if len(re.findall(r'npool\s*(?:\d+)',execPostfix))!=0:
                        execPostfixPrime=re.sub(r'npool\s*(?:\d+)','npool %s'%npool,execPostfix)
                        logging.debug(execPostfixPrime)

            except Exception,e:
                AFLOWpi.run._fancy_error_log(e)



	    
            AFLOWpi.run._oneRun(__submitNodeName__,nscf_calc,nscf_ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='scf',executable=None)
	    
            abortIFRuntimeError(subdir, nscf_ID)
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            oneCalc['__runList__'].append('nscf')
	    nscf_calc
##################################################################################################################
        pdos_calc,pdos_ID = AFLOWpi.scfuj.projwfc(oneCalc,ID)


        if not re.match('northo',execPostfix) or not re.match('no',execPostfix):
            execPostfix+=' -northo 1'


        if 'pdos' not in oneCalc['__runList__']:
            pdosPath = os.path.join(AFLOWpi.prep._ConfigSectionMap('prep','engine_dir'),'projwfc.x')

            AFLOWpi.run._oneRun(__submitNodeName__,pdos_calc,pdos_ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='custom',executable='projwfc.x',execPath=pdosPath)

	    AFLOWpi.prep._form_TB_dir(oneCalc,ID)
#############
            oneCalc['__runList__'].append('pdos')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)
            abortIFRuntimeError(subdir, pdos_ID)

	#when we skip forming TB Ham
	AFLOWpi.prep._form_TB_dir(oneCalc,ID,from_ls=False)

        if 'bands' not in oneCalc['__runList__']:
            want_dict = AFLOWpi.scfuj.WanT_bands(oneCalc,ID)


            for want_ID,want_calc in want_dict.iteritems():
                AFLOWpi.run._oneRun(__submitNodeName__,want_calc,want_ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='custom',execPath='./want_bands.x' )

            

	    AFLOWpi.prep._clean_want_bands(oneCalc,ID)

            for want_ID,want_calc in want_dict.iteritems():
                abortIFRuntimeError(subdir, want_ID)


            oneCalc['__runList__'].append('bands')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)

            '''
            Will need to be filled in for the executable name with whatever wanT executable is called 
            and that executable needs to be moved to the calculation directory tree before this is called
            '''


        if 'want_dos' not in oneCalc['__runList__']:
            want_dos_calc = AFLOWpi.scfuj.WanT_dos(oneCalc,ID,temperature=temperature,num_e=ne,energy_range=en_range,compute_ham=compute_ham,proj_thr=proj_thr,cond_bands=proj_nbnd)

            for want_dos_ID,want_dos in want_dos_calc.iteritems():
                AFLOWpi.run._oneRun(__submitNodeName__,want_dos,want_dos_ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='custom',execPath='./want_dos.x' )

            for want_dos_ID,want_dos in want_dos_calc.iteritems():
                abortIFRuntimeError(subdir, want_dos_ID)


            oneCalc['__runList__'].append('want_dos')
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)



            '''
            Will need to be filled in for the executable name with whatever wanT executable is called 
            and that executable needs to be moved to the calculation directory tree before this is called
            '''

            #save and exit
            oneCalc['__runList__']=[]
            AFLOWpi.prep._saveOneCalc(oneCalc,ID)

            
            return oneCalc,ID


def _run_want_epsilon(__submitNodeName__,oneCalc,ID,en_range=[0.1,8.0],ne=79,compute_ham=False,proj_thr=0.95,proj_sh=5.5,proj_nbnd=None):

	    def abortIFRuntimeError(subdir, ID):
		    outfile = file(os.path.join(subdir, "%s.out"%ID)).read()
		    errorList = re.findall(r'from (.*) : error #.*\n',outfile)
		    if len(errorList) > 0 and 'zmat_hdiag' not in errorList:        
			    logging.error("Error in %s.out -- ABORTING ACBN0 LOOP"%ID)
			    print "Error in %s.out -- ABORTING ACBN0 LOOP"%ID                    
			    raise SystemExit



		    
            want_epsilon_calc = AFLOWpi.scfuj.WanT_epsilon(oneCalc,ID,energy_range=en_range,ne=ne,compute_ham=compute_ham,proj_thr=proj_thr)

	    if AFLOWpi.prep._ConfigSectionMap("run","exec_prefix") != '':
		    execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix")
	    else:
		    execPrefix=''

	    if AFLOWpi.prep._ConfigSectionMap("run","exec_postfix") != '':
		    execPostfix = AFLOWpi.prep._ConfigSectionMap("run","exec_postfix")
	    else:
		    execPostfix=''

            for want_epsilon_ID,want_epsilon in want_epsilon_calc.iteritems():
                AFLOWpi.run._oneRun(__submitNodeName__,want_epsilon,want_epsilon_ID,execPrefix=execPrefix,execPostfix=execPostfix,engine='espresso',calcType='custom',execPath='./want_epsilon.x' )

            for want_epsilon_ID,want_epsilon in want_epsilon_calc.iteritems():
                abortIFRuntimeError(oneCalc["_AFLOWPI_FOLDER_"], want_epsilon_ID)

	    return oneCalc,ID


	    
