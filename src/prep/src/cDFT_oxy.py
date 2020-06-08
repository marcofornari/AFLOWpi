try:
    from .prep import newstepWrapper,_check_lock
except: pass
import AFLOWpi
import os 
import re
import numpy as np

#@newstepWrapper(_check_lock)
def _prep_cDFT_oxy(oneCalc,ID,oxy_dict,initial_step=1.0):


    in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    
    # get list of chemical species in the input
    spec_list = [x.split()[0] for x in in_dict['ATOMIC_SPECIES']['__content__'].split('\n') if len(x)!=0]

    temp_init_alpha  = []
    temp_step_size   = []


    # check if the necessary input parameters are in the pwscf input..if not add them
    for spec in range(len(spec_list)):
        if spec_list[spec] in list(oxy_dict.keys()):

            if not 'lda_plus_u' in list(in_dict['&system'].keys()):
                in_dict['&system']['lda_plus_u']='.true.'
            if not 'hubbard_u(%s)'%(spec+1) in list(in_dict['&system'].keys()):
                in_dict['&system']['hubbard_u(%s)'%(spec+1)]='1.e-8'
            if not 'hubbard_alpha(%s)'%(spec+1) in list(in_dict['&system'].keys()):
                in_dict['&system']['hubbard_alpha(%s)'%(spec+1)]=str(0.0)

            temp_init_alpha.append(initial_step)


    # add the initial value for the independent variables
    # we'll be changing to minimize the "cost" function
    if 'cDFT_grad_hist' not in list(oneCalc.keys()):
        oneCalc['cDFT_grad_hist']=[]
        oneCalc['cDFT_alpha_hist']=[temp_init_alpha]

        new_in_str = AFLOWpi.retr._joinInput(in_dict)
        with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'w') as new_inputfile:
            new_inputfile.write(new_in_str)

        oneCalc['_AFLOWPI_INPUT_'] = new_in_str
        AFLOWpi.prep._saveOneCalc(oneCalc,ID)



    return oneCalc,ID

        



def _check_cDFT_conv(oneCalc,ID,oxy_dict,conv_thr):
    

    
    conv_bool = True

    # get list of chemical species in the input
    in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    spec_list = [x.split()[0] for x in in_dict['ATOMIC_SPECIES']['__content__'].split('\n') if len(x)!=0]

    # grab the output string
    ofs = AFLOWpi.retr._getOutputString(oneCalc,ID)

    # get atomic position labels from ATOMIC_POSITIONS
    atoms = AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])

    val_dict = {}
    #add dummy list of occ value and count
    for spec in range(len(spec_list)):
        atom_label = spec_list[spec]
        if atom_label in list(oxy_dict.keys()):
            val_dict[spec_list[spec]] = [0.0,0]

    # use a regular expression to pull the trace of the occupations from the output data
    for atom_ind in range(len(atoms)):
        atom_label = atoms[atom_ind]

        if atom_label in list(oxy_dict.keys()):

            ais = str((atom_ind+1))
            rgx_str = r"atom\s*"+ais+"\s*Tr\[ns\(na\)\s*\]\s*=\s*([-.\d]*)"
            try:
                sol =  float(re.findall(rgx_str,ofs)[-1])
            except:
                try:
                    rgx_str = r"atom\s*"+ais+"\s*Tr\[ns\(na\)\s*\].*=\s*([-.\d]*)\s+([-.\d]*)\s+([-.\d]*)"
                    sol =  float(re.findall(rgx_str,ofs)[-1][-1])

                except: pass

            # add the entry for the trace of the occupations to a total for that species
            try:
                temp_val = val_dict[atom_label][0]
                temp_cnt = val_dict[atom_label][1]

                val_dict[atom_label][0] = temp_val + sol
                val_dict[atom_label][1] = temp_cnt + 1          
            except Exception as e: print(e) 

    # check to see if the occupations match the target value. 
    # if there is more than one atom for a given species...
    # take the average occupation for that species
    for spec,dat in list(val_dict.items()):
        if abs(oxy_dict[spec]-dat[0]/dat[1])>=conv_thr:
            conv_bool=False


    return conv_bool



def _cDFT_newstep(oneCalc,ID,oxy_dict):

    # get list of chemical species in the input    
    in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    spec_list = [x.split()[0] for x in in_dict['ATOMIC_SPECIES']['__content__'].split('\n') if len(x)!=0]

    # grab the output string
    ofs = AFLOWpi.retr._getOutputString(oneCalc,ID)

    # get atomic position labels from ATOMIC_POSITIONS
    atoms = AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])

    val_dict = {}
    #add dummy list of occ value and count
    for spec in range(len(spec_list)):
        atom_label = spec_list[spec]
        if atom_label in list(oxy_dict.keys()):
            val_dict[spec_list[spec]] = [0.0,0]

    # use a regular expression to pull the trace of the occupations from the output data
    for atom_ind in range(len(atoms)):
        atom_label = atoms[atom_ind]

        if atom_label in list(oxy_dict.keys()):

            ais = str((atom_ind+1))
            rgx_str = r"atom\s*"+ais+"\s*Tr\[ns\(na\)\s*\]\s*=\s*([-.\d]*)"
            try:
                sol =  float(re.findall(rgx_str,ofs)[-1])
            except:
                try:
                    rgx_str = r"atom\s*"+ais+"\s*Tr\[ns\(na\)\s*\].*=\s*([-.\d]*)\s+([-.\d]*)\s+([-.\d]*)"
                    sol =  float(re.findall(rgx_str,ofs)[-1][-1])

                except: pass

            # add the entry for the trace of the occupations to a total for that species
            try:
                temp_val = val_dict[atom_label][0]
                temp_cnt = val_dict[atom_label][1]

                val_dict[atom_label][0] = temp_val + sol
                val_dict[atom_label][1] = temp_cnt + 1          
            except Exception as e: print(e) 


    # check to see if the occupations match the target value. 
    # if there is more than one atom for a given species...
    # take the average occupation for that species
    temp_grad = []
    for spec,dat in list(val_dict.items()):
        temp_grad.append(oxy_dict[spec]-dat[0]/dat[1])
    oneCalc['cDFT_grad_hist'].append(temp_grad)

    # if it's the first iteration just take gradient of the cost  
    # function and the inital step size to calculate the next step
    alpha_n  =  np.array(oneCalc['cDFT_alpha_hist'][-1])
    
    if len(oneCalc['cDFT_grad_hist'])<2:
        new_alpha = alpha_n - alpha_n*np.array(temp_grad)
    # if not the first iteration use the Barzilai-Borwein method
    # to calculate the step size and values for the next step
    else:
        alpha_mo =  np.array(oneCalc['cDFT_alpha_hist'][-2])

        grad_n   =  np.array(oneCalc['cDFT_grad_hist'][-1])
        grad_mo  =  np.array(oneCalc['cDFT_grad_hist'][-2])

        nss = np.dot((alpha_n-alpha_mo).T,(grad_n-grad_mo))/np.sum((grad_n-grad_mo)**2)
        new_alpha = alpha_n - nss*grad_n

    # add the values of the parameters for the next step to the step history
    oneCalc['cDFT_alpha_hist'].append(new_alpha.tolist())

    # change the values of the parameters in the input for the next step
    for i in range(len(new_alpha)):
        in_dict['&system']['hubbard_alpha(%s)'%(i+1)]=str(new_alpha[i])

    new_in_str = AFLOWpi.retr._joinInput(in_dict)

    # save the input file for the next step in the oneCalc
    # dictionary and write the new input to disk
    oneCalc['_AFLOWPI_INPUT_'] = new_in_str

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'w') as new_inputfile:
        new_inputfile.write(new_in_str)


    # reset runtime counters and locks 
    oneCalc['__status__']['Complete'] = False
    oneCalc['__execCounter__'] = 0

    # save oneCalc to disk
    AFLOWpi.prep._saveOneCalc(oneCalc,ID)


    return oneCalc,ID

def _cDFT_cleanup(oneCalc,ID):
    # cleanup any unneeded flags in oneCalc or do any other cleanup needed
    try:
        del oneCalc['cDFT_step_sizes']
    except: pass

    return oneCalc,ID
