from prep import newstepWrapper,_check_lock
import AFLOWpi
import os 
import re
import numpy as np

#@newstepWrapper(_check_lock)
def _prep_cDFT_oxy(oneCalc,ID,oxy_dict,initial_step=1.0):


    in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    
    spec_list = [x.split()[0] for x in in_dict['ATOMIC_SPECIES']['__content__'].split('\n') if len(x)!=0]

    temp_init_alpha  = []
    temp_step_size   = []

    for spec in range(len(spec_list)):
        if spec_list[spec] in oxy_dict.keys():

            if not 'lda_plus_u' in in_dict['&system'].keys():
                in_dict['&system']['lda_plus_u']='.true.'
            if not 'hubbard_u(%s)'%(spec+1) in in_dict['&system'].keys():
                in_dict['&system']['hubbard_u(%s)'%(spec+1)]='1.e-8'
            if not 'hubbard_alpha(%s)'%(spec+1) in in_dict['&system'].keys():
                in_dict['&system']['hubbard_alpha(%s)'%(spec+1)]=str(0.0)

            temp_init_alpha.append(initial_step)


    if 'cDFT_grad_hist' not in oneCalc.keys():
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

    in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    spec_list = [x.split()[0] for x in in_dict['ATOMIC_SPECIES']['__content__'].split('\n') if len(x)!=0]


    ofs = AFLOWpi.retr._getOutputString(oneCalc,ID)

    atoms = AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])

    val_dict = {}
    #add dummy list of occ value and count
    for spec in range(len(spec_list)):
        atom_label = spec_list[spec]
        if atom_label in oxy_dict.keys():
            val_dict[spec_list[spec]] = [0.0,0]

    for atom_ind in range(len(atoms)):
        atom_label = atoms[atom_ind]

        if atom_label in oxy_dict.keys():

            ais = str((atom_ind+1))
            rgx_str = r"atom\s*"+ais+"\s*Tr\[ns\(na\)\s*\]\s*=\s*([-.\d]*)"
            try:
                sol =  float(re.findall(rgx_str,ofs)[-1])
            except:
                try:
                    rgx_str = r"atom\s*"+ais+"\s*Tr\[ns\(na\)\s*\].*=\s*([-.\d]*)\s+([-.\d]*)\s+([-.\d]*)"
                    sol =  float(re.findall(rgx_str,ofs)[-1][-1])

                except: pass

            try:
                temp_val = val_dict[atom_label][0]
                temp_cnt = val_dict[atom_label][1]

                val_dict[atom_label][0] = temp_val + sol
                val_dict[atom_label][1] = temp_cnt + 1          
            except Exception,e: print e 

    
    for spec,dat in val_dict.iteritems():
        if abs(oxy_dict[spec]-dat[0]/dat[1])>=conv_thr:
            conv_bool=False


    return conv_bool



def _cDFT_newstep(oneCalc,ID,oxy_dict):
    
    in_dict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    spec_list = [x.split()[0] for x in in_dict['ATOMIC_SPECIES']['__content__'].split('\n') if len(x)!=0]


    ofs = AFLOWpi.retr._getOutputString(oneCalc,ID)

    atoms = AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])

    val_dict = {}
    #add dummy list of occ value and count
    for spec in range(len(spec_list)):
        atom_label = spec_list[spec]
        if atom_label in oxy_dict.keys():
            val_dict[spec_list[spec]] = [0.0,0]

    for atom_ind in range(len(atoms)):
        atom_label = atoms[atom_ind]

        if atom_label in oxy_dict.keys():

            ais = str((atom_ind+1))
            rgx_str = r"atom\s*"+ais+"\s*Tr\[ns\(na\)\s*\]\s*=\s*([-.\d]*)"
            try:
                sol =  float(re.findall(rgx_str,ofs)[-1])
            except:
                try:
                    rgx_str = r"atom\s*"+ais+"\s*Tr\[ns\(na\)\s*\].*=\s*([-.\d]*)\s+([-.\d]*)\s+([-.\d]*)"
                    sol =  float(re.findall(rgx_str,ofs)[-1][-1])

                except: pass

            try:
                temp_val = val_dict[atom_label][0]
                temp_cnt = val_dict[atom_label][1]

                val_dict[atom_label][0] = temp_val + sol
                val_dict[atom_label][1] = temp_cnt + 1          
            except Exception,e: print e 


    temp_grad = []
    for spec,dat in val_dict.iteritems():
#        temp_grad.append(abs(oxy_dict[spec]-dat[0]/dat[1]))
        temp_grad.append(oxy_dict[spec]-dat[0]/dat[1])
    oneCalc['cDFT_grad_hist'].append(temp_grad)

#    new_alpha = 


    alpha_n  =  np.array(oneCalc['cDFT_alpha_hist'][-1])


    if len(oneCalc['cDFT_grad_hist'])<2:
        new_alpha = alpha_n - alpha_n*np.array(temp_grad)
    else:
        alpha_mo =  np.array(oneCalc['cDFT_alpha_hist'][-2])

        grad_n   =  np.array(oneCalc['cDFT_grad_hist'][-1])
        grad_mo  =  np.array(oneCalc['cDFT_grad_hist'][-2])

        nss = np.dot((alpha_n-alpha_mo).T,(grad_n-grad_mo))/np.sum((grad_n-grad_mo)**2)
        print nss
    #    print nss
    #    print alpha_n
        new_alpha = alpha_n - nss*grad_n

    oneCalc['cDFT_alpha_hist'].append(new_alpha.tolist())

    for i in range(len(new_alpha)):
        in_dict['&system']['hubbard_alpha(%s)'%(i+1)]=str(new_alpha[i])

    new_in_str = AFLOWpi.retr._joinInput(in_dict)

    oneCalc['_AFLOWPI_INPUT_'] = new_in_str

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'.in'),'w') as new_inputfile:
        new_inputfile.write(new_in_str)



    oneCalc['__status__']['Complete'] = False
    oneCalc['__execCounter__'] = 0


    AFLOWpi.prep._saveOneCalc(oneCalc,ID)


    return oneCalc,ID

def _cDFT_cleanup(oneCalc,ID):

    try:
        del oneCalc['cDFT_step_sizes']
    except: pass

    return oneCalc,ID
