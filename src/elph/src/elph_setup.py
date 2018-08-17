import AFLOWpi
import os
import shutil
import glob

def prep_elph(oneCalc,ID):
    phonon_ind = AFLOWpi.prep._get_step(oneCalc,ID,step_type='phonon',last=True)

    phonon_ID = '%s_%02d'%(ID.split('_')[0],phonon_ind[-1])
    H0_file = phonon_ID+'_FD_PHONON/'+phonon_ID+'_FD_PHONON__0001/displaced.0.0.0_01.in'
    H0_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],H0_file)

    elph_folder = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_FD_ELPH')
    
    if not os.path.exists(elph_folder):
        os.mkdir(elph_folder)

    new_H0_file = os.path.join(elph_folder,'H0.in')

    with open(H0_file,'r') as ifo:
        H0_file_str = ifo.read()

    spli = AFLOWpi.retr._splitInput(H0_file_str)

    pos  = AFLOWpi.retr._getPositions(H0_file_str)
    bohr2ang = 0.529177249
    cell = AFLOWpi.retr.getCellMatrixFromInput(H0_file_str)*bohr2ang
    pos  = AFLOWpi.retr._convertFractional(pos,cell)


    labs = AFLOWpi.retr._getPosLabels(H0_file_str)
    spli['ATOMIC_POSITIONS']['__content__'] = AFLOWpi.retr._joinMatrixLabels(labs,pos)
    spli['ATOMIC_POSITIONS']['__modifier__'] = '{crystal}'
    new_H0_file_str = AFLOWpi.retr._joinInput(spli)


    with open(new_H0_file,'w') as ofo:
        ofo.write(new_H0_file_str)



def _gen_elph_in(oneCalc,ID):

    phonon_ind = AFLOWpi.prep._get_step(oneCalc,ID,step_type='phonon',last=True)
    phonon_ID = '%s_%02d'%(ID.split('_')[0],phonon_ind[-1])
    phonon_fd_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],phonon_ID+'_fd.in')    

    with open(phonon_fd_file) as ifo:
        ifs = ifo.read()

    split_FD_in = AFLOWpi.retr._splitInput(ifs)
    nq1 = int(split_FD_in['&inputfd']['nrx1'])
    nq2 = int(split_FD_in['&inputfd']['nrx2'])
    nq3 = int(split_FD_in['&inputfd']['nrx3'])

    phonon_matdyn_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],phonon_ID+'_matdyn_phDOS.in')
    elph_matdyn_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_matdyn_phDOS.in')

    with open(phonon_matdyn_file) as ifo:
        ifs = ifo.read()
    
    split_matdyn_file = AFLOWpi.retr._splitInput(ifs)
    
    kpts_str = AFLOWpi.prep._gen_nosym_kgrid(nq1,nq2,nq3,as_string=True,scale=1.0)          


    del split_matdyn_file['ATOMIC_SPECIES']
    del split_matdyn_file['K_POINTS']

    try:        
        del split_matdyn_file['&input']['l1']
        del split_matdyn_file['&input']['l2']
        del split_matdyn_file['&input']['l3']
    except: pass

    split_matdyn_file['&input']['fldos'] = '%s.phdos'%ID
    split_matdyn_file['&input']['flfrq'] = '%s.phDOS'%ID
    split_matdyn_file['&input']['fleig'] = '%s.phDOS.eig'%ID
    split_matdyn_file['&input']['fldyn'] = '%s.phDOS.dyn'%ID
    split_matdyn_file['&input']['flvec'] = '%s.phDOS.modes'%ID
    

    new_matdyn_file = AFLOWpi.retr._joinInput(split_matdyn_file)
    new_matdyn_file += str(kpts_str)+'\n'
    
    with open(elph_matdyn_file,'w') as ofo:
        ofo.write(new_matdyn_file)


    exec_prefix_serial = AFLOWpi.prep._ConfigSectionMap("run","exec_prefix_serial")

    AFLOWpi.run._oneRun('',oneCalc,'%s_matdyn_phDOS'%ID,execPrefix=exec_prefix_serial,execPostfix='',engine='espresso',calcType='custom',execPath='./matdyn.x' )

    elph_folder = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_FD_ELPH')

    H0_file = os.path.join(elph_folder,'H0.in')

    nat = int(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])['&system']['nat'])

    commensurate_k_modes_fn = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.phDOS.modes'%ID)

    modes = AFLOWpi.elph.read_modes ( commensurate_k_modes_fn, nat )

    #get original cell vectors scaled by alat
    alat,orig_cell = AFLOWpi.retr._getCellParams(oneCalc,ID)
    orig_cell/=alat
    orig_cell = orig_cell.A


    AFLOWpi.elph.write_scf_files ( elph_folder+'/', H0_file, nat, nq1, nq2, nq3, modes, orig_cell )

    elph_input_files = sorted(glob.glob(elph_folder+'/*.in'))


    return elph_input_files


def _elph_pp(oneCalc,ID):
    pass

