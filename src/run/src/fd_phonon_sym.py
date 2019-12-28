import AFLOWpi
import numpy as np
import glob
import os
import re 


#########################################################################################
#########################################################################################
#########################################################################################

def _check(check_dir,force_output_dir):


    pass_list_file= []
    pass_list_mag = []
    fail_list_file= []
    fail_list_mag = []

    for o in sorted(glob.glob("%s/force.*"%force_output_dir)):

        ofn = os.path.basename(o)
        n = os.path.join(check_dir,ofn)
        if not os.path.exists(n):
            print('could not check %s...file not found in checkdir'%ofn)
        odat = np.loadtxt(o)
        ndat = np.loadtxt(n)

        diff = np.abs((odat-ndat)/odat)
        diff[np.isinf(diff)]=0.0
        diff[np.isnan(diff)]=0.0
        diff = np.amax(diff)

        if diff>1.e-5:
            fail_list_file.append(ofn)
            fail_list_mag.append(diff)
        else:
            pass_list_file.append(ofn)
            pass_list_mag.append(diff)

    for o in range(len(pass_list_file)):
        ofn=pass_list_file[o]
        diff=pass_list_mag[o]
        print('pass:','max force percent diff = % 7.2f '%diff,ofn)
    print() 

    for o in range(len(fail_list_file)):
        ofn=fail_list_file[o]
        diff=fail_list_mag[o]
        print('fail:','max force percent diff = % 7.2f '%diff,ofn)
    print() 

#########################################################################################
#########################################################################################
#########################################################################################

def _get_equiv_atom(in_str):
    # get the equivilent atoms for atoms from the output of fd.x


    # list of atoms indices equivilent to the first entry in each list
    equiv_atom = []
    # list of symops indices that transform atoms in equiv_atoms into each other
    equiv_atom_so = []

    # find all atoms and their equivilent atoms in fd.x output
    try:
        equiv_atom_str_list  = re.findall('Atom (\d+) has \d+ equivalent\(s\):([\d+\s]+)\n*',in_str)
        for equiv_atom_str in equiv_atom_str_list:
            equiv_atom_list = list(map(int,equiv_atom_str[1].split()))
            # shift index for list by 1 fort->py convention
            equiv_atom.append([x-1 for x in equiv_atom_list])
    except: pass

    op_str  = re.findall('for symmetry operation\(s\):([\d+\s]+)\n*',in_str)
    for i in range(len(op_str)):
        op_str_list = list(map(int,op_str[i].split()))
            # shift index for list by 1 fort->py convention
        equiv_atom_so.append([x-1 for x in op_str_list])

    # get length of equiv atom array
    counter=0
    for i in range(len(equiv_atom)):
        for j in range(1,len(equiv_atom[i])):
            counter+=1

    equiv_atom_arr       = np.zeros((counter,2),dtype=int)
    equiv_atom_symop_arr = np.zeros((counter),dtype=int)

    # fill in arrays
    counter=0
    for i in range(len(equiv_atom)):
        for j in range(1,len(equiv_atom[i])):            
            equiv_atom_arr[counter][0] = equiv_atom[i][0]
            equiv_atom_arr[counter][1] = equiv_atom[i][j]
            equiv_atom_symop_arr[counter] = equiv_atom_so[i][j]
            counter+=1
            

    return equiv_atom_arr,equiv_atom_symop_arr

#########################################################################################
#########################################################################################
#########################################################################################

def _get_equiv_dir(in_str):
    # get the equivilent directions for atoms from the output of fd.x


    equiv_dir_list = []

    regex_str = 'on atom\s*(\d+)\s*(\S) and (\S) displacements are equivalent for symmetry operation\s*(\d+)'
    equiv_dir_str  = re.findall(regex_str,in_str)
    for i in equiv_dir_str:
        atom_index  = int(i[0])-1
        symop_index = int(i[3])-1
        if i[1]=='x':
            first_dir_ind = 0
        if i[1]=='y':
            first_dir_ind = 1
        if i[1]=='z':
            first_dir_ind = 2

        if i[2]=='x':
            second_dir_ind = 0
        if i[2]=='y':
            second_dir_ind = 1
        if i[2]=='z':
            second_dir_ind = 2

        equiv_dir_list.append([atom_index,symop_index,first_dir_ind,second_dir_ind])


    return equiv_dir_list

#########################################################################################
#########################################################################################
#########################################################################################

def _process_sym_list_str(fd_out_str):
    # process the output of fd.x to find all the equivilent 
    # atoms and associated symmetry operations


    # regex to find all the symmetry operations
    sym_str_list = re.findall(r'isym.*\n.*\n(.*\n.*\n.*\n).*\n(.*\n.*\n.*\n)',fd_out_str)

    nsym=len(sym_str_list)
    rot_crys   = np.zeros((nsym,3,3),dtype=float)
    rot_cart   = np.zeros((nsym,3,3),dtype=float)
    shift_crys = np.zeros((nsym,3,1),dtype=float)
    shift_cart = np.zeros((nsym,3,1),dtype=float)

    # t is index of crystal vs. cart in the string array
    for t in [0,1]:
        for i in range(len(sym_str_list)):
            # regex to split the arrays of the symmetry ops
            split_sym_str = re.findall(r'\((\s*[-\d.]+\s*[-\d.]+\s*[-\d.]+\s*\s+)\)',sym_str_list[i][t])

            # counter for the rows of the sym op
            op_row_counter=0
            # convert lines of text on sym ops into arrays
            for j in range(len(split_sym_str)):
                sym_str_line = [x.strip() for x in split_sym_str[j].split(' ') if len(x.strip())>0]
                sym_arr_line = np.array(list(map(float,sym_str_line)))

                nele = sym_arr_line.shape[0]
                # if length of the line is 3 then it is a rotation
                if nele==3:
                    # modulo so row index is 0,1,2
                    op_row_index = op_row_counter%3                
                    # t=0 is crystal
                    if t==0:
                        rot_crys[i,op_row_index]=sym_arr_line
                    # t=1 is cart
                    if t==1:
                        rot_cart[i,op_row_index]=sym_arr_line
                    # if part of the rotation matrix add to counter
                    op_row_counter+=1
                # if length of the line is 3 then it is a shift
                if nele==1:
                    # t=0 is crystal
                    if t==0:
                        shift_crys[i,op_row_index] = sym_arr_line
                    # t=1 is cart
                    if t==1:
                        shift_cart[i,op_row_index] = sym_arr_line


    return rot_crys,rot_cart,shift_crys,shift_cart

#########################################################################################
#########################################################################################
#########################################################################################

def _get_ss_pos(in_string):
    # load atomic positions from supercell input file 


    # get atomic positions in bohr
    cell = AFLOWpi.retr.getCellMatrixFromInput(in_string)
    # convert cart to crystal
    pos  = AFLOWpi.retr._getPositions(in_string)/0.529177249
    pos = AFLOWpi.retr._convertFractional(pos,cell).getA()

    # make sure all atoms are inside unit cell
    pos[np.where(np.isclose(pos,0.0))]=0.0
    pos=pos%1.0


    return pos

#########################################################################################
#########################################################################################
#########################################################################################

def _load_forces(force_dir,natoms,nrx1,nrx2,nrx3):
    # load the force data


    natoms_red=int(natoms/(nrx1*nrx2*nrx3))
    # load forces and arrange in array
    force_files = sorted(glob.glob('%s/force.*'%force_dir))
    # index 1: +/- direction of displacement
    # index 2: atom that is displaced
    # index 3: forces on atoms from displacement
    # index 4: direction of displacement xyz
    # index 5: direction of forces xyz
    force_arr=np.zeros((2,natoms_red,natoms,3,3),dtype=float)
    # index 1: forces on atoms (no displacement)
    # index 2: direction of forces xyz
    f0=np.zeros((natoms,3),dtype=float)

    for f in force_files:
        fi = np.array(list(map(int,os.path.basename(f).split('.')[1:])))
        # load zero displacement forces
        if np.all(fi == [0,0,0]):
            f0=np.loadtxt(f)
        else:
            # shift index from fortran to python style
            fi-=1
            try:
                force_arr[fi[0],fi[2],:,fi[1]]=np.loadtxt(f)
            except: pass

    return f0,force_arr

#########################################################################################
#########################################################################################
#########################################################################################

def _fill_equiv_directions(f0,force_arr,ed,imap,rot_cart,aisom,shift_crys,rot_crys,pos,innx):
    # fills in parts of force array for equivilent
    # displacement directions on select atoms


    for p in [0,1]:
        for i in range(len(ed)):        
            ai  = ed[i][0]
            soi = ed[i][1]
            d1  = ed[i][2]
            d2  = ed[i][3]

            # mapping of position in supercell
            # rotate atoms of supercell around atom ai
            acmap = _get_symop_ss_map_acenter(pos,shift_crys,rot_crys,ai,soi)

            # symop that rotates forces for equivilent directions
            so=rot_cart[soi]

            ini_force = force_arr[p,ai,:,d1]

            # don't try to fill in forces if starting array hasn't been filled in yet
            if np.all(np.isclose(ini_force+f0,0)):
                continue

            # subtract out zero displacement forces
            rot_force  = ini_force.dot(np.linalg.inv(so))
            # remap forces to atom positions after rotation
            rot_force  = rot_force[acmap]
            # if axes directions change sign from rotation
            pind=aisom[soi,p,d1,1]
            d2 = aisom[soi,p,d1,0]

            # if foward difference deriv then take negative of forces 
            # if rotated displacement is in opposite direction
            if innx==1 and pind!=p:
                pind=0
                rot_force*=-1.0

            force_arr[pind,ai,:,d2] = rot_force


    return force_arr

#########################################################################################
#########################################################################################
#########################################################################################

def _fill_equiv_atoms(f0,force_arr,eai,easoi,rot_cart,imap,aisom,innx):
    # fills in parts of force array for equivilent
    # displacement directions on select atoms


    for i in range(len(eai)):        
        # index of original atom
        ai  = eai[i][0]
        # index of equivilent atom
        bi = eai[i][1]
        # index of sym op that rotates atom a to atom b
        soi = easoi[i]
        # symop that rotates forces for equivilent directions
        so=rot_cart[soi]

        # loop over displacement directions 
        for d1 in range(3):
            # loop over plus/minus displacements
            for p in [0,1]:
                # force in atoms from displacing atom a
                ini_force = force_arr[p,ai,:,d1]

                # don't try to fill in forces if starting array hasn't been filled in yet
                if np.all(np.isclose(ini_force+f0,0)):
                    continue

                # rotate forces
                rot_force  = ini_force.dot(np.linalg.inv(so))
                # remap forces
                amap=imap[i]
                rot_force  = rot_force[amap]
                # get p/m index
                pind=aisom[soi,p,d1,1]
                # get direction index
                d2=aisom[soi,p,d1,0]

                # if foward difference deriv then take negative of forces 
                # if rotated displacement is in opposite direction
                if innx==1 and pind!=p:
                    pind=0
                    rot_force*=-1.0

                force_arr[pind,bi,:,d2] = rot_force


    return force_arr

#########################################################################################
#########################################################################################
#########################################################################################

def _get_axes_trans(rot_cart):
    # find how the axes transform after a rotation
    # e.g. x -> -y from a clockwise rotation of pi/2 around z


    p_axes_trans = np.zeros((rot_cart.shape[0],3),dtype=int)
    n_axes_trans = np.zeros((rot_cart.shape[0],3),dtype=int)
    axes_trans = np.zeros((rot_cart.shape[0],2,3,2),dtype=int)
    dir_ind = np.array([1,2,3])

    #transform pos and neg axes under rotation operation
    for i in range(rot_cart.shape[0]):
        p_axes_trans[i] = np.dot(( 1.0*dir_ind),rot_cart[i],)
        n_axes_trans[i] = np.dot((-1.0*dir_ind),rot_cart[i],)

    # convert sign of rotated axes to sign indice 0 is + 1 is -
    n_signs = np.sign(n_axes_trans)

    n_signs[np.where(n_signs>0.5)]=0
    n_signs[np.where(n_signs<-0.5)]=1

    p_signs = np.sign(p_axes_trans)
    p_signs[np.where(p_signs>0.5)]=0
    p_signs[np.where(p_signs<-0.5)]=1

    # shift axes indices by 1 
    p_axes_trans=np.abs(p_axes_trans)-1
    n_axes_trans=np.abs(n_axes_trans)-1

    # assign to array
    # indice 1: sym op index
    # indice 2: original sign index
    # indice 3: original direction
    # indice 4: first is new direction, second is new sign
    axes_trans[:,0,:,0]=p_axes_trans
    axes_trans[:,1,:,0]=n_axes_trans
    axes_trans[:,0,:,1]=p_signs
    axes_trans[:,1,:,1]=n_signs
    

    return axes_trans

#########################################################################################
#########################################################################################
#########################################################################################

def _get_symop_ss_map(pos,shift_crys,rot_crys,nrx1,nrx2,nrx3,eai,easoi):
    # maps the atomic positions to their equivilent
    # positions after a transformation
            
    index_mapping=np.zeros((easoi.shape[0],pos.shape[0]),dtype=int)

    # scales the translational shift symmetry to supercell
    shift_scale = np.array([1.0/nrx1,1.0/nrx2,1.0/nrx3])
    
    for i in range(index_mapping.shape[0]):
        #symop index
        soi = easoi[i]

        # rotate
        new_pos   = pos.dot(rot_crys[soi])

        # add symop translation scaled to supercell
        new_pos  += shift_crys[soi][:,0]*shift_scale
        
        # fold rotated, shifted supercell back into original supercell
        new_pos[np.where(np.isclose(new_pos,0.0))]=0.0
        new_pos=new_pos%1.0
        new_pos[np.where(np.isclose(new_pos,1.0))]=0.0

        # new atom in original cell
        orig_pos_new_atom = pos[eai[i][1]]

        # original atom in new rotated cell
        orig_pos_old_atom = new_pos[eai[i][0]]

        # amount to shift atoms to shift original
        # atom in rotated cell to first cell
        ut_shift =  orig_pos_old_atom - orig_pos_new_atom 

        # shift old atom in rotated cell into original unit cell
        new_pos-=ut_shift
        
        # fold rotated, shifted supercell back into original supercell
        new_pos[np.where(np.isclose(new_pos,0.0))]=0.0
        new_pos=new_pos%1.0
        new_pos[np.where(np.isclose(new_pos,1.0))]=0.0
        
        # key that maps equivilent positions
        # before and after rotation
        remap_key=np.zeros((pos.shape[0]),dtype=int)
        remap_key[:]=-999

        # roll over atoms rotated out of original
        # supercell back into original supercell
        # and check for equivilent positions
        for j in range(pos.shape[0]):
            for k in range(new_pos.shape[0]):
                if np.all(np.isclose(pos[j],new_pos[k],rtol=1e-5,atol=1e-6,)):
                    remap_key[j] = k
                                       
        if np.unique(remap_key).shape[0]!=pos.shape[0]:
            print("equiv atom remap failed!")
            print('symop #%s:'%(i+1))
            print()
            print('(% 8.6f % 8.6f % 8.6f )  f = (% 8.6f )'%(rot_crys[soi,0,0],rot_crys[soi,0,1],rot_crys[soi,0,2],shift_crys[soi,0,0],))
            print('(% 8.6f % 8.6f % 8.6f )      (% 8.6f )'%(rot_crys[soi,1,0],rot_crys[soi,1,1],rot_crys[soi,1,2],shift_crys[soi,1,0],))
            print('(% 8.6f % 8.6f % 8.6f )      (% 8.6f )'%(rot_crys[soi,2,0],rot_crys[soi,2,1],rot_crys[soi,2,2],shift_crys[soi,2,0],))
            print()
            raise SystemExit

        index_mapping[i]=remap_key


    return index_mapping

#########################################################################################
#########################################################################################
#########################################################################################

def _get_symop_ss_map_acenter(pos,shift_crys,rot_crys,ai,symop):
    # maps the atomic positions to their equivilent positions 
    # after a rotation around a given atom

    # rotate around atom ai
    new_pos   = ((pos-pos[ai]).dot(rot_crys[symop])) + pos[ai]

    # shift atoms at +1 or -1 to 0 in crystal fractional
    tmp_pos=np.copy(new_pos)

    # fold rotated, shifted supercell back into original supercell
    new_pos[np.where(np.isclose(new_pos,0.0))]=0.0
    new_pos=new_pos%1.0
    new_pos[np.where(np.isclose(new_pos,1.0))]=0.0

    # key that maps equivilent positions
    # before and after rotation
    remap_key=np.zeros((pos.shape[0]),dtype=int)
    remap_key[:]=-999

    # roll over atoms rotated out of original
    # supercell back into original supercell
    # and check for equivilent positions
    for j in range(pos.shape[0]):
        for k in range(new_pos.shape[0]):
            if np.all(np.isclose(pos[j],new_pos[k],rtol=1e-5, atol=1e-6)):
                remap_key[j] = k


    if np.unique(remap_key).shape[0]!=pos.shape[0]:
        print("equiv direction remap failed!")
        print('symop #%s:'%(symop+1))
        print()
        print('(% 8.6f % 8.6f % 8.6f )  f = (% 8.6f )'%(rot_crys[symop,0,0],rot_crys[symop,0,1],rot_crys[symop,0,2],shift_crys[symop,0,0],))
        print('(% 8.6f % 8.6f % 8.6f )      (% 8.6f )'%(rot_crys[symop,1,0],rot_crys[symop,1,1],rot_crys[symop,1,2],shift_crys[symop,1,0],))
        print('(% 8.6f % 8.6f % 8.6f )      (% 8.6f )'%(rot_crys[symop,2,0],rot_crys[symop,2,1],rot_crys[symop,2,2],shift_crys[symop,2,0],))
        print()
        raise SystemExit

    
    return remap_key

#########################################################################################
#########################################################################################
#########################################################################################

def _save_forces(f0,force_arr,save_dir):
    # write forces

 
    if not os.path.exists(save_dir):
        try:
            os.mkdir(save_dir)
        except: 
            print('unable to create force output directory')
            print('exiting...')
            raise SystemExit

    # save zero displacement forces
    fn = os.path.join(save_dir,"force.0.0.0")
    np.savetxt(fn,f0,fmt="% 20.16e")            

    # loop over p/m
    for p in range(force_arr.shape[0]):
        # loop over displacement direction
        for d in range(force_arr.shape[3]):
            # loop over atoms
            for a in range(force_arr.shape[1]):
                fn = os.path.join(save_dir,"force.%s.%s.%s"%(p+1,d+1,a+1))
                if not np.isclose(np.sum(np.abs(force_arr[p,a,:,d])),0.0):
                    np.savetxt(fn,force_arr[p,a,:,d],fmt="% 20.16e")


#########################################################################################
#########################################################################################
#########################################################################################

def _symmetrtize_forces(nrx1,nrx2,nrx3,ss_input_file,force_dir,fd_output_file,force_output_dir,check_dir=''):
    # fill in forces of equivilent atoms and displacement directions


    # load undisplaced supercell input
    with open(ss_input_file) as ifo:
        in_string=ifo.read()

    pos = _get_ss_pos(in_string)
    nat = pos.shape[0]

    # get forces from files
    f0,force_arr = _load_forces(force_dir,nat,nrx1,nrx2,nrx3)

    # if the minus displacement direction is all zero then innx=1
    innx=2
    if np.count_nonzero(force_arr[1])==0:
        innx=1

    # remove displacement forces
    force_arr -= f0[None,None,:,None]

    # get symmetry op matrices in crystal and cartesian
    with open(fd_output_file,'r') as ifi:
            fd_out_str = ifi.read()

    rot_crys,rot_cart,shift_crys,shift_cart = _process_sym_list_str(fd_out_str)
    # eai: equivilent atom index
    # eosoi: equivilent atom symop index
    eai,easoi = _get_equiv_atom(fd_out_str)

    # edi: equivilent direction index
    edi = _get_equiv_dir(fd_out_str)

    # find how symmetry operations rotate cartesian axes 
    # aisom: axes index sym op map
    aisom = _get_axes_trans(rot_cart)

    # get mapping from symop transformation for atom indices in supercell
    imap     = _get_symop_ss_map(pos,shift_crys,rot_crys,nrx1,nrx2,nrx3,eai,easoi)

    # fill in equivilent directions in force array
    force_arr = _fill_equiv_directions(f0,force_arr,edi,imap,rot_cart,aisom,shift_crys,rot_crys,pos,innx)

    # fill in forces for atoms equivilent by wyckoff
    force_arr = _fill_equiv_atoms(f0,force_arr,eai,easoi,rot_cart,imap,aisom,innx)

    # reintroduce zero displacement forces
    force_arr += f0[None,None,:,None]

    # write all force files to disk
    _save_forces(f0,force_arr,force_output_dir)





    if check_dir!='':
        _check(check_dir,force_output_dir)
        try:
            _check(check_dir,force_output_dir)
        except: 
            print('check failed...maybe your check_dir is incorrect or not')
            print('all required force files are available for checking')

#########################################################################################
#########################################################################################
#########################################################################################


def _fd_sym(oneCalc,ID,nrx1,nrx2,nrx3):

    force_dir=os.path.join(oneCalc["_AFLOWPI_FOLDER_"],"%s_FD_PHONON"%ID)
    ss_input_file=os.path.join(force_dir,"%s_FD_PHONON__0001/displaced.0.0.0_01.in"%ID)
    fd_output_file=os.path.join(oneCalc["_AFLOWPI_FOLDER_"],"%s_fd.out"%ID)
    force_output_dir=force_dir


    AFLOWpi.run._symmetrtize_forces(nrx1,nrx2,nrx3,ss_input_file,force_dir,fd_output_file,force_output_dir)


