#By Luis Agapito
#April 2014
#%%
#Optimized by Andrew Supka
#September 2016
#%%
import os
import csv
import sys
import numpy as np
from numpy import linalg as la
from scipy import linalg as sla
from Molecule import Molecule
import logging
import integs
import time
import itertools
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
comm=MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

try:
    from cints import contr_coulomb_v3 as ccc
    if rank==0:
        logging.info('Using cints for coulomb integral.')
        print('Using cints for coulomb integral.')
except Exception as e:
    if rank==0:
        logging.warning('cints did not properly import. Switching to pyints.') 
        logging.warning(e)
        print('cints did not properly import. Switching to pyints.')
        print(e)
    from pyints import contr_coulomb_v2 as ccc
#%%


# from mpi4py import MPI

# # -----------------------------------------------------------------------------

# import struct as _struct
# try:
#     from numpy import empty as _empty
#     def _array_new(size, typecode, init=0):
#         a = _empty(size, typecode)
#         a.fill(init)
#         return a
#     def _array_set(ary, value):
#         ary.fill(value)
#     def _array_sum(ary):
#         return ary.sum()
# except ImportError:
#     from array import array as _array
#     def _array_new(size, typecode, init=0):
#         return _array(typecode, [init]) * size
#     def _array_set(ary, value):
#         for i, _ in enumerate(ary):
#             ary[i] = value
#     def _array_sum(ary):
#         return sum(ary, 0)

# # -----------------------------------------------------------------------------

# class Counter(object):

#     def __init__(self, comm, init=0):
#         #
#         size = comm.Get_size()
#         rank = comm.Get_rank()
#         mask = 1
#         while mask < size:
#             mask <<= 1
#         mask >>= 1
#         idx = 0
#         get_idx = []
#         acc_idx = []
#         while mask >= 1:
#             left  = idx + 1
#             right = idx + (mask<<1)
#             if rank < mask:
#                 acc_idx.append( left  )
#                 get_idx.append( right )
#                 idx = left
#             else:
#                 acc_idx.append( right )
#                 get_idx.append( left  )
#                 idx = right
#             rank = rank % mask
#             mask >>= 1
#         #
#         typecode = 'i'
#         datatype = MPI.INT
#         itemsize = datatype.Get_size()
#         #
#         root = 0
#         rank = comm.Get_rank()
#         if rank == root:
#             nlevels = len(get_idx) + 1
#             nentries = (1<<nlevels) - 1
#             self.mem = MPI.Alloc_mem(nentries*itemsize, MPI.INFO_NULL)
#             self.mem[:] = _struct.pack(typecode, init) * nentries
#         else:
#             self.mem = None
#         #
#         self.win = MPI.Win.Create(self.mem, itemsize, MPI.INFO_NULL, comm)
#         self.acc_type = datatype.Create_indexed_block(1, acc_idx).Commit()
#         self.get_type = datatype.Create_indexed_block(1, get_idx).Commit()
#         self.acc_buf = _array_new(len(acc_idx), typecode)
#         self.get_buf = _array_new(len(get_idx), typecode)
#         self.myval = 0

#     def free(self):
#         if self.win:
#             self.win.Free()
#         if self.mem:
#             MPI.Free_mem(self.mem)
#             self.mem = None
#         if self.get_type:
#             self.get_type.Free()
#         if self.acc_type:
#             self.acc_type.Free()

#     def next(self, increment=1):
#         _array_set(self.acc_buf, increment)
#         root = 0
#         self.win.Lock(root)
#         self.win.Get(self.get_buf, root, [0, 1, self.get_type])
#         self.win.Accumulate(self.acc_buf, root, [0, 1, self.acc_type], MPI.SUM)
#         self.win.Unlock(root)
#         nxtval = self.myval + _array_sum(self.get_buf)
#         self.myval += increment
#         return nxtval

# # -----------------------------------------------------------------------------

# class Mutex(object):

#     def __init__(self, comm):
#         self.counter = Counter(comm)

#     def __enter__(self):
#         self.lock()
#         return self

#     def __exit__(self, *exc):
#         self.unlock()
#         return None

#     def free(self):
#         self.counter.free()

#     def lock(self):
#         value = self.counter.next(+1)
#         while value != 0:
#             value = self.counter.next(-1)
#             value = self.counter.next(+1)

#     def unlock(self):
#         self.counter.next(-1)

# # -----------------------------------------------------------------------------


def get_Nmm_spin(Nlm_k,spin_label,Hks,Sks,kpnts_wght):

    #short Fourier transform
    lm_size = Nlm_k.shape[0]
    nbasis  = Nlm_k.shape[1]
    nkpnts  = Nlm_k.shape[2]

    Nlm_aux = np.zeros((lm_size,nbasis),dtype=np.complex128)
    for nk in range(nkpnts):
        #kpnts are in units of 2*pi/alat. alat in Bohrs
        Nlm_aux = Nlm_aux + kpnts_wght[nk]*Nlm_k[:,:,nk]
    Nlm_aux = Nlm_aux/float(np.sum(kpnts_wght))
    Nlm_0 = np.sum(Nlm_aux,axis=1)
    if rank==0:
        print(("get_Nmm_spin: Nlm_0 for spin = %s -->"%spin_label, Nlm_0.real))
    return Nlm_0

def get_hartree_energy_spin(DR_0_up,DR_0_dn,bfs,reduced_basis_2e,fpath):

    
    etemp_U = 0
    etemp_J = 0

    index=np.array(list(itertools.product(reduced_basis_2e,repeat=4)))

    # counter = Counter(comm,init=0)
    index=np.array_split(index,size,0)
    index=comm.scatter(index)
    comm.Barrier()


    for k,l,m,n in index:
    # while True:
    #     c = counter.next()
    #     if c>=index.shape[0]:
    #         break

    #     k,l,m,n = index[c]

        st=time.time()

        myint_U = integs.coulomb(bfs[m],bfs[n],bfs[k],bfs[l],ccc) 
        myint_J = integs.coulomb(bfs[m],bfs[k],bfs[n],bfs[l],ccc) #Pisani

        a_b_0123= DR_0_up[m,n]*DR_0_up[k,l]+DR_0_dn[m,n]*DR_0_dn[k,l]

        etemp_U += (a_b_0123+DR_0_dn[m,n]*DR_0_up[k,l]+DR_0_up[m,n]*DR_0_dn[k,l])*myint_U
        etemp_J += (a_b_0123)*myint_J
#        print(rank,time.time()-st)


    if rank==0:
        etemp_U=comm.reduce(etemp_U)
        etemp_J=comm.reduce(etemp_J)
    else:
        comm.reduce(etemp_U)
        comm.reduce(etemp_J)
    comm.Barrier()


    return etemp_U,etemp_J

def read_basis_unitcell(fpath,latvects,coords,atlabels):
    #nx range of cells in x axis. example nx=range(-2,3)

    #convert from numpy arrays to PyQuante list of tuples
    myatomlist = []
    for i,atomcoords in enumerate(coords):
        myatomlist.append( (atlabels[i].strip(),(atomcoords[0],atomcoords[1],atomcoords[2])) )
        

    atoms=Molecule('unitcell',atomlist = myatomlist,units = 'Angstrom')
    
    #inttol = 1e-6 # Tolerance to which integrals must be equal
    
    basis_file_path = fpath
    bfs = integs.my_getbasis(atoms,basis_file_path)
    if rank==0:
        print("Done generating bfs")
    return bfs


def write_reduced_Dk_spin_v2(fpath,reduced_basis_dm,reduced_basis_2e,spin_label,Hks,Sks):
    #v2 outputs the number of nocc mos
    #spin_label = "up","down","nospin"
    #Similar to write_reduced_Dk_in_k, but the DM is not reduced

    if rank==0:
        print(("write_reduced_Dk_spin_v2: Writing reduced DM(k) for spin=%s"%(spin_label)))
    
    #The size is the full size of the basis
    nbasis  = Hks.shape[0]
    nkpnts  = Hks.shape[2]
    Dk      = np.zeros((nbasis,nbasis,nkpnts),dtype=np.complex128)
    lm_size_2e = reduced_basis_2e.shape[0]
    lm_size_dm = reduced_basis_dm.shape[0]
    Nlm_k   = np.zeros((lm_size_2e,nbasis,nkpnts),dtype=np.complex128)

    #Finding the density matrix at k
    for ik in range(nkpnts):
        #ss = la.inv(sla.sqrtm(Sk)) #S^{-1/2}
        #Hk = Hk.T
        #Sk = Sk.T
        #Mind that Hk has to be in nonorthogonal basis
        Hk = Hks[:,:,ik]
        Sk = Sks[:,:,ik] 
        
        load = False

        w,v =sla.eigh(Hk,Sk) #working with the transposes
        
        #arranging the eigs
        evals     =np.sort(w)
        evecs     =v[:,w.argsort()]

        smearing = 0.0; #change it to +/- 0.0001, or so, if you need .

        indexes  = np.where(evals <=0+smearing)[0] 
        nocc_mo  = indexes.shape[0]
        occ_indexes = indexes[:nocc_mo]

        #Computing the density matrix.
        #n belong to indexes, indexes the occupied MOs
        #D_uv = sum_n c_{un}^{*} . c_vn   
       
        #the basis lm is determined by the input reduced_basis 
        #nocc_mo = indexes.shape[0] #number of occupied MOs

        #lm charge decomposition of each occupied band
        n_lm_dm = np.zeros((lm_size_dm,nocc_mo),dtype=np.complex64) 
        n_lm_2e = np.zeros((lm_size_2e,nocc_mo),dtype=np.complex64) 


        sk_rb_2e=Sk[reduced_basis_2e,:]
        sk_rb_dm=Sk[reduced_basis_dm,:]

        for i_mo in range(nocc_mo):
            cv = evecs[:,i_mo]  #the occupied MOs are in ascending eig order
            n_lm_dm[:,i_mo] = np.conj(cv[reduced_basis_dm]) * (sk_rb_dm.dot(cv))
            n_lm_2e[:,i_mo] = np.conj(cv[reduced_basis_2e]) * (sk_rb_2e.dot(cv))

        Nlm_k[:,:nocc_mo,ik]=n_lm_2e
        uuvv_evecs=evecs[:,occ_indexes]
        n_lm_dm_sum=np.sum(n_lm_dm,0)

        try:
           Dk[:,:,ik] = np.tensordot(np.conj(uuvv_evecs*n_lm_dm_sum),uuvv_evecs,axes=([1],[1])) 
        except:
            for uu in range(nbasis):
                uu_vec=uuvv_evecs[uu]*n_lm_dm_sum 
                for vv in range(nbasis):
                    Dk[uu,vv,ik] = np.vdot(uu_vec,uuvv_evecs[vv]) 


        if ik==0:
           nocc_mo_at_gamma = indexes.shape[0]

    return nocc_mo_at_gamma,Dk,Nlm_k



def read_large_file(fpath,fname):
    fns=fname.split(".")

    bin_file = os.path.join(fpath,fns[0]+".npy")
    
    if os.path.exists(bin_file):
        fin   = open(bin_file,"rb")
        ret=np.load(fin)
        fin.close()
        return ret
    else:
        fin   = open(fpath+'/'+fname,"r")
        ret=np.asarray(list(csv.reader(fin, delimiter=' ',skipinitialspace=True,
                                       quoting=csv.QUOTE_NONNUMERIC)),dtype=np.float32)
        fin.close()
        return ret[:,0]+1j*ret[:,1]



def read_txtdata(fpath,nspin):

    #nspin = 1; non-spin-polarized case
    #nspin = 2; spin-polarized case
    fin   = open(fpath+'/'+'wk.txt',"r")
    kpnts_wght = np.loadtxt(fin)
    fin.close()

    fin   = open(fpath+'/'+'k.txt',"r")

    kpnts = np.loadtxt(fin)
    fin.close()

    if len(kpnts.shape)==1:

        kpnts=kpnts[None]
        kpnts_wght=np.array([kpnts_wght])

    nkpnts  = kpnts.shape[0]
    if rank==0:
        print(("read_txt_data: number of kpoints = %d"%nkpnts))

    kovp_1=read_large_file(fpath,'kovp.txt')

    nbasis  = int(np.sqrt(kovp_1.shape[0]/float(nkpnts)))
    if rank==0:
        print(("read_txt_data: nbasis = %f"%nbasis))

    kovp    = np.reshape(kovp_1,(nbasis,nbasis,nkpnts),order='C')

    for ispin in range(nspin):
        if ispin==0 and nspin==2 : 
           fname = 'kham_up.txt'
        elif ispin==1 and nspin==2 :
           fname = 'kham_dn.txt'
        elif ispin==0 and nspin==1 :
           fname = 'kham.txt'
        else :
            if rank==0:
                print('wrong case 1')

        kham_1 = read_large_file(fpath,fname)

        kham   = np.reshape(kham_1,(nbasis,nbasis,nkpnts),order='C')
        if ispin==0 and nspin==2 : 
           kham_up   = kham
        elif ispin==1 and nspin==2 :
           kham_down = kham
        elif ispin==0 and nspin==1 :
           kham_nospin = kham
        else :
            if rank==0:
                print('wrong case 2')
    
    if nspin == 1: 
       return nkpnts,kpnts,kpnts_wght,kovp,kham_nospin
    elif nspin == 2: 
       return nkpnts,kpnts,kpnts_wght,kovp,kham_up,kham_down
    else:
        if rank==0:
            print("wrong case 3")

def get_DR_0_spin(Dk,spin_label,kpnts_wght):
    if rank==0:
        print(("get_DR_0_spin: Dk reduced spin %s found, shaped %d x %d x %d"%(spin_label,Dk.shape[0],Dk.shape[1],Dk.shape[2])))
    
    nkpnts     =kpnts_wght.shape[0]
    if rank==0:
        print(("get_DR_0_spin: number of kpoints %d"%nkpnts))

    #Overwrite nawf, in case masking of awfc was use
    nawf = Dk.shape[0]
    if rank==0:
        print(("get_DR_0_spin: number of basis %d"%nawf))
        print(("get_DR_0_spin: total kpoints weight %f"%np.sum(kpnts_wght)))
    
    D = np.zeros((nawf,nawf),dtype=np.complex128)
    for nk in range(nkpnts):
        D = D + kpnts_wght[nk]*Dk[:,:,nk]
    D = D/float(np.sum(kpnts_wght))

    return D.real

def test(fpath,reduced_basis_dm,reduced_basis_2e,latvects,coords,atlabels,outfile):

    if rank==0:
        fout = open(fpath+"/"+outfile, "w")
        fout.close()
        fout = open(fpath+"/"+outfile, "r+")
        fout.write("**********************************************************************\n")
        fout.write("* test_dm_solids_spin.py                                             *\n") 
        fout.write("* Computes on-site HF Coulomb + Exchange parameters                  *\n")
        fout.write("* Luis Agapito and Marco Buongiorno-Nardelli, UNT Physics            *\n")
        fout.write("* January 2014                                                       *\n")
        fout.write("**********************************************************************\n")
        fout.write("fpath:        %s\n"%fpath)
        fout.write("outfile:      %s\n"%outfile)
        fout.write("reduced_basis_dm:%s\n"%reduced_basis_dm)
        fout.write("reduced_basis_de:%s\n"%reduced_basis_2e)
        fout.write("latvects:     %s\n"%str(latvects))
        fout.write("coords:       %s\n"%str(coords))
        fout.write("atlabels:     %s\n"%str(atlabels))

    Ha2eV     = 27.211396132 
    Bohr2Angs =  0.529177249

    #%%

    if rank==0:
        print("Generate PyQuante instance of the BFS class")
    bfs     = read_basis_unitcell(fpath,latvects,coords,atlabels)
    nbasis  = len(bfs)
    if rank==0:
        fout.write("PyQuante: Number of basis per prim cell is %d\n"%nbasis)

    ######################################################################
    if rank==0:
        fout.write('Reading the WanT data\n')
        print('Reading the WanT data')
    if nspin == 1: 
       nkpnts,kpnts,kpnts_wght,Sks,Hks_nospin = read_txtdata(fpath,nspin)
    elif nspin == 2: 
       nkpnts,kpnts,kpnts_wght,Sks,Hks_up,Hks_down = read_txtdata(fpath,nspin)
    else:
        if rank==0:
            print('wrong case 1')
    
    if rank==0:
        fout.write('Calculating Nlm_k and reduced D_k''s\n')
        print('Calculating Nlm_k and reduced D_k''s')
    start=time.time()

    if nspin == 1: 
       nocc_mo_gamma,dk,nlm_k = write_reduced_Dk_spin_v2(fpath,reduced_basis_dm,reduced_basis_2e,
                                                         'nospin',Hks_nospin,Sks)
    elif nspin == 2: 
       nocc_mo_gamma,dk_up,nlm_k_up = write_reduced_Dk_spin_v2(fpath,reduced_basis_dm,
                                                               reduced_basis_2e,'up',Hks_up,Sks)
       nocc_mo_gamma,dk_dn,nlm_k_dn = write_reduced_Dk_spin_v2(fpath,reduced_basis_dm,
                                                               reduced_basis_2e,'down',Hks_down,Sks)
                                                
    else:
        if rank==0:
            print('wrong case 2')

    if rank==0:
        fout.write('Calculating Nlm_0\n')
        print('Calculating Nlm_0')
    start=time.time()
    if nspin == 1: 
       Nlm_0_nospin = get_Nmm_spin(nlm_k,'nospin',Hks_nospin,Sks,kpnts_wght)
    elif nspin == 2: 
       Nlm_0_up     = get_Nmm_spin(nlm_k_up,'up'    ,Hks_up    ,Sks,kpnts_wght)
       Nlm_0_down   = get_Nmm_spin(nlm_k_dn,'down'  ,Hks_down  ,Sks,kpnts_wght)
    else:
        if rank==0:
            print('wrong case 3')

    if nspin == 1: 
       Naa=0.0
       lm_size = Nlm_0_nospin.shape[0]
       for m in range(lm_size):
           for mp in range(lm_size):
               if mp == m:
                  continue
               else:
                  Naa = Naa + Nlm_0_nospin[m]*Nlm_0_nospin[mp]
       Nab=0.0
       for m in range(lm_size):
           for mp in range(lm_size):
               Nab = Nab + Nlm_0_nospin[m]*Nlm_0_nospin[mp]
    elif nspin == 2 :
        Naa=0.0
        lm_size = Nlm_0_up.shape[0]
        for m in range(lm_size):
            for mp in range(lm_size):
                if mp == m:
                   continue
                else:
                   Naa = Naa + Nlm_0_up[m]*Nlm_0_up[mp]
        Nbb=0.0
        lm_size = Nlm_0_down.shape[0]
        for m in range(lm_size):
            for mp in range(lm_size):
                if mp == m:
                   continue
                else:
                   Nbb = Nbb + Nlm_0_down[m]*Nlm_0_down[mp]
        Nab=0.0
        for m in range(lm_size):
            for mp in range(lm_size):
                Nab = Nab + Nlm_0_up[m]*Nlm_0_down[mp]
    else:
        if rank==0:
            print('wrong case 4')

    if nspin == 1: 
        if rank==0:
            print(("NaNa + NaNb + NbNa + Nbb = %f"%(2*Nab.real+2*Naa.real)))
        denominator_U = 2*Nab.real+2*Naa.real
        denominator_J = 2*Naa.real
    elif nspin == 2: 
        if rank==0:
            print(("NaNa + NaNb + NbNa + Nbb = %f"%(2*Nab.real+Naa.real+Nbb.real)))
        denominator_U = 2*Nab.real+Naa.real+Nbb.real
        denominator_J = Naa.real+Nbb.real
    else:
        if rank==0:
            print('wrong case')

    if rank==0:
        fout.write("denominator_U = %f\ndenominator_J = %f\n"%(denominator_U,denominator_J))
        print(("denominator_U = %f\ndenominator_J = %f"%(denominator_U,denominator_J)))
        print("Finding the Coulomb and exchange energies")
        fout.write("Started finding the Coulomb and exchange energies at %s\n"%(time.ctime()))


    if  nspin == 1:
        DR_0_up   = get_DR_0_spin(dk,'nospin',kpnts_wght)
        DR_0_down = DR_0_up 
    if  nspin == 2:
        DR_0_up   = get_DR_0_spin(dk_up,'up',kpnts_wght)
        DR_0_down = get_DR_0_spin(dk_dn,'down',kpnts_wght)
    
    if rank==0:
        fout.flush()

    t0   = time.time() 
    U_energy,J_energy = get_hartree_energy_spin(DR_0_up,DR_0_down,bfs,reduced_basis_2e,fpath)
    t1   = time.time() 

    if rank==0:
        print(("Energy Uaa=%+14.10f Ha; Energy Jaa=%+14.10f Ha; %7.3f s"%(U_energy,J_energy,t1-t0)))
        fout.write("Energy Uaa=%+14.10f Ha; Energy Jaa=%+14.10f Ha; %7.3f s\n"%(U_energy,J_energy,t1-t0))

    SI = 0

    U = (U_energy -2*SI)/denominator_U
    J = (J_energy -2*SI)/denominator_J


    if rank==0:
        print(("Parameter U=%f eV"%(U*Ha2eV)))
        fout.write("Parameter U=%f eV\n"%(U*Ha2eV))
        print(("Parameter J=%f eV"%(J*Ha2eV)))
        fout.write("Parameter J=%f eV\n"%(J*Ha2eV))

    if rank==0:    
        if J*Ha2eV == float('Inf'):
            print(("Parameter U_eff = %f eV"%(U*Ha2eV)))
            fout.write("Parameter U_eff = %f eV\n"%(U*Ha2eV))
        else:
            print(("Parameter U_eff = %f eV"%((U-J)*Ha2eV)))
            fout.write("Parameter U_eff = %f eV\n"%((U-J)*Ha2eV))

    tb = time.time()

    if rank==0:
        fout.write("Finished finding the Coulomb energy at %s, elapsed %f s\n"%(time.ctime(),tb-ta))
        fout.close() 

if __name__ == '__main__':
    Bohr2Angs =  0.529177249
    inputfile = sys.argv[1]
    do_sk=True
    try:
        if sys.argv[2]=="skip_sk":
            do_sk=False
    except:
        pass
    input_data = {}
    f = open(inputfile)
    data = f.readlines()
    ta = time.time()
    for line in data:
        line = line.strip()
        if line and not line.startswith("#"):
           line = line.strip()
           # parse input, assign val=es to variables

           key, value = line.split("=")
           input_data[key.strip()] = value.strip()
    f.close()
    
    
    fpath         = input_data['fpath']
    outfile       = input_data['outfile']
    nspin         = int(input_data['nspin'])
    reduced_basis_dm = np.fromstring(input_data['reduced_basis_dm'], dtype=int, sep=',' ) 
    reduced_basis_2e = np.fromstring(input_data['reduced_basis_2e'], dtype=int, sep=',' ) 
    latvects      = np.fromstring(input_data['latvects'], dtype=float, sep=',' ) 
    latvects      = np.reshape(latvects,(3,3))*Bohr2Angs
    atlabels      = input_data['atlabels']
    atlabels      = atlabels.strip(",").split(",")
    coords        = np.fromstring(input_data['coords'], dtype=float, sep=',' ) 
    coords        = np.reshape(coords,(-1,3))*Bohr2Angs

    test(fpath,reduced_basis_dm,reduced_basis_2e,latvects,coords,atlabels,outfile)



