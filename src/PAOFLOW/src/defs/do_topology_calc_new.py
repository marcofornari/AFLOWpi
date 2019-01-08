

# PAOFLOW
#
# Utility to construct and operate on Hamiltonians from the Projections of DFT wfc on Atomic Orbital bases (PAO)
#
# Copyright (C) 2016-2018 ERMES group (http://ermes.unt.edu, mbn@unt.edu)
#
# Reference:
# M. Buongiorno Nardelli, F. T. Cerasoli, M. Costa, S Curtarolo,R. De Gennaro, M. Fornari, L. Liyanage, A. Supka and H. Wang,
# PAOFLOW: A utility to construct and operate on ab initio Hamiltonians from the Projections of electronic wavefunctions on
# Atomic Orbital bases, including characterization of topological materials, Comp. Mat. Sci. vol. 143, 462 (2018).
#
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

from scipy import fftpack as FFT
import numpy as np
import cmath
import sys,os

from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

from kpnts_interpolation_mesh import *
from do_non_ortho import *
from do_momentum import *
from load_balancing import *
from constants import *
from clebsch_gordan import *
from do_eigh_calc import *

import pfaffian as pf

from get_R_grid_fft import *

import matplotlib.pyplot as plt
from get_degeneracies import *
from do_perturb_split import *

# initialize parallel execution
comm=MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def do_topology_calc(HRs,SRs,non_ortho,kq,E_k,v_kp,R,Rfft,R_wght,idx,alat,b_vectors,nelec,bnd,Berry,ipol,jpol,spin_Hall,spol,spin_orbit,sh,nl,eff_mass,inputpath,npool,a_vectors):


    np.set_printoptions(formatter={'float':lambda x: " % 12.3f"%x})
    degen = get_degeneracies(E_k,bnd)

    # Compute Z2 invariant and topological properties on a selected path in the BZ

    nkpi=kq.shape[1]
    nawf,nawf,nk1,nk2,nk3,nspin = HRs.shape



    # Compute Z2 according to Fu, Kane and Mele (2007)
    # Define TRIM points in 2(0-3)/3D(0-7)
    if nspin == 1 and spin_Hall:
        nktrim = 16
        ktrim = np.zeros((nktrim,3),dtype=float)
        ktrim[0] = np.zeros(3,dtype=float)                                  #0 0 0 0
        ktrim[1] = b_vectors[0,:]/2.0                                       #1 1 0 0
        ktrim[2] = b_vectors[1,:]/2.0                                       #2 0 1 0
        ktrim[3] = b_vectors[0,:]/2.0+b_vectors[1,:]/2.0                    #3 1 1 0
        ktrim[4] = b_vectors[2,:]/2.0                                       #4 0 0 1
        ktrim[5] = b_vectors[1,:]/2.0+b_vectors[2,:]/2.0                    #5 0 1 1
        ktrim[6] = b_vectors[2,:]/2.0+b_vectors[0,:]/2.0                    #6 1 0 1
        ktrim[7] = b_vectors[0,:]/2.0+b_vectors[1,:]/2.0+b_vectors[2,:]/2.0 #7 1 1 1
        ktrim[8:16] = -ktrim[:8]
        # Compute eigenfunctions at the TRIM points
        E_ktrim,v_ktrim = do_eigh_calc(HRs,SRs,ktrim,R_wght,R,idx,non_ortho)
        # Define time reversal operator
        theta = -1.0j*clebsch_gordan(nawf,sh,nl,1)
        wl = np.zeros((nktrim/2,nawf,nawf),dtype=complex)
        for ik in range(nktrim/2):
            wl[ik,:,:] = np.conj(v_ktrim[ik,:,:,0].T).dot(theta).dot(np.conj(v_ktrim[ik+nktrim/2,:,:,0]))
            wl[ik,:,:] = wl[ik,:,:]-wl[ik,:,:].T  # enforce skew symmetry
        delta_ik = np.zeros(nktrim/2,dtype=complex)
        for ik in range(nktrim/2):
            delta_ik[ik] = pf.pfaffian(wl[ik,:nelec,:nelec])/np.sqrt(LAN.det(wl[ik,:nelec,:nelec]))

        f=open(os.path.join(inputpath,'Z2'+'.dat'),'w')
        p2D = np.real(np.prod(delta_ik[:4]))
        if p2D+1.0 < 1.e-5:
            v0 = 1
        elif p2D-1.0 < 1.e-5:
            v0 = 0
        f.write('2D case: v0 = %1d \n' %(v0))
        p3D = np.real(np.prod(delta_ik))
        if p3D+1.0 < 1.e-5:
            v0 = 1
        elif p3D-1.0 < 1.e-5:
            v0 = 0
        p3D = delta_ik[1]*delta_ik[3]*delta_ik[6]*delta_ik[7]
        if p3D+1.0 < 1.e-5:
            v1 = 1
        elif p3D-1.0 < 1.e-5:
            v1 = 0
        p3D = delta_ik[2]*delta_ik[3]*delta_ik[5]*delta_ik[7]
        if p3D+1.0 < 1.e-5:
            v2 = 1
        elif p3D-1.0 < 1.e-5:
            v2 = 0
        p3D = delta_ik[4]*delta_ik[6]*delta_ik[5]*delta_ik[7]
        if p3D+1.0 < 1.e-5:
            v3 = 1
        elif p3D-1.0 < 1.e-5:
            v3 = 0
        f.write('3D case: v0;v1,v2,v3 = %1d;%1d,%1d,%1d \n' %(v0,v1,v2,v3))
        f.close()

    # Compute momenta and kinetic energy

    kq_aux = scatter_full(kq.T,npool)
    kq_aux = kq_aux.T
    # Compute R*H(R)
    Rfft = np.reshape(Rfft,(nk1*nk2*nk3,3),order='C')


    RP=np.copy(R)

    RP2=get_R_grid_fft_d(nk1,nk2,nk3,a_vectors)

    HRs  = np.reshape(HRs,(nawf*nawf,nk1*nk2*nk3,nspin),order='C')
    HRs = np.swapaxes(HRs,0,1)
    HRs = np.reshape(HRs,(nk1*nk2*nk3,nawf,nawf,nspin),order='C')

    Rfft_aux = scatter_full(RP,npool)
    HRs_aux = scatter_full(HRs,npool)

################################################################################################################################
################################################################################################################################
################################################################################################################################

    if spin_Hall:
        # Compute spin current matrix elements
        # Pauli matrices (x,y,z)
        sP=0.5*np.array([[[0.0,1.0],[1.0,0.0]],[[0.0,-1.0j],[1.0j,0.0]],[[1.0,0.0],[0.0,-1.0]]])
        if spin_orbit:
            # Spin operator matrix  in the basis of |l,m,s,s_z> (TB SO)
            Sj = np.zeros((nawf,nawf),dtype=complex)
            for i in range(nawf/2):
                Sj[i,i] = sP[spol][0,0]
                Sj[i,i+1] = sP[spol][0,1]
            for i in range(nawf/2,nawf):
                Sj[i,i-1] = sP[spol][1,0]
                Sj[i,i] = sP[spol][1,1]
        else:
            # Spin operator matrix  in the basis of |j,m_j,l,s> (full SO)
            Sj = clebsch_gordan(nawf,sh,nl,spol)

        #jdHks = np.zeros((3,nawf,nawf,nkpi,nspin),dtype=complex)
        jks = np.zeros((kq_aux.shape[1],3,bnd,bnd,nspin),dtype=complex)

################################################################################################################################
################################################################################################################################
################################################################################################################################


    diag_ind = np.diag_indices(bnd)
    dHks = np.zeros((kq_aux.shape[1],3,nawf,nawf,nspin),dtype=complex)
    velk = np.zeros((kq_aux.shape[1],3,bnd,nspin),dtype=float)

    for l in range(3):
        if l==0:
            sg=np.abs(np.sign(RP2[:,1])*np.sign(RP2[:,2]))
        if l==1:
            sg=np.abs(np.sign(RP2[:,2])*np.sign(RP2[:,0]))
        if l==2:
            sg=np.abs(np.sign(RP2[:,0])*np.sign(RP2[:,1]))


        dHRs  = np.zeros((HRs_aux.shape[0],nawf,nawf,nspin),dtype=complex)
        for ispin in range(nspin):
            for n in range(HRs_aux.shape[1]):
                for m in range(HRs_aux.shape[2]):
                    dHRs[:,n,m,ispin] = -1.0j*alat*HRs_aux[:,n,m,ispin]*Rfft_aux[:,l]*2*np.pi#*sg



        # Compute dH(k)/dk on the path

        # Load balancing
        dHRs = gather_full(dHRs,npool)    

        if rank!=0:
            dHRs = np.zeros((nk1*nk2*nk3,nawf,nawf,nspin),dtype=complex)            
        comm.Bcast(dHRs)

        dHRs = np.reshape(dHRs,(nk1*nk2*nk3,nawf*nawf,nspin),order='C')            
        dHRs = np.swapaxes(dHRs,0,1)
        dHRs = np.reshape(dHRs,(nawf,nawf,nk1*nk2*nk3,nspin),order='C')            

        dHks_aux = np.zeros((kq_aux.shape[1],nawf,nawf,nspin),dtype=complex) # read data arrays from tasks

        dHks_aux[:,:,:,:] = band_loop_H(nspin,nawf,dHRs,kq_aux,R)

        dHRs = None


        # Compute momenta
        dHks[:,l] = dHks_aux

        for ik in range(dHks_aux.shape[0]):
            for ispin in range(nspin):
                velk[ik,l,:,ispin] = (np.conj(v_kp[ik,:,:,ispin].T).dot \
                    (dHks_aux[ik,:,:,ispin]).dot(v_kp[ik,:,:,ispin])[diag_ind[0],diag_ind[1]]).real


        if spin_Hall:
            for ik in range(dHks_aux.shape[0]):
                for ispin in range(nspin):
                    jks[ik,l,:,:,ispin] = (np.conj(v_kp[ik,:,:,ispin].T).dot \
                        (0.5*(np.dot(Sj,dHks_aux[ik,:,:,ispin])+np.dot(dHks_aux[ik,:,:,ispin],Sj))).dot(v_kp[ik,:,:,ispin]))[:bnd,:bnd]
        if spin_Hall:
            Omj_znk = np.zeros((pks.shape[0],bnd),dtype=float)
            Omj_zk = np.zeros((pks.shape[0],1),dtype=float)


################################################################################################################################
################################################################################################################################
################################################################################################################################

    # RP2=np.reshape(RP2,(nk1,nk2,nk3,3),order="C")
    # HRs=np.reshape(HRs,(nk1,nk2,nk3,nawf,nawf,nspin),order="C")

    # flip=np.ones((nk1,nk2,nk2))
    # flip[(nk1/2):]*=-1
    # flip[:,(nk1/2):]*=-1
    # np.save("HRs.npy",HRs)
    # print HRs[1:,1:,1,0,0,0].real*RP2[1:,1:,1,0]*RP2[1:,1:,1,0]#*np.sign(RP2[1:,1:,2,0])
    # print
    # print HRs[1:,1:,-1,0,0,0].real*RP2[1:,1:,-1,0]*RP2[1:,1:,-1,0]

    # raise SystemExit
    #        print gr[ik]
    eff_mass = True
    if eff_mass == True: 
        tks = np.zeros((kq_aux.shape[1],6,bnd,nspin))
        d2HRs = np.zeros((nawf,nawf,nk1*nk2*nk3,nspin),order="C",dtype=complex)

        d2Hks_aux = np.zeros((3,3,kq_aux.shape[1],nawf,nawf,nspin),order="C",dtype=complex)
        for l in range(3):
            for lp in range(3):
                if (l==0 and lp==1) or (l==1 and lp==0):
                    sg=np.abs(np.sign(RP2[:,2]))
                if (l==0 and lp==2) or (l==2 and lp==0):
                    sg=np.abs(np.sign(RP2[:,1]))
                if (l==2 and lp==1) or (l==1 and lp==2):
                    sg=np.abs(np.sign(RP2[:,0]))
                if l==0 and lp==0:
                    sg=np.abs(np.sign(RP2[:,1])*np.sign(RP2[:,2]))
                if l==1 and lp==1:
                    sg=np.abs(np.sign(RP2[:,2])*np.sign(RP2[:,0]))
                if l==2 and lp==2:
                    sg=np.abs(np.sign(RP2[:,0])*np.sign(RP2[:,1]))

                
                for ispin in range(nspin):
                    for m in range(HRs.shape[1]):
                        for n in range(HRs.shape[2]):
                            d2HRs[m,n,:,ispin] = -1.0*HRs[:,m,n,ispin]*RP2[:,l]*RP2[:,lp]*(2*np.pi*alat)**2#*sg


                # Compute d2H(k)/dk*dkp on the path


                d2Hks_aux[l,lp,:,:,:,:] = band_loop_H(nspin,nawf,d2HRs[:,:,:,:],kq_aux,R)

        d2HRs = None



        ij_ind = np.array([[0,0],[1,1],[2,2],[0,1],[1,2],[0,2,]],dtype=int)

        for ik in range(d2Hks_aux.shape[2]):
            for ij in range(6):
                for ispin in range(nspin):
                    l   = ij_ind[ij,0]
                    lp  = ij_ind[ij,1]

                    d2t = d2Hks_aux[l,lp,ik,:,:,ispin]                                                    
                    vp = v_kp[ik,:,:,ispin]

                    rt = do_perturb_split(d2t,vp,degen[ispin][ik])

#                    tks[ik,ij,:,ispin] = (np.dot(np.conj(vp.T),np.dot(d2t,vp))[diag_ind[0],diag_ind[1]]).real
                    tks[ik,ij,:,ispin] = (rt[diag_ind[0],diag_ind[1]]).real





        d2Hks_aux=None
        mkm1 = np.zeros((tks.shape[0],bnd,3,3,nspin))        
        E_temp = np.zeros((bnd,nawf),order="C",dtype=complex)

        # Compute effective mass
        for ispin in range(mkm1.shape[4]):
            for ik in range(mkm1.shape[0]):
                E_temp = ((E_k[ik,:,ispin]-E_k[ik,:,ispin][:,None])[:,:]).T

                E_temp[np.where(np.abs(E_temp)<1.e-4)]=np.inf

                for ij in range(ij_ind.shape[0]):
                    ipol = ij_ind[ij,0]
                    jpol = ij_ind[ij,1]

                    mkm1[ik,:,ipol,jpol,ispin] = tks[ik,ij,:,ispin]

                    pksp_i,pksp_j = do_perturb_split_twoop(dHks[ik,ipol,:,:,ispin],
                                                           dHks[ik,jpol,:,:,ispin],
                                                           v_kp[ik,:,:,ispin],
                                                           degen[ispin][ik])




                    mkm1[ik,:,ipol,jpol,ispin] += np.sum((((pksp_i*pksp_j.T +\
                                                            pksp_j*pksp_i.T) / E_temp).real),axis=1)[:bnd]
#

#                    print np.sum((((pksp_i*pksp_j.T + pksp_j*pksp_i.T) / E_temp).real)[:8,:8],axis=1)/ \
#                        np.sum((((pksp_i*pksp_j.T + pksp_j*pksp_i.T) / E_temp).real)[:8,8:],axis=1)5H
#                    print


#                raise SystemExit
        mkm1[:,:,1,0] = mkm1[:,:,0,1]
        mkm1[:,:,2,1] = mkm1[:,:,1,2]
        mkm1[:,:,2,0] = mkm1[:,:,0,2]
        

        tks=None

        mkm1*=0.003324201



        dem = np.zeros_like(E_k)
        # print mkm1[0,0,:,:,0]
        # print np.sum(mkm1[0,1:4,:,:,0],axis=0)
        # print np.sum(mkm1[0,4:7,:,:,0],axis=0)
        # print mkm1[0,7,:,:,0]


        for ispin in range(mkm1.shape[4]):
            for ik in range(mkm1.shape[0]):

                for n in range(mkm1.shape[1]):
                    if ik==0:
                        print mkm1[ik,n,:,:,0]
                        print
                    # if ik==0 and n==0:
                    #     for a in degen[ispin][ik]:
                    #         try:
                    #             print np.sum(mkm1[0,a,:,:,0],axis=0)
                    #             print
                    #         except: pass

                    try:
                        effm = LAN.eigvalsh(mkm1[ik,n,:,:,ispin]).real


                        if np.prod(effm)<0:                                                                      
                            dos_em = -np.prod(1/np.abs(effm))**(1.0/3.0)                               
                        else:                                                                          
                            dos_em =  np.prod(1/np.abs(effm))**(1.0/3.0)                                       

                        dem[ik,n,ispin] = dos_em

                    except:
                        mkm1[ik,n,:,:,ispin]=0.0

        mkm1 = gather_full(mkm1,npool)
        dem = gather_full(dem,npool)
        
        if rank == 0:
            ij_ind = np.array([[0,0],[1,1],[2,2],[0,1],[1,2],[0,2]],dtype=int)
            for ij in range(ij_ind.shape[0]):
                ipol1 = ij_ind[ij,0]
                jpol1 = ij_ind[ij,1]
            
                for ispin in range(nspin):
                    f=open(os.path.join(inputpath,'effmass'+'_'+str(LL[ipol1])+str(LL[jpol1])+'_'+str(ispin)+'.dat'),'w')
                    for ik in range(nkpi):
                        s="%d\t"%ik
                        for  j in np.real(mkm1[ik,:bnd,ipol1,jpol1,ispin]):s += "% 3.5f "%j
                        s+="\n"
                        f.write(s)
                    f.close()

            for ispin in range(nspin):
                f=open(os.path.join(inputpath,'dos_effmass'+'_'+str(ispin)+'.dat'),'w')
                for ik in range(nkpi):
                    s="%d\t"%ik
                    for  j in dem[ik,:bnd,ispin]:s += "% 3.5f "%j
                    s+="\n"
                    f.write(s)
                f.close()

        mkm1=None
        raise SystemExit
################################################################################################################################
################################################################################################################################
################################################################################################################################    
    HRs_aux = None
    HRs = None
    




        #if rank == 0:
        #    plt.matshow(abs(tks[0,ipol,jpol,:,:,0]))
        #    plt.colorbar()
        #    plt.show()



    # Compute Berry curvature
    if Berry or spin_Hall:
        deltab = 0.05
        mu = -0.2 # chemical potential in eV)
        Om_znk = np.zeros((pks.shape[0],bnd),dtype=float)
        Om_zk = np.zeros((pks.shape[0],1),dtype=float)
        for ik in range(pks.shape[0]):
            for n in range(bnd):
                for m in range(bnd):
                    if m!= n:
                        if Berry:
                            Om_znk[ik,n] += -1.0*np.imag(pks[ik,jpol,n,m,0]*pks[ik,ipol,m,n,0]-pks[ik,ipol,n,m,0]*pks[ik,jpol,m,n,0]) / \
                            ((E_k[ik,m,0] - E_k[ik,n,0])**2 + deltab**2)
                        if spin_Hall:
                            Omj_znk[ik,n] += -2.0*np.imag(jks[ik,ipol,n,m,0]*pks[ik,jpol,m,n,0]) / \
                            ((E_k[ik,m,0] - E_k[ik,n,0])**2 + deltab**2)
            Om_zk[ik] = np.sum(Om_znk[ik,:]*(0.5 * (-np.sign(E_k[ik,:bnd,0]) + 1)))  # T=0.0K
            if spin_Hall: Omj_zk[ik] = np.sum(Omj_znk[ik,:]*(0.5 * (-np.sign(E_k[ik,:bnd,0]-mu) + 1)))  # T=0.0K

    if Berry:
        Om_zk = gather_full(Om_zk,npool)
        if rank == 0:
            f=open(os.path.join(inputpath,'Omega_'+str(LL[spol])+'_'+str(LL[ipol])+str(LL[jpol])+'.dat'),'w')
            for ik in range(nkpi):
                f.write('%3d  %.5f \n' %(ik,-Om_zk[ik,0]))
            f.close()
    if spin_Hall:
        Omj_zk = gather_full(Omj_zk,npool)
        if rank == 0:
            f=open(os.path.join(inputpath,'Omegaj_'+str(LL[spol])+'_'+str(LL[ipol])+str(LL[jpol])+'.dat'),'w')
            for ik in range(nkpi):
                f.write('%3d  %.5f \n' %(ik,Omj_zk[ik,0]))
            f.close()


    velk = gather_full(velk,npool)
    if rank==0:
        for ispin in range(nspin):
            for l in range(3):
                f=open(os.path.join(inputpath,'velocity_'+str(l)+'_'+str(ispin)+'.dat'),'w')
                for ik in range(nkpi):
                    s="%d\t"%ik
                    for  j in velk[ik,l,:bnd,ispin]:s += "% 3.5f "%j
                    s+="\n"
                    f.write(s)
                f.close()

    return()

def band_loop_H(nspin,nawf,HRaux,kq,R):

    kdot = np.zeros((kq.shape[1],R.shape[0]),dtype=complex,order="C")
    kdot = np.tensordot(R,2.0j*np.pi*kq,axes=([1],[0]))
    np.exp(kdot,kdot)

    auxh = np.zeros((nawf,nawf,kq.shape[1],nspin),dtype=complex,order="C")

    for ispin in range(nspin):
        auxh[:,:,:,ispin]=np.tensordot(HRaux[:,:,:,ispin],kdot,axes=([2],[0]))

    kdot  = None
    auxh = np.transpose(auxh,(2,0,1,3))
    return auxh


def get_R_grid_fft_test(nk1,nk2,nk3,a_vectors):
    nrtot = nk1*nk2*nk3
    R = np.zeros((nrtot,3),dtype=float)
    Rfft = np.zeros((3,nk1,nk2,nk3),dtype=float,order="C")
    R_wght = np.ones((nrtot),dtype=float)
    idx = np.zeros((nk1,nk2,nk3),dtype=int)

    for i in range(nk1):
        for j in range(nk2):
            for k in range(nk3):
                n = k + j*nk3 + i*nk2*nk3
                Rx = float(i)/float(nk1)
                Ry = float(j)/float(nk2)
                Rz = float(k)/float(nk3)
                Rx -= 0.5
                Ry -= 0.5
                Rz -= 0.5
                if Rx < 0.5: Rx=Rx+1.0
                if Ry < 0.5: Ry=Ry+1.0
                if Rz < 0.5: Rz=Rz+1.0
                Rx -= int(Rx)
                Ry -= int(Ry)
                Rz -= int(Rz)

                R[n,:] = Rx*nk1,Ry*nk2,Rz*nk3
                Rfft[:,i,j,k] = R[n,:]
                idx[i,j,k]=n

#    Rfft = FFT.ifftshift(Rfft,axes=(1,2,3))
    R = np.ascontiguousarray(np.swapaxes(np.reshape(Rfft,(3,nrtot),order="C"),0,1))
    print Rfft
    raise SystemExit
#    for ik in range(nrtot):
#        print R[ik]
#        R[ik] = np.dot(R[ik],a_vectors)

    return(R,Rfft,R_wght,nrtot,idx)
