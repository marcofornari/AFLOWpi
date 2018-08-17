#
# PAOFLOW
#
# Utility to construct and operate on Hamiltonians from the Projections of DFT wfc on Atomic Orbital bases (PAO)
#
# Copyright (C) 2016,2017 ERMES group (http://ermes.unt.edu, mbn@unt.edu)
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
from scipy import fftpack as FFT
import numpy as np
try:
    import psutil
except: pass

import cmath
import sys,os

from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
from write_PAO_eigs import *
from kpnts_interpolation_mesh import *
from do_non_ortho import *
from load_balancing import *
from communication import *
from get_K_grid_fft import *
import scipy.linalg as LAN
# initialize parallel execution
comm=MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()




def do_bands_calc(HRaux,SRaux,kq,R_wght,R,idx,read_S,npool,nk1,nk2,nk3,b_vectors):
    # Load balancing
    nawf,nawf,nktot,nspin = HRaux.shape            
    kq_aux = scatter_full(kq.T,npool)
    kq_aux = kq_aux.T
 
    if read_S:
        Sks_aux = band_loop_S(nspin,nk1,nk2,nk3,nawf,SRaux,R_wght,kq_aux,R,idx)
    else: Sks_aux = None

    print R.shape
    print kq_aux.shape
#    raise SystemExit
    HRaux = band_loop_H_inv(nspin,nktot,nawf,HRaux,R_wght,kq_aux,R,idx)
    


    kq,kq_wght,_,idk = get_K_grid_fft(nk1,nk2,nk3,b_vectors)

    print HRaux.shape
    R=kq_aux.T
    kq_aux = scatter_full(kq.T,npool)
    kq_aux = kq_aux.T    
#    raise SystemExit
    Hks_aux = band_loop_H(nspin,HRaux.shape[2],nawf,HRaux,R_wght,kq_aux,R,idx)
    print Hks_aux.shape
    return Hks_aux



def band_loop_H(nspin,nktot,nawf,HRaux,R_wght,kq,R,idx):


    HRaux = np.reshape(HRaux,(nawf,nawf,nktot,nspin),order='C')
    kdot = np.zeros((kq.shape[1],R.shape[0]),dtype=complex,order="C")
    kdot = np.tensordot(R,2.0j*np.pi*kq,axes=([1],[0]))
    np.exp(kdot,kdot)

    auxh = np.zeros((nawf,nawf,kq.shape[1],nspin),dtype=complex,order="C")

    for ispin in xrange(nspin):
        auxh[:,:,:,ispin]=np.tensordot(HRaux[:,:,:,ispin],kdot,axes=([2],[0]))

    kdot  = None
    return auxh*1.0/nktot

def band_loop_H_inv(nspin,nktot,nawf,HRaux,R_wght,kq,R,idx):



    HRaux = np.reshape(HRaux,(nawf,nawf,nktot,nspin),order='C')
    kdot = np.zeros((kq.shape[1],R.shape[0]),dtype=complex,order="C")
    kdot = np.tensordot(R,-2.0j*np.pi*kq,axes=([1],[0]))
    np.exp(kdot,kdot)

    auxh = np.zeros((nawf,nawf,kq.shape[1],nspin),dtype=complex,order="C")

    for ispin in xrange(nspin):
        auxh[:,:,:,ispin]=np.tensordot(HRaux[:,:,:,ispin],kdot,axes=([2],[0]))

    kdot  = None
    return auxh


def band_loop_S(nspin,nk1,nk2,nk3,nawf,SRaux,R_wght,kq,R,idx):

    nsize = kq.shape[1]
    auxs  = np.zeros((nawf,nawf,nsize),dtype=complex)

    for ik in xrange(kq.shape[1]):
        for i in xrange(nk1):
            for j in xrange(nk2):
                for k in xrange(nk3):
                    phase=R_wght[idx[i,j,k]]*cmath.exp(2.0*np.pi*kq[:,ik].dot(R[idx[i,j,k],:])*1j)
                    auxs[:,:,ik] += SRaux[:,:,i,j,k]*phase

    return(auxs)
