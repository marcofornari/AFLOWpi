import do_tetra
import numpy as np


E_k = np.load('32.npy')

print E_k.shape



print 'init done'


dos = np.zeros((1801,2))
en_range = np.linspace(-12,6,1801,endpoint=True)
for i in range(en_range.shape[0]):
    dos[i,0] += en_range[i]



bnd=8
nspin=1

tet = tetra_int(E_k)


nk1=32
nk2=32
nk3=32
nktot=nk1*nk2*nk3
nbnd = E_k.shape[1]

tetra = np.zeros((nbnd,nk1*nk2*nk3,nspin,6,4),order="C")

for ispin in range(0,nspin):
    for n in range(0,bnd):
        tetra[n,:,ispin] = tet.setup_lin_tet_opt(n,ispin)
        print tetra.shape
raise SystemExit
all_weights = []
all_inds    = []
for i in range(en_range.shape[0]):
        inds,wbu = tet.get_weights_mu(en_range[i],tetra)
        inds_mu.append(inds)
        wb_mu.append(wbu)

        all_inds.append(inds_mu)
        all_weights.append(wb_mu)

    
np.savetxt('dos.txt',dos)

import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import glob
import scipy.signal as si
import AFLOWpi
import os



dat = np.loadtxt('dos.txt')
plt.plot(dat[:,0],dat[:,1])


plt.savefig("./CARRIER_CONC.pdf",layout="tight",bbox="tight")
plt.close()

