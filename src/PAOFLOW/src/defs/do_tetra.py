
import numpy as np
from communication import *
#################################################################################
#################################################################################
#################################################################################



def get_dos_mu(mu,edge_en):
     if np.all(edge_en>mu) or np.all(edge_en<mu):
          return (),np.array([])

     nktot   =  edge_en.shape[0]

     weights = np.zeros(edge_en.shape[0])

     for i in range(6):
          case1 = np.where(np.logical_and(edge_en[:,i,0]<mu,edge_en[:,i,1]>mu))[0]
          case2 = np.where(np.logical_and(edge_en[:,i,1]<mu,edge_en[:,i,2]>mu))[0]
          case3 = np.where(np.logical_and(edge_en[:,i,2]<mu,edge_en[:,i,3]>mu))[0]


          ee_c1 = np.ascontiguousarray(edge_en[case1,i].T)
          ee_c2 = np.ascontiguousarray(edge_en[case2,i].T)
          ee_c3 = np.ascontiguousarray(edge_en[case3,i].T)


          try:
               weights[case1] += 3.0*(mu-ee_c1[0])**2/((ee_c1[1]-ee_c1[0])*(ee_c1[2]-ee_c1[0])*(ee_c1[3]-ee_c1[0]))
          except Exception as e: pass

          try:
               E21 = ee_c2[1]-ee_c2[0]
               E31 = ee_c2[2]-ee_c2[0]
               E32 = ee_c2[2]-ee_c2[1]
               E41 = ee_c2[3]-ee_c2[0]
               E42 = ee_c2[3]-ee_c2[1]
               E43 = ee_c2[3]-ee_c2[2]
               EMN = mu-ee_c2[1]

               weights[case2] += (3.0*E21 + 6.0*EMN - 3.0*(E31+E42)*EMN**2/(E32*E42))/(E31*E41)
          except Exception as e: pass

          try:
               weights[case3] += 3.0*(mu-ee_c3[3])**2/((ee_c3[3]-ee_c3[0])*(ee_c3[3]-ee_c3[1])*(ee_c3[3]-ee_c3[2]))
          except Exception as e: pass

     inds = np.where(weights!=0.0)[0]


     return inds,weights[inds]/6.0


def get_weights_tetra(tetra_dat,emin,emax,ne):

     tetra_en = tetra_dat[0]

     en_arr = np.linspace(emin,emax,ne,endpoint=True)


     by_en_wgts = []
     by_en_inds = []
     for en in range(en_arr.shape[0]):
          by_spin_wgts = []
          by_spin_inds = []

          for ispin in range(tetra_en.shape[2]):
               by_n_wgts = []
               by_n_inds = []
               for n in range(tetra_en.shape[1]):
#                    inds,weights = get_weights_mu(en_arr[en],tetra_en[:,n,ispin])
                    inds,weights = get_dos_mu(en_arr[en],tetra_en[:,n,ispin])
                    by_n_wgts.append(weights)
                    by_n_inds.append(inds)
               by_spin_wgts.append(by_n_wgts)
               by_spin_inds.append(by_n_inds)

          by_en_wgts.append(by_spin_wgts)
          by_en_inds.append(by_spin_inds)


     tetra_dat.append(en_arr)
     tetra_dat.append(by_en_inds)
     tetra_dat.append(by_en_wgts)


     return tetra_dat
#################################################################################
#################################################################################
#################################################################################

def get_tetra_ind(nk1,nk2,nk3):

     edge_ind = np.array([[0, 1, 2, 5],
                          [0, 2, 4, 5],
                          [2, 4, 5, 6],
                          [2, 5, 6, 7],
                          [2, 3, 5, 7],
                          [1, 2, 3, 5]],dtype=int)

     corner_ind = np.zeros((8,3),dtype=int)

     tet_ind = np.zeros((nk1,nk2,nk3,6,4,3),dtype=int)

     for k1 in range(nk1):
          for k2 in range(nk2):
               for k3 in range(nk3):
                    corner_ind[0] = [k1,k2,k3] # 0,0,0
                    corner_ind[1] = [(k1+1)%nk1,k2,k3] # 1,0,0
                    corner_ind[2] = [k1,(k2+1)%nk2,k3] # 0,1,0
                    corner_ind[4] = [k1,k2,(k3+1)%nk3] # 0,0,1
                    corner_ind[3] = [(k1+1)%nk1,(k2+1)%nk2,k3] # 1,1,0
                    corner_ind[5] = [(k1+1)%nk1,k2,(k3+1)%nk3] # 1,0,1
                    corner_ind[6] = [k1,(k2+1)%nk2,(k3+1)%nk3] # 0,1,1
                    corner_ind[7] = [(k1+1)%nk1,(k2+1)%nk2,(k3+1)%nk3] # 1,1,1

                    for t in range(6):
                         tet_ind[k1,k2,k3,t] = corner_ind[edge_ind[t]]



     return tet_ind

#################################################################################
#################################################################################
#################################################################################

def get_opt_tet_inds(tifo):
     nk1 = tifo.shape[0]
     nk2 = tifo.shape[1]
     nk3 = tifo.shape[2]

     tiff = np.zeros((nk1,nk2,nk3,6,20,3),dtype=int,order="C")



     mda=np.array([nk1,nk2,nk3],dtype=int)

     tif = np.zeros((20,3),dtype=int,order="C")
     # for k1 in range(nk1):
     #      for k2 in range(nk2):
     #           for k3 in range(nk3):
     #                for t in range(6):
     #                     tif[0]=tifo[k1,k2,k3,t,0]
     #                     tif[1]=tifo[k1,k2,k3,t,1]
     #                     tif[2]=tifo[k1,k2,k3,t,2]
     #                     tif[3]=tifo[k1,k2,k3,t,3]

     #                     tif_p1=tif[:4]+1

     #                     tif[4] =((2*tif_p1[0]-tif_p1[1])%mda)-1
     #                     tif[5] =((2*tif_p1[1]-tif_p1[2])%mda)-1
     #                     tif[6] =((2*tif_p1[2]-tif_p1[3])%mda)-1
     #                     tif[7] =((2*tif_p1[3]-tif_p1[0])%mda)-1
     #                     tif[8] =((2*tif_p1[0]-tif_p1[2])%mda)-1
     #                     tif[9] =((2*tif_p1[1]-tif_p1[3])%mda)-1
     #                     tif[10]=((2*tif_p1[2]-tif_p1[0])%mda)-1
     #                     tif[11]=((2*tif_p1[3]-tif_p1[1])%mda)-1
     #                     tif[12]=((2*tif_p1[0]-tif_p1[3])%mda)-1
     #                     tif[13]=((2*tif_p1[1]-tif_p1[0])%mda)-1
     #                     tif[14]=((2*tif_p1[2]-tif_p1[1])%mda)-1
     #                     tif[15]=((2*tif_p1[3]-tif_p1[2])%mda)-1

     #                     tif[16]=((tif_p1[3]-(tif_p1[0])+tif_p1[1])%mda)-1
     #                     tif[17]=((tif_p1[0]-(tif_p1[1])+tif_p1[2])%mda)-1
     #                     tif[18]=((tif_p1[1]-(tif_p1[2])+tif_p1[3])%mda)-1
     #                     tif[19]=((tif_p1[2]-(tif_p1[3])+tif_p1[0])%mda)-1

     #                     tiff[k1,k2,k3,t] = tif
     # return tiff

     tifc  = np.copy(tifo)
     tifc += 1

     tiff[:,:,:,:,0] = tifc[:,:,:,:,0]
     tiff[:,:,:,:,1] = tifc[:,:,:,:,1]
     tiff[:,:,:,:,2] = tifc[:,:,:,:,2]
     tiff[:,:,:,:,3] = tifc[:,:,:,:,3]

     tiff[:,:,:,:,4] = 2*tifc[:,:,:,:,0]-tifc[:,:,:,:,1]
     tiff[:,:,:,:,5] = 2*tifc[:,:,:,:,1]-tifc[:,:,:,:,2]
     tiff[:,:,:,:,6] = 2*tifc[:,:,:,:,2]-tifc[:,:,:,:,3]
     tiff[:,:,:,:,7] = 2*tifc[:,:,:,:,3]-tifc[:,:,:,:,0]
     tiff[:,:,:,:,8] = 2*tifc[:,:,:,:,0]-tifc[:,:,:,:,2]
     tiff[:,:,:,:,9] = 2*tifc[:,:,:,:,1]-tifc[:,:,:,:,3]
     tiff[:,:,:,:,10]= 2*tifc[:,:,:,:,2]-tifc[:,:,:,:,0]
     tiff[:,:,:,:,11]= 2*tifc[:,:,:,:,3]-tifc[:,:,:,:,1]
     tiff[:,:,:,:,12]= 2*tifc[:,:,:,:,0]-tifc[:,:,:,:,3]
     tiff[:,:,:,:,13]= 2*tifc[:,:,:,:,1]-tifc[:,:,:,:,0]
     tiff[:,:,:,:,14]= 2*tifc[:,:,:,:,2]-tifc[:,:,:,:,1]
     tiff[:,:,:,:,15]= 2*tifc[:,:,:,:,3]-tifc[:,:,:,:,2]

     tiff[:,:,:,:,16]= tifc[:,:,:,:,3]-tifc[:,:,:,:,0]+tifc[:,:,:,:,1]
     tiff[:,:,:,:,17]= tifc[:,:,:,:,0]-tifc[:,:,:,:,1]+tifc[:,:,:,:,2]
     tiff[:,:,:,:,18]= tifc[:,:,:,:,1]-tifc[:,:,:,:,2]+tifc[:,:,:,:,3]
     tiff[:,:,:,:,19]= tifc[:,:,:,:,2]-tifc[:,:,:,:,3]+tifc[:,:,:,:,0]

     tiff[...,0] %= nk1
     tiff[...,1] %= nk2
     tiff[...,2] %= nk3

     tiff -= 1          

     tiff[...,0] %= nk1
     tiff[...,1] %= nk2
     tiff[...,2] %= nk3




     return tiff

#################################################################################
#################################################################################
#################################################################################

def sort_tetra_ind(E_k,tetra_ind):

     nk1,nk2,nk3 = E_k.shape
     tet_E_temp  = np.zeros(4)

     for k1 in range(nk1):
          for k2 in range(nk2):
               for k3 in range(nk3):
                    for t in range(6):
                         oti = tetra_ind[k1,k2,k3,t] # set of k1,k2,k3 for one tet

                         tet_E_temp[0] = E_k[oti[0,0],oti[0,1],oti[0,2]]
                         tet_E_temp[1] = E_k[oti[1,0],oti[1,1],oti[1,2]]
                         tet_E_temp[2] = E_k[oti[2,0],oti[2,1],oti[2,2]]
                         tet_E_temp[3] = E_k[oti[3,0],oti[3,1],oti[3,2]]



                         tetra_ind[k1,k2,k3,t] = oti[np.argsort(tet_E_temp)]

     oti = tet_E_temp = None 

     return tetra_ind

#################################################################################
#################################################################################
#################################################################################

def get_opt_tet_adj_E(E_k,tet_opt_ind):

     nk1,nk2,nk3 = E_k.shape

     wlsm = np.zeros((4,20),dtype=float,order="C")

     wlsm[0, 0: 4] = [1440,    0,   30,    0,] 
     wlsm[1, 0: 4] = [   0, 1440,    0,   30,] 
     wlsm[2, 0: 4] = [  30,    0, 1440,    0,] 
     wlsm[3, 0: 4] = [   0,   30,    0, 1440,] 

     wlsm[0, 4: 8] = [ -38,    7,   17,  -28,] 
     wlsm[1, 4: 8] = [ -28,  -38,    7,   17,] 
     wlsm[2, 4: 8] = [  17,  -28,  -38,    7,] 
     wlsm[3, 4: 8] = [   7,   17,  -28,  -38,] 

     wlsm[0, 8:12] = [ -56,    9,  -46,    9,] 
     wlsm[1, 8:12] = [   9,  -56,    9,  -46,] 
     wlsm[2, 8:12] = [ -46,    9,  -56,    9,] 
     wlsm[3, 8:12] = [   9,  -46,    9,  -56,] 

     wlsm[0,12:16] = [ -38,  -28,   17,    7,] 
     wlsm[1,12:16] = [   7,  -38,  -28,   17,] 
     wlsm[2,12:16] = [  17,    7,  -38,  -28,] 
     wlsm[3,12:16] = [ -28,   17,    7,  -38,] 

     wlsm[0,16:20] = [ -18,  -18,   12,  -18,] 
     wlsm[1,16:20] = [ -18,  -18,  -18,   12,] 
     wlsm[2,16:20] = [  12,  -18,  -18,  -18,] 
     wlsm[3,16:20] = [ -18,   12,  -18,  -18,] 

     wlsm /= 1260.0

     E_k_ot = np.zeros((nk1,nk2,nk3,6,4),dtype=float,order="C")
     E_k_20 = np.zeros((nk1,nk2,nk3,6,20),dtype=float,order="C")

     for k1 in range(nk1):
          for k2 in range(nk2):
               for k3 in range(nk3):
                    for t in range(6):
                         oti = tet_opt_ind[k1,k2,k3,t] # set of k1,k2,k3 for one tet

                         E_k_ot[k1,k2,k3,t] = np.dot(wlsm,E_k[oti[:,0],oti[:,1],oti[:,2]])

     oti = None

     return E_k_ot

#################################################################################
#################################################################################
#################################################################################

# def setup_lin_tet_opt(E_k,nk1,nk2,nk3,npool):

#      E_k_copy = gather_scatter(np.ascontiguousarray(E_k),1,npool)

#      _,nbnd,nspin = E_k_copy.shape

#      E_k_copy = np.reshape(E_k_copy,(nk1,nk2,nk3,nbnd,nspin),order="C")

#      tet_ind = get_tetra_ind(nk1,nk2,nk3)

#      edge_en = np.zeros((nbnd,nk1*nk2*nk3,nspin,6,4),order="C")

#      for ispin in range(nspin):
#           for n in range(nbnd):
#                tet_ind_temp = sort_tetra_ind(np.ascontiguousarray(E_k_copy[:,:,:,n,ispin]),tet_ind)
#                tet_opt_ind = get_opt_tet_inds(tet_ind_temp)
#                edge_en_temp  = get_opt_tet_adj_E(np.ascontiguousarray(E_k_copy[:,:,:,n,ispin]),tet_opt_ind)

#                edge_en[n,:,ispin] = np.reshape(edge_en_temp,(nk1*nk2*nk3,6,4),order="C")

     
#      edge_en = gather_scatter(edge_en,1,npool)

#      return np.ascontiguousarray(np.swapaxes(edge_en,0,1))
     
#################################################################################
#################################################################################
#################################################################################
def setup_lin_tet_opt(E_k,nk1,nk2,nk3,npool):

     E_k_copy = gather_scatter(np.ascontiguousarray(E_k),1,npool)

     _,nbnd,nspin = E_k_copy.shape

     E_k_copy = np.reshape(E_k_copy,(nk1,nk2,nk3,nbnd,nspin),order="C")

     tet_ind = get_tetra_ind(nk1,nk2,nk3)

     edge_en = np.zeros((nbnd,nk1*nk2*nk3,nspin,6,4),order="C")

     for ispin in range(nspin):
          for n in range(nbnd):
               tet_ind_temp  = sort_tetra_ind(np.ascontiguousarray(E_k_copy[:,:,:,n,ispin]),tet_ind)
               tet_opt_ind   = get_opt_tet_inds(tet_ind_temp)
               edge_en_temp  = get_opt_tet_adj_E(np.ascontiguousarray(E_k_copy[:,:,:,n,ispin]),tet_opt_ind)

               edge_en[n,:,ispin] = np.reshape(edge_en_temp,(nk1*nk2*nk3,6,4),order="C")





     edge_en = gather_scatter(edge_en,1,npool)

     edge_en = np.swapaxes(edge_en,0,1)

     edge_en_ind = np.argsort(edge_en,axis=4)


     edge_en = np.sort(edge_en,axis=4)

     edge_en = np.ascontiguousarray(edge_en)


     return [edge_en,nk1,nk2,nk3,npool,edge_en_ind]
     
#################################################################################
#################################################################################
#################################################################################

def get_weights_mu(mu,edge_en):
     if np.all(edge_en>mu) or np.all(edge_en<mu):
          return (),np.array([])




     wlsm = np.zeros((20,4),dtype=float,order="C")

     wlsm[ 0: 4, 0] = [1440,    0,   30,    0,] 
     wlsm[ 0: 4, 1] = [   0, 1440,    0,   30,] 
     wlsm[ 0: 4, 2] = [  30,    0, 1440,    0,] 
     wlsm[ 0: 4, 3] = [   0,   30,    0, 1440,] 
                 
     wlsm[ 4: 8, 0] = [ -38,    7,   17,  -28,] 
     wlsm[ 4: 8, 1] = [ -28,  -38,    7,   17,] 
     wlsm[ 4: 8, 2] = [  17,  -28,  -38,    7,] 
     wlsm[ 4: 8, 3] = [   7,   17,  -28,  -38,] 
                 
     wlsm[ 8:12, 0] = [ -56,    9,  -46,    9,] 
     wlsm[ 8:12, 1] = [   9,  -56,    9,  -46,] 
     wlsm[ 8:12, 2] = [ -46,    9,  -56,    9,] 
     wlsm[ 8:12, 3] = [   9,  -46,    9,  -56,] 
                 
     wlsm[12:16, 0] = [ -38,  -28,   17,    7,] 
     wlsm[12:16, 1] = [   7,  -38,  -28,   17,] 
     wlsm[12:16, 2] = [  17,    7,  -38,  -28,] 
     wlsm[12:16, 3] = [ -28,   17,    7,  -38,] 
                 
     wlsm[16:20, 0] = [ -18,  -18,   12,  -18,] 
     wlsm[16:20, 1] = [ -18,  -18,  -18,   12,] 
     wlsm[16:20, 2] = [  12,  -18,  -18,  -18,] 
     wlsm[16:20, 3] = [ -18,   12,  -18,  -18,] 

     wlsm /= 1260.0

     nktot   =  edge_en.shape[0]

     weights = np.zeros((edge_en.shape[0],6,4),order="C")

     weights_aux = np.zeros((4,edge_en.shape[0]))

     for i in range(6):
          case1 = np.where(np.logical_and(edge_en[:,i,0]<mu,edge_en[:,i,1]>mu))[0]
          case2 = np.where(np.logical_and(edge_en[:,i,1]<mu,edge_en[:,i,2]>mu))[0]
          case3 = np.where(np.logical_and(edge_en[:,i,2]<mu,edge_en[:,i,3]>mu))[0]


          ee_c1 = np.ascontiguousarray(edge_en[case1,i].T)
          ee_c2 = np.ascontiguousarray(edge_en[case2,i].T)
          ee_c3 = np.ascontiguousarray(edge_en[case3,i].T)



          try:
              E1  = mu - ee_c1[0]
              E21 = ee_c1[1] - ee_c1[0]
              E31 = ee_c1[2] - ee_c1[0]
              E41 = ee_c1[3] - ee_c1[0]
              C   = E1**3/(24.0*E21*E31*E41)
              weights_aux[0,case1] = C*(4.0-E1*(1.0/E21 + 1.0/E31 + 1.0/E41))
              weights_aux[1,case1] = C*E1/E21
              weights_aux[2,case1] = C*E1/E31
              weights_aux[3,case1] = C*E1/E41

              # dC  = E1**2/(8.0*E21*E31*E41)
              # weights_aux[0,case1] = dC*(4 - E1*(1.0/E21 + 1.0/E31 + 1.0/E41)) - C*(1.0/E21 + 1.0/E31 + 1.0/E41)
              # weights_aux[1,case1] = dC*E1/E21 + C/E21
              # weights_aux[2,case1] = dC*E1/E31 + C/E31
              # weights_aux[3,case1] = dC*E1/E41 + C/E41

          except Exception as e: print(e)

          try:
              E1  =  mu - ee_c2[0]
              E2  =  mu - ee_c2[1]
              E3  =  mu - ee_c2[2]
              E4  =  mu - ee_c2[3]
              E31 =  ee_c2[2] - ee_c2[0]
              E32 =  ee_c2[2] - ee_c2[1]
              E41 =  ee_c2[3] - ee_c2[0]
              E42 =  ee_c2[3] - ee_c2[1]
              C1  =  E1**2/(24.0*E41*E31)
              C2  = -(E1*E2*E3)/(24.0*E41*E32*E31)
              C3  = -(E2**2*E4)/(24.0*E42*E32*E41)
              weights_aux[0,case2] = C1-E3/E31*(C1+C2)-E4/E41*(C1+C2+C3) 
              weights_aux[1,case2] = (C1+C2+C3)-E3/E31*(C2+C3)-E4/E42*C3
              weights_aux[2,case2] = E1/E31*(C1+C2) + E2/E32*(C2+C3)
              weights_aux[3,case2] = E1/E41*(C1+C2+C3) + E2/E42*C3
                                     

              # dC1 =  E1/(12.0*E41*E31)
              # dC2 = -(E1*E2 + E1*E3 + E2*E3)/(24.0*E41*E32*E31)
              # dC3 = -(E2*(2.0*E4 + E2))/(24.0*E42*E32*E41)
              # weights_aux[0,case2] = dC1 - (dC1 + dC2)*E3/E31 - (C1 + C2)/E31 - \
              #                        (dC1 + dC2 + dC3)*E4/E41 - (C1 + C2 + C3)/E41
              # weights_aux[1,case2] = dC1 + dC2 + dC3 - (dC2 + dC3)*E3/E32 - \
              #                        (C2 + C3)/E32 - dC3*E4/E42 - C3/E42
              # weights_aux[2,case2] = (dC1 + dC2)*E1/E31 + (C1 + C2)/E31 + \
              #                        (dC2 + dC3)*E2/E32 + (C2 + C3)/E32
              # weights_aux[3,case2] = E1/E41*(dC1 + dC2 + dC3) + (C1 + C2 + C3)/E41 + \
              #                        dC3*E2/E42 + C3/E42





          except Exception as e: print(e)

          try:
              E4  = mu - ee_c3[3]
              E41 = ee_c3[3] - ee_c3[0]
              E42 = ee_c3[3] - ee_c3[1]
              E43 = ee_c3[3] - ee_c3[2]
              C   = -1.0/24.0 * E4**3 / E41 / E42 / E43
              weights_aux[0,case3] = 1.0/24.0 + C*E4/E41
              weights_aux[1,case3] = 1.0/24.0 + C*E4/E42
              weights_aux[2,case3] = 1.0/24.0 + C*E4/E43
              weights_aux[3,case3] = 1.0/24.0

              # dC  = -E4**2/(8.0*E41*E42*E43)
              # weights_aux[0,case3] = dC*E4/E41 + C/E41
              # weights_aux[1,case3] = dC*E4/E42 + C/E42
              # weights_aux[2,case3] = dC*E4/E43 + C/E43
              # weights_aux[3,case3] = -dC*(4.0 + E4*(1.0/E41 + 1.0/E42 + 1.0/E43)) - C*(1.0/E41 + 1.0/E42 + 1.0/E43)



          except Exception as e: print(e)

          # try:
              
          #      wgt = np.zeros((4,ee_c1.shape[1]),order="C")

          #      "f12 = (e0-e2)/(e1-e2)"
          #      f12 = (mu-ee_c1[1])/(ee_c1[0]-ee_c1[1])
          #      "f13 = (e0-e3)/(e1-e3)"        
          #      f13 = (mu-ee_c1[2])/(ee_c1[0]-ee_c1[2])
          #      "f14 = (e0-e4)/(e1-e4)"
          #      f14 = (mu-ee_c1[3])/(ee_c1[0]-ee_c1[3])

          #      f21 = 1.0 - f12
          #      f31 = 1.0 - f13
          #      f41 = 1.0 - f14

          #      "G  =  3.0 * f21 * f31 * f41 / (e0-e1)"
          #      G  =  3.0 * f21 * f31 * f41 / (mu-ee_c1[0])

          #      weights_aux[0,case1] =  (f12 + f13 + f14) * G
          #      weights_aux[1,case1] =  f21 * G
          #      weights_aux[2,case1] =  f31 * G
          #      weights_aux[3,case1] =  f41 * G

          #      weights_aux[:,case1] = wgt 

          # except Exception,e: print e


          # try:

          #      wgt = np.zeros((4,ee_c2.shape[1]),order="C")

          #      "f13 = (e0-e3)/(e1-e3)"
          #      f13 = (mu-ee_c2[2])/(ee_c2[0]-ee_c2[2])
          #      "f14 = (e0-e4)/(e1-e4)"
          #      f14 = (mu-ee_c2[3])/(ee_c2[0]-ee_c2[3])
          #      "f23 = (e0-e3)/(e2-e3)"
          #      f23 = (mu-ee_c2[2])/(ee_c2[1]-ee_c2[2])
          #      "f24 = (e0-e4)/(e2-e4)"
          #      f24 = (mu-ee_c2[3])/(ee_c2[1]-ee_c2[3])

          #      f31 = 1.0 - f13
          #      f41 = 1.0 - f14
          #      f32 = 1.0 - f23
          #      f42 = 1.0 - f24                 
               


          #      G = 1.0/(ee_c2[3]-ee_c2[0])

          #      weights_aux[0,case2] =  f14 / 3.0 + f13*f31*f23 * G
          #      weights_aux[1,case2] =  f23 / 3.0 + f24*f24*f32 * G 
          #      weights_aux[2,case2] =  f32 / 3.0 + f31*f31*f23 * G
          #      weights_aux[3,case2] =  f41 / 3.0 + f42*f24*f32 * G

          # except Exception,e: print e


          # try:                 

          #      wgt = np.zeros((4,ee_c3.shape[1]),order="C")

          #      "f14 = (e0-e4)/(e1-e4)"
          #      f14 = (mu-ee_c3[3])/(ee_c3[0]-ee_c3[3])
          #      "f24 = (e0-e4)/(e2-e4)"
          #      f24 = (mu-ee_c3[3])/(ee_c3[1]-ee_c3[3])
          #      "f34 = (e0-e4)/(e3-e4)"
          #      f34 = (mu-ee_c3[3])/(ee_c3[2]-ee_c3[3])
          #      "G  =  3.0_dp * f14 * f24 * f34 / (e4-e0)"                 
          #      G  =  f14 * f24 * f34 / (ee_c3[3]-mu)

          #      weights_aux[0,case3] =  f14 * G                     
          #      weights_aux[1,case3] =  f24 * G                     
          #      weights_aux[2,case3] =  f34 * G                     
          #      weights_aux[3,case3] =  (3.0 - f14 - f24 - f34 ) * G

          # except Exception,e: print e

          
#          weights[:,i] += weights_aux.T/6
          weights[case1,i] = (np.dot(wlsm[:4],weights_aux[:,case1])/6).T
          weights[case2,i] = (np.dot(wlsm[:4],weights_aux[:,case2])/6).T
          weights[case3,i] = (np.dot(wlsm[:4],weights_aux[:,case3])/6).T
#          for ik in range(weights_aux.shape[1]):
#               if not np.all(weights_aux[:,ik]==0):
#                    weights[ik] += np.sum(np.dot(wlsm,weights_aux[:,ik]))



     inds = np.where(np.any(weights!=0.0,axis=(1,2)))[0]
#     print np.sum(weights[inds],axis=(1,2))
#     raise SystemExit


     return inds,weights[inds]
   
              


#################################################################################
#################################################################################
#################################################################################                 

def get_tet_E(E_k):

     nk1,nk2,nk3 = E_k.shape

     E_k_ot = np.zeros((nk1,nk2,nk3,6,4),dtype=float,order="C")

     # 0,0,0 
     E_k_ot[:,:,:,0,0] = np.copy(E_k)
     E_k_ot[:,:,:,1,0] = np.copy(E_k)

     # 0,1,0 
     en_temp = np.ascontiguousarray(np.roll(E_k,-1,axis=1))
     E_k_ot[:,:,:,0,2] = en_temp
     E_k_ot[:,:,:,1,1] = en_temp
     E_k_ot[:,:,:,2,0] = en_temp
     E_k_ot[:,:,:,3,0] = en_temp
     E_k_ot[:,:,:,4,0] = en_temp
     E_k_ot[:,:,:,5,1] = en_temp

     # 1,1,0 
     en_temp = np.ascontiguousarray(np.roll(en_temp,-1,axis=0))
     E_k_ot[:,:,:,4,1] = en_temp
     E_k_ot[:,:,:,5,2] = en_temp

     # 1,0,0 
     en_temp = np.ascontiguousarray(np.roll(E_k,-1,axis=0))
     E_k_ot[:,:,:,0,1] = en_temp
     E_k_ot[:,:,:,5,0] = en_temp

     # 1,0,1 
     en_temp = np.ascontiguousarray(np.roll(en_temp,-1,axis=2))
     E_k_ot[:,:,:,0,3] = en_temp
     E_k_ot[:,:,:,1,3] = en_temp
     E_k_ot[:,:,:,2,2] = en_temp
     E_k_ot[:,:,:,3,1] = en_temp
     E_k_ot[:,:,:,4,2] = en_temp
     E_k_ot[:,:,:,5,3] = en_temp

     # 0,0,1 
     en_temp = np.ascontiguousarray(np.roll(E_k,-1,axis=2))
     E_k_ot[:,:,:,1,2] = en_temp
     E_k_ot[:,:,:,2,1] = en_temp

     # 0,1,1 
     en_temp = np.ascontiguousarray(np.roll(en_temp,-1,axis=1))
     E_k_ot[:,:,:,2,3] = en_temp
     E_k_ot[:,:,:,3,2] = en_temp

     # 1,1,1 
     en_temp = np.ascontiguousarray(np.roll(en_temp,-1,axis=0))
     E_k_ot[:,:,:,3,3] = en_temp
     E_k_ot[:,:,:,4,3] = en_temp

     return E_k_ot

#################################################################################
#################################################################################
#################################################################################

def setup_lin_tet(E_k,nk1,nk2,nk3,npool):

     E_k_copy = gather_scatter(np.ascontiguousarray(E_k),1,npool)

     _,nbnd,nspin = E_k_copy.shape

     E_k_copy = np.reshape(E_k_copy,(nk1,nk2,nk3,nbnd,nspin),order="C")

     edge_en = np.zeros((nbnd,nk1*nk2*nk3,nspin,6,4),order="C")     

     for ispin in range(nspin):
          for n in range(nbnd):
               edge_en_temp = get_tet_E(np.ascontiguousarray(E_k_copy[:,:,:,n,ispin]))
               edge_en[n,:,ispin] = np.reshape(edge_en_temp,(nk1*nk2*nk3,6,4),order="C")

     edge_en = gather_scatter(edge_en,1,npool)

     edge_en = np.sort(edge_en,axis=4)

     return [np.ascontiguousarray(np.swapaxes(edge_en,0,1)),nk1,nk2,nk3,npool,None]




def get_tet_I(integrand,ee_ind,npool,nk1,nk2,nk3):

     I_copy = gather_scatter(np.ascontiguousarray(integrand),1,npool)

     _,nbnd = I_copy.shape

     I_copy = np.reshape(I_copy,(nk1,nk2,nk3,nbnd),order="C")


     I_tet = np.zeros((nbnd,nk1*nk2*nk3,6,4),order="C")     
     for n in range(nbnd):
          I_tet_temp = get_tet_E(I_copy[:,:,:,n])
          I_tet[n] = np.reshape(I_tet_temp,(nk1*nk2*nk3,6,4),order="C")


     I_tet = gather_scatter(I_tet,1,npool)

     I_tet = np.swapaxes(I_tet,0,1)

     for k in range(I_tet.shape[0]):
          for nb in range(I_tet.shape[1]):
               for tet in range(I_tet.shape[2]):
                    ind = ee_ind[k,nb,tet]
                    I_tet[k,nb,tet,:] = I_tet[k,nb,tet,ind]
     return I_tet

     
