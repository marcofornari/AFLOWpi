import AFLOWpi
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import colors
from matplotlib import cm
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit


def radial_grid(origin=[0.0,0.0,0.0],nk_r=3,nk_theta=30,nk_phi=30,phi_range=[0.0,2.0*np.pi],theta_range=[0.0,np.pi],r_max=0.05):


    grouped_radial=False
    #add one because origin doesn't count
    nk_r_plus_one=nk_r+1
    r_range=[0.0,r_max]

    nk_theta_arr = np.linspace(theta_range[0],theta_range[1],num=nk_theta,     endpoint=False)[1:]
    nk_phi_arr   = np.linspace(  phi_range[0],  phi_range[1],num=nk_phi+1,     endpoint=True)[:-1]
    r_grid   = np.linspace(    r_range[0],    r_range[1],num=nk_r_plus_one,endpoint=True)
    r_grid=r_grid
    nk_r_arr = np.linspace(  r_range[0],    r_range[1],num=nk_r_plus_one,endpoint=True)

    nk_r_arr=nk_r_arr[1:]

    len_ar = (nk_r*(nk_theta-1)*(nk_phi))+(nk_r+1)*2-nk_r-1
    nk_str = np.zeros((len_ar,3),)
    temp_rad  = np.array(np.meshgrid(nk_theta_arr,nk_r_arr,nk_phi_arr,)).T.reshape(-1,3)

    #swap axes so it goes phi,theta,R
    temp_rad[:,[2, 1]]=temp_rad[:,[1, 2]]
    temp_rad[:,[0, 1]]=temp_rad[:,[1, 0]]

    #the rest is the radial points
    nk_str[:temp_rad.shape[0]]=temp_rad
    theta_zero=np.zeros(((nk_r)+1,3))

    theta_zero[:,2] = r_grid

    #do R=0 as last entry
    nk_str[-theta_zero.shape[0]:] =  theta_zero 

    #return in units cartesian 2pi/alat

    return radsphere_to_cart(nk_str,origin=origin)


def define_vec_direction(E_k,nk_r,nk_phi,nk_theta):

    #reorganize the flattened eigenvalues
    #into rays coming from the origin to 
    #do the fit for a given theta,phi along R
    
    xy_zero=E_k[-((nk_r+1)):]


    xy_zero_center_index=0
    origin=xy_zero[xy_zero_center_index]

    theta_zero = xy_zero
    no_xy_zero=E_k[:-(2*(nk_r+1)-1)]

    nbnd=E_k.shape[1]
    ray =np.zeros(( ((nk_phi)*(nk_theta-1)) ,(nk_r+1),nbnd  ),dtype=E_k.dtype)

    for n in range(no_xy_zero.shape[0]/nk_r):
        ray[n+1][1:] = no_xy_zero[n*nk_r:(n+1)*nk_r]
        ray[n+1][0]  = origin
        
    ray[0]=theta_zero
#    ray[1:]=np.roll(ray[1:],nk_r/2,axis=1)

    return ray


def radsphere_to_cart(in_radial,origin=np.array([0.0,0.0,0.0]),invert_R=False):
    if in_radial.shape[1]==3:
        R=in_radial[:,2]
        in_cart=np.zeros(in_radial.shape)
    else:
        R=np.ones(in_radial.shape[0],dtype=np.float64)
        if invert_R==True:
            R*=-1.0
        in_cart=np.zeros((in_radial.shape[0],3))

    #X = r*sin(theta)*cos(phi)
    in_cart[:,0] = R*np.sin(in_radial[:,1])*np.cos(in_radial[:,0])
    #Y = r*sin(theta)*sin(phi)
    in_cart[:,1] = R*np.sin(in_radial[:,1])*np.sin(in_radial[:,0])
    #Z = r*cos(theta)
    in_cart[:,2] = R*np.cos(in_radial[:,1])

    return in_cart+origin

def cart_to_radsphere(in_cart,origin=np.array([0.0,0.0,0.0])):
    in_cart[:,:3]-=origin
    in_radial=np.zeros(in_cart.shape)
    eps=1.e-12
    try:
        #r     = x^2+y^2+z^2)^(0.5)
        in_radial[:,2] = np.sqrt(in_cart[:,0]**2+in_cart[:,1]**2+in_cart[:,2]**2) 
        #theta = arccos(z/r)
        in_radial[:,1] = np.arccos((in_cart[:,2]+eps)/(in_radial[:,2]+eps))%(np.pi)
        #phi   = arctan(y/x)
        in_radial[:,0] = np.arctan2(in_cart[:,1],in_cart[:,0])%(2.0*np.pi)
    except:
        pass

    in_cart[:,:3]+=origin

    return in_radial


def _read_matdyn_out(oneCalc,ID):


    inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])




    return  np.loadtxt("%s.phVEL.gp"%ID,usecols=(1,2,3))

def _write_matdyn_in(oneCalc,ID,in_cart_flat):


    inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

    masses=[]
    species_string = inputDict['ATOMIC_SPECIES']['__content__']
    
    split_species_string = species_string.split('\n')
    for i in split_species_string:
        try:
            masses.append(float(i.split()[1]))
        except:
            pass

    amass_str=''
    for i in range(len(masses)):
        amass_str+='   amass(%d) = %f,\n'% (i+1,masses[i])



    matdyn_in = open("%s_matdyn_phVEL.in"%ID,"w")
    matdyn_template='''
 &input
   asr='crystal',
%s

   flfrc='%s.fc',
   flfrq='%s.phVEL', 


   nosym=.TRUE.
   q_in_cryst_coord = .false. 
   fd=.true.
   na_ifc=.true.
   q_in_band_form=.false.,
   eigen_similarity=.false.
 /
%s
'''%(amass_str,ID,ID,in_cart_flat.shape[0])


    matdyn_in.write(matdyn_template)
    np.savetxt(matdyn_in,in_cart_flat)
    matdyn_in.close()


def do_sound_velocity(__submitNodeName__,oneCalc,ID,dk_theta=0.1,dk_phi=0.2,dk_r=0.0125,r_max=0.05,theta_range=[-np.pi/2.0,np.pi/2.0],phi_range=[0.0,2.0*np.pi],origin=[0.0,0.0,0.0],nspin=1,kpi=0,read_S=False,shift=0.0,run_matdyn=True):



    theta_range= [0.0,np.pi]
    phi_range  = [0.0,np.pi]

    nk_theta=int(np.abs(np.diff(theta_range))[0]/dk_theta)
    nk_phi=int(np.abs(np.diff(phi_range))[0]/dk_phi)
    nk_r=int(r_max/dk_r)


    nk_theta =  30
    nk_phi   =  30
    nk_r     =  8


    # Define k-point mesh for radial grid (flattened)
    in_cart_flat=AFLOWpi.run.radial_grid(origin=origin,nk_r=nk_r,nk_theta=nk_theta,nk_phi=nk_phi,
                                         theta_range=theta_range,phi_range=phi_range,r_max=r_max,)


    
    AFLOWpi.run._write_matdyn_in(oneCalc,ID,in_cart_flat)

    execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix_serial")
    if run_matdyn:
        AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s_matdyn_phVEL'%ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='custom',execPath='./matdyn.x' )



    nkpi=in_cart_flat.shape[0]
    kq=np.zeros((3,nkpi))
    kq=np.copy(in_cart_flat.T)



    alat = float(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])["&system"]["celldm(1)"])
#    cell = AFLOWpi.retr._cellMatrixToString(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])["CELL_PARAMETERS"]["__content__"]))
#    cell*=alat 
    

    #2pi/alat units
    rad_points = np.linspace(0.0,r_max,nk_r+1,endpoint=True)

    #alat to angstrom     
    alat *= 0.529177249
    #angstrom to m
    alat /= 1.e10

    conv = 1.0
    #convert cm-1 to THz
    conv *= 0.0299792458
    #convert THz to Hz
    conv *= 1.e12

    # 1/s/m = m/s
    conv *= alat
#    print conv
    #load freq    
    W_qp = _read_matdyn_out(oneCalc,ID)
    #reshape to angular directions
    W_qp_radial = AFLOWpi.run.define_vec_direction(W_qp,nk_r,nk_phi,nk_theta)





    fit_func = lambda x,a: a*x

    sol = np.zeros(W_qp_radial.shape[0])
    res=[]
    print((W_qp_radial[:,:,0]))


    for branch in range(3):



        sol = np.polyfit(rad_points,W_qp_radial[:,:,branch].T,4)

        print((sol.shape))
        res.append(np.mean(sol[3])*conv)




    return res

    in_cart_flat=radial_grid(origin=origin,nk_r=1,nk_theta=nk_theta,nk_phi=nk_phi,
                             theta_range=theta_range,phi_range=phi_range,r_max=r_max,)







    pos_R_TA = cart_to_radsphere(in_cart_flat[:-2])
    pos_R_TA[:,-1] = np.ravel(sol[0,:,-2])
    neg_R_TA = np.copy(pos_R_TA)
    neg_R_TA[:,-1] *= -1.0
    pos_R_TA = radsphere_to_cart(pos_R_TA)
    neg_R_TA = radsphere_to_cart(neg_R_TA)
    pos_R_TA= pos_R_TA.reshape(nk_phi,nk_theta-1,3)
    neg_R_TA= neg_R_TA.reshape(nk_phi,nk_theta-1,3)
    pos_R_TA = np.concatenate([pos_R_TA,neg_R_TA],axis=0)

    pos_R_TA1 = cart_to_radsphere(in_cart_flat[:-2])
    pos_R_TA1[:,-1] = np.ravel(sol[1,:,-2])
    neg_R_TA1 = np.copy(pos_R_TA1)
    neg_R_TA1[:,-1] *= -1.0
    pos_R_TA1 = radsphere_to_cart(pos_R_TA1)
    neg_R_TA1 = radsphere_to_cart(neg_R_TA1)
    pos_R_TA1= pos_R_TA1.reshape(nk_phi,nk_theta-1,3)
    neg_R_TA1= neg_R_TA1.reshape(nk_phi,nk_theta-1,3)
    pos_R_TA1 = np.concatenate([pos_R_TA1,neg_R_TA1],axis=0)

    pos_R_LA = cart_to_radsphere(in_cart_flat[:-2])
    pos_R_LA[:,-1] = np.ravel(sol[2,:,-2])
    neg_R_LA = np.copy(pos_R_LA)
    neg_R_LA[:,-1] *= -1.0
    pos_R_LA = radsphere_to_cart(pos_R_LA)
    neg_R_LA = radsphere_to_cart(neg_R_LA)
    pos_R_LA= pos_R_LA.reshape(nk_phi,nk_theta-1,3)
    neg_R_LA= neg_R_LA.reshape(nk_phi,nk_theta-1,3)
    pos_R_LA = np.concatenate([pos_R_LA,neg_R_LA],axis=0)

    # """reorganize eigenvalues for K point in radial grid"""
    # #into rays from the origin of the radial grid

    # plot_arr = np.zeros((nk_phi*nk_theta,3))
    # plot_arr[nk_phi:] = W_qp_radial[:-2]
    # plot_arr[:nk_phi] = W_qp_radial[-2]

    # plot_arr2 = np.zeros((nk_phi*nk_theta*2,3))

    # plot_arr2[(nk_phi*nk_theta):] = plot_arr
    # plot_arr2[:(nk_phi*nk_theta)] =-plot_arr
    # plot_arr2 = plot_arr2.reshape(nk_phi*2,nk_theta,3)




#    return res

    band_min1=band_min2=band_min3=0.0



    fig  = plt.figure(        figsize=(12,4)) 
    gs = gridspec.GridSpec(22,3)


    ax1 = plt.subplot(gs[:20, 0])
    ax2 = plt.subplot(gs[:20, 1])
    ax3 = plt.subplot(gs[:20, 2])        


# #######
    ax1 = plt.subplot(gs[:20,0], projection='3d')
    ax1.view_init(0,0)
    ax2 = plt.subplot(gs[:20,1], projection='3d')
    ax2.view_init(0,0)
    ax3 = plt.subplot(gs[:20,2], projection='3d')
    ax3.view_init(0,0)



#     ax1.grid(b=False)
#     ax1.grid(b=False)
#     ax1.grid(b=False)
#     ax2.grid(b=False)
#     ax2.grid(b=False)
#     ax2.grid(b=False)
#     ax3.grid(b=False)
#     ax3.grid(b=False)
#     ax3.grid(b=False)

#     phi   = np.zeros((nk_phi*2*nk_theta))
#     theta = np.zeros((nk_phi*2*nk_theta))

#     angles=cart_to_radsphere(in_cart_flat)
#     print angles
#     plot_arr[nk_phi:] = angles[:-2]
#     plot_arr[:nk_phi] = angles[-2]

#     plot_arr2 = np.zeros((nk_phi*nk_theta*2,3))

#     plot_arr2[(nk_phi*nk_theta):] = plot_arr
#     plot_arr2[:(nk_phi*nk_theta)] =-plot_arr


#     plot_arr2[(nk_phi*nk_theta):] = plot_arr
#     plot_arr2[:(nk_phi*nk_theta)] =-plot_arr
#     print plot_arr2
#     plot_arr2 = plot_arr2.reshape(nk_phi*2,nk_theta,3)

#     plot_arr2 = plot_arr2.reshape(nk_phi*2,nk_theta,3)

#     phi   = np.reshape(phi  ,(nk_phi*2,nk_theta))
#     theta = np.reshape(theta,(nk_phi*2,nk_theta))


#     in_cart = radsphere_to_cart(in_radial,origin=np.array([0.0,0.0,0.0]),invert_R=False)

#     Z1 = radsphere_to_cart(in_radial,origin=np.array([0.0,0.0,0.0]),invert_R=False)
#     Z2 = radsphere_to_cart(in_radial,origin=np.array([0.0,0.0,0.0]),invert_R=False)
#     Z3 = radsphere_to_cart(in_radial,origin=np.array([0.0,0.0,0.0]),invert_R=False)
# ########
#    plot_arr = plot_arr.reshape(nk_phi,nk_theta,3)

    # ax1 = plt.subplot(gs[:20, 0])
    # ax2 = plt.subplot(gs[:20, 1])
    # ax3 = plt.subplot(gs[:20, 2])        

    # ax1.grid(b=False)
    # ax1.grid(b=False)
    # ax1.grid(b=False)
    # ax2.grid(b=False)
    # ax2.grid(b=False)
    # ax2.grid(b=False)
    # ax3.grid(b=False)
    # ax3.grid(b=False)
    # ax3.grid(b=False)

    cmap0=cm.Purples
    cmap1=cm.Purples
    cmap2=cm.Purples


    ax1.set_title(r'Ta')
    ax2.set_title(r"Ta'")
    ax3.set_title(r'La')
    ax1.set_xticks([])
    ax2.set_xticks([])
    ax3.set_xticks([])
    ax1.set_yticks([])
    ax2.set_yticks([])
    ax3.set_yticks([])
    ax1.set_xticks([])
    ax2.set_xticks([])
    ax3.set_xticks([])



    rta = np.sqrt(pos_R_TA[:,:,0]**2+pos_R_TA[:,:,1]**2+pos_R_TA[:,:,2]**2) 
    cm_af = plt.get_cmap('afmhot')
#    print rta
    rrs=np.int32(np.ceil(255*(rta+np.abs(rta).max()) / (2.0*np.abs(rta).max())))
    face_col = cm.coolwarm(rrs)                                                                                                              
    norm = matplotlib.colors.Normalize(vmin=np.amin(sol[0,:,-2]), vmax=np.amax(sol[0,:,-2]))                          
    num=int(np.ceil((np.amax(rta)-np.amin(rta))/(2.0*np.amax(sol[0,:,-2])*255)))                 
    # dt = np.linspace(np.amin(R_f_tp),np.amax(R_f_tp),                                
    #                  num=num,endpoint=True)                                         
    # if dt[0]<0.0 and dt[-1]>0.0:                                                                                                           
    #      dt -= dt[np.argmin(np.abs(dt))]                                           




    #cmap=cm_af,alpha=1.0,facecolors=face_col,
                                 
    surf1 = ax1.plot_surface(pos_R_TA[:,:,0],pos_R_TA[:,:,1] , pos_R_TA[:,:,2], rstride=1, cstride=1, cmap=cm_af,alpha=1.0,facecolors=face_col,
                             linewidth=0, antialiased=True,)
    surf1 = ax2.plot_surface(pos_R_TA1[:,:,0],pos_R_TA1[:,:,1] ,pos_R_TA1[:,:,2],rstride=1,cstride=1,cmap=cm_af,alpha=1.0,facecolors=face_col,
                             linewidth=0, antialiased=True,)
    surf1 = ax3.plot_surface(pos_R_LA[:,:,0],pos_R_LA[:,:,1] , pos_R_LA[:,:,2], rstride=1, cstride=1,cmap=cm_af,alpha=1.0,facecolors=face_col,
                             linewidth=0, antialiased=True,)



    ax1.set_xlim(-np.amax(pos_R_TA),np.amax(pos_R_TA))
    ax2.set_xlim(-np.amax(pos_R_TA1),np.amax(pos_R_TA1))
    ax3.set_xlim(-np.amax(pos_R_LA),np.amax(pos_R_LA))
    ax1.set_ylim(-np.amax(pos_R_TA),np.amax(pos_R_TA))
    ax2.set_ylim(-np.amax(pos_R_TA1),np.amax(pos_R_TA1))
    ax3.set_ylim(-np.amax(pos_R_LA),np.amax(pos_R_LA))
    ax1.set_zlim(-np.amax(pos_R_TA),np.amax(pos_R_TA))
    ax2.set_zlim(-np.amax(pos_R_TA1),np.amax(pos_R_TA1))
    ax3.set_zlim(-np.amax(pos_R_LA),np.amax(pos_R_LA))

    ax1.set_aspect('equal')        
    ax2.set_aspect('equal')        
    ax3.set_aspect('equal')     
   
    plt.savefig('sound_vel.pdf')
    plt.close()    

#    return np.mean(W_qp_radial,axis=0)

#    surf1 = ax2.plot_surface(neg_R[:,:,0],neg_R[:,:,1] , neg_R[:,:,2], rstride=1, cstride=1, cmap=cm.RdGy,
#                             linewidth=0, antialiased=True,)

    band_min1=np.amin(W_qp_radial[:,0])
    band_max1=np.amax(W_qp_radial[:,0])
#    surf1 = ax1.imshow(plot_arr[:,:,0],cmap=cmap1,interpolation="None",
#                       vmin=band_min1,vmax=band_max1)

    band_min2=np.amin(W_qp_radial[:,1])
    band_max2=np.amax(W_qp_radial[:,1])
#    surf1 = ax2.imshow(plot_arr[:,:,1],cmap=cmap1,interpolation="None",
#                       vmin=band_min2,vmax=band_max2)

    band_min3=np.amin(W_qp_radial[:,2])
    band_max3=np.amax(W_qp_radial[:,2])
#    surf1 = ax3.imshow(plot_arr[:,:,2],cmap=cmap1,interpolation="None",
#                       vmin=band_min3,vmax=band_max3)


 #    ax4 = plt.subplot(gs[21, 0])
#     ax5 = plt.subplot(gs[21, 1])
#     ax6 = plt.subplot(gs[21, 2])        
# /
#    cb1 = matplotlib.colorbar.ColorbarBase(ax4, cmap=cmap1,norm=plt.Normalize(band_min1,band_max1), orientation='horizontal')
#    cb2 = matplotlib.colorbar.ColorbarBase(ax5, cmap=cmap1,norm=plt.Normalize(band_min2,band_max2), orientation='horizontal')
#    cb3 = matplotlib.colorbar.ColorbarBase(ax6, cmap=cmap1,norm=plt.Normalize(band_min3,band_max3), orientation='horizontal')
   # cb2.ax.set_title(r'$\frac{d^{2}E(k_{r})}{dk_{r}^{2}}$',fontsize=22,va="bottom")




#        for i in xrange(W_qp_radial.shape[0]):
#            fit[i] = curve_fit(fit_func,rad_points,W_qp_radial[i,:,branch],maxfev=10000)
#        res.append(np.mean(fit))
#        print popt



