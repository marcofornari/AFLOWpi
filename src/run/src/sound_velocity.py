import AFLOWpi
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import colors
from matplotlib import cm
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



def radial_grid(origin=[0.0,0.0,0.0],nk_r=3,nk_theta=30,nk_phi=30,phi_range=[0.0,2.0*np.pi],theta_range=[0.0,np.pi],r_max=0.05):


    grouped_radial=False
    #add one because origin doesn't count
    nk_r_plus_one=nk_r+1
    r_range=[-r_max,r_max]

    nk_theta_arr = np.linspace(theta_range[0],theta_range[1],num=nk_theta,     endpoint=False)[1:]
    nk_phi_arr   = np.linspace(  phi_range[0],  phi_range[1],num=nk_phi+1,     endpoint=True)[:-1]
    if grouped_radial==True:
        pos = np.linspace(    0.003125,  r_max,num=nk_r/2,endpoint=True)
        neg=np.flipud(-1.0*pos)
        r_grid=np.concatenate((neg,np.array([0.0]),pos))

        nk_r_arr = np.fft.ifftshift(r_grid)
    else:
        r_grid   = np.linspace(    r_range[0],    r_range[1],num=nk_r_plus_one,endpoint=True)
        nk_r_arr = np.fft.ifftshift(np.linspace(  r_range[0],    r_range[1],num=nk_r_plus_one,endpoint=True))

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
    xy_zero_center_index=xy_zero.shape[0]/2
    origin=xy_zero[xy_zero_center_index]

    theta_zero = xy_zero
    no_xy_zero=E_k[:-(2*(nk_r+1)-1)]

    nbnd=E_k.shape[1]
    ray =np.zeros(( ((nk_phi)*(nk_theta-1)) ,(nk_r+1),nbnd  ),dtype=E_k.dtype)

    for n in xrange(no_xy_zero.shape[0]/nk_r):
        ray[n+1][1:] = no_xy_zero[n*nk_r:(n+1)*nk_r]
        ray[n+1][0]  = origin

    ray[0]=theta_zero
    ray[1:]=np.roll(ray[1:],nk_r/2,axis=1)

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


    eigs = np.loadtxt("%s.phVEL.gp"%ID)

    return eigs[:,1:4]

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

    nk_theta =  100
    nk_phi   =  100
    nk_r     =  1

    # Define k-point mesh for radial grid (flattened)
    in_cart_flat=radial_grid(origin=origin,nk_r=nk_r,nk_theta=nk_theta,nk_phi=nk_phi,
                             theta_range=theta_range,phi_range=phi_range,r_max=r_max,)



    _write_matdyn_in(oneCalc,ID,in_cart_flat)

    execPrefix=AFLOWpi.prep._ConfigSectionMap("run","exec_prefix_serial")
    if run_matdyn:
        AFLOWpi.run._oneRun(__submitNodeName__,oneCalc,'%s_matdyn_phVEL'%ID,execPrefix=execPrefix,execPostfix='',engine='espresso',calcType='custom',execPath='./matdyn.x' )



    nkpi=in_cart_flat.shape[0]
    kq=np.zeros((3,nkpi))
    kq=np.copy(in_cart_flat.T)

    #get it in the right format for the fft
    
    E_kp = _read_matdyn_out(oneCalc,ID)

    #filter only the bands we want


    alat = float(AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])["&system"]["celldm(1)"])
    #alat to angstrom
    alat *= 0.529177249
    #convert to Hz
    E_kp*=0.0299792458*1.e12*alat

    E_kp_radial = E_kp/(r_max*1.e10)

    """reorganize eigenvalues for K point in radial grid"""
    #into rays from the origin of the radial grid

    plot_arr = np.zeros((nk_phi*nk_theta,3))
    plot_arr[nk_phi:] = E_kp_radial[:-2]
    plot_arr[:nk_phi] = E_kp_radial[-2]

    plot_arr2 = np.zeros((nk_phi*nk_theta*2,3))

    plot_arr2[(nk_phi*nk_theta):] = plot_arr
    plot_arr2[:(nk_phi*nk_theta)] =-plot_arr
    plot_arr2 = plot_arr2.reshape(nk_phi*2,nk_theta,3)



    


    band_min1=band_min2=band_min3=0.0



    fig  = plt.figure(        figsize=(12,4)) 
    gs = gridspec.GridSpec(22,3)


    ax1 = plt.subplot(gs[:20, 0])
    ax2 = plt.subplot(gs[:20, 1])
    ax3 = plt.subplot(gs[:20, 2])        


# #######
    ax1 = plt.subplot(gs[:20,0], projection='3d')
#     ax1.view_init(0, 0)
    ax2 = plt.subplot(gs[:20,1], projection='3d')
#     ax2.view_init(45,45)
    ax3 = plt.subplot(gs[:20,2], projection='3d')
#     ax3.view_init(0, 90)



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
    plot_arr = plot_arr.reshape(nk_phi,nk_theta,3)

    ax1 = plt.subplot(gs[:20, 0])
    ax2 = plt.subplot(gs[:20, 1])
    ax3 = plt.subplot(gs[:20, 2])        

    ax1.grid(b=False)
    ax1.grid(b=False)
    ax1.grid(b=False)
    ax2.grid(b=False)
    ax2.grid(b=False)
    ax2.grid(b=False)
    ax3.grid(b=False)
    ax3.grid(b=False)
    ax3.grid(b=False)

    cmap3=cm.Purples
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



    band_min1=np.amin(E_kp_radial[:,0])
    band_max1=np.amax(E_kp_radial[:,0])
#    surf1 = ax1.imshow(plot_arr[:,:,0],cmap=cmap1,interpolation="None",
#                       vmin=band_min1,vmax=band_max1)

    band_min2=np.amin(E_kp_radial[:,1])
    band_max2=np.amax(E_kp_radial[:,1])
#    surf1 = ax2.imshow(plot_arr[:,:,1],cmap=cmap1,interpolation="None",
#                       vmin=band_min2,vmax=band_max2)

    band_min3=np.amin(E_kp_radial[:,2])
    band_max3=np.amax(E_kp_radial[:,2])
#    surf1 = ax3.imshow(plot_arr[:,:,2],cmap=cmap1,interpolation="None",
#                       vmin=band_min3,vmax=band_max3)

    ax4 = plt.subplot(gs[21, 0])
    ax5 = plt.subplot(gs[21, 1])
    ax6 = plt.subplot(gs[21, 2])        


    cb1 = matplotlib.colorbar.ColorbarBase(ax4, cmap=cmap1,norm=plt.Normalize(band_min1,band_max1), orientation='horizontal')
    cb2 = matplotlib.colorbar.ColorbarBase(ax5, cmap=cmap1,norm=plt.Normalize(band_min2,band_max2), orientation='horizontal')
    cb3 = matplotlib.colorbar.ColorbarBase(ax6, cmap=cmap1,norm=plt.Normalize(band_min3,band_max3), orientation='horizontal')
   # cb2.ax.set_title(r'$\frac{d^{2}E(k_{r})}{dk_{r}^{2}}$',fontsize=22,va="bottom")

    plt.savefig('sound_vel.pdf')
    plt.close()    
    return np.mean(E_kp_radial,axis=0)
