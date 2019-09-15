import re
import numpy as np
import os
import matplotlib
matplotlib.use('pdf')
from matplotlib import pylab
from matplotlib import pyplot as plt
from cycler import cycler


def _get_raman_spectrum(oneCalc,ID,e1=None,e2=None,axis="z",p=0.0,temp=300,w_range=None,freq_lines=False,data=False,plot=True):


    atpolre = re.compile(r"\d+\n\s*([\.\-\dE\+]+)\s*([\.\-\dE\+]+)\s*([\.\-\dE\+]+)\s*\n\s*([\.\-\dE\+]+)\s*([\.\-\dE\+]+)\s*([\.\-\dE\+]+)\s*\n\s*([\.\-\dE\+]+)\s*([\.\-\dE\+]+)\s*([\.\-\dE\+]+)\s*")

    with open(os.path.join(oneCalc["_AFLOWPI_FOLDER_"],"%s_raman.out"%ID)) as ifo:
        instr=ifo.read()



    instr = re.split("Raman alpha tensor",instr)[1]
    print(instr)
    instr=re.sub("-"," -",instr)
    alpha_str = atpolre.findall(instr)




    # print alpha_str
    
    # ###in the rare case the raman iutput is messed up###
    # for i in range(len(alpha_str)):
    #     alpha_str[i]=list(alpha_str[i])
    #     for j in range(len(alpha_str[i])):

    #         alpha_str[i][j]=re.sub("-"," -",alpha_str[i][j])
    #         print alpha_str[i][j]
    alpha =  np.array([list(map(float,i)) for i in alpha_str]).reshape((len(alpha_str),3,3))



    fn=os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"%s.phBAND.eig"%ID)
    with open(fn,"r") as fo:
        fs=fo.read()


    re_qp=re.compile("freq.*=.*=\s*([-.0-9]+).*\n((?:\s*\(\s*[-.0-9]+\s*[-.0-9]+\s*[-.0-9]+\s*[-.0-9]+\s*[-.0-9]+\s*[-.0-9]+\s*\)\s*\n)+)")
    #
    qs=re_qp.findall(fs)


    
    q_list = np.array([list(map(float,k.strip().split())) for k in re.findall("q\s*=(.*)\n",fs)])
    num_q  = len(q_list)


    by_eig=[]
    for i in qs:
        eig=[j for j in i[1].split('\n') if len(j.strip())!=0]
        qs=[list(map(float,k.split()[1:-1])) for k in eig]
        eig_val=[float(i[0])]

        by_eig.append(eig_val)

    gamma_w = np.ravel(np.array(by_eig[:alpha.shape[0]]))


    min_w = np.amin(gamma_w[3:])
    max_w = np.amax(gamma_w[3:])

    nsample=2000
    if w_range==None:
        w_sample = np.linspace(min_w*0.98+1,max_w*1.02,nsample,endpoint=True)
    else:
        w_sample = np.linspace(w_range[0],w_range[1],nsample,endpoint=True)



    
    angles = np.linspace(0,np.pi,19,endpoint=True)
    e_i = 1.0/((1.0+p)**0.5)*np.ones((angles.shape[0],3))
    e_r = 1.0/((1.0+p)**0.5)*np.ones((angles.shape[0],3))
    if axis=="x":
        ordering=[2,0,1]
    if axis=="y":
        ordering=[1,2,0]
    if axis=="z":
        ordering=[0,1,2]
    e_r[:,ordering[0]] *= -np.sin(angles)
    e_r[:,ordering[1]] *= np.cos(angles)
    e_r[:,ordering[2]] *= p
    e_i[:,ordering[0]] *= np.cos(angles)
    e_i[:,ordering[1]] *= np.sin(angles)
    e_i[:,ordering[2]] *= p



    res=np.zeros((w_sample.shape[0],angles.shape[0],3,3))

    edotAdote = np.zeros((3,3))

    temp *= 8.621738e-5
    eps=1.e-12

    for ang in range(angles.shape[0]):
        for n in range(3,alpha.shape[0]):
            edotAdote = (np.tensordot(e_i[ang][:,None],np.tensordot(alpha[n],e_r[ang],axes=([0],[0]))[:,None],axes=([1],[1])))**2
            for w in range(w_sample.shape[0]):                
                res[w,ang] += edotAdote*np.nan_to_num(((1.0/(np.exp((np.abs(gamma_w[n]-w_sample[w])*1.2398e-4)/temp)-1.0))+1.0)/(gamma_w[n]*1.2398e-4)) 

    rs1=res.shape[0]
    rs2=res.shape[1]

    res  = res.reshape((rs1*rs2,3,3))
    res /= np.amax(res,axis=0)+eps
    res  = res.reshape((rs1,rs2,3,3))

    angles*=180.0/np.pi

    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 15,
        }


    colors=["#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
            "#9467bd",#
            "#17becf",
            "#bcbd22",
            "#7f7f7f",
            "#e377c2",
            "#8c564b",
            "#d62728",
            "#2ca02c",
            "#ff7f0e",
            "#1f77b4",
]


    colors=list(map(matplotlib.colors.to_rgb,colors))

#    delta=0.1
#    for w in xrange(w_sample.shape[0]):                
#        (np.exp(-((ene-eig)/delta)**2)/delta)/np.sqrt(np.pi)
    directions = ["x","y","z"]
    for l in range(3):
        for m in range(3):
            if data:
                fn = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"RAMAN_%s%s_%s.dat"%(directions[l],directions[m],ID))
                np.savetxt(fn,np.concatenate((w_sample[:,None],res[:,:,l,m]),axis=1),fmt="% 12.6f")
            if plot:
                plt.figure(figsize=(9,12))
                if freq_lines:
                    for f in range(3,gamma_w.shape[0]):
                        plt.axvline(gamma_w[f],color="k",linestyle="--")
                for ang in range(angles.shape[0]):
                    plt.plot(w_sample,res[:,ang,l,m]+ang,color=colors[int(float(ang)%float(len(colors)))])
                    plt.text(w_sample[-1]-0.1*(w_sample[-1]-w_sample[0]), ang+0.1, r'%5s$^{\circ}$'%angles[ang], fontdict=font)
                plt.ylim(-0.5,angles.shape[0]+0.5)
                plt.xlim(np.amin(w_sample),np.amax(w_sample))
                plt.xlabel("raman shift (cm$^{-1}$)",fontsize=16)
                plt.ylabel("Intensity (a.u.)",fontsize=16)
    #            plt.title(r"$y(%s%s)\bar{y}$"%(directions[l],directions[m]),fontsize=18)
                plt.yticks([])
                fn = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],"RAMAN_%s%s_%s.pdf"%(directions[l],directions[m],ID))
                plt.savefig(fn)
                plt.close()

