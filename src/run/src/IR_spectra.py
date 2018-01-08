import re
import numpy as np
import AFLOWpi


def _get_IR_spectra(oneCalc,ID):
    np.set_printoptions(precision=6,suppress=True)
    dynmat_re = re.compile("\s*Dynamical\s*Matrix\s*in.*\n.*\n.*\n(?=([\s\d\-\.]+\n)+)\s*")

    with open ("%s.phBAND.dyn"%ID) as ifo:
        instr= ifo.read()
    dynstr = dynmat_re.findall(instr)[0]


    dynmat =  np.array([map(float,x.split()) for x in dynstr.split("\n") if (len(x.strip())!=0)*(len(x.split())>2)])[:,[0,2,4]]

    nfreq = int(np.sqrt(float(dynmat.shape[0]/3)))*3
    print nfreq
    dynmat = dynmat.reshape(nfreq,nfreq)


    natom=nfreq/3
    dynmat = dynmat.reshape(natom,nfreq,3)
    print dynmat



    born_str = AFLOWpi.run._pull_born_out(oneCalc,ID)
    born =  np.array([map(float,x.split()) for x in born_str.split("\n") if (len(x.strip())!=0)*(len(x.split())>2)])
    
    

    born = born.reshape((natom,3,3))
    IR =np.zeros((nfreq))
    print natom

    for nu in xrange(nfreq):
        temp=np.zeros(3)
        for na in xrange(natom):
            for jpol in xrange(3):
                temp += born[na,:,jpol]*dynmat[na,nu,jpol]
        IR[nu]+=2.0*np.sum(temp**2)
            

#            print IR[alpha]
#        print born[alpha].dot(dynmat[alpha,nu])


    print IR
