import AFLOWpi
import numpy as np
from collections import OrderedDict

def read_modes( fn, natoms ):

  q_pnts = OrderedDict()
  q_eigs = OrderedDict()
  f_cnt  = OrderedDict()

  # Open File
  with open(fn, 'r') as f:
    q_indx = 0
    l = f.readline()

    # Categorize Modes for Each q
    while l != '':
      ls = l.split()
      if len(ls) > 0:
        if ls[0] == 'q':
          qs = [float(ls[2]), float(ls[3]), float(ls[4])]
          f.readline()

          # Temporary Storage Arrays
          qk_tmp = []
          qp_tmp = []
          qe_tmp = []
          freqs = np.zeros((3*natoms), dtype=float)

          # Read the 3N Modes
          for i in xrange(3*natoms):
            ls = f.readline().split()

            freqs[i] = np.round(float(ls[-2]),4)
            eigs = np.zeros((natoms,3), dtype=complex)

            # Parse Complex Displacement for Each Atom
            for j in xrange(natoms):
              ls = f.readline().split()
              eigs[j,0] = complex(ls[1])+complex(ls[2]+'j')
              eigs[j,1] = complex(ls[3])+complex(ls[4]+'j')
              eigs[j,2] = complex(ls[5])+complex(ls[6]+'j')

            qk_tmp.append(str(q_indx)+'_'+str(i))
            qp_tmp.append(qs)
            qe_tmp.append(eigs)

          # Check for Degenerate Eigenvalues
          fkey = '_'.join([str(fq) for fq in freqs])
          if fkey not in f_cnt:
            f_cnt[fkey] = [1, q_indx]
            # Save Eigenvalues for his q Point
            for j in xrange(len(qk_tmp)):
              q_pnts[qk_tmp[j]] = qp_tmp[j]
              q_eigs[qk_tmp[j]] = qe_tmp[j]
          else:
            f_cnt[fkey][0] += 1

          q_indx += 1
      l = f.readline()

  # Sort Keys and Return
  spnts = sorted(q_pnts.keys())
  seigs = sorted(q_eigs.keys())
  assert(spnts == seigs)

  return [q_pnts, q_eigs, f_cnt]


def get_phase( R, q ):
  return np.exp(np.dot(2.0j*np.pi*q,R))

def displace_atomic_position( scdm, line, disp, q, ephase, d_frac=0.013 ):

  ls = line.split()

  new_pos = map(complex, ls[1:])

  phase = np.dot(q, new_pos)
  eiqr = complex(np.cos(phase), np.sin(phase))

  eiqr = get_phase

  for i in range(len(new_pos)):
    disp[i] = np.real(disp[i]*eiqr)
  disp_norm = np.linalg.norm(disp)
  if abs(disp_norm) > 1e-3:
    disp = disp/disp_norm

  for i in range(len(new_pos)):
    new_pos[i] += scdm[i] + d_frac*disp[i]

  # Build Array
  new_pos = [ls[0]] + [str('%.6f'%np.round(np.real(i),6)) for i in new_pos]

  return('  '.join(new_pos)+'\n')


def write_scf_files( odir, fscf, natoms, nr1, nr2, nr3, modes , cell ,dfrac = 0.01 ):

  #read in supercell input file
  with open(fscf) as ifo:
    ifs = ifo.read()

  #split supercell file input
  spli = AFLOWpi.retr._splitInput(ifs)

  pos  = AFLOWpi.retr._getPositions(ifs).A
  labs = AFLOWpi.retr._getPosLabels(ifs)

  pos_cart = np.copy(pos)
  pos_cart[:,0] *= nr1
  pos_cart[:,1] *= nr2
  pos_cart[:,2] *= nr3
  pos_cart = pos_cart.dot(cell)

  super_cell = np.copy(cell)
  super_cell[0] *= nr1
  super_cell[1] *= nr2
  super_cell[2] *= nr3

  nat_orig = pos_cart.shape[0]/(nr1*nr2*nr3)


  # Write New scf Files
  for key in modes[0].keys():

    ephase = []
    qp = np.asarray(modes[0][key])
    # calculate phase factor for atoms in supercell
    for i in range(pos_cart.shape[0]):
      ephase.append(get_phase(pos_cart[i], qp))

    ephase = np.asarray(ephase,dtype=complex)
    
    ephase /= np.tile(ephase[:nat_orig],nr1*nr2*nr3).T
#    print ephase.real
    eig_disp = np.repeat(np.reshape(modes[1][key],(nat_orig,3),order='C'),nr1*nr2*nr3,axis=0)


    eig_disp[:,0] *= ephase
    eig_disp[:,1] *= ephase
    eig_disp[:,2] *= ephase

#    print eig_disp



    # convert eigdisplacements in 2pi/alat to crystal cords
    

    eig_disp_real = (eig_disp.real*dfrac).dot(np.linalg.inv(super_cell))
    pos += eig_disp_real
    


    #add new displaced positions to input
    spli['ATOMIC_POSITIONS']['__content__'] = AFLOWpi.retr._joinMatrixLabels(labs,pos)
    #join split input
    disp_input_str = AFLOWpi.retr._joinInput(spli)

    #write displaced supercell input
    with open(odir+key+'.in', 'w') as f:
      f.write(disp_input_str)


  # Write a Legend with KeyValue Pairs - { FilePrefix : qArray }
  with open(odir+'q_legend.txt', 'w') as f:
    f.writelines(['%s : %s\n'%(k,v) for k,v in modes[0].iteritems()])

  # Write a Legend with KeyValue Pairs - { Frequencies : [Degeneracy, qIndex] }
  with open(odir+'degeneracy_legend.txt', 'w') as f:
    f.writelines(['%s %s\n'%(v[1], v[0]) for k,v in modes[2].iteritems()])

