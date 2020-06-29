import AFLOWpi
import numpy as np
from tempfile import NamedTemporaryFile
import xml.etree.cElementTree as ET
import os

def resolve_old_qe_l(sd,onamd,slist):


    for sp in sd.keys(): 
        if sd[sp]==0:
            target="S"
        if sd[sp]==1:
            target="P"
        if sd[sp]==2:
            target="D"
        if sd[sp]==3:
            target="F"
        
        spec_orb=onamd[np.where(sp==slist)[0]]
        for i in spec_orb:
            if i[1]==target:
                sd[sp]=i
                break


    return sd

############################################################################################
############################################################################################
############################################################################################


def read_old_QE_xml( fn):

    spec_dict={}
    pp_dict={}

    for event,elem in ET.iterparse(fn,events=('start','end')):
        if event == 'end':

            # number of atoms
            if elem.tag == 'IONS':
                species_list=[]
                nspecs = int(elem.findall("NUMBER_OF_SPECIES")[0].text.split()[0])
                for n in range(nspecs):
                    string = "SPECIE."+str(n+1)
                    species = elem.findall(string+"/ATOM_TYPE")[0].text.split()[0]
                    pseudos = elem.findall(string+"/PSEUDO")[0].text.split()[0]
                    pp_dict[species]=pseudos
                    species_list.append(species)


                atoms = []
                natoms = int(elem.findall("NUMBER_OF_ATOMS")[0].text.split()[0])
                tau = np.zeros((natoms,3), dtype=float)
                for n in range(natoms):
                    string = "ATOM."+str(n+1)
                    aux = elem.findall(string)[0]
                    atoms.append(aux.attrib['SPECIES'][:-1])


            if elem.tag == 'EXCHANGE_CORRELATION':
                HL = elem.findall("HUBBARD_L")[0].text.split()
                for i in range(len(HL)):   
                    spec_dict[species_list[i]]=int(HL[i])


    atoms=np.array(atoms)

    return spec_dict,pp_dict,atoms



############################################################################################
############################################################################################
############################################################################################


def read_QE_xml( fn):

    spec_dict={}
    pp_dict={}

    for event,elem in ET.iterparse(fn,events=('start','end')):
        if event == 'end':
            try:
                for ohu in elem.findall("dft/dftU/Hubbard_U"):
                
                    lab=ohu.attrib['label'].upper()
                    spe=ohu.attrib['specie']
                    spec_dict[spe]=lab
            except: pass



            lspecies = elem.findall("atomic_species/species")
            for n in lspecies:
                pp_dict[n.attrib['name']]=n.findall('pseudo_file')[0].text


            if elem.tag == "output":
                atoms = []
                natoms = int(elem.findall("atomic_structure")[0].attrib['nat'])
                tau = np.zeros((natoms,3), dtype=float)
                latoms = elem.findall("atomic_structure/atomic_positions/atom")
                for n in range(natoms):
                    atoms.append(latoms[n].attrib['name'])


    return spec_dict,pp_dict,atoms



############################################################################################
############################################################################################
############################################################################################

def read_shell ( workpath,savedir,species,atoms,    spin_orb=False):
    # reads in shelks from pseudo files





    # Get Shells for each species
    sdict = {}
    jchid = {}
    onamd = {}
    jchia = None

    for k,v in species.items():
      sdict[k],jchid[k],onamd[k] = read_pseudopotential(os.path.join(workpath,savedir,v))

    #double the l=0 if spin orbit
    if spin_orb:
        for s,p in sdict.items():
            tmp_list=[]
            tmp_list_l=[]
            for o in range(len(p)):
                tmp_list.append(p[o])
                tmp_list_l.append(onamd[s][o])
                # if l=0 include it twice
                if p[o]==0 or len(jchid[s])==0:
                    tmp_list.append(p[o])
                    tmp_list_l.append(onamd[s][o])

            sdict[s] = np.array(tmp_list)
            onamd[s] = np.array(tmp_list_l)

            # when using scalar rel pseido with spin orb..
            if len(jchid[s])==0:
                tmp=[]
                for o in sdict[s][::2]:
                    if o==0:
                        tmp.extend([0.5])
                    if o==1:
                        tmp.extend([0.5,1.5])
                    if o==2:
                        tmp.extend([1.5,2.5])
                    if o==3:
                        tmp.extend([2.5,3.5])
                jchid[s]=np.array(tmp)


        jchia = np.hstack([jchid[a] for a in atoms])


    # index of which orbitals belong to which atom in the basis
    a_index = np.hstack([[a]*np.sum((2*sdict[atoms[a]])+1) for a in range(len(atoms))])

    # value of l
    shell   = np.hstack([sdict[a] for a in atoms])


    by_lab=[]
    for a in range(len(atoms)):
        orbs=sdict[atoms[a]]
        for o in range(len(orbs)):
            by_lab.extend(((onamd[atoms[a]][o]+" ")*(2*orbs[o]+1)).split()            )


    by_lab=np.array(by_lab)

    return shell,a_index,jchia,by_lab

############################################################################################
############################################################################################
############################################################################################

def read_pseudopotential ( fpp ):
  '''
  Reads a psuedopotential file to determine the included shells and occupations.

  Arguments:
      fnscf (string) - Filename of the pseudopotential, copied to the .save directory

  Returns:
      sh, nl (lists) - sh is a list of orbitals (s-0, p-1, d-2, etc)
                       nl is a list of occupations at each site
      sh and nl are representative of one atom only
  '''

  import numpy as np
  import xml.etree.cElementTree as ET
  import re

  sh = []
  on = []
  # fully rel case
  jchi=[]

  # clean xnl before reading
  with open(fpp) as ifo:
      temp_str=ifo.read()

  temp_str = re.sub('&',' ',temp_str)
  f = NamedTemporaryFile(mode='w',delete=True)
  f.write(temp_str)

  try:
      iterator_obj = ET.iterparse(f.name,events=('start','end'))
      iterator     = iter(iterator_obj)
      event,root   = next(iterator)

      for event,elem in iterator:        
          try:
              for i in elem.findall("PP_PSWFC/"):
                  sh.append(int(i.attrib['l']))
                  on.append(i.attrib['label'])
          except Exception as e:
              pass
          for i in elem.findall("PP_SPIN_ORB/"):
              try:
                  jchi.append(float(i.attrib["jchi"]))
              except: pass
              
      jchi = np.array(jchi)
      sh   = np.array(sh)

  except Exception as e:

      with open(fpp) as ifo:
          ifs=ifo.read()
      res=re.findall("(.*)\s*Wavefunction",ifs)[1:]      
      sh=np.array(list(map(int,list([x.split()[1] for x in res]))))
      on=np.array(list([x.split()[0] for x in res]))

  return sh,jchi,on

############################################################################################
############################################################################################
############################################################################################


def read_U_orbs(workpath,savedir):

    old_qe=False
    try:
        dfn=os.path.join(workpath,savedir,"data-file-schema.xml")
        sd,ppd,atoms=read_QE_xml(dfn)
    except:
        old_qe=True
        dfn=os.path.join(workpath,savedir,"data-file.xml")
        sd,ppd,atoms=read_old_QE_xml(dfn)

    atoms=np.array(atoms)

    shell,a_index,_,onamd = read_shell ( workpath,savedir,ppd,atoms,spin_orb=False)
    slist=atoms[a_index]

    if old_qe:
        sd=resolve_old_qe_l(sd,onamd,slist)

    orb_dict={}
    orb_dict_red={}
    for sp in ppd.keys(): 
        allorb=np.where(np.logical_and(sp==slist,sd[sp]==onamd))[0].tolist()

        orb_dict[sp]=allorb
        acount=len(np.where(atoms==sp)[0])
        rac=int(len(allorb)/acount)

        orb_dict_red[sp]=allorb[:rac]

    return orb_dict,orb_dict_red
