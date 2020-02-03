import AFLOWpi
import os

def _write_qe_eps_in(oneCalc,ID,intersmear=0.136,wmin=0.0,wmax=30.0,nw=600,smeartype='gauss',intrasmear=0.0,metalcalc=False,jdos=False,offdiag=False,occ=False):
    

    calctype='eps'
    prefix = oneCalc['_AFLOWPI_PREFIX_']
    exten='eps'
    AFLOWpi.run._gen_eps_in(oneCalc,ID,calctype,exten,intersmear=intersmear,
                              wmin=wmin,wmax=wmax,nw=nw,smeartype=smeartype,
                              intrasmear=intrasmear,metalcalc=metalcalc)
    if jdos:
        calctype='jdos'
        exten='eps_jdos'
        AFLOWpi.run._gen_eps_in(oneCalc,ID,calctype,exten,intersmear=intersmear,
                                wmin=wmin,wmax=wmax,nw=nw,smeartype=smeartype,
                                intrasmear=intrasmear,metalcalc=metalcalc)                            
    if occ:
        calctype='occ'
        exten='eps_occ'
        AFLOWpi.run._gen_eps_in(oneCalc,ID,calctype,exten,intersmear=intersmear,
                                wmin=wmin,wmax=wmax,nw=nw,smeartype=smeartype,
                                intrasmear=intrasmear,metalcalc=metalcalc)        
    if offdiag:
        calctype='offdiag'
        exten='eps_offdiag'
        AFLOWpi.run._gen_eps_in(oneCalc,ID,calctype,exten,intersmear=intersmear,
                                wmin=wmin,wmax=wmax,nw=nw,smeartype=smeartype,
                                intrasmear=intrasmear,metalcalc=metalcalc)
        

def _gen_eps_in(oneCalc,ID,calculation,exten,intersmear=0.136,wmin=0.0,wmax=30.0,nw=600,smeartype='gauss',intrasmear=0.0,metalcalc=False):
 
    prefix=oneCalc['_AFLOWPI_PREFIX_']
    hdir= oneCalc['_AFLOWPI_FOLDER_']

    in_str="""&inputpp
  outdir="./"
  prefix=%s
  calculation=%s
/
&energy_grid
  smeartype=%s
  intersmear=%s
  intrasmear=%s
  wmax=%s
  wmin=%s
  nw=%s
!  shift=0.0d0
!  metalcalc=%s
/
"""%(prefix,calculation,smeartype,intersmear,intrasmear,wmax,wmin,nw,metalcalc)

    with open(os.path.join(hdir,'%s_%s.in'%(ID,exten)),'w') as ofs:
        ofs.write(in_str)

def _rename_qe_eps(oneCalc,ID):
    prefix=oneCalc['_AFLOWPI_PREFIX_']
    hdir= oneCalc['_AFLOWPI_FOLDER_']

    for dat_type in ["epsr","epsi","ieps","eels"]:
        old=os.path.join(hdir,dat_type+"_"+prefix+".dat")
        new=os.path.join(hdir,ID+"_"+dat_type+".dat")
        try:
            os.rename(old,new)
        except Exception as e: print(e)
            
