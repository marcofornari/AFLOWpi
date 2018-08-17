import AFLOWpi
import numpy as np
import re
import subprocess
import os 

def _getHighSymPoints(oneCalc,ID=None):
    '''
    Searching for the ibrav number in the input file for the calculation
    to determine the path for the band structure calculation

    Arguments:
          oneCalc (dict): a dictionary containing properties about the AFLOWpi calculation
                
    Keyword Arguments:
          ID (str): ID string for the particular calculation and step          

    Returns:
          special_points (list): list of the HSP names
          band_path (str): path in string form

    '''


    
    ibrav = 0
    ibravRegex = re.compile('ibrav[\s]*=[\s]*([\d]+)\s*[,\n]*')

    ibrav=int(ibravRegex.findall(oneCalc['_AFLOWPI_INPUT_'])[-1])




    alat,cellOld = AFLOWpi.retr._getCellParams(oneCalc,ID)        
    orig_ibrav = oneCalc['_AFLOWPI_ORIG_IBRAV_']
    if orig_ibrav!=0:
        ibrav = orig_ibrav
    else:
        return AFLOWpi.retr._getHighSymPoints_aflow(oneCalc,ID=ID)

    cellOld=cellOld.getA()

#    raise SystemExit
    #get a,b,c of QE convention conventional cell from primitive lattice vecs
    if ibrav == 1:
        a=np.abs(cellOld[0][0])
    if ibrav == 2:
        a=np.abs(cellOld[0][0])*2.0
    if ibrav == 3:
        a=np.abs(cellOld[0][0])*2.0
    if ibrav == 4:
        a=np.abs(cellOld[0][0])
        b=a
        c=np.abs(cellOld[2][2])
    if ibrav == 5:
        a=np.sqrt(cellOld[0].dot(cellOld[0].T))
        alpha=np.arccos(cellOld[0].dot(cellOld[1].T)/(a**2))
        alpha_deg = alpha*180.0/np.pi
    if ibrav == 6:
        a=np.abs(cellOld[0][0])
        c=np.abs(cellOld[2][2])
    if ibrav == 7:
        a=np.abs(cellOld[0][0])*2.0
        c=np.abs(cellOld[2][2])*2.0
    if ibrav == 8:
        a=np.abs(cellOld[0][0])
        b=np.abs(cellOld[1][1])
        c=np.abs(cellOld[2][2])

    if ibrav == 9:
        a=np.abs(cellOld[0][0])*2.0
        b=np.abs(cellOld[1][1])*2.0
        c=np.abs(cellOld[2][2])

    if ibrav == 10:
        a=np.abs(cellOld[0][0])*2.0
        b=np.abs(cellOld[1][1])*2.0
        c=np.abs(cellOld[2][2])*2.0

    if ibrav == 11:
        a=np.abs(cellOld[0][0])*2.0
        b=np.abs(cellOld[1][1])*2.0
        c=np.abs(cellOld[2][2])*2.0

    if ibrav == 12:
        a=np.sqrt(cellOld[0].dot(cellOld[0].T))
        b=np.sqrt(cellOld[1].dot(cellOld[1].T))
        c=np.sqrt(cellOld[2].dot(cellOld[2].T))
        alpha=np.arccos(cellOld[1].dot(cellOld[2].T)/(c*b))
        beta =np.arccos(cellOld[0].dot(cellOld[2].T)/(a*c))
        gamma=np.arccos(cellOld[0].dot(cellOld[1].T)/(a*b))

    if ibrav==14:
        cellOld=np.linalg.inv(cellOld.T).T
        a=np.sqrt(cellOld[0].dot(cellOld[0].T))
        b=np.sqrt(cellOld[1].dot(cellOld[1].T))
        c=np.sqrt(cellOld[2].dot(cellOld[2].T))
        alpha=np.arccos(cellOld[1].dot(cellOld[2].T)/(c*b))*(180.0/np.pi)
        beta =np.arccos(cellOld[0].dot(cellOld[2].T)/(a*c))*(180.0/np.pi)
        gamma=np.arccos(cellOld[0].dot(cellOld[1].T)/(a*b))*(180.0/np.pi)
        alat*=  0.529177249


    if   ibrav==1:  ibrav_var =  'CUB'
    elif ibrav==2:  ibrav_var =  'FCC'
    elif ibrav==3:  ibrav_var =  'BCC'
    elif ibrav==4:  ibrav_var =  'HEX'
    elif ibrav==6:  ibrav_var =  'TET'
    elif ibrav==8:  ibrav_var =  'ORC'
    elif ibrav==9:  ibrav_var =  'ORCC'
    elif ibrav==11: ibrav_var =  'ORCI'
    elif ibrav==12: ibrav_var =  'MCL'

    elif ibrav==5:
        if   alpha_deg < 90.0: ibrav_var = 'RHL1'
        elif alpha_deg > 90.0: ibrav_var = 'RHL2'
    elif ibrav==7:
        if(c < a):   ibrav_var =  'BCT1'
        elif(c > a): ibrav_var =  'BCT2'
        else:        ibrav_var =  'BCC'
    elif ibrav==10:

        if    (1.0/a**2 > 1.0/b**2+1.0/c**2): ibrav_var =  'ORCF1'
        elif  np.isclose(1.0/a**2, 1.0/b**2+1.0/c**2,1.e-2): ibrav_var =  'ORCF3'
        elif  (1.0/a**2 < 1.0/b**2+1.0/c**2): ibrav_var =  'ORCF2'

    elif(int(ibrav)==14):
        print alpha,beta,gamma
        minAngle = np.amin([alpha,beta,gamma])
        maxAngle = np.amax([alpha,beta,gamma])
        if alpha==90.0 or beta==90.0 or gamma==90.0:
            if alpha>=90.0 or beta>=90.0 or gamma>=90.0: ibrav_var =  'TRI2A'
            if alpha<=90.0 or beta<=90.0 or gamma<=90.0: ibrav_var =  'TRI2B'
        elif minAngle>90.0:                              ibrav_var =  'TRI1A'
        elif maxAngle<90:                                ibrav_var =  'TRI1B'
        else: ibrav_var =  'TRI1A'
###############################################################################
###############################################################################
    if ibrav_var=='CUB':
        band_path = 'gG-X-M-gG-R-X|M-R'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                           'M'   : (0.5, 0.5, 0.0),
                           'R'   : (0.5, 0.5, 0.5),
                           'X'   : (0.0, 0.5, 0.0)}
                           
    if ibrav_var=='FCC':
        band_path = 'gG-X-W-K-gG-L-U-W-L-K|U-X'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                          'K'    : (0.375, 0.375, 0.750),
                          'L'    : (0.5, 0.5, 0.5),
                          'U'    : (0.625, 0.250, 0.625),
                          'W'    : (0.5, 0.25, 0.75),
                          'X'    : (0.5, 0.0, 0.5)}
                          
    if ibrav_var=='BCC':
        band_path = 'gG-H-N-gG-P-H|P-N'
        special_points = {'gG'   : (0, 0, 0),
                          'H'    : (0.5, -0.5, 0.5),
                          'P'    : (0.25, 0.25, 0.25,), 
                          'N'    : (0.0, 0.0, 0.5)}
            
    if ibrav_var=='HEX':
        band_path = 'gG-M-K-gG-A-L-H-A|L-M|K-H'
        special_points = {'gG'   : (0, 0, 0),
                          'A'    : (0.0, 0.0, 0.5),
                          'H'    : (1.0/3.0, 1.0/3.0, 0.5),
                          'K'    : (1.0/3.0, 1.0/3.0, 0.0),
                          'L'    : (0.5, 0.0, 0.5),
                          'M'    : (0.5, 0.0, 0.0)}
        
    if ibrav_var=='RHL1':
        eta = (1.0 + 4.0*np.cos(alpha))/(2.0 + 4.0*np.cos(alpha))
        nu =0.75-eta/2.0
        band_path = 'gG-L-B1|B-Z-gG-X|Q-F-P1-Z|L-P'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                          'B'    : (eta, 0.5, 1.0-eta),
                          'B1'   : (0.5, 1.0-eta, eta-1.0),
                          'F'    : (0.5, 0.5, 0.0),
                          'L'    : (0.5, 0.0, 0.0),
                          'L1'   : (0.0, 0.0, -0.5),
                          'P'    : (eta, nu, nu),
                          'P1'   : (1.0-nu, 1.0-nu, 1.0-eta),
                          'P2'   : (nu, nu, eta-1.0),
                          'Q'    : (1.0-nu, nu, 0.0),
                          'X'    : (nu, 0.0, -nu),
                          'Z'    : (0.5, 0.5, 0.5)}
         
    if ibrav_var=='RHL2':
        eta=1.0/(2.0*np.tan(alpha/2.0)**2)
        nu =0.75-eta/2.0
        band_path = 'gG-P-Z-Q-gG-F-P1-Q1-L-Z'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                          'F'    : (0.5, -0.5, 0.0),
                          'L'    : (0.5, 0.0, 0.0),
                          'P'    : (1.0-nu, -nu, 1.0-nu),
                          'P1'   : (nu, nu-1.0, nu-1.0),
                          'Q'    : (eta, eta, eta),
                          'Q1'   : (1.0-eta, -eta, -eta),
                          'Z'    : (0.5, -0.5, 0.5)} 

    if ibrav_var=='TET':
        band_path = 'gG-X-M-gG-Z-R-A-Z|X-R|M-A'
        special_points = {'gG'   : (0.0, 0.0, 0.0),
                          'A'    : (0.5, 0.5, 0.5),
                          'M'    : (0.5, 0.5, 0.0),
                          'R'    : (0.0, 0.5, 0.5),
                          'X'    : (0.0, 0.5, 0.0),
                          'Z'    : (0.0, 0.0, 0.5)}

    if ibrav_var=='BCT1':
       eta = (1.0+(c/a)**2)/4.0
       band_path = 'gG-X-M-gG-Z-P-N-Z1-M|X-P'
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'M'     : (-0.5, 0.5, 0.5),
                         'N'     : (0.0, 0.5, 0.0),
                         'P'     : (0.25, 0.25, 0.25),
                         'X'     : (0.0, 0.0, 0.5),
                         'Z'     : (eta, eta, -eta),
                         'Z1'    : (-eta, 1.0-eta, eta)}
         
    if ibrav_var=='BCT2':
       band_path = 'gG-X-Y-gS-gG-Z-gS1-N-P-Y1-Z|X-P'
       eta = (1.0+(a/c)**2)/4.0
       zeta = 0.5*(a/c)**2

       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'N'     : (0.0, 0.5, 0.0),
                         'P'     : (0.25, 0.25, 0.25),
                         'gS'    : (-eta, eta, eta),
                         'gS1'   : (eta, 1-eta, -eta),
                         'X'     : (0.0, 0.0, 0.5),
                         'Y'     : (-zeta, zeta, 0.5),
                         'Y1'    : (0.5, 0.5, -zeta),
                         'Z'     : (0.5, 0.5, -0.5)}
         
    if ibrav_var=='ORC':
         band_path = 'gG-X-S-Y-gG-Z-U-R-T-Z|Y-T|U-X|S-R'
         special_points = {'gG'  : (0.0, 0.0, 0.0),
                           'R'   : (0.5, 0.5, 0.5),
                           'S'   : (0.5, 0.5, 0.0),
                           'T'   : (0.0, 0.5, 0.5),
                           'U'   : (0.5, 0.0, 0.5),
                           'X'   : (0.5, 0.0, 0.0),
                           'Y'   : (0.0, 0.5, 0.0),
                           'Z'   : (0.0, 0.0, 0.5)}

    if ibrav_var=='ORCC':
       band_path = 'gG-X-S-R-A-Z-gG-Y-X1-A1-T-Y|Z-T'
       zeta=(1.0+((a/b)**2))/4.0
       special_points = {'gG'    : (  0.0, 0.0     , 0.0),
                         'A'     : ( zeta, zeta    , 0.5),
                         'A1'    : (-zeta, 1.0-zeta, 0.5),
                         'R'     : (  0.0, 0.5     , 0.5),
                         'S'     : (  0.0, 0.5     , 0.0),
                         'T'     : ( -0.5, 0.5     , 0.5),
                         'X'     : ( zeta, zeta    , 0.0),
                         'X1'    : (-zeta, 1.0-zeta, 0.0),
                         'Y'     : ( -0.5, 0.5     , 0.0),
                         'Z'     : (  0.0, 0.0     , 0.5)}


    if ibrav_var=='ORCF1':
       band_path = 'gG-Y-T-Z-gG-X-A1-Y|T-X1|X-A-Z|L-gG'
       eta =(1.0+(a/b)**2+(a/c)**2)/4.0
       zeta=(1.0+(a/b)**2-(a/c)**2)/4.0
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'A'     : (0.5, 0.5 + zeta, zeta),
                         'A1'    : (0.5, 0.5-zeta, 1.0-zeta),
                         'L'     : (0.5, 0.5, 0.5),
                         'T'     : (1.0, 0.5, 0.5),
                         'X'     : (0.0, eta, eta),
                         'X1'    : (1.0, 1.0-eta, 1.0-eta),
                         'Y'     : (0.5, 0.0, 0.5),
                         'Z'     : (0.5, 0.5, 0.0)}

    if ibrav_var=='ORCF2':
       band_path = 'gG-Y-C-D-X-gG-Z-D1-H-C|C1-Z|X-H1|H-Y|L-gG'
       eta =(1.0+(a/b)**2-(a/c)**2)/4.0
       phi =(1.0+(c/b)**2-(c/a)**2)/4.0
       delta =(1.0+(b/a)**2-(b/c)**2)/4.0

       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'C'     : (0.5, 0.5-eta, 1.0-eta),
                         'C1'    : (0.5, 0.5+eta, eta),
                         'D'     : (0.5-delta, 0.5, 1.0-delta),
                         'D1'    : (0.5+delta, 0.5, delta),
                         'L'     : (0.5, 0.5, 0.5),
                         'H'     : (1.0-phi, 0.5-phi, 0.5),
                         'H1'    : (phi, 0.5+phi, 0.5),
                         'X'     : (0.0, 0.5, 0.5),
                         'Y'     : (0.5, 0.0, 0.5),
                         'Z'     : (0.5, 0.5, 0.0),}

    if ibrav_var=='ORCF3':
       band_path = 'gG-Y-T-Z-gG-X-A1-Y|X-A-Z|L-R'
       eta =(1.0+(a/b)**2+(a/c)**2)/4.0
       zeta=(1.0+(a/b)**2+(a/c)**2)/4.0
       special_points = {'gG'    : (0.0, 0.0, 0.0),
                         'A'     : (0.5, 0.5 + zeta, zeta),
                         'A1'    : (0.5, 0.5-zeta, 1.0-zeta),
                         'L'     : (0.5, 0.5, 0.5),
                         'T'     : (1.0, 0.5, 0.5),
                         'X'     : (0.0, eta, eta),
                         'X1'    : (1.0, 1.0-eta, 1.0-eta),
                         'Y'     : (0.5, 0.0, 0.5),
                         'Z'     : (0.5, 0.5, 0.0)}

         
    if ibrav_var=='ORCI':
         band_path = 'gG-X-L-T-W-R-X1-Z-gG-Y-S-W|L1-Y|Y1-Z'
         chi   = (1.0  + (a/c)**2)/(4.0)
         eta   = (1.0  + (b/c)**2)/(4.0)
         delta = ( (b/c)**2 - (a/c)**2    )/(4.0)
         mu    = ( (b/c)**2 + (a/c)**2    )/(4.0)
         special_points = {'gG'   : (0, 0, 0),
                           'L'    : (-mu, mu, 0.5-delta),
                           'L1'   : (mu, -mu, 0.5+delta),
                           'L2'   : (0.5-delta, 0.5+delta, -mu),
                           'R'    : (0.0, 0.5, 0.0),
                           'S'    : (0.5, 0.0, 0.0),
                           'T'    : (0.0, 0.0, 0.5),
                           'W'    : (0.25,0.25,0.25),
                           'X'    : (-chi, chi, chi),
                           'X1'   : (chi, 1.0-chi, -chi),
                           'Y'    : (eta, -eta, eta),
                           'Y1'   : (1.0-eta, eta, -eta),
                           'Z'    : (0.5, 0.5, -0.5)}

    if ibrav_var=='MCL':
         #abc->cba
         eta =  (1.0 - (b/a)*np.cos(np.pi-gamma))/(2.0*np.sin(np.pi-gamma)**2)
         nu =   0.5  - eta*(a/b)*np.cos(np.pi-gamma)
         band_path = 'gG-Y-H-C-E-M1-A-X-gG-Z-D-M|Z-A|D-Y|X-H1'
         special_points = {
                           'gG'    : (0.0, 0.0    , 0.0    ),
                           'A'     : (0.5, 0.5    , 0.0    ),
                           'C'     : (0.0, 0.5    , 0.5    ),
                           'D'     : (0.5, 0.0    , 0.5    ),
                           'D1'    : (0.5, 0.0    ,-0.5    ),
                           'E'     : (0.5, 0.5    , 0.5    ),
                           'H'     : (0.0, eta    , 1.0-nu ),
                           'H1'    : (0.0, 1.0-eta, nu     ),
                           'H2'    : (0.0, eta    ,-nu     ),
                           'M'     : (0.5, eta    , 1.0-nu ),
                           'M1'    : (0.5, 1.0-eta, nu     ),
                           'M2'    : (0.5, eta    ,-nu     ),
                           'X'     : (0.0, 0.5    , 0.0    ),
                           'Y'     : (0.0, 0.0    , 0.5    ),
                           'Y1'    : (0.0, 0.0    ,-0.5    ),
                           'Z'     : (0.5, 0.0    , 0.0    )}
         
   
    if ibrav_var=='TRI1A':         
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG' 
        special_points = {'gG'    : (0.0,0.0,0.0),
                          'L'     : (0.5,0.5,0.0),
                          'M'     : (0.0,0.5,0.5),
                          'N'     : (0.5,0.0,0.5),
                          'R'     : (0.5,0.5,0.5),
                          'X'     : (0.5,0.0,0.0),
                          'Y'     : (0.0,0.5,0.0),
                          'Z'     : (0.0,0.0,0.5),}
        
    if ibrav_var=='TRI2A':        
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG'
        special_points = {'gG'    : (0.0,0.0,0.0),
                          'L'     : (0.5,0.5,0.0),
                          'M'     : (0.0,0.5,0.5),
                          'N'     : (0.5,0.0,0.5),
                          'R'     : (0.5,0.5,0.5),
                          'X'     : (0.5,0.0,0.0),
                          'Y'     : (0.0,0.5,0.0),
                          'Z'     : (0.0,0.0,0.5),}
 
    if ibrav_var=='TRI1B':        
        band_path = "X-gG-Y|L-gG-Z|N-gG-M|R-gG"
        special_points = {'gG'    : ( 0.0, 0.0,0.0),
                          'L'     : ( 0.5,-0.5,0.0),
                          'M'     : ( 0.0, 0.0,0.5),
                          'N'     : (-0.5,-0.5,0.5),
                          'R'     : ( 0.0,-0.5,0.5),
                          'X'     : ( 0.0,-0.5,0.0),
                          'Y'     : ( 0.5, 0.0,0.0),
                          'Z'     : (-0.5, 0.0,0.5),}

    if ibrav_var=='TRI2B':        
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG'
        special_points = {'gG'    : ( 0.0, 0.0,0.0),
                          'L'     : ( 0.5,-0.5,0.0),
                          'M'     : ( 0.0, 0.0,0.5),
                          'N'     : (-0.5,-0.5,0.5),
                          'R'     : ( 0.0,-0.5,0.5),
                          'X'     : ( 0.0,-0.5,0.0),
                          'Y'     : ( 0.5, 0.0,0.0),
                          'Z'     : (-0.5, 0.0,0.5),}


    aflow_conv = np.identity(3)
    qe_conv    = np.identity(3)

    if ibrav==2:
        aflow_conv = np.asarray([[ 0.0, 1.0, 1.0],[ 1.0, 0.0, 1.0],[ 1.0, 1.0, 0.0]])/2.0     
        qe_conv    = np.asarray([[-1.0, 0.0, 1.0],[ 0.0, 1.0, 1.0],[-1.0, 1.0, 0.0]])/2.0
    if ibrav==3:
        aflow_conv = np.asarray([[-1.0, 1.0, 1.0],[ 1.0,-1.0, 1.0],[ 1.0, 1.0,-1.0]])/2.0     
        qe_conv    = np.asarray([[ 1.0, 1.0, 1.0],[-1.0, 1.0, 1.0],[-1.0,-1.0, 1.0]])/2.0     
#    if ibrav==4:
#        s32 = np.sqrt(3)/2.0
#        aflow_conv = np.asarray([[ 0.5,-s32, 0.0],[ 0.5, s32, 0.0],[ 0.0, 0.0, 1.0]])     
#        qe_conv    = np.asarray([[ 1.0, 0.0, 0.0],[-0.5, s32, 0.0],[ 0.0, 0.0, 1.0]])     
    if ibrav==7:
        aflow_conv = np.asarray([[-1.0, 1.0, 1.0],[ 1.0,-1.0, 1.0],[ 1.0, 1.0,-1.0]])/2.0
        qe_conv    = np.asarray([[ 1.0,-1.0, 1.0],[ 1.0, 1.0, 1.0],[-1.0,-1.0, 1.0]])/2.0
    if ibrav==9:
        aflow_conv = np.asarray([[ 1.0,-1.0, 0.0],[ 1.0, 1.0, 0.0],[ 0.0, 0.0, 2.0]])/2.0
        qe_conv    = np.asarray([[ 1.0, 1.0, 0.0],[-1.0, 1.0, 0.0],[ 0.0, 0.0, 2.0]])/2.0
    if ibrav==10:
        aflow_conv = np.asarray([[ 0.0, 1.0, 1.0],[ 1.0, 0.0, 1.0],[ 1.0, 1.0, 0.0]])/2.0
        qe_conv    = np.asarray([[ 1.0, 0.0, 1.0],[ 1.0, 1.0, 0.0],[ 0.0, 1.0, 1.0]])/2.0  
    if ibrav==11:
        aflow_conv = np.asarray([[-1.0, 1.0, 1.0],[ 1.0,-1.0, 1.0],[ 1.0, 1.0,-1.0]])/2.0
        qe_conv    = np.asarray([[ 1.0, 1.0, 1.0],[-1.0, 1.0, 1.0],[-1.0,-1.0, 1.0]])/2.0
    if ibrav==12:
        aflow_conv = np.asarray([[ 0.0, 0.0, 1.0],[ 0.0, 1.0, 0.0],[ 1.0, 0.0, 0.0]])
        qe_conv    = np.asarray([[ 1.0, 0.0, 0.0],[ 0.0, 1.0, 0.0],[ 0.0, 0.0, 1.0]])



                                   

    for k,v in special_points.iteritems():
        first  = np.array(v).dot(np.linalg.inv(aflow_conv))
        if ibrav == 9:
            second = qe_conv.T.dot(first)
        else:
            second = qe_conv.dot(first)
        special_points[k]=tuple((second).tolist())

    return special_points, band_path


def _path_by_lattice_variation(ibrav_var):
    if ibrav_var=='CUB':
        band_path = 'gG-X-M-gG-R-X|M-R'
    if ibrav_var=='FCC':
        band_path = 'gG-X-W-K-gG-L-U-W-L-K|U-X'
    if ibrav_var=='BCC':
        band_path = 'gG-H-N-gG-P-H|P-N'
    if ibrav_var=='HEX':
        band_path = 'gG-M-K-gG-A-L-H-A|L-M|K-H'
    if ibrav_var=='RHL1':
        band_path = 'gG-L-B1|B-Z-gG-X|Q-F-P1-Z|L-P'
    if ibrav_var=='RHL2':
        band_path = 'gG-P-Z-Q-gG-F-P1-Q1-L-Z'
    if ibrav_var=='TET':
        band_path = 'gG-X-M-gG-Z-R-A-Z|X-R|M-A'
    if ibrav_var=='BCT1':
        band_path = 'gG-X-M-gG-Z-P-N-Z1-M|X-P'
    if ibrav_var=='BCT2':
        band_path = 'gG-X-Y-gS-gG-Z-gS1-N-P-Y1-Z|X-P'
    if ibrav_var=='ORC':
        band_path = 'gG-X-S-Y-gG-Z-U-R-T-Z|Y-T|U-X|S-R'
    if ibrav_var=='ORCC':
        band_path = 'gG-X-S-R-A-Z-gG-Y-X1-A1-T-Y|Z-T'
    if ibrav_var=='ORCF1':
        band_path = 'gG-Y-T-Z-gG-X-A1-Y|T-X1|X-A-Z|L-gG'
    if ibrav_var=='ORCF2':
        band_path = 'gG-Y-C-D-X-gG-Z-D1-H-C|C1-Z|X-H1|H-Y|L-gG'
    if ibrav_var=='ORCF3':
        band_path = 'gG-Y-T-Z-gG-X-A1-Y|X-A-Z|L-R'
    if ibrav_var=='ORCI':
        band_path = 'gG-X-L-T-W-R-X1-Z-gG-Y-S-W|L1-Y|Y1-Z'
    if ibrav_var=='MCL':
        band_path = 'gG-Y-H-C-E-M1-A-X-gG-Z-D-M|Z-A|D-Y|X-H1'
    if ibrav_var=='MCLC1':
        band_path = 'gG-Y-F-L-I|I1-Z-F1|Y-X1|X-gG-N|M-gG'
    if ibrav_var=='MCLC2':
        band_path = 'gG-Y-F-L-I|I1-Z-F1|N-gG-M'
    if ibrav_var=='MCLC3':
        band_path = 'gG-Y-F-H-Z-I-F1|H1-Y1-X-gG-N|M-gG'
    if ibrav_var=='MCLC4':
        band_path = 'gG-Y-F-H-Z-I|H1-Y1-X-gG-N|M-gG'
    if ibrav_var=='MCLC5':
        band_path = 'gG-Y-F-L-I|I1-Z-H-F1|H1-Y1-X-gG-N|M-gG'
    if ibrav_var=='TRI1A':         
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG' 
    if ibrav_var=='TRI2A':        
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG'
    if ibrav_var=='TRI1B':        
        band_path = "X-gG-Y|L-gG-Z|N-gG-M|R-gG"
    if ibrav_var=='TRI2B':        
        band_path = 'X-gG-Y|L-gG-Z|N-gG-M|R-gG'

    return band_path

def _getHighSymPoints_aflow(oneCalc,ID=None):

    in_str = oneCalc['_AFLOWPI_INPUT_']

    AFLOWSYM_LOC = os.path.join(AFLOWpi.__path__[0],'AFLOWSYM')
    AFLOW_EXE    = os.path.join(AFLOWSYM_LOC,'aflow')

    find_sym_process = subprocess.Popen('%s --kpath --grid=1'%AFLOW_EXE,stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,shell=True)
    output = find_sym_process.communicate(input=in_str)[0]

    ibrav_var = re.findall('K_POINTS\s*crystal\s*!\s*(\S*)\s*',output)[0]


    

    band_path = _path_by_lattice_variation(ibrav_var)

    SP_string = re.findall('([-\d\.]+)\s*([-\d\.]+)\s*([-\d\.]+)\s*.*\!\s*(.*)\s*\/\/ nk.*',output)

    label_map = {'\Gamma':'gG',
                 '\Sigma':'gS',
                 '\Sigma_1':'gS1',
                 'Z_1':'Z1',
                 'Y_1':'Y1',
                 'Y_2':'Y2',
                 'Y_3':'Y3',
                 'X_1':'X1',
                 'X_2':'X2',
                 'A_1':'A1',
                 'C_1':'C1',
                 'D_1':'D1',
                 'H_1':'H1',
                 'H_2':'H2',
                 'L_1':'L1',
                 'L_2':'L2',
                 'B_1':'B1',
                 'P_1':'P1',
                 'P_2':'P2',
                 'N_1':'N1',
                 'M_1':'M1',
                 'M_2':'M2',
                 'Q_1':'Q1',
                 'F_1':'F1',
                 'F_2':'F2',
                 'F_3':'P3',
                 'I_1':'I1',}




    kp_dict={}
    for i in SP_string:
        kx = float(i[0])
        ky = float(i[1])
        kz = float(i[2])
        point_label = i[3].strip()

        try:
            point_label = label_map[point_label]
        except:
            pass

        kp_dict[point_label]=(kx,ky,kz)

    return kp_dict,band_path

