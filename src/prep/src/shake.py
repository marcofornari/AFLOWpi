import AFLOWpi
import numpy as np
import os

def _shake_atoms(oneCalc,ID,dist=0.1,weight_by_mass=False,atoms=[]):


    pos_str = AFLOWpi.retr._getPositions(oneCalc['_AFLOWPI_INPUT_'],matrix=False)
    pos_str,flags = AFLOWpi.retr.detachPosFlags(pos_str)
    labels = AFLOWpi.retr._getPosLabels(oneCalc['_AFLOWPI_INPUT_'])

    alat,cell = AFLOWpi.retr._getCellParams(oneCalc,ID)

    if np.abs(cell.getA()[0][0])<=2.0:
        cell*=alat

    pos = AFLOWpi.retr._getPositions(oneCalc['_AFLOWPI_INPUT_'],matrix=True)
    cart_pos=pos.dot(cell)
    
    phi   = np.random.random((cart_pos.shape[0]))*2.0*np.pi
    theta = np.random.random((cart_pos.shape[0]))*1.0*np.pi

    shift = np.zeros((cart_pos.shape[0],3))

    weights = np.ones(cart_pos.shape[0])
    if weight_by_mass:
        for i in range(len(labels)):
            sl = labels[i].strip("0123456789")
            weights[i] = AFLOWpi.prep._getAMass(sl)[0]
        weights = np.amin(weights)/weights

    if len(atoms)!=0:
        for i in range(len(labels)):
            if labels[i] not in shake_atoms:
                weights[i]=0.0


    dist*=weights

    shift[:,0] = dist*np.sin(theta)*np.cos(phi)
    shift[:,1] = dist*np.sin(theta)*np.sin(phi)
    shift[:,2] = dist*np.cos(theta)
    print(cell)

    BohrToAngstrom = 0.529177249
    cart_pos+=shift/BohrToAngstrom


    pos = cart_pos.dot(np.linalg.inv(cell))

    pos_str = AFLOWpi.retr._joinMatrixLabels(labels,pos)
    pos_str = AFLOWpi.retr.attachPosFlags(pos_str,flags)

    split_input = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    split_input['ATOMIC_POSITIONS']['__content__']=pos_str
    oneCalc['_AFLOWPI_INPUT_'] = AFLOWpi.retr._joinInput(split_input)

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID),'w') as ofo:
        ofo.write(oneCalc['_AFLOWPI_INPUT_'])

    return oneCalc,ID


    
