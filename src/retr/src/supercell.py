# ***************************************************************************
# *                                                                         *
# *          AFLOWpi - Central Michigan University University, 2017         *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2017 - Andrew Supka and Marco Fornari - AFLOW.ORG consortium
#
#  This file is part of AFLOWpi software.
#
#  AFLOWpi is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ***************************************************************************

import numpy as np
import AFLOWpi
import copy
import decimal


def supercell(inputString,numX=1,numY=1,numZ=1):
    inputString = AFLOWpi.prep._transformInput(inputString)    
    return AFLOWpi.retr._constructSupercell(inputString,numX=numX,numY=numY,numZ=numZ)




# def _expandBoundaries(labels,symMatrix,numX,numY,numZ,beginX=0,beginY=0,beginZ=0):
#     superList=[]
#     try:
#         symMatrix=symMatrix.getA()
#     except:
#         pass

#     for entry in range(len(symMatrix)):
#         for x in range(0,numX):
#             first  = (symMatrix[entry][0]+x)/float(abs(numX))
#             second = symMatrix[entry][1]
#             third  = symMatrix[entry][2]
#             superList.append([labels[entry],first,second,third])

#     newPos=[x[1:] for x in superList]
#     labels=[x[0] for x in superList]
#     superList=[]
#     #################################################################
#     #################################################################
#     for entry in range(len(newPos)):
#         for y in range(0,numY):
#             first  = newPos[entry][0]
#             second = (newPos[entry][1]+y)/float(abs(numY))
#             third  = newPos[entry][2]
#             superList.append([labels[entry],first,second,third])

#     labels=[x[0] for x in superList]
#     newPos=[x[1:] for x in superList]
#     superList=[]
#     #################################################################
#     #################################################################
#     for entry in range(len(newPos)):
#         for z in range(0,numZ):
#             first  = newPos[entry][0]
#             second = newPos[entry][1]
#             third  = (newPos[entry][2]+z)/float(abs(numZ))
#             superList.append([labels[entry],first,second,third])

#     labels=[x[0] for x in superList]
#     newPos=[x[1:] for x in superList]

#     orig_list=[]
#     orig_atom_ss_index = [(x-1)*numX*numY*numZ for x in range(1,len(symMatrix)+1)]


#     for i in range(len(orig_atom_ss_index)):
        
#         popped=superList.pop(orig_atom_ss_index[i])
#         superList.insert(i,popped)

#     symMatrix= np.array([x[1:] for x in superList])
#     labels = np.array([x[0] for x in superList])

#     return labels,symMatrix

def _expandBoundaries(labels,symMatrix,numX,numY,numZ,beginX=0,beginY=0,beginZ=0):

    try:
        symMatrix=symMatrix.getA()
    except:
        pass
    
    nat = symMatrix.shape[0]
    ss_coords=np.zeros((numX,numY,numZ,nat,3),order="C")
    ss_coords[:,:,:]=symMatrix

    labs=np.tile(np.array(labels),numX*numY*numZ)

    for x in range(numX):
        ss_coords[x,:,:,:,0]+=1.0*x
    for y in range(numY):
        ss_coords[:,y,:,:,1]+=1.0*y
    for z in range(numZ):
        ss_coords[:,:,z,:,2]+=1.0*z

    ss_coords[:,:,:,:,0]/=numX
    ss_coords[:,:,:,:,1]/=numY
    ss_coords[:,:,:,:,2]/=numZ

    ss_coords = np.reshape(ss_coords,(numX*numY*numZ*nat,3),order="C")

    return labs,ss_coords


def _constructSupercell(inputString,numX=1,numY=1,numZ=1,stringOrMatrix='String',newVectors=True):
    splitInput =  AFLOWpi.retr._splitInput(inputString)
    
    if '{crystal}' != splitInput['ATOMIC_POSITIONS']['__modifier__']:
        logging.error('unit in AFLOWpi.retr._constructSupercell not for ATOMIC_POSITIONS MUST BE {crystal}')
        return inputString

    cellParamMatrix = AFLOWpi.retr.getCellMatrixFromInput(inputString)

    labels =  AFLOWpi.retr._getPosLabels(inputString)
    symMatrix = AFLOWpi.retr._getPositions(inputString)

    symMatrix_orig=copy.deepcopy(symMatrix)
    labels_orig=copy.deepcopy(labels)

    coordold,flags = AFLOWpi.retr.detachPosFlags(AFLOWpi.qe.regex.atomic_positions(inputString))

    outputString=''
    superList=[]

    if newVectors==True:
        labels,symMatrix=AFLOWpi.retr._expandBoundaries(labels,symMatrix,numX,numY,numZ)
    else:
        labels,symMatrix=AFLOWpi.retr._expandBoundariesNoScale(labels,symMatrix,numX,numY,numZ)

    for entry in range(len(symMatrix)):
        posLineStr = ' '.join(['%20.14f' % (decimal.Decimal(str(np.around(i,9)))) for i in symMatrix[entry]])+'\n'
        outputString+='%4s %8s' % (labels[entry],posLineStr)

    
    splitInput['&system']['nat']=str(len(labels))




    splitInput['&system']['celldm(1)']=str(float(splitInput['&system']['celldm(1)']))
    alat=float(splitInput['&system']['celldm(1)'])

    cellParamMatrix[0]*=float(numX)    
    cellParamMatrix[1]*=float(numY)
    cellParamMatrix[2]*=float(numZ)

    cellParamMatrix/=alat
    splitInput['CELL_PARAMETERS']['__content__']=AFLOWpi.retr._cellMatrixToString(cellParamMatrix)
    splitInput['CELL_PARAMETERS']['__modifier__']="{alat}"

    splitInput['ATOMIC_POSITIONS']['__content__']=outputString

    returnString = AFLOWpi.retr._joinInput(splitInput)

    returnString = AFLOWpi.run.reduce_kpoints(returnString,[numX,numY,numZ])

    returnString = AFLOWpi.prep._transformInput(returnString)    

    return returnString


