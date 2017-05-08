import numpy
import AFLOWpi
import copy
import decimal


def supercell(inputString,numX=1,numY=1,numZ=1):
    inputString = AFLOWpi.prep._transformInput(inputString)    
    return AFLOWpi.retr._constructSupercell(inputString,numX=numX,numY=numY,numZ=numZ)




def _expandBoundaries(labels,symMatrix,numX,numY,numZ,beginX=0,beginY=0,beginZ=0):
    superList=[]
    try:
        symMatrix=symMatrix.getA()
    except:
        pass

    for entry in range(len(symMatrix)):
        for x in range(0,numX):
            first  = (symMatrix[entry][0]+x)/float(abs(numX))
            second = symMatrix[entry][1]
            third  = symMatrix[entry][2]
            superList.append([labels[entry],first,second,third])

    newPos=[x[1:] for x in superList]
    labels=[x[0] for x in superList]
    superList=[]
    #################################################################
    #################################################################
    for entry in range(len(newPos)):
        for y in range(0,numY):
            first  = newPos[entry][0]
            second = (newPos[entry][1]+y)/float(abs(numY))
            third  = newPos[entry][2]
            superList.append([labels[entry],first,second,third])

    labels=[x[0] for x in superList]
    newPos=[x[1:] for x in superList]
    superList=[]
    #################################################################
    #################################################################
    for entry in range(len(newPos)):
        for z in range(0,numZ):
            first  = newPos[entry][0]
            second = newPos[entry][1]
            third  = (newPos[entry][2]+z)/float(abs(numZ))
            superList.append([labels[entry],first,second,third])

    labels=[x[0] for x in superList]
    newPos=[x[1:] for x in superList]

    orig_list=[]
    orig_atom_ss_index = [(x-1)*numX*numY*numZ for x in range(1,len(symMatrix)+1)]


    for i in range(len(orig_atom_ss_index)):
        
        popped=superList.pop(orig_atom_ss_index[i])
        superList.insert(i,popped)

    symMatrix= numpy.array([x[1:] for x in superList])
    labels = numpy.array([x[0] for x in superList])

    return labels,symMatrix

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
        posLineStr = ' '.join(['%20.14f' % (decimal.Decimal(str(numpy.around(i,9)))) for i in symMatrix[entry]])+'\n'
        outputString+='%4s %8s' % (labels[entry],posLineStr)

    
    splitInput['&system']['nat']=str(len(labels))




    splitInput['&system']['celldm(1)']=str(float(splitInput['&system']['celldm(1)'])*numX)


    scaleY = float(numY)/float(numX)
    if 'celldm(2)' in splitInput['&system'].keys():
        splitInput['&system']['celldm(2)']=str(float(splitInput['&system']['celldm(2)'])*scaleY)
    else:
        splitInput['&system']['celldm(2)']=str(float(scaleY))

    scaleZ = float(numZ)/float(numX)
    if 'celldm(3)' in splitInput['&system'].keys():
        splitInput['&system']['celldm(3)']=str(float(splitInput['&system']['celldm(3)'])*scaleZ)
    else:
        splitInput['&system']['celldm(3)']=str(float(scaleZ))



    splitInput['ATOMIC_POSITIONS']['__content__']=outputString

    returnString=AFLOWpi.retr._joinInput(splitInput)
    iso = AFLOWpi.prep.isotropy()
    iso.qe_input(returnString,accuracy=0.001)
    print iso.convert()
    return returnString


