
import AFLOWpi
import re
import subprocess
import os 

def _standardize_alat(in_str):

    si = AFLOWpi.retr._splitInput(in_str)
    if "celldm(1)" in [x.lower() for x in si["&system"].keys()]:
        alat=float(si["&system"]["celldm(1)"])*0.529177
        return alat
    elif "a" in [x.lower() for x in si["&system"].keys()]:
        alat=float(si["&system"]["a"])
        return alat
    try:

        AFLOWSYM_LOC = os.path.join(AFLOWpi.__path__[0],'AFLOWSYM')
        AFLOW_EXE    = os.path.join(AFLOWSYM_LOC,'aflow')
        find_sym_process = subprocess.Popen('%s --edata=1.e-9'%AFLOW_EXE,stdin=subprocess.PIPE,stdout=subprocess.PIPE,shell=True)
        output = find_sym_process.communicate(input=in_str)[0]

        standard_input = re.findall('SCONV.*\n((?:.*\n)+)',output)[0]

        split_stand_input = AFLOWpi.retr._splitInput(standard_input)
        cell = AFLOWpi.retr._cellStringToMatrix(split_stand_input['CELL_PARAMETERS']['__content__']).A
        # print cell
        alat=cell[0][0]

        if alat==0:
            alat=cell[1][0]
            if alat==0:
                alat=cell[2][0]

    except Exception,e:
        print e
        alat=0.529177

    return alat
