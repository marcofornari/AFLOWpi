import AFLOWpi
import os 


def _prep_berry(oneCalc,ID,gdir,kp_mult):

    inputDict=AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    inputDict["&control"]["lberry"] = ".true."
    inputDict["&control"]["gdir"] = str(gdir)
    inputDict["&control"]["calculation"] = "'nscf'"

    grid = map(int,inputDict["K_POINTS"]["__content__"].split())

    grid[gdir-1] = int(kp_mult*grid[gdir-1])
    inputDict["K_POINTS"]["__content__"] = " ".join(map(str,grid))

    inputDict["&control"]["nppstr"] = "%s"%grid[gdir-1]

    inputString = AFLOWpi.retr._joinInput(inputDict)
    ID_gdir = ID+"_gdir%s"%gdir

    with open(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%ID_gdir),'w') as newIn:
        newIn.write(inputString)

    return oneCalc,ID_gdir
