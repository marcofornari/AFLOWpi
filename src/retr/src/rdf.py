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

import AFLOWpi
import os
import shutil
import subprocess

def atomicDistances(calcs,runlocal=False,inpt=False,outp=True):
    engineDir  = AFLOWpi.prep._ConfigSectionMap("prep",'engine_dir')	
    if os.path.isabs(engineDir) == False:
        configFileLocation = AFLOWpi.prep._getConfigFile()
        engineDir =  os.path.join(configFileLocation, enginedir)
    distXPath = os.path.join(engineDir,'dist.x')
    try:
        if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()!='false':
            AFLOWpi.prep.totree(distXPath,calcs)
        for ID,oneCalc in list(calcs.items()):
            if runlocal:
                if outp==True:
                    AFLOWpi.retr._getDist(oneCalc,ID,outp=True)
                if inpt==True:
                    AFLOWpi.retr._getDist(oneCalc,ID,outp=False)
            else:
                if outp==True:
                    AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','AFLOWpi.retr._getDist(oneCalc,ID,outp=True)\n')
                if inpt==True:
                    AFLOWpi.prep._addToBlock(oneCalc,ID,'RUN','AFLOWpi.retr._getDist(oneCalc,ID,outp=False)\n')

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)


def _getDist(oneCalc,ID,outp=True):
    try:
        if outp==True:
            AFLOWpi.retr._writeInputFromOutput(oneCalc,ID)
            os.rename(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_new.in'),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_dist.in')) 
        else:
            shutil.copyfile('%s.in'%ID,'%s_dist.in'%ID)

    except Exception as e:
        AFLOWpi.run._fancy_error_log(e)

    #fix to clear restart_mode in input file to avoid error with dist.x
    dist_ID='%s_dist'%ID        
    infil = AFLOWpi.retr._getInputFileString(oneCalc,dist_ID)
    inDict=AFLOWpi.retr._splitInput(infil)
    inDict['&control']['restart_mode']='"from_scratch"'
    infil = AFLOWpi.retr._joinInput(inDict)
    dest_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_dist.in')

    with open(dest_file,'w') as fo:
        fo.write(infil)

    if AFLOWpi.prep._ConfigSectionMap('prep','copy_execs').lower()=='false':
        engineDir = AFLOWpi.prep._ConfigSectionMap('prep','engine_dir')
        if os.path.isabs(enginedir) == False:
            configFileLocation = AFLOWpi.prep._getConfigFile()
            enginedir =  os.path.join(configFileLocation, enginedir)

        execPath=os.path.join(engineDir,'dist.x')
    else:
        execPath='./dist.x'
    try:
        subprocess.Popen('%s < %s_dist.in > /dev/null' % (execPath,ID),cwd=oneCalc['_AFLOWPI_FOLDER_'],shell=True).wait()
        if outp==True:
            os.rename(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'dist.out'),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_output_dist.out'))
        else:
            os.rename(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'dist.out'),os.path.join(oneCalc['_AFLOWPI_FOLDER_'],ID+'_input_dist.out'))

    except Exception as e:
        logging.error('dist.x did not run properly')
        print('ERROR: dist.x did not run properly')

