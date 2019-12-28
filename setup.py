# ***************************************************************************
# *                                                                         *
# *          AFLOWpi - Central Michigan University University, 2017         *
# *                                                                         *
# ***************************************************************************
#
#  Copyright 2017 - Andrew Supka and Marco Fornari - AFLOW.ORG consortium
#
#  This file is part of AFLOW software.
#
#  AFLOW is free software: you can redistribute it and/or modify
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


import sys
import os
from distutils.core import setup, Command,Extension
from distutils.util import change_root, convert_path
from distutils.command.install import install
import shutil
sys.path.append(os.path.curdir)
import site

#os.environ["CC"] = "gcc"
#os.environ["CXX"] = "g++"




ext_modules=[Extension("cints",sources=["src/scfuj/extensions/cints.c"]) ]


import glob
import re
import os
import inspect

orig_dir=os.path.abspath(os.curdir)
def __init__gen(src_folder):


   os.chdir('%s'%src_folder)
   files = glob.glob('./*.py')

   new_file=''

   for i in files:
       try:
           if '__init__'!=i[2:-3]:
               with open(i,'r') as pypy:
                   py_string  = pypy.read()

               file_name=os.path.basename(i)[:-3]

               functions = re.findall(r'^def\s*(\w+)\(.*\):',py_string,re.M)
               for j in functions:
                  new_file+='from  .%s import %s\n'%(file_name,j)

               classes = re.findall(r'^class (\w+)(?:\(.*?\))?:',py_string,re.M)
               for j in classes:
                  new_file+='from  .%s import %s\n'%(file_name,j)

       except Exception as e:
           print('CRITICAL ERROR DURING GENERATION OF __init__.py') 
           print(e)
           raise SystemExit

   with open('./__init__.py','w') as initFile:
       initFile.write(new_file)


   os.chdir(orig_dir)

modules=['prep','run','pseudo','plot','retr','plot','aflowlib','db','scfuj','elph',"environ"]
for j in modules:
   try:
      os.remove('./src/%s/src/__init__.py'%j)
   except:
      pass

for j in modules:
   try:
      print(('generating ./src/%s/src/__init__.py'%j))
      __init__gen('./src/%s/src/'%j)

   except Exception as e:
      print(e)


for j in modules:
   try:
      os.remove('./src/%s/tests/__init__.py'%j)
   except:
      pass

for j in modules:
   try:

      print(('generating ./src/%s/tests/__init__.py'%j))
      __init__gen('./src/%s/tests/'%j)      
   except Exception as e:
      pass


try:
   


   PAOPY_SRC = glob.glob('src/PAOFLOW/src/*.py')
   PAOPY_DEF = glob.glob('src/PAOFLOW/src/defs/*.py')

   ACBN0=['scfuj/acbn0_support/integs.py',
          'scfuj/acbn0_support/acbn0.py',
          'scfuj/acbn0_support/Molecule.py',
          'scfuj/acbn0_support/pyints.py',]

   setup(name = "AFLOWpi",
         version = "0.9.9",
         description = "Medium Throughput Framework for Quantum Espresso",
         author = "Andrew Supka,Marco Fornari",
         author_email = "supka1ar@cmich.edu",



         packages = ['AFLOWpi',
                   'AFLOWpi.qe',
                   'AFLOWpi.scfuj',
                   'AFLOWpi.pseudo',
                   'AFLOWpi.run',
                   'AFLOWpi.retr',
                   'AFLOWpi.plot',
                   'AFLOWpi.db',
                   'AFLOWpi.aflowlib',
                   'AFLOWpi.prep',
                   'AFLOWpi.environ',
                   'AFLOWpi.elph',
                   "PAOFLOW",
                   "PAOFLOW.defs",
                     
],
         package_dir = {'AFLOWpi'       :'src',
                      'AFLOWpi.qe'      :'src/qe',
                      'AFLOWpi.prep'    :'src/prep/src/',
                      'AFLOWpi.run'     :'src/run/src/',
                      'AFLOWpi.scfuj'   :'src/scfuj/src/',
                      'AFLOWpi.retr'    :'src/retr/src/',
                      'AFLOWpi.plot'    :'src/plot/src/',
                      'AFLOWpi.pseudo'  :'src/pseudo/src/',
                      'AFLOWpi.db'      :'src/db/src/',
                      'AFLOWpi.elph'      :'src/elph/src/',
                      'AFLOWpi.environ'      :'src/environ/src/',
                      'AFLOWpi.aflowlib':'src/aflowlib/src/',
                      'PAOFLOW'         :'src/PAOFLOW/',   
                      'PAOFLOW.defs'         :'src/PAOFLOW/src/'   
},


                      
         package_data = {
                                  'AFLOWpi':['PAOFLOW/src/*/*','PAOFLOW/src/*.py','PAOFLOW/examples/*.py',
                                             'scfuj/acbn0_support/*','AFLOWSYM/*'],

                                  },


       ext_modules = ext_modules,
       long_description = """Install Script for AFLOWpi""",) 
   



except Exception as e:
   print(e)
   print('Something went wrong...exiting')
   exit



def binaries_directory():
      """Return the installation directory, or None"""
      # taken from stackoverflow
      if '--user' in sys.argv:
         paths = (site.getusersitepackages(),)
      else:
         py_version = '%s.%s' % (sys.version_info[0], sys.version_info[1])
         paths = (s % (py_version) for s in (
            sys.prefix + '/lib/python%s/dist-packages/',
            sys.prefix + '/lib/python%s/site-packages/',
            sys.prefix + '/local/lib/python%s/dist-packages/',
            sys.prefix + '/local/lib/python%s/site-packages/',
            '/Library/Python/%s/site-packages/',
         ))
         

      for path in paths:
         if os.path.exists(path):
            return path
      print('no installation path found', file=sys.stderr)
      return None

try:
   inst_dir=binaries_directory()


   AFLOW_EXEC = os.path.join(inst_dir,'AFLOWpi','AFLOWSYM','aflow')
   if not os.access(AFLOW_EXEC,3):
      os.chmod(AFLOW_EXEC,733)      

except Exception as e:
    print(('Could not install AFLOW binary:',e))
   



