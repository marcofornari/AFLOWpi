import glob
import os 
import shutil
import AFLOWpi
import sys
#modules=['prep','run','pseudo','plot','retr','plot','aflowlib','db','scfuj',]
#don't include aflowlib and db modules for now
modules=['run','pseudo','plot','retr','plot','scfuj','prep',]

if os.path.exists('./build/'):
    shutil.rmtree('./build/')
os.mkdir('./build/')
os.mkdir('./build/AFLOWpi/')
sys.path.insert(0, os.path.abspath('../src'))


for mod in sorted(modules):
    mod_files=sorted(glob.glob('../src/%s/src/*.py'%mod))
    mod_files=[x for x in mod_files if '__init__.py' not in x]
    print(mod_files)
    single_mod_str=''
    for mod_file in mod_files:
        with open(mod_file,'r') as mod_file_obj:
            mod_file_str=mod_file_obj.read()
        single_mod_str+=mod_file_str+'\n'
    with open('./build/AFLOWpi/%s.py'%mod,'w') as mod_file_obj:
        mod_file_obj.write(single_mod_str)


shutil.copy('../src/__init__.py','./build/')
os.chdir('./build/')
os.system('sphinx-apidoc ./AFLOWpi/ -f -F  -H AFLOWpi -A "Andrew Supka" -V 1.1.0 -o ./')
shutil.copy('../conf_nap.py','./conf.py')
os.system('make latex')
os.system('make html')
#os.chdir('./_build/latex/')
#os.system('pdflatex AFLOWpi.tex')

shutil.copy('./AFLOWpi.pdf','../../../AFLOWpi_UserGuide.pdf')
