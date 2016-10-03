import re
import subprocess
import AFLOWpi

def __grab__hostname(oneCalc,ID):
    status = subprocess.Popen("uname -n",stdout=subprocess.PIPE,shell=True)
    hostname = status.communicate()[0]
    return hostname.strip('\n').strip()

def __grab__kernel_version(oneCalc,ID):
    status = subprocess.Popen("uname -r",stdout=subprocess.PIPE,shell=True)
    kernel_version = status.communicate()[0]
    return kernel_version.strip('\n').strip()

def __grab__architecture(oneCalc,ID):
    status = subprocess.Popen("uname -m",stdout=subprocess.PIPE,shell=True)
    arch = status.communicate()[0]
    return arch.strip('\n').strip()

def __grab__os(oneCalc,ID):
    status = subprocess.Popen("lsb_release -a",stdout=subprocess.PIPE,shell=True)
    output = status.communicate()[0]
    try:
        flavor = re.findall('.*Distributor ID:\s*(\w*)\s*',output)[0]
    except:
        flavor=''
    return flavor.strip('\n').strip()

def __grab__os_version(oneCalc,ID):
    status = subprocess.Popen("lsb_release -a",stdout=subprocess.PIPE,shell=True)
    output = status.communicate()[0]
    try:
        version = re.findall('.*Release:\s*([\w.-_]*)\s*',output)[0]
    except:
        version=''
    return version.strip('\n').strip()
