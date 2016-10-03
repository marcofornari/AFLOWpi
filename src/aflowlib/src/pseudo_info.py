import AFLOWpi
import re
import os
import logging
import shutil
import cPickle
import glob
import fnmatch
import copy 
import contextlib
import sys
import cStringIO
import subprocess


def __getPPFileString(oneCalc,ID):
        pseudoPathList = []
	for k,v in oneCalc.iteritems():
		pseudodir = AFLOWpi.prep._ConfigSectionMap('prep','pseudo_dir')
		if os.path.isabs(pseudodir) == False:
			configFileLocation = AFLOWpi.prep._getConfigFile()
			configFileLocation = os.path.dirname(configFileLocation)
			pseudodir =  os.path.join(configFileLocation, pseudodir)
	
		if len(re.findall(ur"_AFLOWPI_[A-Za-z0-9]{1,4}PSEUDO_",k)) !=0:
			pseudoPath =  os.path.join(pseudodir,v)

			if os.path.exists(pseudoPath):               
				pseudoPathList.append(pseudoPath)

	try:
		with open(pseudoPathList[0],'r') as ppFileObj:
			ppFileString = ppFileObj.read()
	except:
		logging.error('Could not find PP file %s')
	return ppFileString

def __grab__dft_type(oneCalc,ID):
	outString = ''
	ppFileString = __getPPFileString(oneCalc,ID)
	outString = __DFTType(ppFileString)

	return outString


def __DFTType(ppFileString):
	funPP = re.compile(ur'[Ff]unctional\s*[:=]*\s*(\w*)[\s\n]',)
	typePP = re.compile(ur'Pseudopotential type:\s*(\w*)[\s\n]',)

	if len(re.findall(ur'Generated using "atomic" code',ppFileString))!=0:
		logging.info("Using 'atomic' Generated PP")
#		print "Using 'atomic' Generated PP"
		funPP = re.compile(ur'[Ff]unctional\s*[:=]*\s*(\w*)[\s\n]',)
		typePP = re.compile(ur'Pseudopotential type:\s*(\w*)[\s\n]',)
	if len(re.findall(ur'Generated using Vanderbilt code',ppFileString))!=0:
		logging.info('Using Vanderbilt Generated PP')
#		print "Using Vanderbilt Generated PP"
		typePP = re.compile(ur'\s(\w*)\s*Exchange-Correlation functional')
		funPP = re.compile(ur'^\s*(\w*).*Exchange-Correlation functional',re.M)
	if len(re.findall(ur'converted with fhi2upf.x',ppFileString))!=0 or len(re.findall(ur'Generated using FHI98PP',ppFileString))!=0:       
		funPP = re.compile(ur'[Ff]unctional\s*[:=]*\s*([\w_-]*)[\s\n]',re.M)

	if len(re.findall(ur'Generated in analytical, separable form',ppFileString))!=0:
		funPP = re.compile(ur'[Ff]unctional\s*[:=]*\s*(\w*)-*',re.M)		

	if len(re.findall(ur'metaGGA',ppFileString))!=0:
		funPP = re.compile(ur'^\s*(\w*).*Exchange-Correlation functional',re.M)
		typePP = re.compile(ur'\s(\w*)\s*Exchange-Correlation functional')
		
	if len(re.findall(ur'Exchange-Correlation functional',ppFileString))!=0:
		funPP = re.compile(ur'^\s*(\w*).*Exchange-Correlation functional',re.M)
		typePP = re.compile(ur'\s(\w*)\s*Exchange-Correlation functional')
		
	funPPList= funPP.findall(ppFileString)
	typePPList = typePP.findall(ppFileString)

	try:
		if typePPList[0]=='USPP':
			typePPList[0]='US'
		if typePPList[0]=='META':
			typePPList[0]='GGA'

	except Exception,e:

		pass
	try:
		if len(re.findall(ur'SLA',funPPList[0]))!=0:
			funPPList[0]='LDA'
	except Exception,e:

		pass

	return ','.join(funPPList)+','+','.join(typePPList)
