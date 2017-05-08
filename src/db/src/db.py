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

import sqlite3
import AFLOWpi
import copy
from collections import OrderedDict
import os
import re


def _cleanKeyDict(oneCalc):
    cleanRegex=re.compile(r'[!@#$%^&*~{}\[\]\'?><(),.;"]+')
    keyCleanDict=OrderedDict()
    keyClean=''
    for key in oneCalc.keys():
        if len(cleanRegex.findall(key)):
            keyCleanSplit=cleanRegex.split(key)
            keyClean=''.join(keyCleanSplit)
            keyCleanDict[key]=keyClean
        else:
            keyCleanDict[key]=key

    return keyCleanDict


def read_db(db_path):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    columns = c.execute('select * from main')
    column_names = [d[0] for d in columns.description]

    calc_data = {}

    for row in columns:

            # build dict
            info = dict(zip(column_names, row))
            try:
                calc_data[row[0]] =info
            except Exception,e:
                print e


    for k,v in calc_data.iteritems():
        for j,l in v.iteritems():
            print '%s %-20s : %s'%(k,j,l)
    
def _dictToTable(inputDict,dbPath,tableName):
    
    conn = sqlite3.connect(dbPath)
    c = conn.cursor()
    # Create table
    c.execute('''CREATE TABLE IF NOT EXISTS %s
           (ID TEXT NOT NULL PRIMARY KEY)''' % tableName)
    sample = inputDict.items()[0][1]
    cleanedKey = AFLOWpi.db._cleanKeyDict(sample)
    cleanedKey['ID']='ID'

    for key in sample.keys():
        try:
            c.execute('SELECT %s FROM %s' % (cleanedKey[key],tableName))
        except Exception,f:
#            print 1            
            try:
                c.execute('''ALTER TABLE %s ADD COLUMN %s TEXT;''' % (tableName,cleanedKey[key])) 
                
            except sqlite3.OperationalError,e:
#                print 2
                print f
                print e
                
 

    for ID,oneCalc in inputDict.iteritems():
        inputDict[ID]['ID']=ID
        for key,value in oneCalc.iteritems():
            try:
                exeStr = "INSERT OR REPLACE INTO %s" % tableName +" {}".format(tuple('%s' % cleanedKey[key] for key in oneCalc.keys()))+"VALUES {}".format(tuple('%s' % value for value in oneCalc.values()))
#                print exeStr
                c.execute(exeStr)

            except Exception,e:
                print 3
                print e
                pass



        del oneCalc['ID']
        
        
    conn.commit()
    conn.close()


def export(project,set_name='',config='',author='',affiliation='',exclude_steps=[],exclude_properties=[],discard_incomplete=True,db_path='./'):
    '''
    Parses the espresso input,espresso output,the AFLOWpi dictionary for 
    each calculation (the AFLOWpi metadata), and the aflowkeys and populates
    a database in the AFLOWpi folder of the project and set (if there is one)
    with data regarding each calculation with the primary keys of the
    database being the calculation ID hash checksum of the espresso input 
    file. All entries in the database are strings (for consistency) 
    Fields that have no value from the regular espression search are
    entered as blank strings. (for example in a non variable cell scf there
    would be no final cell or atomic position parameters so they would be 
    left as empty strings in the database)


    Arguments:
     - calcs -- dictionary of dictionaries of calculations
     - aflowkeys -- the aflowkeys
    '''


    print "Extracting output from calculations and saving to database"
   
    # aflowcopy = copy.deepcopy(aflowkeys)
    # aflowcopy['calc_set']=aflowcopy['set']
    # del aflowcopy['set']
 
    aflowDict=OrderedDict()
    # inputDict=OrderedDict()
    # outputDict=OrderedDict()
    """
    build a dictionary from the parsed input files for each of 
    the calculation with the key being the ID of each calculation
    and the value being the parsed input file dictFlag=True so that
    the function __parseRef returns a dictionary instead of an input
    string (if it were the default False)
    """
#    sampleDict = calcs[random.choice(calcs.keys())]          
#    baseDir = os.path.dirname(sampleDict['_AFLOWPI_FOLDER_'])
    if set_name=="":
        dbName = '%s_%s.db' % (project,set_name)
    else:
        dbName = '%s_%s.db' % (project,set_name)

    db_path = os.path.join(db_path,dbName)

#    for ID,oneCalc in calcs.iteritems():
#        inputDict[ID] = AFLOWpi.prep._parseRef(oneCalc['_AFLOWPI_INPUT_'],dictFlag=True)

#    for ID,oneCalc in calcs.iteritems():        
#        outputDict[ID] = AFLOWpi.retr.getCellOutput(ID,oneCalc['_AFLOWPI_FOLDER_'])

#    for ID,oneCalc in calcs.iteritems():
#        aflowDict[ID] = aflowcopy

    aflowDict = AFLOWpi.aflowlib.export(project,set_name=set_name,config=config,author=author,affiliation=affiliation,exclude_steps=[],exclude_properties=[],discard_incomplete=discard_incomplete,write=False)

    

    AFLOWpi.db._dictToTable(aflowDict,db_path,'main')
#    _dictToTable(calcs,dbPath,'meta')

#    try:
#        _dictToTable(aflowDict,dbPath,'aflowkeys')
#    except:
#        print 'Error in parsing aflowkeys'

#    try:
#        _dictToTable(outputDict,dbPath,'output')
#        print "Successfully saved calculation results to database"
#    except:
#        print 'Error in parsing calculation output. One or more calculations are missing proper output'
    
        

