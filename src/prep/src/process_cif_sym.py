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
import numpy as np
import re



############################################################
############################################################
def _process_cif_sym_shift(symops):
    re_shift = re.compile('[+-]*\d\/\d')
    addition = re_shift.findall(symops)[0].strip(' ')
    sign=1.0
    try:
        minus = addition.index('-')
        sign=-1.0
    except:
        pass
    numer=float(addition[addition.index('/')-1])
    denom=float(addition[addition.index('/')+1])


    return sign*numer/denom
############################################################
############################################################

def _process_cif_sym(symops):

        ident = np.identity(3,dtype=np.float64)
        operations = np.zeros((len(symops),3,3))
#        operations[:,]=ident
        shift = np.zeros((len(symops),3))

        re_frac         = re.compile('([\+\-]*\d+/\d+)')
        re_frac_replace = re.compile('[\+\-]*\d+/\d+')

        for i in range(len(symops)):

            if symops[i][0][0]!='-' and symops[i][0][0]!='+':
                symops[i][0]='+'+symops[i][0]

            if symops[i][1][0]!='-' and symops[i][1][0]!='+':
                symops[i][1]='+'+symops[i][1]

            if symops[i][2][0]!='-' and symops[i][2][0]!='+':
                symops[i][2]='+'+symops[i][2]

            symops[i][0]=symops[i][0].upper().strip()
            symops[i][1]=symops[i][1].upper().strip()
            symops[i][2]=symops[i][2].upper().strip()


            try:
                shift[i][0]=AFLOWpi.prep._process_cif_sym_shift(re_frac.findall(symops[i][0])[0])
                symops[i][0] = re_frac.sub('',symops[i][0])
            except:pass
            try:
                shift[i][1]=AFLOWpi.prep._process_cif_sym_shift(re_frac.findall(symops[i][1])[0])
                symops[i][1]=re_frac.sub('',symops[i][1])
            except: pass
            try:
                shift[i][2]=AFLOWpi.prep._process_cif_sym_shift(re_frac.findall(symops[i][2])[0])
                symops[i][2]=re_frac.sub('',symops[i][2])
            except: pass
                




                
            symops[i][0] = symops[i][0].replace('X','np.array([[1.0],[0.0],[0.0]])')
            symops[i][0] = symops[i][0].replace('Y','np.array([[0.0],[1.0],[0.0]])')
            symops[i][0] = symops[i][0].replace('Z','np.array([[0.0],[0.0],[1.0]])')

            symops[i][1] = symops[i][1].replace('X','np.array([[1.0],[0.0],[0.0]])')
            symops[i][1] = symops[i][1].replace('Y','np.array([[0.0],[1.0],[0.0]])')
            symops[i][1] = symops[i][1].replace('Z','np.array([[0.0],[0.0],[1.0]])')

            symops[i][2] = symops[i][2].replace('X','np.array([[1.0],[0.0],[0.0]])')
            symops[i][2] = symops[i][2].replace('Y','np.array([[0.0],[1.0],[0.0]])')
            symops[i][2] = symops[i][2].replace('Z','np.array([[0.0],[0.0],[1.0]])')


            operations[i,:,0] =  eval(symops[i][0])[:,0]
            operations[i,:,1] =  eval(symops[i][1])[:,0]
            operations[i,:,2] =  eval(symops[i][2])[:,0]



            
        return shift,operations
        
        
