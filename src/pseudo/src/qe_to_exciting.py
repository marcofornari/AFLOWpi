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


def __qe_to_exiting_scf_input(qe_input):
    cell = AFLOWpi.retr.getCellMatrixFromInput(qe_input).getA()
    cell_list=[]
    
    for i in range(len(cell)):
        cell_list.append(' '.join([str(i) for i in cell[i]])) 

    ID='test'

    exciting_input='''<title>%s</title>''' %ID
    
    species_pot='../../../species/'

    exciting_input+='''
      <structure speciespath="%s">'''%species_pot
    exciting_input+='''
        <crystal scale="1.0" >
          <basevect>%s</basevect>
          <basevect>%s</basevect>
          <basevect>%s</basevect>
    </crystal>''' % (cell_list[0],cell_list[1],cell_list[2])


    pos = AFLOWpi.retr._getPositions(qe_input,matrix=True).getA()
    lab = AFLOWpi.retr._getPosLabels(qe_input)

    pos_dict={}

    for i in range(len(lab)):
        pos_str=' '.join([str(j) for j in pos[i]])
        if lab[i] not in list(pos_dict.keys()):
          
            pos_dict[lab[i]]=[pos_str]
        else:
            pos_dict[lab[i]].append(pos_str)

    for species,positions in list(pos_dict.items()):
        exciting_input+='''</species>
        <species speciesfile="%s.xml">'''%species
        for i in range(len(positions)):
            exciting_input+='''
            <atom coord= bfcmt="%s"></atom>'''%positions[i]

        exciting_input+='''
            </species>
        '''
    exciting_input+='''</structure>
        '''

    input_dict = AFLOWpi.retr._splitInput(qe_input)

    k_list = input_dict['K_POINTS']['__content__'].split()
    
    g_dense = k_list[:3]
    g_shift = [float(i)/2.0 for i in k_list[3:]]

    exciting_input+='''<groundstate vkloff="%s  %s  %s" ngridk="%s %s %s"'
    mixer='msec' nosource="true" tforce="true"></groundstate>
    ''' % (g_shift[0],g_shift[1],g_shift[2],g_dense[0],g_dense[1],g_dense[2])


    print(exciting_input)
