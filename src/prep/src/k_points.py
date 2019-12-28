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
import numpy


def _gen_nosym_kgrid(nk1,nk2,nk3,sk1=0,sk2=0,sk3=0,as_string=False,scale=1.0):

    shift1 = 1.0/nk1*float(sk1)/2.0
    shift2 = 1.0/nk2*float(sk2)/2.0
    shift3 = 1.0/nk3*float(sk3)/2.0

    kdist1 = 1.0/nk1
    kdist2 = 1.0/nk2
    kdist3 = 1.0/nk3

    tot = nk1*nk2*nk3
    nk_str = '%s'%tot
#    if shift_gamma==True:
    shift1-=0.5
    shift2-=0.5
    shift3-=0.5
    if as_string==True:
        for i in range(nk1):
            for j in range(nk2):
                for k in range(nk3):
                    nk_str+='\n%8.8f %8.8f %8.8f'%(((float(i)*kdist1+shift1)*scale),
                                                   ((float(j)*kdist2+shift2)*scale),
                                                   ((float(k)*kdist3+shift3)*scale))
    else:
        nk_str=[]#numpy.zeros([tot,3])
        for i in range(nk1):
            for j in range(nk2):
                for k in range(nk3):
                    nk_str.append([(float(i)*kdist1+shift1)*scale,
                                   (float(j)*kdist2+shift2)*scale,
                                   (float(k)*kdist3+shift3)*scale])
    nk_str=numpy.asarray(nk_str)

    return nk_str








