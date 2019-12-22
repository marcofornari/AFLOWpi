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
import scipy
import os
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot


def _gen_smear_conv_calcs(oneCalc,ID,num_points=4,smear_type='mp',smear_variance=0.3,calc_type='scf'):

    input_dict = AFLOWpi.retr._splitInput(oneCalc["_AFLOWPI_INPUT_"])

    try:
        initial_degauss = input_dict['&system']['degauss']
    except:
        initial_degauss = '0.01'

        input_dict['&system']['occupations']='"smearing"'
        if 'smearing' not in list(input_dict['&system'].keys()):
            input_dict['&system']['smearing']='"%s"'%smear_type

    input_dict['&system']['degauss']=initial_degauss

    fl_degauss = float(initial_degauss)
    ub_degauss = fl_degauss*(1.0+smear_variance)
    lb_degauss = fl_degauss*(1.0-smear_variance)
    
    ub_degauss = fl_degauss*(5.0)
    lb_degauss = fl_degauss*(0.5)

    ub_degauss = fl_degauss+0.005
    lb_degauss = fl_degauss-0.005
    degauss_vals = numpy.linspace(lb_degauss,ub_degauss,num_points)


    input_files = []

    for val in degauss_vals:
        input_dict['&system']['degauss'] = val
        input_dict['&control']['calculation'] = "'%s'"%calc_type
        one_input = AFLOWpi.retr._joinInput(input_dict)

        input_files.append(one_input)


    return input_files



def _extrapolate_smearing(oneCalc,ID):

    folder = oneCalc['_AFLOWPI_FOLDER_']
    calc_log = os.path.join(folder,'SMEARING','AFLOWpi','calclogs','step_01.log')

    calcs=AFLOWpi.prep._load_log_from_filename(calc_log)

    de_set=[]
    smear_vals=[]

    calcs = AFLOWpi.retr.grabEnergyOut(calcs)
    for k,v in list(calcs.items()):

        oc_in_dict = AFLOWpi.retr._splitInput(v['_AFLOWPI_INPUT_'])
        smear_type = oc_in_dict['&system']['smearing']
        smear_val = float(oc_in_dict['&system']['degauss'])
        en         = float(v['Energy'])
#        force=AFLOWpi.retr.getForce(v,k,string=False)
        if en != 0.0:
#            de_set.append(force)
            de_set.append(en)
            smear_vals.append(smear_val)


    de_set = numpy.asarray(de_set)
    #CHECK THIS
    smear_vals = numpy.asarray(smear_vals)
#    (x_0,x_1) = 
#    de_set=numpy.abs(de_set)
    max_de = numpy.amax(smear_vals)
    min_de = numpy.amin(smear_vals)

    if min_de<=0:
        min_end=2.0*min_de
    else:
        min_end=-2.0*min_de

    if max_de<=0:
        max_end=-2.0*max_de
    else:
        max_end= max_de*2.0

    extrap_smear= numpy.linspace(0.0,max_end,200)
#    interp_en = interp_en.tolist()
#    interp_en.append(0.0)
#    interp_en = numpy.asarray(interp_en)

    params=[smear_vals,de_set,]

#    extrap_smear = numpy.polyval(numpy.polyfit(en_set,smear_vals,2),interp_en)

#    en_data_diff=numpy.abs(numpy.gradient(params[1]))
#    en_data_diff=numpy.gradient(params[1])
    fit = numpy.polyfit(params[0],params[1],2)
    fit_en =numpy.polyval(fit,extrap_smear)
#    en_deriv = numpy.polyfit(params[0],en_data_diff,1)
#    print fit
    en_deriv=numpy.polyder(fit)
#    en_deriv=fit

    print(en_deriv)
#    en_deriv=fit.deriv()
#    fir1d= numpy.polyfit(params[0],params[1],1)
    interp_en =numpy.polyval(en_deriv,extrap_smear)

#    interp_en = en_deriv(extrap_smear)
#    at_zero = en_deriv[0]
    at_zero=en_deriv[0]/en_deriv[1]
    print(at_zero)
    print(de_set)
#    print en_data_diff
    inputDict = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])

    if at_zero<0.001:
        #if we don't need smearing get rid of it.
        try:
            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','degauss',at_zero,del_value=True)
        except:
            pass
        try:
            #CHECK THIS!!!!
            smu = eval(inputDict['&system']['occupations'])
            if smu == 'smearing':
                oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','occupations','"smearing"',del_value=True)

        except:
            pass
        try:
            oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','smearing',smear_type,del_value=True)
        except: 
            pass
        

    else:
        #if we need smearing
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','degauss',at_zero,del_value=False)
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','occupations','"smearing"',del_value=False)
        oneCalc,ID = AFLOWpi.prep._modifyNamelistPW(oneCalc,ID,'&system','smearing',smear_type,del_value=False)


    return_input = AFLOWpi.retr._joinInput(inputDict)
    oneCalc['_AFLOWPI_INPUT_']=return_input


    matplotlib.rc("font", family="serif")      #to set the font type
    matplotlib.rc("font", size=10)             #to set the font size

    matplotlib.pyplot.plot(smear_vals,de_set,"b.")
    matplotlib.pyplot.plot(extrap_smear,fit_en,"r")
#    matplotlib.pyplot.plot(smear_vals,en_data_diff,"b.")
#    matplotlib.pyplot.plot(extrap_smear,interp_en,"r")
    matplotlib.pyplot.ylabel("|$\Delta$E|")
    matplotlib.pyplot.xlabel("degauss (Ry)")
    subdir=oneCalc['_AFLOWPI_FOLDER_']
    fileplot = os.path.join(subdir,'SMEARING_%s_%s.pdf' % (AFLOWpi.retr._getStoicName(oneCalc,strip=True),ID))

    matplotlib.pyplot.savefig(fileplot)

    return oneCalc,ID


