import AFLOWpi
import os
import glob
import shutil
import re
import numpy 


def _grab_elastic_generated_inputs(oneCalc,ID):
    glob_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'Structures_ESPRESSO')
    input_files = glob.glob(glob_path+'/*.in')

    return input_files

def _prep_elastic(oneCalc,ID,eta_max=0.005,num_dist=49,use_stress = True,order=2):

    if ID in oneCalc['prev']:
        return

    try:
        os.rm(os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'ElaStic_PW.in'))
    except:
        pass

    split_input = AFLOWpi.retr._splitInput(oneCalc['_AFLOWPI_INPUT_'])
    pos = AFLOWpi.retr.detachPosFlags(split_input['ATOMIC_POSITIONS']['__content__'])[0]
    split_input['ATOMIC_POSITIONS']['__content__']=pos

    oneCalc['_AFLOWPI_INPUT_'] = AFLOWpi.retr._joinInput(split_input)
    
    write_in = oneCalc['_AFLOWPI_INPUT_']
    elastic_qe_in = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'ElaStic_PW.in')
    with open(elastic_qe_in,'w') as eqfo:
        eqfo.write(write_in)

    elastic_ID = ID+'_ElaStic'

    if use_stress:
        calc_meth=2
    else:
        calc_meth=1

    elastic_in_string='''%s
%s
%s
%s'''%(calc_meth,order,eta_max,num_dist)

    elastic_in_file = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.in'%elastic_ID)
    with open(elastic_in_file,'w') as eifo:
        eifo.write(elastic_in_string)

    os.system('ElaStic_Setup_ESPRESSO<%s'%elastic_in_file)



def _copy_qe_out_file(oneCalc,ID):
    stripped_step_ID='_'.join(ID.split('_')[:-1])
    dist_num = ID.split('_')[0]

    orig_path = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s.out'%ID)
    new_path  = os.path.join('../../%s/%s/'%(dist_num,stripped_step_ID),'%s.out'%(stripped_step_ID))

    try:
        shutil.copy(orig_path,new_path)
    except Exception,e:
        AFLOWpi.run._fancy_error_log(e)

def _grab_el_plot_data(file_path):
    fit_data_by_order={}

    with open(file_path,'r') as dfo:
        dfs = dfo.read()

    per_fit_order = re.findall('#.*(\d+).*\n((?:[-0-9.\s]*\n)*)',dfs)
    for order,fit_data in per_fit_order:
        fit_array=[]
        for i in fit_data.split('\n'):
            to_float = map(float,i.split())
            if len(to_float)!=0:
                fit_array.append(to_float)

        fit_data_by_order[order]=fit_array

    return fit_data_by_order

def _find_flat_region(data):
    eta = []
    val = []
    #Sort the data into two arrays. one for eta and one for the values
    for k,v in data:
        eta.append(k)
        val.append(v)

    #set the max variance from one step to the next for the plateau region.
    #if the variance from one step to the next is greater than the value
    #below then that is considered the end of the plateau.
    variance=10.0

    plateau_index=0
    #reverse the order or eta and values to start with lowest eta first
    val = [x for x in reversed(val)]
    eta = [x for x in reversed(eta)]

    #find difference between one point and the next
    grad = numpy.gradient(val).tolist()
    found=False
    for x in range(len(grad)):
        #if we find a values between one eta and the next is less than
        #our designated variance threshold then start looking for the end
        #if it. If we haven't found a plateau and the variance is greater
        #than the threshold keep looking for the plateau. If we have found
        #the plateau and we find the end of it then record the index of the
        #end of the plateau and use the eta at that index as our max eta
        if grad[x]<variance:
            if found==False:
                found=True

            plateau_index=x
        else:
            if found==False:
                continue
            else:
                break

    print 'plateau ends at index:',plateau_index

    return eta[plateau_index]


def _pp_elastic(oneCalc,ID,order=2,use_stress=True):

    res_analyze=''
    result_input_file = ''
    result_output_file = ''
    if use_stress:
        os.system('ElaStic_Analyze_Stress')
        if order==2:
            res_analyze = 'ElaStic_Result_Stress_2nd'
            result_input_file = 'ElaStic_2nd.in'
            result_output_file = 'ElaStic_2nd.out'
        if order==3:
            res_analyze = 'ElaStic_Result_Stress_3rd'
            result_input_file = 'ElaStic_3rd.in'
            result_output_file = 'ElaStic_3rd.out'
    else:
        os.system('ElaStic_Analyze_Energy')

        if order==2:
            res_analyze = 'ElaStic_Result_Energy_2nd'
            result_input_file = 'ElaStic_2nd.in'
            result_output_file = 'ElaStic_2nd.out'
        if order==3:
            res_analyze = 'ElaStic_Result_Energy_3rd'
            result_input_file = 'ElaStic_3rd.in'
            result_output_file = 'ElaStic_3rd.out'



    svs_folder = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'Stress-vs-Strain')

    fit_files = glob.glob(svs_folder+'/*_d1S.dat')


    files_sep_dist={}
    results_input_str=''
    for data_file in fit_files:
        distortion_ID = os.path.basename(data_file).split('_')[0]
        try:
            files_sep_dist[distortion_ID].append(data_file)
        except:
            files_sep_dist[distortion_ID]=[data_file]
    
    
    for distortion_ID,fit_files in files_sep_dist.iteritems():
        order_list=[]
        m_eta_list=[]
        for data_file in fit_files:
            distortion_ID = os.path.basename(data_file).split('_')[0]
            fit_data = AFLOWpi.run._grab_el_plot_data(data_file)
            #set a global max eta for all fittings
            global_max_eta=0.0
            #set the order of the polynomial to fit to.
            poly_fit_choice = 0

            if use_stress:
                fit_index=1
            else:
                fit_index=2
            #iterate through all orders of poly fit data
            for poly_order,data in fit_data.iteritems():
                #find max eta for one order of poly fit
                max_eta = AFLOWpi.run._find_flat_region(data)
                #check if the max eta is greater or equal to the current global
                #max eta. If the value is equal then prefer a higher order
                #polynomial fit to a lower one.
                if max_eta>=global_max_eta:
                    global_max_eta=max_eta
                    fit_index=poly_order
            #append the values for max eta and poly fit order 
            #to a list to later be used to create the the ElaStic_2nd.in
            #file.
            order_list.append(fit_index)
            m_eta_list.append(global_max_eta)
        #write the entry in the file
        results_input_str+=distortion_ID+' '+' '.join(map(str,m_eta_list))+'\n'
        results_input_str+=' '*len(distortion_ID)+' '.join(map(str,order_list))+'\n'

    #write the input for ElaStic_Result_Stress_2nd
    res_input = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],result_input_file)
    with open(res_input,'w') as rifo:
        rifo.write(results_input_str)

    #run ElaStic_Result_Stress_2nd or whatever result analyzer code
    #that is to be run.
    os.system(res_analyze)
    #move results to file with unique name
    res_output = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],result_output_file)
    new_res_output = os.path.join(oneCalc['_AFLOWPI_FOLDER_'],'%s_elastic.out'%ID)

    shutil.move(res_output,new_res_output)
