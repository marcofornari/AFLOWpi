import AFLOWpi
import shutil
import os
import logging

def _copy_to_fig_dir(oneCalc,filename):
    proj_name=oneCalc["PROJECT"]
    set_name=oneCalc["SET"]
    
    fig_dir = AFLOWpi.prep._ConfigSectionMap('prep','fig_dir')
    calc_path=os.path.basename(oneCalc["_AFLOWPI_FOLDER_"])


    if fig_dir=="":
        return
    if set_name!="":
        save_path=os.path.join(fig_dir,proj_name,set_name)
    else:
        save_path=os.path.join(fig_dir,proj_name)


    if not os.path.exists(fig_dir):
        try:
            os.mkdir(fig_dir)
        except:
            AFLOWpi.run._fancy_error_log(e)
            return

    save_proj_path=os.path.join(fig_dir,proj_name)    
    if not os.path.exists(save_proj_path):
        try:
            os.mkdir(save_proj_path)
        except:
            return

    if not os.path.exists(save_path):
        try:
            os.mkdir(save_path)
        except:
            return

    save_path=os.path.join(save_path,calc_path)

    if not os.path.exists(save_path):
        try:
            os.mkdir(save_path)
        except:
            return

    shutil.copy(filename,save_path)
