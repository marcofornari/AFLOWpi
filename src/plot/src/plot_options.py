import AFLOWpi

def _get_plot_ext_type():
    plot_ext_type='.png'
    plot_ext_type = AFLOWpi.prep._ConfigSectionMap('plot','plot_file_type')
    if plot_ext_type == "":
            plot_ext_type = "png"   
    if plot_ext_type.lower() not in ["png","pdf"]:
            plot_ext_type = "png"   

    return plot_ext_type


def _get_title_option():
    plot_title_op=False
    plot_title_op = AFLOWpi.prep._ConfigSectionMap('plot','plot_file_type')
    if plot_title_op not in [False,True]:
        plot_title_op = False
        
    return plot_title_op
