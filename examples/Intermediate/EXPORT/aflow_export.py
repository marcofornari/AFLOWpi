import AFLOWpi

# export command loads each calculation in the set
# from the last step it successfully completed.
# and bundles data from all steps in a tar file 
# for easy export to the aflowlib database
AFLOWpi.aflowlib.export('export_example',
                        set_name='zincblende_mHT',
                        config='./do_calcs.config')
