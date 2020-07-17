import AFLOWpi

# export command loads each calculation in the set
# from the last step it successfully completed.
# and saves the data from all steps as a entry 
# in a sqlite database with each calculation 
# having a unique key ID
AFLOWpi.db.export('export_example',
                  set_name='zincblende_mHT',
                  config='./do_calcs.config')

# when we want to read the db and construct a dictionary  
# of dictionaries object that we can use for data mining
# over all the entries in the .db file
AFLOWpi.db.read_db('./export_example_zincblende_mHT.db')
