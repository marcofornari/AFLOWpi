import AFLOWpi


def _setup_iterative(calcs,conv_func,false_func,true_func):

    add = 'conv_bool = ' + conv_func
    AFLOWpi.prep.addToAll_(calcs,block='ITERATIVE',addition=add)

    add = 'if not conv_bool:'
    AFLOWpi.prep.addToAll_(calcs,block='ITERATIVE',addition=add)

    add = '    '+false_func
    AFLOWpi.prep.addToAll_(calcs,block='ITERATIVE',addition=add)

    for ID,oneCalc in calcs.iteritems():
        add = '''    AFLOWpi.run._submitJob("%s",oneCalc,__submitNodeName__,forceOneJob=False)'''%ID
        AFLOWpi.prep._addToBlock(oneCalc,ID,'ITERATIVE',add)


    add = 'else:'
    AFLOWpi.prep.addToAll_(calcs,block='ITERATIVE',addition=add)

    add = '    '+true_func
    AFLOWpi.prep.addToAll_(calcs,block='ITERATIVE',addition=add)


