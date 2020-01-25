def invert_atom_pos_map(pos,equiv_atom,inv_flag):
    # finds mapping of atomic positions for inversion
    
    for i in range(inv_flag.shape[0]):

        if not inv_flag[i]:
            new_pos=pos
        else:
            new_pos=(-pos)%1.0

        remap_key=np.zeros(pos.shape[0],dtype=int)
        for j in range(pos.shape[0]):
            for k in range(new_pos.shape[0]):
                if np.all(np.isclose(pos[j],new_pos[k],rtol=1e-4,atol=1e-4,)):
                    remap_key[j] = k

        equiv_atom[i] = equiv_atom[i][remap_key]

    return equiv_atom
