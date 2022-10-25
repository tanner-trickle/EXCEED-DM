# Prints information about the bloch, STO basis states

import numpy as np

def print_bloch_STO_basis_summary(data, type = 'init'):
    
    if f'elec_states/{type}/bloch/STO_basis' in data:
        
        STO_basis_data = data['elec_states'][type]['bloch']['STO_basis']
        
        # config
        n_r_vec_grid = STO_basis_data['config']['n_r_vec_grid'][...]
        n_x_grid = STO_basis_data['config']['n_x_grid'][...]
        
        print("\tSTO basis:")
        print()
        
        print("\t\tConfig:")
        print()
        print(f"\t\t\tNumber of lattice vectors in sum: {n_r_vec_grid}")
        print(f"\t\t\tNumber of x vectors in unit cell: {n_x_grid}")
        print()
        
        # states
        
        n_states = len(STO_basis_data['state_info']['energy_list'][...])
        
        n_k = np.amax(STO_basis_data['state_info']['k_id_list'][...])
        n_bands = int(n_states/n_k)

        i_min = np.amin(STO_basis_data['state_info']['i_list'][...])
        i_max = np.amax(STO_basis_data['state_info']['i_list'][...])
        
        E_min = np.amin(STO_basis_data['state_info']['energy_list'][...])
        E_max = np.amax(STO_basis_data['state_info']['energy_list'][...])
        
        n_min = np.amin(np.transpose(STO_basis_data['state_info']['nlm_list'][...])[:, 0])
        n_max = np.amax(np.transpose(STO_basis_data['state_info']['nlm_list'][...])[:, 0])
        
        print("\t\tState Info:")
        print()
        print(f"\t\t\tNumber of states: {n_states}")
        print(f"\t\t\t\tNumber of k points: {n_k}")
        print(f"\t\t\t\tNumber of bands: {n_bands}")
        print(f"\t\t\t\tMinimum principal quantum number, n: {n_min}")
        print(f"\t\t\t\tMaximum principal quantum number, n: {n_max}")
        if type == 'init':
            print(f"\t\t\t\t\tIndex range: {i_min} -> {i_max}")
        print()
        print(f"\t\t\tMinimum energy state: {E_min}")
        print(f"\t\t\tMaximum energy state: {E_max}")
        print()