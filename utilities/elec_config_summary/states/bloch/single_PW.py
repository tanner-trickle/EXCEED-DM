# Prints information about the bloch, single PW states

import numpy as np

def print_bloch_single_PW_summary(data, type = 'init'):
    
    if f'elec_states/{type}/bloch/single_PW' in data:
        
        single_PW_data = data['elec_states'][type]['bloch']['single_PW']
        
        # config
        n_x_grid = single_PW_data['config']['n_x_grid'][...]
        
        print("\tSingle PW:")
        print()
        
        print("\t\tConfig:")
        print()
        print(f"\t\t\tNumber of x vectors in unit cell: {n_x_grid}")
        print()
        
        # states
        
        n_states = len(single_PW_data['state_info']['energy_list'][...])
        
        n_k = np.amax(single_PW_data['state_info']['k_id_list'][...])
        n_bands = int(n_states/n_k)

        i_min = np.amin(single_PW_data['state_info']['i_list'][...])
        i_max = np.amax(single_PW_data['state_info']['i_list'][...])
        
        E_min = np.amin(single_PW_data['state_info']['energy_list'][...])
        E_max = np.amax(single_PW_data['state_info']['energy_list'][...])
        
        print("\t\tState Info:")
        print()
        print(f"\t\t\tNumber of states: {n_states}")
        print(f"\t\t\t\tNumber of k points: {n_k}")
        print(f"\t\t\t\tNumber of bands: {n_bands}")
        if type == 'init':
            print(f"\t\t\t\t\tIndex range: {i_min} -> {i_max}")
        print()
        print(f"\t\t\tMinimum energy state: {E_min}")
        print(f"\t\t\tMaximum energy state: {E_max}")
        print()