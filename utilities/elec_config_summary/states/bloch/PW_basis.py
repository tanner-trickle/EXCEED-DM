# Prints information about the bloch, PW basis states

import numpy as np

def print_bloch_PW_basis_summary(data, type = 'init'):
    
    if f'elec_states/{type}/bloch/PW_basis' in data:
        
        PW_basis_data = data['elec_states'][type]['bloch']['PW_basis']
        
        # config
        G_list_red = np.transpose(PW_basis_data['config']['G_list_red'][...])
        
        n_G = len(G_list_red)
        
        G_grid_min_x = np.amin(G_list_red[:, 0])
        G_grid_min_y = np.amin(G_list_red[:, 1])
        G_grid_min_z = np.amin(G_list_red[:, 2])
        
        G_grid_max_x = np.amax(G_list_red[:, 0])
        G_grid_max_y = np.amax(G_list_red[:, 1])
        G_grid_max_z = np.amax(G_list_red[:, 2])
        
        print("\tPW basis:")
        print()
        
        print("\t\tConfig:")
        print()
        print(f"\t\t\tNumber of G vectors: {n_G}")
        print(f"\t\t\t\tG grid min: {G_grid_min_x}, {G_grid_min_y}, {G_grid_min_z}")
        print(f"\t\t\t\tG grid max: {G_grid_max_x}, {G_grid_max_y}, {G_grid_max_z}")
        print()
        
        # states
        
        n_states = len(PW_basis_data['state_info']['energy_list'][...])
        
        n_k = np.amax(PW_basis_data['state_info']['k_id_list'][...])
        n_bands = int(n_states/n_k)
        
        i_min = np.amin(PW_basis_data['state_info']['i_list'][...])
        i_max = np.amax(PW_basis_data['state_info']['i_list'][...])
        
        E_min = np.amin(PW_basis_data['state_info']['energy_list'][...])
        E_max = np.amax(PW_basis_data['state_info']['energy_list'][...])
        
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