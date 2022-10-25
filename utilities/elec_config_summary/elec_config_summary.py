# Prints a summary of the states in the electronic configuration file.

import sys
import h5py

# bloch
from states.bloch.PW_basis import print_bloch_PW_basis_summary
from states.bloch.STO_basis import print_bloch_STO_basis_summary
from states.bloch.single_PW import print_bloch_single_PW_summary

def print_elec_config_summary(filename):
    
    print("-- Electronic Configuration File Summary --")
    print()
    print(f"Filename: {filename}")
    print()
    
    with h5py.File(filename, 'r') as data:
        
        print("Initial States:")
        print()
        
        # bloch
        print_bloch_PW_basis_summary(data, type = 'init')
        print_bloch_STO_basis_summary(data, type = 'init')
        print_bloch_single_PW_summary(data, type = 'init')
        
        print("Final States:")
        print()
        
        # bloch
        print_bloch_PW_basis_summary(data, type = 'fin')
        print_bloch_STO_basis_summary(data, type = 'fin')
        print_bloch_single_PW_summary(data, type = 'fin')

if __name__ == "__main__":
    
    elec_config_filename = sys.argv[1]
    
    print_elec_config_summary(elec_config_filename)