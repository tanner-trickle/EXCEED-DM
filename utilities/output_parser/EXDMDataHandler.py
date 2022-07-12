# Useful functions for parsing the output of an EXCEED-DMv1.0.0 calculation.

import h5py
import numpy as np

from itertools import groupby

def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)

def get_index(ele, li, eps = 1e-3, offset = 1):
    """
        Returns the index of the 'closest' element in list to ele.
    """
    index = np.abs(li - ele).argmin()
    diff = np.abs(li[index] - ele)
    if diff > eps:
        print('-- WARNING --')
        print('  Searching ', str(li), ' for ', str(ele), '. The closest value is',
              str(li[index]) )
    return index + offset

def get_index_2d(ele, li, eps = 1e-3, offset = 1):
    
    found_match = False
    
    for index, li_ele in enumerate(li):
        
        match = 0
        
        for e, ele_part in enumerate(ele):
            
            if np.abs( ele_part - li_ele[e] ) < eps:
                match += 1
                
        if match == len(ele):
            
            match_index = index + offset
            found_match = True
            
            break
            
    if found_match:
            
        return match_index
    
    else:
        
        print('-- ERROR --')
        print('  ', ele, 'is not in', li)
    
def rebin_1d(array, x_width, new_x_width):
    """
        Rebins a 1d array with previous width x_width and new width, new_x_width.
    """
    
    if new_x_width < x_width:
        
        print('--- WARNING ---')
        print('    New x_width <= old x_width and therefore cannot rebin. Returning original array.')
        
        return [ x_width, array ]
        
    elif new_x_width == x_width:
        
        return [ x_width, array ]
    
    else:
        
        # make sure that the new_x_width is an integer multiple of the x_width
        div_bin_width = int(new_x_width/x_width)

        rebinned_arr = [ array[i:i+div_bin_width] for i in range(0, len(array), div_bin_width) ]

        for i, ele in enumerate(rebinned_arr):
            if len(ele) != div_bin_width:
                rebinned_arr[i] = np.append(ele, np.zeros(div_bin_width - len(ele)))

        return [ new_x_width, np.sum(rebinned_arr, axis = 1) ]

class EXDMData:
    
    def __init__(self, filename = ''):
        
        self.filename = filename
        
        self.hdf5_data = h5py.File(self.filename, 'r')
        
    def get_masses_eV(self):
        return self.hdf5_data['dm_model/mX'][...]
        
    def get_masses_MeV(self):
        return 10.**(-6)*self.get_masses_eV()
    
    def get_med_FF(self):
        return self.hdf5_data['dm_model/med_FF'][...]
    
    def get_expt_M(self):
        return self.hdf5_data['experiment/M'][...]
        
    def get_expt_T(self):
        return self.hdf5_data['experiment/T'][...]

    def get_material_band_gap(self):
        return self.hdf5_data['material/band_gap'][...]
        
    def get_numerics_scatter_binned_rate_E_bin_width(self):
        return self.hdf5_data['numerics_binned_scatter_rate/E_bin_width'][...]
    
    def get_numerics_absorption_rate_widths(self):
        return np.transpose(self.hdf5_data['numerics_absorption_rate/widths'][...])

    def get_binned_scatter_rate_qE(self, 
                          mass_MeV = 1, 
                          med_FF = 2., 
                          vE = [0, 0, 240], 
                          i_list = [0]):
        
        mass_idx = get_index(mass_MeV, self.get_masses_MeV())
        med_FF_idx = get_index(med_FF, self.get_med_FF())
        
        if i_list == [0]:
            if len(self.get_med_FF()) > 1:

                binned_scatter_rate_qE = self.hdf5_data['binned_scatter_rate'][f'model_{med_FF_idx}'][f'mass_{mass_idx}']['total_binned_scatter_rate'][...]
                
            else:
                binned_scatter_rate_qE = self.hdf5_data['binned_scatter_rate'][f'mass_{mass_idx}']['total_binned_scatter_rate'][...]

        return binned_scatter_rate_qE + 10**(-100)
    
    def get_absorption_rates(self, g2 = 10**(-26), 
                             expt_M_kg = 1,
                             expt_T_year = 1, 
                             width = [0.2, 1, 0.2]):
        
        masses_eV = self.get_masses_eV()
        
        # normalization
        norm = g2*(expt_M_kg/self.get_expt_M())*(expt_T_year/self.get_expt_T())
        
        m_abs_rate = []
        for mass in masses_eV:
            
            mass_idx = get_index(mass, masses_eV)
            
            if len(self.get_numerics_absorption_rate_widths()) == 1:
                
                abs_rate = self.hdf5_data['absorption_rate'][f'mass_{mass_idx}']['absorption_rate'][...]
                
            else:
                
                width_idx = get_index_2d(width,self.get_numerics_absorption_rate_widths())
                
                abs_rate = self.hdf5_data['absorption_rate'][f'width_{width_idx}'][f'mass_{mass_idx}']['absorption_rate'][...]
                
            m_abs_rate.append([ mass, norm*np.abs(abs_rate) + 10**(-100) ])
            
        return np.transpose(np.array(sorted(m_abs_rate , key=lambda ele: ele[0])))
    
    def get_abs_g_constraint(self, 
                             n_cut = 3,
                             expt_M_kg = 1,
                             expt_T_year = 1, 
                             width = [0.2, 1, 0.2]):
        
        [ masses_eV, abs_rate ] = self.get_absorption_rates(g2 = 1, 
                             expt_M_kg = expt_M_kg,
                             expt_T_year = expt_T_year, 
                             width = width)
        
        return [ masses_eV, np.sqrt( n_cut / abs_rate ) ]
        
    def get_binned_scatter_rate_E(self, 
                          mass_MeV = 1, 
                          med_FF = 2., 
                          sigma_cm2 = 10**(-40), 
                          expt_M_kg = 1,
                          expt_T_year = 1, 
                          E_bin_width = 1):
        """
            Returns the (dimensionless) binned (in E) scatter rate for a given mass. 
        """
        
        # normalization
        norm = sigma_cm2*(expt_M_kg/self.get_expt_M())*(expt_T_year/self.get_expt_T())
        
        # sum over the q bins
        binned_scatter_rate_E = np.sum(self.get_binned_scatter_rate_qE(
                                            mass_MeV = mass_MeV, 
                                            med_FF = med_FF), 
                                axis = 1)
        
        # rebin to specified width
        [ new_E_bin_width, rebin_scatter_rate_E ] = rebin_1d(binned_scatter_rate_E,
                                                             self.get_numerics_scatter_binned_rate_E_bin_width(),
                                                             E_bin_width)
        
        E_bin_LHS = [ self.get_material_band_gap() + n*new_E_bin_width for n in range(len(rebin_scatter_rate_E)) ]
              
        return [ E_bin_LHS, norm*rebin_scatter_rate_E ]
    
    def get_cs_reach(self, 
                     med_FF = 2., 
                     n_cut = 3,
                     expt_M_kg = 1,
                     expt_T_year = 1,
                     expt_E_threshold = 0., 
                     i_list = [0]):
        
        band_gap = self.get_material_band_gap()
        E_width = self.get_numerics_scatter_binned_rate_E_bin_width()
        
        threshold_idx = int(max(0., np.floor( 
            ( expt_E_threshold - band_gap ) / E_width )))
        
        cs_reach = []
        
        for m, mass in enumerate(self.get_masses_MeV()):
            
            total_rate = (expt_M_kg/self.get_expt_M())*(expt_T_year/self.get_expt_T())*np.sum(self.get_binned_scatter_rate_qE( 
                          mass_MeV = mass, 
                          med_FF = med_FF,
                          i_list = i_list)[threshold_idx:, :])
            
            cs_reach.append([ mass, n_cut / total_rate ])
            
        return np.transpose(np.array(sorted(cs_reach , key=lambda ele: ele[0])))
            
            
        
        
        
        
        
        
        
        
        
    