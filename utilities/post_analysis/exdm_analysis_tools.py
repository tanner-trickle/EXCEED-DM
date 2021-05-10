# Class to aid in the post-processing of the output of EXCEED-DM

import h5py
import numpy as np

cmet_to_inv_eV = 5.0677*10.**4
kg_year = 2.6876*10**58

class DMEOutput:
    
    def __init__(self, filename):
            
        # read the hdf5 file
        self.data = h5py.File(filename, 'r')
        
    def get_E_bin_width(self):
        """
            Returns the bin width in energy.
        """
        
        return self.data['numerics']['E_bin_width'][...][0]
    
    def get_band_gap(self):
        """
            Returns the band gap.
        """
        
        return self.data['material']['band_gap'][...]
    
    def get_core_init_states(self, n_list, l_list):
        """
            Given a core type, returns the list of init_states
        """

        core_elec_config = self.data['core']['core_elec_config'][...].T

        init_states = []

        for s, state in enumerate(core_elec_config):

            if state[1] in n_list and state[2] in l_list:

                init_states.append(s + 1)

        return init_states
        
#     def get_bin_omega_min(self):
#         """
        
#             Each bin in omega space extends from [ omega_min, omega_max ]. This function
#             returns the list of omega_mins that were comptued for.
        
#         """
        
#         band_gap = self.data['material']['band_gap'][...]
        
#         E_bin_width = self.data['numerics']['E_bin_width'][...]
        
#         n_E_bins = int(self.data['numerics']['n_E_bins'][...])
        
#         omega_min_list = np.zeros(n_E_bins)
        
#         for n in range(n_E_bins):
            
#             omega_min_list[n] = band_gap + n*E_bin_width
            
#         return omega_min_list
        
    def get_masses(self):
        """
            
            Returns the list of masses calculations were run for.
            
            Units : eV
        
        """
        
        mass_list = self.data['particle_physics']['mX'][...]
        
        return mass_list
    
    def get_FDMs(self):
        """
            
            Returns the list of DM mediator form factors calculations were run for.
            
            Units : Nnoe
        
        """
        
        return self.data['particle_physics']['FDM_list'][...]
    
    def get_times(self):
        """
            
            Returns the list of times calculations were run for.
            
            Units : Nnoe
        
        """
        
        return self.data['particle_physics']['t_list'][...]
    
    
    def get_mass_key(self, mass):
        """
            
            Finds the key for the mass point 
            
        """
        
        mass_list = self.get_masses()
        
        ind = (np.abs(mass_list - mass)).argmin()
    
        mass_diff = np.abs(mass_list[ind] - mass)

        if mass_diff > 1:

            print('-- WARNING --')
            print('  Searching data for mass = '+str(mass)
                  +' and the closest value = '+str(mass_list[ind])
                  +' is not that close.')

        return 'm_'+str(ind + 1)
    
    def get_fdm_key(self, f_pow):
        """
            
            Finds the key of the mediator type
            
            e.g.
                Light : f_pow = 2
                Heavy : f_pow = 0
            
        """
        
        f_list = self.get_FDMs()
        
        ind = (np.abs(f_list - f_pow)).argmin()
    
        f_diff = np.abs(f_list[ind] - f_pow)
        if f_diff > 10**(-2):

            print('-- WARNING --')
            print('  Searching data for FDM = '+str(f_pow)
                  +' and the closest value = '+str(f_list[ind])
                  +' is not that close.')

        return 'f_'+str(ind + 1)
    
    def get_time_key(self, time):
    
        t_list = self.get_times()
        ind = (np.abs(t_list - time)).argmin()

        t_diff = np.abs(t_list[ind] - time)

        if t_diff > 10**(-2):

            print('-- WARNING --')
            print('  Searching data for t = '+str(time)
                  +' and the closest value = '+str(t_list[ind])
                  +' is not that close.')

        return 't_'+str(ind + 1)
    
    def get_2d_binned_rate_i(self, init, mass, f_pow, time = 0, 
                            ref_cs_cmet2 = 10**(-40)):
        """
            
            Returns the 2d, binned transition rate from state i -> { all final } for a given
            mass, f_pow, and time.
            
            Units : (kg yr)^(-1)
            
            Dim : [ n_E_bins = 1, n_q_bins + 1 ]
            
            Note : rate calculations are run taking sigma_e = 1, and therefore to get
            the rate per kg-year one must multipy by a reference cross section with units
            of eV^(-2)
            
        """
        # avoid 0 entries
        eps = 10.**(-100)
        
        ref_cs = ref_cs_cmet2*cmet_to_inv_eV**2
        sigma_exp = ref_cs*kg_year
        
        init_key = 'init_'+str(int(init))
        
        mass_key = self.get_mass_key(mass)
        f_key = self.get_fdm_key(f_pow)
        t_key = self.get_time_key(time)

        dR = self.data['rates'][t_key][f_key][mass_key][init_key]['binned_i'][...]
        
        return sigma_exp*dR + eps
    
    def get_binned_rate_E(self, mass, f_pow, time = 0, 
                          init_states = [], 
                          ref_cs_cmet2 = 10**(-40), 
                          bin_width = 0.):
        """
            
            Returns the transition rate, binned in E from states init_states -> { all final } for a given
            mass, f_pow, and time.
            
            Units : (kg yr)^(-1)
            
            Dim : [ 2, n_E_bins + 1]
            
            First element is the left bin position in eV, second element is the transition rate.
            
            Note : rate calculations are run taking sigma_e = 1, and therefore to get
            the rate per kg-year one must multipy by a reference cross section with units
            of eV^(-2)
            
        """
        
        if bin_width == 0:
            bin_width = self.get_E_bin_width()
        
        # avoid 0 entries
        eps = 10.**(-100)
        
        # this already includes the + 1
        n_E_bins = int(self.data['numerics']['n_E_bins'][...])
        
        if init_states == []:
            # just get tht total

            ref_cs = ref_cs_cmet2*cmet_to_inv_eV**2
            sigma_exp = ref_cs*kg_year

            mass_key = self.get_mass_key(mass)
            f_key = self.get_fdm_key(f_pow)
            t_key = self.get_time_key(time)

            dR_E = sigma_exp*np.sum(self.data['rates'][t_key][f_key][mass_key]['total_binned'][...], axis = 1) + eps

        else:
            # sum inidividual elements

            dR_E = np.zeros(n_E_bins)

            for init in init_states:

                dR_E += np.sum(self.get_2d_binned_rate_i(init, mass, f_pow, time, 
                                                        ref_cs_cmet2 = ref_cs_cmet2), axis = 1)
                
        dR_E_out = simple_rebin(dR_E, self.get_E_bin_width(), bin_width) + eps
        bin_mins = bin_min(bin_width, len(dR_E_out), self.get_band_gap())
            
        output = np.zeros((2, len(dR_E_out)))
        output[0, :] = bin_mins
        output[1, :] = dR_E_out
        
        return output

            
    def get_dR_E_list(self, f_pow, time = 0, init_states = []):
        """

            Returns the binned transition rate from init_states -> { all final } for each mass

            Dim : [n_masses, n_E_bins + 1]

        """

        dR_E_list = []
        
        mass_list = self.get_masses()
        
        mass_list.sort()

        for mass in mass_list:

            dR_E_list.append(self.get_binned_rate_E(mass, f_pow, 
                                                    time = time,
                                                   init_states = init_states, 
                                                   ref_cs_cmet2 = 1))

        return dR_E_list 

    
def cross_section(dR_E_list, threshold = 0, bin_width = 1, n_cut = 3, exp = 1):
    
    t_bin = int(threshold/bin_width)
    rate = [ np.sum(dR_E[t_bin:]) for dR_E in dR_E_list ]
    
    cs = []
    
    for r in rate:
        
        if r > 0:
            cs.append(n_cut/(r*exp))
        else:
            cs.append(10**(-10))
        
    return cs

def simple_rebin(li, old_bin_width, new_bin_width):
    """
        Takes in a list, li, binned with old_width and returns a rebinned list
    """
    if new_bin_width < old_bin_width:
        print('-- WARNING --')
        print('  Cannot rebin to a smaller size. Returning original with bin width = '
              +str(old_bin_width)
             )
        return li
    
    if old_bin_width == new_bin_width:
        return li
    
    # make sure that the new_bin_width is an integer multiple of the old bin width
    div_bin_width = int(new_bin_width/old_bin_width)
    
    rebinned_li = [ li[i:i+div_bin_width] for i in range(0, len(li), div_bin_width) ]
    
    for i, ele in enumerate(rebinned_li):
        if len(ele) != div_bin_width:
            rebinned_li[i] = np.append(ele, np.zeros(div_bin_width - len(ele)))
                
    return np.sum(rebinned_li, axis = 1)
            
    
def bin_min(width, n_ele, offset = 0):
    """
        Minimum edge of a bin given the number of bins, their width, and an offset. (Usually taken to be the band gap here.)
    """
    
    return width*np.arange(n_ele) + offset
        
        
        
        
        
        
        
        
        
        