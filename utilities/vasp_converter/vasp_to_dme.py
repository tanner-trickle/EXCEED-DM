# Python script to convert VASP output {vasprun.xml, WAVECAR} to 
# a compact version to be read by dm-electron.

# If all-electron reconstruction is requested (default) POTCAR
# is also needed.

import click

def get_data(vasprun_file):
    import numpy as np
    from pymatgen.io.vasp import Vasprun
    from pymatgen.electronic_structure.core import Spin

    vr = Vasprun(vasprun_file)
    
    structure = vr.final_structure
    a_vecs = structure.lattice.matrix
    b_vecs = structure.lattice.reciprocal_lattice.matrix
    pc_vol_A = structure.lattice.volume
    k_red = np.array(vr.actual_kpoints)
    k_weight = np.array(vr.actual_kpoints_weights)*2
    n_k = len(k_red)

    energy = vr.eigenvalues[Spin.up][:,:,0] # Eigenenergies taken from first spin channel only
    n_elec = int(vr.parameters['NELECT'])
    n_val = int(n_elec/2)
    n_k = len(vr.actual_kpoints)
    n_bands = energy.shape[1]
    n_cond = n_bands - n_val
    
    return {
        'n_elec': n_elec,
        'n_k': n_k,
        'n_bands': n_bands,
        'n_val': n_val,
        'n_cond': n_cond,
        'energy': energy,
        'a_vecs_A': a_vecs,
        'b_vecs_A': b_vecs,
        'k_red': k_red,
        'k_weight': k_weight,
        'pc_vol_A': pc_vol_A
    }


def get_wavefunction_data(wavecar_filename):
    import numpy as np
    from tqdm.auto import tqdm
    from pymatgen.io.vasp.outputs import Wavecar

    wc = Wavecar(wavecar_filename)
    nkpts = wc.nk
    nbands = wc.nb

    all_gpoints = np.concatenate(wc.Gpoints).astype(int)
    unique_gpoints = np.unique(all_gpoints, axis=0)
    ngpoints = len(unique_gpoints)
    
    # create a dictionary mapping of gpoint: index
    # the gpoints are converted to tuples of ints, 
    # as lists cannot be used as dictionary keys
    g_index_map = dict(zip(map(tuple, unique_gpoints), np.arange(ngpoints)))
    
    # create empty wavefunction data array
    wfc_data = np.zeros((nbands, nkpts, ngpoints), dtype=complex)
    
    for kpt in tqdm(range(nkpts), total=nkpts):
        gpoints = wc.Gpoints[kpt].astype(int)
        
        # put coeffs with shape (ngpoints, nbands)
        coeffs = np.column_stack(wc.coeffs[kpt])
        
        # get the position of each gpoint in the unique_gpoints array
        # using the index mapping
        order = np.array([g_index_map[g] for g in map(tuple, gpoints)])
        
        wfc_data[:, kpt, order] = coeffs.T
    
    return unique_gpoints, ngpoints, wfc_data


def get_wavefunction_data_pawpy(
    wavecar_file, potcar_file, vasprun_file, cutoff=1000,
):
    import numpy as np
    from tqdm.auto import tqdm
    from pymatgen.io.vasp import Vasprun
    from pymatgen.io.vasp.inputs import Potcar
    from pymatgen.io.vasp.outputs import Wavecar
    from pawpyseed.core.wavefunction import Wavefunction, CoreRegion
    from pawpyseed.core.momentum import MomentumMatrix
    from pawpyseed.core import pawpyc

    # Set up Wavefunction and MomentumMatrix objects
    vr = Vasprun(vasprun_file)
    structure = vr.final_structure
    dim = np.array([vr.parameters["NGX"], vr.parameters["NGY"], vr.parameters["NGZ"]])
    symprec = vr.parameters["SYMPREC"]
    potcar = Potcar.from_file(potcar_file)
    pwf = pawpyc.PWFPointer(wavecar_file, vr)
    wf = Wavefunction(structure, pwf, CoreRegion(potcar), dim, symprec, True)
    mm = MomentumMatrix(wf, cutoff)
    
    # g-point grid and number of k-points and bands
    gpoints = mm.momentum_grid
    nkpts = wf.nwk
    nbands = wf.nband
    ngpoints = len(gpoints)
 
    wfc_data = np.zeros((nbands, nkpts, ngpoints), dtype=complex)

    for kpt_idx in tqdm(range(nkpts), total=nkpts):
        for band_idx in range(nbands):
            wfc_data[band_idx, kpt_idx, :] = mm.get_reciprocal_fullfw(band_idx, kpt_idx, 0)
    
    return gpoints, ngpoints, wfc_data

@click.command(
    context_settings=dict(help_option_names=["-h", "--help"])
)
@click.argument('mat_name')
@click.argument('band_gap', type=float)
@click.option('-v', '--vasprun', 'vasprun_filename', default='vasprun.xml', help='vasprun.xml file')
@click.option('-w', '--wavecar', 'wavecar_filename', default='WAVECAR', help='WAVECAR file')
@click.option('-p', '--potcar', 'potcar_filename', default='POTCAR', help='POTCAR file')
@click.option('-o', '--output', 'out_filename', default='doctor_in.h5', help='Output filename')
@click.option('--no-pawpyseed', 'use_pawpy', is_flag=True, default=True, help="Don't use pawpyseed")
@click.option('-c', '--cutoff', default=1000, type=float, help='pawpyseed cutoff [eV]', show_default=True)
def create_input(
    mat_name, 
    band_gap, 
    vasprun_filename,
    wavecar_filename,
    potcar_filename,  
    out_filename,
    use_pawpy,
    cutoff
):
    """Create input for DOCTOR using VASP outputs"""
    import h5py
    import warnings
    
    warnings.filterwarnings("ignore",  module="pymatgen")

    data = get_data(vasprun_filename)
    
    if use_pawpy:
        unique_gpoints, ngpoints, wfc_data = get_wavefunction_data_pawpy(
            wavecar_filename, potcar_filename, vasprun_filename, cutoff=cutoff
        )
    else:
        warnings.warn(
            "WARNING: Projector augmented wavefunctions have been disabled.\n"
            "All electron wavefunction reconstruction has not been done and cutoff will be ignored."
        )
        unique_gpoints, ngpoints, wfc_data = get_wavefunction_data(wavecar_filename)
    
    with h5py.File(out_filename, "w") as f:
    
        f.create_dataset('mat_name', data=mat_name)
        f.create_dataset('band_gap', data=band_gap)
        f.create_dataset('n_k', data=data['n_k'])
        f.create_dataset('n_val', data=data['n_val'])
        f.create_dataset('n_cond', data=data['n_cond'])
        f.create_dataset('pc_vol_A', data=data['pc_vol_A'])
        f.create_dataset('k_weight', data=data['k_weight'])
        f.create_dataset('k_red', data=data['k_red'].T)
        f.create_dataset('a_vecs_A', data=data['a_vecs_A'].T)
        f.create_dataset('b_vecs_A', data=data['b_vecs_A'].T)
        f.create_dataset('energy_bands_raw', data=data['energy'].T)
        
        f.create_dataset('n_in_G', data=ngpoints)
        f.create_dataset('in_G_grid_red', data=unique_gpoints.T)
        
        for band_idx, band_data in enumerate(wfc_data):   
            f.create_dataset('in_wfc_FT_r/{}'.format(band_idx + 1), data=band_data.real.T)
            f.create_dataset('in_wfc_FT_c/{}'.format(band_idx + 1), data=band_data.imag.T)



if __name__ == "__main__":
    create_input()
