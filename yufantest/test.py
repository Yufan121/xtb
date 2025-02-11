nnxtb = '/scratch/pawsey0799/yx7184/githubs/NNxTB/scratchfolder_par_6_separate_ele_mace'
import sys
sys.path.append(nnxtb) 
import subprocess
import argparse
from test_gfn2_derivative import *
from data_processing.DataSets import parse_xyz
# from training_final_val import 



def call_xtb(xtb_exec, xtbpath, xyzfile, target):
    # Call the xTB program with the specified parameters
    cmd = ['sh', 'test_xtb_fortran.sh', xyzfile, xtbpath, '--output']
    cmd.append(target)
    cmd.append(xtb_exec)
    # join to a string and print
    print(' '.join(cmd))
    

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"xTB execution failed: {result.stderr}")
    
    return result

def obtain_vanilla_xtb(xtb_vanilla=None, xtbpath=None, xyzfile=None):
    # Call vanilla xTB
    # Parse energies and forces from the output
    result = call_xtb(xtb_vanilla, xtbpath, xyzfile, '--fe')
    xyz_name = os.path.basename(xyzfile)
    xyz_name = os.path.splitext(xyz_name)[0]
    print(f'xyz_name: {xyz_name}')
    energies, forces = parse_fe(f'{xyz_name}.engrad')
    
    result = call_xtb(xtb_vanilla, xtbpath, xyzfile, '--vib')
    frequencies = parse_vib_freq('vibspectrum')

    # remove f'{xyz_name}.engrad' and vibspectrum
    os.remove(f'{xyz_name}.engrad')
    os.remove('vibspectrum')
    
    return energies, forces, frequencies

def obtain_modified_xtb(xtb_modified=None, xtbpath=None, xyzfile=None, param_file=None, aid=None, key=None, delta=None, zero_params=False):
    # Modify the parameter file
    param_file.prep_per_atom_param()

    if zero_params:
        # Set all parameters to zero
        param_file.set_param([{key: 0 for key in param_file.ele_param_enum}])
    else:
        # Perturb parameter
        param_file.perturb_param(aid=aid, index=0, key=key, offset=-delta)  # Subtract value from element parameter
        param_file.perturb_param(aid=aid, index=0, key=key, offset=delta)  # Add value to atom parameter

    # Write the modified parameter file
    param_file.write_parameters()

    # Call modified xTB
    result = call_xtb(xtb_modified, xtbpath, xyzfile, '--fe')
    xyz_name = os.path.basename(xyzfile)
    xyz_name = os.path.splitext(xyz_name)[0]
    print(f'xyz_name: {xyz_name}')
    energies, forces = parse_fe(f'{xyz_name}.engrad')
    
    result = call_xtb(xtb_modified, xtbpath, xyzfile, '--vib')
    frequencies = parse_vib_freq('vibspectrum')

    # remove f'{xyz_name}.engrad' and vibspectrum
    os.remove(f'{xyz_name}.engrad')
    os.remove('vibspectrum')

    return energies, forces, frequencies

def compare_results(vanilla, modified, name):
    if not np.allclose(vanilla, modified, atol=0):
        raise ValueError(f"{name} mismatch: Vanilla vs Modified")
    print(f"{name} match.")
    
    
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Process xTB calculations.')
    parser.add_argument('--xtb_modified', type=str, default='xtb', help='Path to the modified xTB executable.')
    parser.add_argument('--xtb_vanilla', type=str, default='/scratch/pawsey0799/yx7184/xtb_install_vanilla/bin/xtb', help='Path to the vanilla xTB executable.')
    parser.add_argument('--xtbpath', type=str, default='/path/to/xtb/share/xtb', help='Path to the xTB share directory.')
    parser.add_argument('--xyzfile', type=str, default='molecule.xyz', help='Path to the input XYZ file.')    
    
    args = parser.parse_args()

    xtb_modified = args.xtb_modified
    xtb_vanilla = args.xtb_vanilla
    xtbpath = args.xtbpath
    xyzfile = args.xyzfile
    
    glob_param_enum = ['ks', 'kp', 'kd', 'ksd', 'kpd', 'kdiff', 'enscale', 'ipeashift', 'gam3s', 'gam3p', 'gam3d1', 'gam3d2', 'aesshift', 'aesexp', 'aesrmax', 'alphaj', 'a1', 'a2', 's8', 's9', 'aesdmp3', 'aesdmp5', 'kexp', 'kexplight']
    ele_param_enum = ['lev', 'exp', 'EN','GAM', 'GAM3', 'KCNS', 'KCNP', 'KCND', 'DPOL', 'QPOL', 'REPA', 'REPB', 'POLYS', 'POLYP', 'POLYD', 'LPARP', 'LPARD', 'mpvcn', 'mprad']
    len_ele_param_enum = [3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

    
    elem_param = f'{xtbpath}/param_gfn2-xtb.txt'
    atom_param = f'{xtbpath}/param_gfn2-xtb-per-atom.txt'
    ref_dict = parse_xyz(xyzfile) 
    print(f'ref_dict: {ref_dict}')
    
    at_list = ref_dict['atomic_numbers']
    print(f'at_list: {at_list}')
    
    param_file = ParameterFile(elem_param, atom_param, ele_param_enum, glob_param_enum, at_list) # Update with correct parameters

    # Vanilla xTB
    energies_vanilla, forces_vanilla, frequencies_vanilla = obtain_vanilla_xtb(xtb_vanilla, xtbpath, xyzfile)

    # Modified xTB with all zero parameters
    energies_zero, forces_zero, frequencies_zero = obtain_modified_xtb(xtb_modified, xtbpath, xyzfile, param_file, zero_params=True)

    # Modified xTB with perturbed parameters
    # energies_perturbed, forces_perturbed, frequencies_perturbed = obtain_modified_xtb(xtb_modified, xtbpath, xyzfile, param_file, aid=0, key='EN', delta=1, zero_params=False)


    # Compare the results
    try:
        compare_results(energies_vanilla, energies_zero, "Energies")
        compare_results(forces_vanilla, forces_zero, "Forces")
        compare_results(frequencies_vanilla, frequencies_zero, "Frequencies")
        print("All comparisons passed.")
    except ValueError as e:
        print(f"Error: {e}")


    # # Compare the results
    # print("Vanilla xTB Energies:", energies_vanilla)
    # print("Zero Params xTB Energies:", energies_zero)
    # # print("Perturbed Params xTB Energies:", energies_perturbed)

    # print("Vanilla xTB Forces:", forces_vanilla)
    # print("Zero Params xTB Forces:", forces_zero)
    # # print("Perturbed Params xTB Forces:", forces_perturbed)

    # print("Vanilla xTB Frequencies:", frequencies_vanilla)
    # print("Zero Params xTB Frequencies:", frequencies_zero)
    # # print("Perturbed Params xTB Frequencies:", frequencies_perturbed)
