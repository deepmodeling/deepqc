#These 3 functions are copied from dpgen to generate ABACUS INPUT , KPT and STRU file.

bohr2ang = 0.52917721067

def make_abacus_scf_kpt(fp_params):
    # Make KPT file for abacus pw scf calculation.
    # KPT file is the file containing k points infomation in ABACUS scf calculation.
    k_points = [1, 1, 1, 0, 0, 0]
    if "k_points" in fp_params:
        k_points = fp_params["k_points"]
        if len(k_points) != 6:
            raise RuntimeError("k_points has to be a list containig 6 integers specifying MP k points generation.")
    ret = "K_POINTS\n0\nGamma\n"
    for i in range(6):
        ret += str(k_points[i]) + " "
    return ret

def make_abacus_scf_input(fp_params):
    # Make INPUT file for abacus pw scf calculation.
    ret = "INPUT_PARAMETERS\n"
    ret += "calculation scf\n"
    # assert(fp_params['ntype'] >= 0 and type(fp_params["ntype"]) == int),  "'ntype' should be a positive integer."
    # ret += "ntype %d\n" % fp_params['ntype']
    #ret += "pseudo_dir ./\n"
    if "ecutwfc" in fp_params:
        assert(fp_params["ecutwfc"] >= 0) ,  "'ntype' should be non-negative."
        ret += "ecutwfc %f\n" % fp_params["ecutwfc"]
    if "scf_thr" in fp_params:
        ret += "scf_thr %e\n" % fp_params["scf_thr"]
    if "scf_nmax" in fp_params:
        assert(fp_params['scf_nmax'] >= 0 and type(fp_params["scf_nmax"])== int), "'scf_nmax' should be a positive integer."
        ret += "scf_nmax %d\n" % fp_params["scf_nmax"]    
    if "basis_type" in fp_params:
        assert(fp_params["basis_type"] in ["pw", "lcao", "lcao_in_pw"]) , "'basis_type' must in 'pw', 'lcao' or 'lcao_in_pw'."
        ret+= "basis_type %s\n" % fp_params["basis_type"]
    if "dft_functional" in fp_params:
        ret += "dft_functional %s\n" % fp_params["dft_functional"]
    if "gamma_only" in fp_params:
        assert(fp_params["gamma_only"] ==0 or fp_params["gamma_only"] ==1 ) , "'gamma_only' should be 0 or 1."
        ret+= "gamma_only %d\n" % fp_params["gamma_only"]  
    if "mixing_type" in fp_params:
        assert(fp_params["mixing_type"] in ["plain", "kerker", "pulay", "pulay-kerker", "broyden"])
        ret += "mixing_type %s\n" % fp_params["mixing_type"]
    if "mixing_beta" in fp_params:
        assert(fp_params["mixing_beta"] >= 0 and fp_params["mixing_beta"] < 1), "'mixing_beta' should between 0 and 1."
        ret += "mixing_beta %f\n" % fp_params["mixing_beta"]
    if "symmetry" in fp_params:
        #assert(fp_params["symmetry"] == 0 or fp_params["symmetry"] == 1), "'symmetry' should be either 0 or 1."
        ret += "symmetry %d\n" % fp_params["symmetry"]
    if "nbands" in fp_params:
        if(type(fp_params["nbands"]) == int and fp_params["nbands"] > 0):
            ret += "nbands %d\n" % fp_params["nbands"]
        else:
            print("warnning: Parameter [nbands] given is not a positive integer, the default value of [nbands] in ABACUS will be used. ")
    if "nspin" in fp_params:
        assert(fp_params["nspin"] == 1 or fp_params["nspin"] == 2 or fp_params["nspin"] == 4), "'nspin' can anly take 1, 2 or 4"
        ret += "nspin %d\n" % fp_params["nspin"]
    if "ks_solver" in fp_params:
        assert(fp_params["ks_solver"] in ["cg", "dav", "lapack", "genelpa", "hpseps", "scalapack_gvx"]), "'ks_sover' should in 'cgx', 'dav', 'lapack', 'genelpa', 'hpseps', 'scalapack_gvx'."
        ret += "ks_solver %s\n" % fp_params["ks_solver"]
    if "smearing_method" in fp_params:
        assert(fp_params["smearing_method"] in ["gaussian", "fd", "fixed", "mp", "mp2", "mv"]), "'smearing' should in 'gaussian', 'fd', 'fixed', 'mp', 'mp2', 'mv'. "
        ret += "smearing_method %s\n" % fp_params["smearing_method"]
    if "smearing_sigma" in fp_params:
        assert(fp_params["smearing_sigma"] >= 0), "'smearing_sigma' should be non-negative."
        ret += "smearing_sigma %f\n" % fp_params["smearing_sigma"]
    if (("kspacing" in fp_params) and (fp_params["k_points"] is None) and (fp_params["gamma_only"] == 0)):
        assert(fp_params["kspacing"] > 0), "'kspacing' should be positive."
        ret += "kspacing %f\n" % fp_params["kspacing"]
    if "cal_force" in fp_params:
        assert(fp_params["cal_force"] == 0  or fp_params["cal_force"] == 1), "'cal_force' should be either 0 or 1."
        ret += "cal_force %d\n" % fp_params["cal_force"]
    if "cal_stress" in fp_params:
        assert(fp_params["cal_stress"] == 0  or fp_params["cal_stress"] == 1), "'cal_stress' should be either 0 or 1."
        ret += "cal_stress %d\n" % fp_params["cal_stress"]    
    #paras for deepks
    if "deepks_out_labels" in fp_params:
        assert(fp_params["deepks_out_labels"] == 0 or fp_params["deepks_out_labels"] == 1), "'deepks_out_labels' should be either 0 or 1."
        ret += "deepks_out_labels %d\n" % fp_params["deepks_out_labels"]
    if "deepks_scf" in fp_params:
        assert(fp_params["deepks_scf"] == 0  or fp_params["deepks_scf"] == 1), "'deepks_scf' should be either 0 or 1."
        ret += "deepks_scf %d\n" % fp_params["deepks_scf"]
    if "deepks_bandgap" in fp_params:
        assert(fp_params["deepks_bandgap"] == 0  or fp_params["deepks_bandgap"] == 1), "'deepks_bandgap' should be either 0 or 1."
        ret += "deepks_bandgap %d\n" % fp_params["deepks_bandgap"]
    if "deepks_v_delta" in fp_params:
        assert(fp_params["deepks_v_delta"] == 0  or fp_params["deepks_v_delta"] == 1 or fp_params["deepks_v_delta"] == 2), "'deepks_v_delta' should be either 0/1/2."
        ret += "deepks_v_delta %d\n" % fp_params["deepks_v_delta"]
    if "model_file" in fp_params:
        ret += "deepks_model %s\n" % fp_params["model_file"]
    if "out_wfc_lcao" in fp_params:
        ret += "out_wfc_lcao %s\n" % fp_params["out_wfc_lcao"]
    if fp_params["dft_functional"] == "hse":
        ret += "exx_pca_threshold 1e-4\n"
        ret += "exx_c_threshold 1e-4\n"
        ret += "exx_dm_threshold 1e-4\n"
        ret += "exx_schwarz_threshold 1e-5\n"
        ret += "exx_cauchy_threshold 1e-7\n"
        ret += "exx_ccp_rmesh_times 1\n"
    return ret

def make_abacus_scf_stru(sys_data, fp_pp_files, fp_params):
    atom_names = sys_data['atom_names']  # Get the list of atom names, e.g., ['Cs', 'Pb', 'I']
    atom_numbs = sys_data['atom_numbs']  # Get the number of each atom type, e.g., [4, 4, 12]
    
    # Create a dictionary to store the pseudopotential and orbital files for each atom in the system
    valid_pp_files = {}
    valid_orb_files = {}
    
    # Iterate over the provided pp_files, validate and select the pseudopotential files that match atoms in the current system
    for pp_file in fp_pp_files:
        # Manually get the file name, assuming the path is separated by '/'
        filename = pp_file.split('/')[-1]  # Extracts the file name from '../../../Pb_ONCV_PBE-1.0.upf' -> 'Pb_ONCV_PBE-1.0.upf'
        element_name = filename.split('_')[0]  # Extract the element name, e.g., 'Pb_ONCV_PBE-1.0.upf' -> 'Pb'
        
        # If the element exists in the current system's atom_names, check if a pseudopotential file has already been assigned
        if element_name in atom_names:
            # If a pseudopotential file for this element already exists, raise an error
            assert element_name not in valid_pp_files, f"Pseudopotential file for element {element_name} already exists, cannot assign multiple files for the same element."
            # Otherwise, add the pseudopotential file to the valid_pp_files dictionary
            valid_pp_files[element_name] = pp_file

    # Ensure that all elements have corresponding pseudopotential files
    for atom in atom_names:
        assert atom in valid_pp_files, f"No pseudopotential file found for element {atom}, please check the provided pp_files."
    
    # Continue building the STRU file
    ret = "ATOMIC_SPECIES\n"
    for iatom in range(len(atom_names)):
        atom = atom_names[iatom]
        # Now select the correct pseudopotential file from valid_pp_files
        ret += f"{atom} 1.00 {valid_pp_files[atom]}\n"  # Write atom name, mass (1.00), and corresponding pseudopotential file
    
    # Continue generating other parts, such as LATTICE_CONSTANT, LATTICE_VECTORS, etc.
    if "lattice_constant" in fp_params:
        ret += "\nLATTICE_CONSTANT\n"
        ret += f"{fp_params['lattice_constant']}\n\n"  # in Bohr
    else:
        ret += "\nLATTICE_CONSTANT\n"
        ret += f"{1 / bohr2ang}\n\n"  # Default value is 1/bohr2ang
    
    ret += "LATTICE_VECTORS\n"
    cell = sys_data["cells"][0].reshape([3, 3])
    for ix in range(3):
        for iy in range(3):
            ret += f"{cell[ix][iy]} "
        ret += "\n"
    
    # Continue processing atomic coordinates
    ret += "\nATOMIC_POSITIONS\n"
    ret += f"{fp_params['coord_type']}\n\n"
    
    natom_tot = 0
    coord = sys_data['coords'][0]
    for iele in range(len(atom_names)):
        ret += f"{atom_names[iele]}\n"
        ret += "0.0\n"  # Reference energy
        ret += f"{atom_numbs[iele]}\n"
        for iatom in range(atom_numbs[iele]):
            ret += f"{coord[natom_tot, 0]:.12f} {coord[natom_tot, 1]:.12f} {coord[natom_tot, 2]:.12f} 0 0 0\n"
            natom_tot += 1
    assert natom_tot == sum(atom_numbs), "The total number of atoms does not match."
    
    # If the basis type is localized orbitals (lcao)
    if "basis_type" in fp_params and fp_params["basis_type"] == "lcao":
        ret += "\nNUMERICAL_ORBITAL\n"
        # Iterate over the provided orb_files, validate and select the orbital files that match atoms in the current system
        for orb_file in fp_params["orb_files"]:
            # Manually get the file name, assuming the path is separated by '/'
            filename = orb_file.split('/')[-1]  # Extracts the file name
            element_name = filename.split('_')[0]  # Extract the element name, e.g., 'Pb_gga_7au_100Ry_2s2p1d.orb' -> 'Pb'
            
            # If the element exists in the current system's atom_names, check if an orbital file has already been assigned
            if element_name in atom_names:
                # If an orbital file for this element already exists, raise an error
                assert element_name not in valid_orb_files, f"Orbital file for element {element_name} already exists, cannot assign multiple files for the same element."
                # Otherwise, add the orbital file to the valid_orb_files dictionary
                valid_orb_files[element_name] = orb_file
        
        # Ensure that all elements have corresponding orbital files
        for atom in atom_names:
            assert atom in valid_orb_files, f"No orbital file found for element {atom}, please check the provided orb_files."
        
        # Write the valid orbital files into ret
        for atom in atom_names:
            ret += f"{valid_orb_files[atom]}\n"

    # If DeepKS calculation is enabled
    if "deepks_scf" in fp_params and fp_params["deepks_out_labels"] == 1:
        ret += "\nNUMERICAL_DESCRIPTOR\n"
        ret += f"{fp_params['proj_file'][0]}\n"
    
    return ret