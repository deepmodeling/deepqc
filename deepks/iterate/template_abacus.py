import os
import numpy as np
from glob import glob
from deepks.utils import flat_file_list, load_dirs
from deepks.utils import get_sys_name, load_sys_paths
from deepks.task.task import PythonTask
from deepks.task.task import BatchTask, GroupBatchTask, DPDispatcherTask
from deepks.task.workflow import Sequence
from deepks.iterate.template import check_system_names, make_cleanup
from deepks.iterate.generator_abacus import make_abacus_scf_kpt, make_abacus_scf_input, make_abacus_scf_stru

MODEL_FILE = "model.pth"
CMODEL_FILE = "model.ptg"

NAME_TYPE = {   'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
            'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
        'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
        'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
        'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
        'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
        'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
        'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
        'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
        'Ba': 56, #'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
            ## 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
            ## 'Er': 68, 'Tm': 69, 'Yb': 70, 
            ## 'Lu': 71, 
        'Hf': 72, 'Ta': 73,
        'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
        'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 
            ## 'Po': 84, #'At': 85,
            ## 'Rn': 86, #'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
            ## 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
            ## 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
            ## 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
            ## 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Uut': 113,
            ## 'Fl': 114, 'Uup': 115, 'Lv': 116, 'Uus': 117, 'Uuo': 118
        } #dict
TYPE_NAME ={v:k for k, v in NAME_TYPE.items()}
TYPE_INDEX = {k:v for k, v in NAME_TYPE.items()}

ABACUS_CMD="bash run_abacus.sh"

DEFAULT_SCF_ARGS_ABACUS={
    "orb_files": ["orb"],  #atomic number order
    "pp_files": ["upf"],  #atomic number order
    "proj_file": ["orb"], 
    "ntype": 0,
    "nspin": 1,
    "symmetry": 0,
    "nbands": None,
    "ecutwfc": 50,
    "scf_thr": 1e-7,
    "scf_nmax": 50,
    "dft_functional": "pbe", 
    "basis_type": "lcao",
    "gamma_only": 1,
    "k_points": None,
    "kspacing": None,
    "smearing_method":"fixed",
    "smearing_sigma":0.001,
    "mixing_type": "pulay",
    "mixing_beta": 0.4,
    "cal_force": 0,
    "cal_stress": 0,
    "deepks_bandgap": 0,
    "deepks_out_labels":1,
    "deepks_scf":0,
    "lattice_constant": 1,
    "lattice_vector": np.eye(3,dtype=int),
    "coord_type": "Cartesian",
    "run_cmd": "mpirun",
    "sub_size": 1,
    "abacus_path": "/usr/local/bin/ABACUS.mpi",
}

def coord_to_atom(path):
    coords = np.load(f"{path}/coord.npy")
    nframes = coords.shape[0]
    # get type_map.raw and type.raw, use it
    with open(f"{path}/type_map.raw") as fp:
        my_type_map =[TYPE_INDEX[i] for i in fp.read().split()]
    atom_types = np.loadtxt(f"{path}/type.raw", ndmin=1).astype(int)
    atom_types = np.array([int(my_type_map[i-1]) for i in atom_types])\
        .reshape(1,-1).repeat(nframes,axis=0)
    atom_data = np.insert(coords, 0, values=atom_types, axis=2)
    return atom_data


def make_scf_abacus(systems_train, systems_test=None, *,
             train_dump="data_train", test_dump="data_test", cleanup=None, 
             dispatcher=None, resources =None, no_model=True, group_size=1,
             workdir='00.scf', share_folder='share', model_file=None,
             orb_files=[], pp_files=[], proj_file=[],  **scf_abacus):
    #share orb_files and pp_files
    from deepks.iterate.iterate import check_share_folder
    for i in range (len(orb_files)):
        orb_files[i] = check_share_folder(orb_files[i], orb_files[i], share_folder)
    for i in range (len(pp_files)):
        pp_files[i] = check_share_folder(pp_files[i], pp_files[i], share_folder)
        #share the traced model file
    for i in range (len(proj_file)):
        proj_file[i] = check_share_folder(proj_file[i], proj_file[i], share_folder)
   # if(no_model is False):
        #model_file=os.path.abspath(model_file)
        #model_file = check_share_folder(model_file, model_file, share_folder)
    orb_files=[os.path.abspath(s) for s in flat_file_list(orb_files, sort=False)]
    pp_files=[os.path.abspath(s) for s in flat_file_list(pp_files, sort=False)]
    proj_file=[os.path.abspath(s) for s in flat_file_list(proj_file, sort=False)]
    forward_files=orb_files+pp_files+proj_file
    pre_scf_abacus = make_convert_scf_abacus(
        systems_train=systems_train, systems_test=systems_test,
        no_model=no_model, workdir='.', share_folder=share_folder, 
        model_file=model_file, resources=resources,
        dispatcher=dispatcher, orb_files=orb_files, pp_files=pp_files, 
        proj_file=proj_file, **scf_abacus)
    run_scf_abacus = make_run_scf_abacus(systems_train, systems_test,
        no_model=no_model, model_file=model_file, group_data=False,
        workdir='.', outlog="log.scf", share_folder=share_folder, 
        dispatcher=dispatcher, resources=resources, group_size=group_size,
        forward_files=forward_files, 
        **scf_abacus)
    post_scf_abacus = make_stat_scf_abacus(
        systems_train, systems_test,
        train_dump=train_dump, test_dump=test_dump, workdir=".", 
        **scf_abacus)
    # concat
    seq = [pre_scf_abacus, run_scf_abacus, post_scf_abacus]
    #seq = [post_scf_abacus]
    #seq = [pre_scf_abacus]
    if cleanup:
        clean_scf = make_cleanup(
            ["slurm-*.out", "task.*/err", "fin.record"],
            workdir=".")
        seq.append(clean_scf)
    #make sequence
    return Sequence(
        seq,
        workdir=workdir
    )


### need parameters: orb_files, pp_files, proj_file
def convert_data(systems_train, systems_test=None, *, 
                no_model=True, model_file=None, pp_files=[], 
                lattice_vector=np.eye(3, dtype=int), dispatcher=None,**pre_args):
    #trace a model (if necessary)
    if not no_model:
        if model_file is not None:
            from deepks.model import CorrNet
            model = CorrNet.load(model_file)
            model.compile_save(CMODEL_FILE)
            #set 'deepks_scf' to 1, and give abacus the path of traced model file
            pre_args.update(deepks_scf=1, model_file=os.path.abspath(CMODEL_FILE))
        else:
            raise FileNotFoundError(f"No required model file in {os.getcwd()}")
    # split systems into groups
    nsys_trn = len(systems_train)
    nsys_tst = len(systems_test)
    #ntask_trn = int(np.ceil(nsys_trn / sub_size))
    #ntask_tst = int(np.ceil(nsys_tst / sub_size))
    train_sets = [systems_train[i::nsys_trn] for i in range(nsys_trn)]
    test_sets = [systems_test[i::nsys_tst] for i in range(nsys_tst)]
    systems=systems_train+systems_test
    sys_paths = [os.path.abspath(s) for s in load_sys_paths(systems)]
    from pathlib import Path
    if dispatcher=="dpdispatcher" and \
        pre_args["dpdispatcher_machine"]["context_type"].upper().find("LOCAL")==-1:
        #write relative path into INPUT and STRU
        orb_files=pre_args["orb_files"]
        proj_file=pre_args["proj_file"]
        orb_files=["../../../"+str(os.path.basename(s)) for s in orb_files]
        pp_files=["../../../"+str(os.path.basename(s)) for s in pp_files]
        proj_file=["../../../"+str(os.path.basename(s)) for s in proj_file]
        pre_args["orb_files"]=orb_files
        pre_args["proj_file"]=proj_file
        if not no_model:
            pre_args["model_file"]="../../../"+CMODEL_FILE
    #init sys_data (dpdata)
    for i, sset in enumerate(train_sets+test_sets):
        try:
            atom_data = np.load(f"{sys_paths[i]}/atom.npy")
        except FileNotFoundError:
            atom_data = coord_to_atom(sys_paths[i])
        if os.path.isfile(f"{sys_paths[i]}/box.npy"):
            cell_data = np.load(f"{sys_paths[i]}/box.npy")
        nframes = atom_data.shape[0]
        natoms = atom_data.shape[1]
        atoms = atom_data[0,:,0]
        #atoms.sort() # type order
        types = np.unique(atoms) #index in type list
        ntype = types.size
        from collections import Counter
        nta = Counter(atoms) #dict {itype: nta}, natom in each type
        if not os.path.exists(f"{sys_paths[i]}/ABACUS"):
            os.mkdir(f"{sys_paths[i]}/ABACUS")
        #pre_args.update({"lattice_vector":lattice_vector})
        #if "stru_abacus.yaml" exists, update STRU args in pre_args:
        pre_args_new=dict(zip(pre_args.keys(),pre_args.values()))
        if os.path.exists(f"{sys_paths[i]}/stru_abacus.yaml"):
            from deepks.utils import load_yaml
            stru_abacus = load_yaml(f"{sys_paths[i]}/stru_abacus.yaml")
            for k,v in stru_abacus.items():
                pre_args_new[k]=v
        for f in range(nframes):
            if not os.path.exists(f"{sys_paths[i]}/ABACUS/{f}"):
                os.mkdir(f"{sys_paths[i]}/ABACUS/{f}")
            ###create STRU file
            if not os.path.isfile(f"{sys_paths[i]}/ABACUS/{f}/STRU"):
                Path(f"{sys_paths[i]}/ABACUS/{f}/STRU").touch()
            #create sys_data for each frame
            frame_data=atom_data[f]
            #frame_sorted=frame_data[np.lexsort(frame_data[:,::-1].T)] #sort cord by type
            sys_data={'atom_names':[TYPE_NAME[it] for it in nta.keys()], 'atom_numbs': list(nta.values()), 
                        #'cells': np.array([lattice_vector]), 'coords': [frame_sorted[:,1:]]}
                        'cells': np.array([lattice_vector]), 'coords': [frame_data[:,1:]]}
            if os.path.isfile(f"{sys_paths[i]}/box.npy"):
                sys_data={'atom_names':[TYPE_NAME[it] for it in nta.keys()], 'atom_numbs': list(nta.values()),
                        'cells': [cell_data[f]], 'coords': [frame_data[:,1:]]}
            #write STRU file
            with open(f"{sys_paths[i]}/ABACUS/{f}/STRU", "w") as stru_file:
                stru_file.write(make_abacus_scf_stru(sys_data, pp_files, pre_args_new))
            #write INPUT file
            with open(f"{sys_paths[i]}/ABACUS/{f}/INPUT", "w") as input_file:
                input_file.write(make_abacus_scf_input(pre_args))


            #write KPT file if k_points is explicitly specified or for gamma_only case
            if pre_args["k_points"] is not None or pre_args["gamma_only"] is True:
                with open(f"{sys_paths[i]}/ABACUS/{f}/KPT","w") as kpt_file:
                    kpt_file.write(make_abacus_scf_kpt(pre_args_new))


def make_convert_scf_abacus(systems_train, systems_test=None,
                no_model=True, model_file=None, resources=None, **pre_args):
    # if no test systems, use last one in train systems
    systems_train = [os.path.abspath(s) for s in load_sys_paths(systems_train)]
    systems_test = [os.path.abspath(s) for s in load_sys_paths(systems_test)]
    #share model file if needed
    link_prev = pre_args.pop("link_prev_files", [])
    if not systems_test:
        systems_test.append(systems_train[-1])
        # if len(systems_train) > 1:
        #     del systems_train[-1]
    check_system_names(systems_train)
    check_system_names(systems_test)
    #update pre_args
    if not no_model:
        assert model_file is not None
        link_prev.append((model_file, "model.pth"))
    if resources is not None and "task_per_node" in resources:
        task_per_node = resources["task_per_node"]
    pre_args.update(
        systems_train=systems_train, 
        systems_test=systems_test,
        model_file=model_file,
        no_model=no_model, 
        task_per_node = task_per_node, 
        **pre_args)
    return PythonTask(
        convert_data, 
        call_kwargs=pre_args,
        outlog="convert.log",
        errlog="err",
        workdir='.', 
        link_prev_files=link_prev
    )


def make_run_scf_abacus(systems_train, systems_test=None,  
                outlog="out.log",  errlog="err.log", group_size=1,
                resources=None, dispatcher=None, 
                share_folder="share", workdir=".", link_systems=True, 
                dpdispatcher_machine=None, dpdispatcher_resources=None,
                no_model=True, **task_args):
    #basic args
    link_share = task_args.pop("link_share_files", [])
    link_prev = task_args.pop("link_prev_files", [])
    link_abs = task_args.pop("link_abs_files", [])
    forward_files = task_args.pop("forward_files", [])
    backward_files = task_args.pop("backward_files", [])
    if not no_model:
        forward_files.append("../"+CMODEL_FILE) #relative to work_base: system
    #get systems
    systems_train = [os.path.abspath(s) for s in load_sys_paths(systems_train)]
    systems_test = [os.path.abspath(s) for s in load_sys_paths(systems_test)]
    if not systems_test:
        systems_test.append(systems_train[-1])
        # if len(systems_train) > 1:
        #     del systems_train[-1]
    check_system_names(systems_train)
    check_system_names(systems_test)
    #systems=systems_train+systems_test
    sys_train_paths = [os.path.abspath(s) for s in load_sys_paths(systems_train)]
    sys_train_base = [get_sys_name(s) for s in sys_train_paths]
    sys_train_name = [os.path.basename(s) for s in sys_train_base]
    sys_test_paths = [os.path.abspath(s) for s in load_sys_paths(systems_test)]
    sys_test_base = [get_sys_name(s) for s in sys_test_paths]
    sys_test_name = [os.path.basename(s) for s in sys_test_base]
    sys_paths=sys_train_paths + sys_test_paths
    sys_base=sys_train_base+sys_test_base
    sys_name=sys_train_name+sys_test_name
    if link_systems:
        target_dir="systems"
        src_files = sum((glob(f"{base}*") for base in sys_base), [])
        for fl in src_files:
            dst = os.path.join(target_dir, os.path.basename(fl))
            link_abs.append((fl, dst)) 
    #set parameters
    if resources is not None and "task_per_node" in resources:
        task_per_node = resources["task_per_node"]
    run_cmd = task_args.pop("run_cmd", "mpirun")
    abacus_path = task_args.pop("abacus_path", None)
    assert abacus_path is not None
    #make task
    task_list=[]
    if dispatcher=="dpdispatcher":
        if dpdispatcher_resources is not None and "cpu_per_node" in dpdispatcher_resources:
            assert task_per_node <= dpdispatcher_resources["cpu_per_node"]
        #make task_list
        from dpdispatcher import Task
        singletask={
            "command": None, 
            "task_work_path": "./",
            "forward_files":[],
            "backward_files": [], 
            "outlog": outlog,
            "errlog": errlog
        }
        for i, pth in enumerate(sys_paths):
            try:
                atom_data = np.load(f"{str(pth)}/atom.npy")
            except FileNotFoundError:
                atom_data = coord_to_atom(str(pth))
            nframes = atom_data.shape[0]
            for f in range(nframes):
                singletask["command"]=str(f"cd {sys_name[i]}/ABACUS/{f}/ &&  \
                    {run_cmd} -n {task_per_node} {abacus_path} > {outlog} 2>{errlog}  &&  \
                    echo {f}`grep convergence ./OUT.ABACUS/running_scf.log` > conv  &&  \
                    echo {f}`grep convergence ./OUT.ABACUS/running_scf.log`")
                singletask["task_work_path"]="."
                singletask["forward_files"]=[str(f"./{sys_name[i]}/ABACUS/{f}/")]
                singletask["backward_files"]=[str(f"./{sys_name[i]}/ABACUS/{f}/")]
                task_list.append(Task.load_from_dict(singletask))
        return DPDispatcherTask(
            task_list,
            work_base="systems",
            outlog=outlog,
            share_folder=share_folder,
            link_share_files=link_share,
            link_prev_files=link_prev,
            link_abs_files=link_abs,
            machine=dpdispatcher_machine,
            resources=dpdispatcher_resources,
            forward_files=forward_files,
            backward_files=backward_files
        )
    else:
        batch_tasks=[]
        for i, pth in enumerate(sys_paths):
            try:
                atom_data = np.load(f"{str(pth)}/atom.npy")
            except FileNotFoundError:
                atom_data = coord_to_atom(str(pth))
            nframes = atom_data.shape[0]
            for f in range(nframes):
                batch_tasks.append(BatchTask(
                    cmds=str(f"cd {sys_name[i]}/ABACUS/{f}/ &&  \
                    {run_cmd} -n {task_per_node} {abacus_path} > {outlog} 2>{errlog}  &&  \
                    echo {f}`grep convergence ./OUT.ABACUS/running_scf.log` > conv  &&  \
                    echo {f}`grep convergence ./OUT.ABACUS/running_scf.log`"),
                    workdir="systems",
                    forward_files=[str(f"./{sys_name[i]}/ABACUS/{f}/")],
                    backward_files=[str(f"./{sys_name[i]}/ABACUS/{f}/")]
                )) 
        return GroupBatchTask(
            batch_tasks,
            group_size=group_size, 
            workdir="./",
            dispatcher=dispatcher,
            resources=resources,
            outlog=outlog,
            share_folder=share_folder,
            link_share_files=link_share,
            link_prev_files=link_prev,
            link_abs_files=link_abs,
            forward_files=forward_files,
            backward_files=backward_files
        )
    



def gather_stats_abacus(systems_train, systems_test, 
                train_dump, test_dump, cal_force=0, cal_stress=0, deepks_bandgap=0, **stat_args):
    sys_train_paths = [os.path.abspath(s) for s in load_sys_paths(systems_train)]
    sys_test_paths = [os.path.abspath(s) for s in load_sys_paths(systems_test)]
    sys_train_paths = [get_sys_name(s) for s in sys_train_paths]
    sys_test_paths = [get_sys_name(s) for s in sys_test_paths]
    sys_train_names = [os.path.basename(s) for s in sys_train_paths]
    sys_test_names = [os.path.basename(s) for s in sys_test_paths]
    if train_dump is None:
        train_dump = "."
    if test_dump is None:
        test_dump = "."
    #concatenate data (train)
    if not os.path.exists(train_dump):
        os.mkdir(train_dump)
    for i in range(len(systems_train)):
        if not os.path.exists(train_dump + '/' + sys_train_names[i]):
            os.mkdir(train_dump + '/' + sys_train_names[i])
        try:
            atom_data = np.load(f"{sys_train_paths[i]}/atom.npy")
        except FileNotFoundError:
            atom_data = coord_to_atom(sys_train_paths[i])
        nframes = atom_data.shape[0]
        c_list=np.full((nframes,1), False)
        d_list=[]
        e0_list=[]
        f0_list=[]
        s0_list=[]
        o0_list=[]
        e_list=[]
        f_list=[]
        s_list=[]
        o_list=[]
        op_list=[]
        gvx_list=[]
        gvepsl_list=[]
        for f in range(nframes):
            with open(f"{sys_train_paths[i]}/ABACUS/{f}/conv","r") as conv_file:
                ic=conv_file.read().split()
                if "achieved" in ic and "not" not in ic:
                    c_list[(int)(ic[0])]=True
            des = np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_dm_eig.npy")
            d_list.append(des)
            ene = np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_ebase.npy")
            e0_list.append(ene/2)    #Ry to Hartree
            ene = np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_etot.npy")
            e_list.append(ene/2)
            if(cal_force):
                fcs=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_fbase.npy")
                f0_list.append(fcs/2)    #Ry to Hartree
                fcs=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_ftot.npy")
                f_list.append(fcs/2)
                if os.path.exists(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_gradvx.npy"):
                    gvx=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_gradvx.npy")
                    gvx_list.append(gvx)
            if(cal_stress):
                scs=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_sbase.npy")
                s0_list.append(scs/2)    #Ry to Hartree
                scs=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_stot.npy")
                s_list.append(scs/2)
                if os.path.exists(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_gvepsl.npy"):
                    gvepsl=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_gvepsl.npy")
                    gvepsl_list.append(gvepsl)
            if(deepks_bandgap):
                ocs=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_obase.npy")
                o0_list.append(ocs/2)      
                ocs=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_otot.npy")
                o_list.append(ocs/2)
                if os.path.exists(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_orbpre.npy"):
                    orbital_precalc=np.load(f"{sys_train_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_orbpre.npy")
                    op_list.append(orbital_precalc)             
        np.save(f"{train_dump}/{sys_train_names[i]}/conv.npy", c_list)
        dm_eig=np.array(d_list)   #concatenate
        np.save(f"{train_dump}/{sys_train_names[i]}/dm_eig.npy", dm_eig)
        e_base=np.array(e0_list)
        np.save(f"{train_dump}/{sys_train_names[i]}/e_base.npy", e_base)    #Ry to Hartree
        e_ref=np.load(f"{sys_train_paths[i]}/energy.npy")
        np.save(f"{train_dump}/{sys_train_names[i]}/energy.npy", e_ref)
        np.save(f"{train_dump}/{sys_train_names[i]}/atom.npy", atom_data)
        np.save(f"{train_dump}/{sys_train_names[i]}/l_e_delta.npy", e_ref-e_base)
        np.save(f"{train_dump}/{sys_train_names[i]}/e_tot.npy", np.array(e_list))
        if(cal_force): 
            f_base=np.array(f0_list)
            np.save(f"{train_dump}/{sys_train_names[i]}/f_base.npy", f_base)
            f_ref=np.load(f"{sys_train_paths[i]}/force.npy")
            np.save(f"{train_dump}/{sys_train_names[i]}/force.npy", f_ref)
            np.save(f"{train_dump}/{sys_train_names[i]}/l_f_delta.npy", f_ref-f_base)
            np.save(f"{train_dump}/{sys_train_names[i]}/f_tot.npy", np.array(f_list))
            if len(gvx_list) > 0:
                np.save(f"{train_dump}/{sys_train_names[i]}/grad_vx.npy", np.array(gvx_list))
        if(cal_stress): 
            s_base=np.array(s0_list)
            np.save(f"{train_dump}/{sys_train_names[i]}/s_base.npy", s_base)
            s_ref=np.load(f"{sys_train_paths[i]}/stress.npy")
            s_ref=s_ref[:,[0,1,2,4,5,8]] #only train the upper-triangle part
            np.save(f"{train_dump}/{sys_train_names[i]}/stress.npy", s_ref)
            np.save(f"{train_dump}/{sys_train_names[i]}/l_s_delta.npy", s_ref-s_base)
            np.save(f"{train_dump}/{sys_train_names[i]}/s_tot.npy", np.array(s_list))
            if len(gvepsl_list) > 0:
                np.save(f"{train_dump}/{sys_train_names[i]}/grad_vepsl.npy", np.array(gvepsl_list))
        if(deepks_bandgap): 
            o_base=np.array(o0_list)
            np.save(f"{train_dump}/{sys_train_names[i]}/o_base.npy", o_base)
            o_ref=np.load(f"{sys_train_paths[i]}/orbital.npy")
            np.save(f"{train_dump}/{sys_train_names[i]}/orbital.npy", o_ref)
            np.save(f"{train_dump}/{sys_train_names[i]}/l_o_delta.npy", o_ref-o_base)
            np.save(f"{train_dump}/{sys_train_names[i]}/o_tot.npy", np.array(o_list))
            if len(op_list) > 0:
                np.save(f"{train_dump}/{sys_train_names[i]}/orbital_precalc.npy", np.array(op_list))
    #concatenate data (test)
    if not os.path.exists(test_dump):
            os.mkdir(test_dump)
    for i in range(len(systems_test)):
        if not os.path.exists(test_dump + '/' + sys_test_names[i]):
            os.mkdir(test_dump + '/' + sys_test_names[i])
        try:
            atom_data = np.load(f"{sys_test_paths[i]}/atom.npy")
        except FileNotFoundError:
            atom_data = coord_to_atom(sys_test_paths[i])
        nframes = atom_data.shape[0]
        c_list=np.full((nframes,1), False)
        d_list=[]
        e0_list=[]
        f0_list=[]
        s0_list=[]
        o0_list=[]
        e_list=[]
        f_list=[]
        s_list=[]
        o_list=[]
        op_list=[]
        gvx_list=[]
        gvepsl_list=[]
        for f in range(nframes):
            with open(f"{sys_test_paths[i]}/ABACUS/{f}/conv","r") as conv_file:
                ic=conv_file.read().split()
                if "achieved" in ic and "not" not in ic:
                    c_list[(int)(ic[0])]=True
            des = np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_dm_eig.npy")
            d_list.append(des)
            ene = np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_ebase.npy")
            e0_list.append(ene/2)    #Ry to Hartree
            ene = np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_etot.npy")
            e_list.append(ene/2)
            if(cal_force):
                fcs=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_fbase.npy")
                f0_list.append(fcs/2)    #Ry to Hartree
                fcs=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_ftot.npy")
                f_list.append(fcs/2)
                if os.path.exists(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_gradvx.npy"):
                    gvx=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_gradvx.npy")
                    gvx_list.append(gvx)
            if(cal_stress):
                scs=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_sbase.npy")
                s0_list.append(scs/2)    #Ry to Hartree
                scs=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_stot.npy")
                s_list.append(scs/2)
                if os.path.exists(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_gvepsl.npy"):
                    gvepsl=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_gvepsl.npy")
                    gvepsl_list.append(gvepsl)
            if(deepks_bandgap):
                ocs=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_obase.npy")
                o0_list.append(ocs/2)      
                ocs=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_otot.npy")
                o_list.append(ocs/2)
                if os.path.exists(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_orbpre.npy"):
                    orbital_precalc=np.load(f"{sys_test_paths[i]}/ABACUS/{f}/OUT.ABACUS/deepks_orbpre.npy")
                    op_list.append(orbital_precalc)
        dm_eig=np.array(d_list)   #concatenate
        np.save(f"{test_dump}/{sys_test_names[i]}/dm_eig.npy", dm_eig)
        e_base=np.array(e0_list)
        np.save(f"{test_dump}/{sys_test_names[i]}/e_base.npy", e_base)
        e_ref=np.load(f"{sys_test_paths[i]}/energy.npy")
        np.save(f"{test_dump}/{sys_test_names[i]}/energy.npy", e_ref)
        np.save(f"{test_dump}/{sys_test_names[i]}/atom.npy", atom_data)
        np.save(f"{test_dump}/{sys_test_names[i]}/l_e_delta.npy", e_ref-e_base)
        np.save(f"{test_dump}/{sys_test_names[i]}/e_tot.npy", np.array(e_list))
        if(cal_force): 
            f_base=np.array(f0_list)
            np.save(f"{test_dump}/{sys_test_names[i]}/f_base.npy", f_base)
            f_ref=np.load(f"{sys_test_paths[i]}/force.npy")
            np.save(f"{test_dump}/{sys_test_names[i]}/force.npy", f_ref)
            np.save(f"{test_dump}/{sys_test_names[i]}/l_f_delta.npy", f_ref-f_base)
            np.save(f"{test_dump}/{sys_test_names[i]}/f_tot.npy", np.array(f_list))
            if len(gvx_list)>0:
                np.save(f"{test_dump}/{sys_test_names[i]}/grad_vx.npy", np.array(gvx_list))
        if(cal_stress): 
            s_base=np.array(s0_list)
            np.save(f"{test_dump}/{sys_test_names[i]}/s_base.npy", s_base)
            s_ref=np.load(f"{sys_test_paths[i]}/stress.npy")
            s_ref=s_ref[:,[0,1,2,4,5,8]] #only train the upper-triangle part
            np.save(f"{test_dump}/{sys_test_names[i]}/stress.npy", s_ref)
            np.save(f"{test_dump}/{sys_test_names[i]}/l_s_delta.npy", s_ref-s_base)
            np.save(f"{test_dump}/{sys_test_names[i]}/s_tot.npy", np.array(s_list))
            if len(gvepsl_list)>0:
                np.save(f"{test_dump}/{sys_test_names[i]}/grad_vepsl.npy", np.array(gvepsl_list))
        if(deepks_bandgap): 
            o_base=np.array(o0_list)
            np.save(f"{test_dump}/{sys_test_names[i]}/o_base.npy", o_base)
            o_ref=np.load(f"{sys_test_paths[i]}/orbital.npy")
            np.save(f"{test_dump}/{sys_test_names[i]}/orbital.npy", o_ref)
            np.save(f"{test_dump}/{sys_test_names[i]}/l_o_delta.npy", o_ref-o_base)
            np.save(f"{test_dump}/{sys_test_names[i]}/o_tot.npy", np.array(o_list))
            if len(op_list) > 0:
                np.save(f"{test_dump}/{sys_test_names[i]}/orbital_precalc.npy", np.array(op_list))
        np.save(f"{test_dump}/{sys_test_names[i]}/conv.npy",c_list)
    #check convergence and print in log
    from deepks.scf.stats import print_stats
    print_stats(systems=systems_train, test_sys=systems_test,
            dump_dir=train_dump, test_dump=test_dump, group=False, 
            with_conv=True, with_e=True, e_name="e_tot", 
               with_f=True, f_name="f_tot")
    return


def make_stat_scf_abacus(systems_train, systems_test=None, *, 
                  train_dump="data_train", test_dump="data_test", cal_force=0, cal_stress=0, deepks_bandgap=0,
                  workdir='.', outlog="log.data", **stat_args):
    # follow same convention for systems as run_scf
    systems_train = [os.path.abspath(s) for s in load_sys_paths(systems_train)]
    systems_test = [os.path.abspath(s) for s in load_sys_paths(systems_test)]
    if not systems_test:
        systems_test.append(systems_train[-1])
        # if len(systems_train) > 1:
        #     del systems_train[-1]
    # load stats function
    stat_args.update(
        systems_train=systems_train,
        systems_test=systems_test,
        train_dump=train_dump,
        test_dump=test_dump,
        cal_force=cal_force,
        cal_stress=cal_stress,
        deepks_bandgap=deepks_bandgap)
    # make task
    return PythonTask(
        gather_stats_abacus,
        call_kwargs=stat_args,
        outlog=outlog,
        errlog="err",
        workdir=workdir
    )



