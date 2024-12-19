import os
import sys
import numpy as np
from glob import glob
from deepks.utils import check_list
from deepks.utils import flat_file_list
from deepks.utils import get_sys_name, load_sys_paths
from deepks.task.task import PythonTask, ShellTask
from deepks.task.task import BatchTask, GroupBatchTask, DPDispatcherTask
from deepks.task.workflow import Sequence
from deepks.utils import QCDIR


SCF_CMD = " ".join([
    "{python} -u",
    "-m deepks.scf.run"
    # os.path.join(QCDIR, "scf/run.py") # this is the backup choice
])

TRN_CMD = " ".join([
    "{python} -u",
    "-m deepks.model.train"
    # os.path.join(QCDIR, "train/train.py") # this is the backup choice
])


DEFAULT_SCF_RES = {
    "time_limit": "24:00:00",
    "cpus_per_task": 8,
    "mem_limit": 8,
    "envs": {
        "PYSCF_MAX_MEMORY": 8000
    }
}

DEFAULT_SCF_SUB_RES = {
    "numb_node": 1,
    "task_per_node": 1,
    "cpus_per_task": 8,
    "exclusive": True
}

DEFAULT_TRN_RES = {
    "time_limit": "24:00:00",
    "cpus_per_task": 8,
    # "numb_gpu": 1, # do not use gpu by default
    "mem_limit": 8
}

DEFAULT_DPDISPATCHER_RES = {
    "number_node": 1,
    "cpu_per_node": 8,
    "group_size": 1
}

def check_system_names(systems):
    sys_names = [get_sys_name(os.path.basename(s)) for s in systems]
    if len(set(sys_names)) != len(systems):
        raise ValueError("Systems have duplicated base names. Not supported yet.")


def make_cleanup(pattern="slurm-*.out", workdir=".", **task_args):
    pattern = check_list(pattern)
    pattern = " ".join(pattern)
    assert pattern
    return ShellTask(
        f"rm -r {pattern}",
        workdir=workdir,
        **task_args
    )


def make_scf_task(*, workdir=".",
                  arg_file="scf_input.yaml", source_arg=None,
                  model_file="model.pth", source_model=None,
                  proj_basis=None, source_pbasis=None,
                  systems="systems.raw", link_systems=True, 
                  dump_dir="results", share_folder="share", 
                  outlog="log.scf", group_data=None,
                  dispatcher=None, resources=None, 
                  python="python", **task_args):
    # set up basic args
    command = SCF_CMD.format(python=python)
    link_share = task_args.pop("link_share_files", [])
    link_prev = task_args.pop("link_prev_files", [])
    link_abs = task_args.pop("link_abs_files", [])
    forward_files = task_args.pop("forward_files", [])
    backward_files = task_args.pop("backward_files", [])
    sys_name = None
    #set up optional args
    if arg_file:
        command += f" {arg_file}"
        if source_arg is not None:
            link_share.append((source_arg, arg_file))
        forward_files.append(arg_file)
    if model_file:
        command += f" -m {model_file}"
        if model_file.upper() != "NONE":
            if source_model is not None:
                link_prev.append((source_model, model_file))
            forward_files.append(model_file)
    if proj_basis:
        command += f" -P {proj_basis}"
        if source_pbasis is not None:
            link_share.append((source_pbasis, proj_basis))
        forward_files.append(proj_basis)
    if systems:
        # check system paths and make forward files
        sys_paths = [os.path.abspath(s) for s in load_sys_paths(systems)]
        sys_base = [get_sys_name(s) for s in sys_paths]
        sys_name = [os.path.basename(s) for s in sys_base]
        if link_systems:
            target_dir = "systems"
            src_files = sum((glob(f"{base}*") for base in sys_base), [])
            for fl in src_files:
                dst = os.path.join(target_dir, os.path.basename(fl))
                link_abs.append((fl, dst))
            forward_files.append(target_dir)
            sys_str= os.path.join(target_dir, "*")
        else: # cannot forward files here
            sys_str = " ".join(sys_paths)
        command += f" -s {sys_str}"
    if dump_dir:
        command += f" -d {dump_dir}"
        if sys_name:
            for nm in sys_name:
                backward_files.append(os.path.join(dump_dir, nm))
        else:  # backward whole folder, may cause problem
            backward_files.append(dump_dir)
    if group_data is not None:
        command += " -G" if group_data else " -NG"
    # make task
    return BatchTask(
        command, 
        workdir=workdir,
        dispatcher=dispatcher,
        resources=resources,
        outlog=outlog,
        share_folder=share_folder,
        link_share_files=link_share,
        link_prev_files=link_prev,
        link_abs_files=link_abs,
        forward_files=forward_files,
        backward_files=backward_files,
        **task_args
    )


def make_run_scf(systems_train, systems_test=None, *,
                 train_dump="data_train", test_dump="data_test", 
                 no_model=False, group_data=None,
                 workdir='.', share_folder='share', outlog="log.scf",
                 source_arg="scf_input.yaml", source_model="model.pth",
                 source_pbasis=None, dispatcher=None, resources=None, 
                 sub_size=1, group_size=1, ingroup_parallel=1, 
                 sub_res=None, python='python', **task_args):
    # if no test systems, use last one in train systems
    systems_train = [os.path.abspath(s) for s in load_sys_paths(systems_train)]
    systems_test = [os.path.abspath(s) for s in load_sys_paths(systems_test)]
    if not systems_test:
        systems_test.append(systems_train[-1])
        # if len(systems_train) > 1:
        #     del systems_train[-1]
    check_system_names(systems_train)
    check_system_names(systems_test)
    # split systems into groups
    nsys_trn = len(systems_train)
    nsys_tst = len(systems_test)
    ntask_trn = int(np.ceil(nsys_trn / sub_size))
    ntask_tst = int(np.ceil(nsys_tst / sub_size))
    train_sets = [systems_train[i::ntask_trn] for i in range(ntask_trn)]
    test_sets = [systems_test[i::ntask_tst] for i in range(ntask_tst)]
    # make subtasks
    model_file = "../model.pth" if not no_model else "NONE"
    proj_basis = "../proj_basis.npz" if source_pbasis else None
    nd = max(len(str(ntask_trn+ntask_tst)), 2)
    if sub_res is None:
        sub_res = {}
    sub_res = {**DEFAULT_SCF_SUB_RES, **sub_res}
    trn_tasks = [
        make_scf_task(systems=sset, workdir=f"task.trn.{i:0{nd}}", 
                      arg_file="../scf_input.yaml", source_arg=None,
                      model_file=model_file, source_model=None,
                      proj_basis=proj_basis, source_pbasis=None,
                      dump_dir=f"../{train_dump}", group_data=group_data,
                      link_systems=True, resources=sub_res, python=python)
        for i, sset in enumerate(train_sets)
    ]
    tst_tasks = [
        make_scf_task(systems=sset, workdir=f"task.tst.{i:0{nd}}", 
                      arg_file="../scf_input.yaml", source_arg=None,
                      model_file=model_file, source_model=None,
                      proj_basis=proj_basis, source_pbasis=None,
                      dump_dir=f"../{test_dump}", group_data=group_data, 
                      link_systems=True, resources=sub_res, python=python)
        for i, sset in enumerate(test_sets)
    ]
    # set up optional args
    link_share = task_args.pop("link_share_files", [])
    link_share.append((source_arg, "scf_input.yaml"))
    if source_pbasis:
        link_share.append((source_pbasis, "proj_basis.npz"))
    link_prev = task_args.pop("link_prev_files", [])
    if not no_model:
        link_prev.append((source_model, "model.pth"))
    if resources is None:
        resources = {}
    resources = {**DEFAULT_SCF_RES, "numb_node": ingroup_parallel, **resources}
    # make task
    return GroupBatchTask(
        trn_tasks + tst_tasks,
        workdir=workdir,
        group_size=group_size,
        ingroup_parallel=ingroup_parallel,
        dispatcher=dispatcher,
        resources=resources,
        outlog=outlog,
        errlog="err",
        share_folder=share_folder,
        link_share_files=link_share,
        link_prev_files=link_prev
    )


def make_stat_scf(systems_train, systems_test=None, *, 
                  train_dump="data_train", test_dump="data_test", group_data=False,
                  workdir='.', outlog="log.data", **stat_args):
    # follow same convention for systems as run_scf
    systems_train = [os.path.abspath(s) for s in load_sys_paths(systems_train)]
    systems_test = [os.path.abspath(s) for s in load_sys_paths(systems_test)]
    if not systems_test:
        systems_test.append(systems_train[-1])
        # if len(systems_train) > 1:
        #     del systems_train[-1]
    # load stats function
    from deepks.scf.stats import print_stats
    stat_args.update(
        systems=systems_train,
        test_sys=systems_test,
        dump_dir=train_dump,
        test_dump=test_dump,
        group=group_data)
    # make task
    return PythonTask(
        print_stats,
        call_kwargs=stat_args,
        outlog=outlog,
        errlog="err",
        workdir=workdir
    )


def make_scf(systems_train, systems_test=None, *,
             train_dump="data_train", test_dump="data_test",
             no_model=False, workdir='00.scf', share_folder='share',
             source_arg="scf_input.yaml", source_model="model.pth",
             source_pbasis=None, dispatcher=None, resources=None, 
             dpdispatcher_machine=None, dpdispatcher_resources=None, 
             sub_size=1, group_size=1, ingroup_parallel=1, 
             sub_res=None, python='python', 
             cleanup=False, **task_args):
    run_scf = make_run_scf(
        systems_train, systems_test,
        train_dump=train_dump, test_dump=test_dump, 
        no_model=no_model, group_data=False,
        workdir=".", outlog="log.scf", share_folder=share_folder, 
        source_arg=source_arg, source_model=source_model, source_pbasis=source_pbasis,
        dispatcher=dispatcher, resources=resources, 
        dpdispatcher_machine=dpdispatcher_machine, 
        dpdispatcher_resources=dpdispatcher_resources,
        group_size=group_size, ingroup_parallel=ingroup_parallel,
        sub_size=sub_size, sub_res=sub_res, python=python, **task_args
    )
    post_scf = make_stat_scf(
        systems_train=systems_train, systems_test=systems_test,
        train_dump=train_dump, test_dump=test_dump, workdir=".", 
        outlog="log.data", group_data=False
    )
    # concat
    seq = [run_scf, post_scf]
    if cleanup:
        clean_scf = make_cleanup(
            ["slurm-*.out", "task.*/err", "fin.record"],
            workdir=".")
        seq.append(clean_scf)
    # make sequence
    return Sequence(
        seq,
        workdir=workdir
    )


def make_train_task(*, workdir=".",
                    arg_file="train_input.yaml", source_arg=None,
                    restart_model=None, source_model=None, 
                    proj_basis=None, source_pbasis=None,
                    save_model="model.pth", group_data=False,
                    data_train="data_train", source_train=None,
                    data_test="data_test", source_test=None,
                    share_folder="share", outlog="log.train",
                    dispatcher=None, resources=None, 
                    dpdispatcher_machine=None, 
                    dpdispatcher_resources=None,
                    python="python", **task_args):
    # set up basic args
    command = TRN_CMD.format(python=python)
    link_share = task_args.pop("link_share_files", [])
    link_prev = task_args.pop("link_prev_files", [])
    forward_files = task_args.pop("forward_files", [])
    backward_files = task_args.pop("backward_files", [])
    # set up optional args
    if arg_file:
        command += f" {arg_file}"
        if source_arg is not None:
            link_share.append((source_arg, arg_file))
        forward_files.append(arg_file)
    if restart_model:
        command += f" -r {restart_model}"
        if source_model is not None:
            link_prev.append((source_model, restart_model))
        forward_files.append(restart_model)
    if proj_basis:
        command += f" -P {proj_basis}"
        if source_pbasis is not None:
            link_share.append((source_pbasis, proj_basis))
        forward_files.append(proj_basis)
    if data_train:
        command += f" -d {data_train}" + ("" if group_data else "/*")
        if source_train is not None:
            link_prev.append((source_train, data_train))
        forward_files.append(data_train)
    if data_test:
        command += f" -t {data_test}" + ("" if group_data else "/*")
        if source_test is not None:
            link_prev.append((source_test, data_test))
        forward_files.append(data_test)
    if save_model:
        command += f" -o {save_model}"
        backward_files.append(save_model)
    if resources is None:
        resources = {}
    resources = {**DEFAULT_TRN_RES, **resources}
    # make task
    task_list=[]
    if dispatcher=="dpdispatcher":
        from dpdispatcher import Task
        task_list.append(Task(command=command, task_work_path='./', 
                        forward_files=[], backward_files=[], 
                        outlog=outlog, errlog="err.train"))
        if dpdispatcher_resources is None:
            dpdispatcher_resources = {}
        dpdispatcher_resources = {**DEFAULT_DPDISPATCHER_RES, **dpdispatcher_resources}
        return DPDispatcherTask(
            task_list, 
            work_base=workdir,
            machine=dpdispatcher_machine, 
            resources=dpdispatcher_resources, 
            outlog=outlog,
            errlog='err',
            share_folder=share_folder,
            link_share_files=link_share,
            link_prev_files=link_prev,
            forward_files=forward_files,
            backward_files=backward_files,
            **task_args
        )
    else:
        return BatchTask(
            command,
            workdir=workdir,
            dispatcher=dispatcher,
            resources=resources,
            outlog=outlog,
            errlog='err',
            share_folder=share_folder,
            link_share_files=link_share,
            link_prev_files=link_prev,
            forward_files=forward_files,
            backward_files=backward_files,
            **task_args
        )


def make_run_train(source_train="data_train", source_test="data_test", *,
                   restart=True, source_model="model.pth", save_model="model.pth", 
                   source_pbasis=None, source_arg="train_input.yaml", 
                   workdir=".", share_folder="share", outlog="log.train",
                   dispatcher=None, resources=None, 
                   dpdispatcher_machine=None, 
                   python="python", **task_args):
    # just add some presetted arguments of make_train_task
    # have not implement parrallel training for multiple models
    restart_model = "old_model.pth" if restart else None
    proj_basis = "proj_basis.npz" if source_pbasis else None
    return make_train_task(
        workdir=workdir, 
        arg_file="train_input.yaml", source_arg=source_arg,
        restart_model=restart_model, source_model=source_model, 
        proj_basis=proj_basis, source_pbasis=source_pbasis,
        save_model=save_model, group_data=False,
        data_train="data_train", source_train=source_train,
        data_test="data_test", source_test=source_test,
        share_folder=share_folder, outlog=outlog,
        dispatcher=dispatcher, resources=resources,
        dpdispatcher_machine=dpdispatcher_machine, 
        python=python, **task_args
    )


def make_test_train(data_paths, model_file="model.pth", *,
                    output_prefix="test", group_results=True, 
                    workdir='.', outlog="log.test", **test_args):
    from deepks.model.test import main as test_func
    test_args.update(
        data_paths=data_paths,
        model_file=model_file,
        output_prefix=output_prefix,
        group=group_results)
    # make task
    return PythonTask(
        test_func,
        call_kwargs=test_args,
        outlog=outlog,
        errlog="err",
        workdir=workdir
    )


def make_train(source_train="data_train", source_test="data_test", *,
               restart=True, source_model="model.pth", save_model="model.pth",
               source_pbasis=None, source_arg="train_input.yaml", 
               workdir="01.train", share_folder="share",
               dispatcher=None, resources=None, 
               dpdispatcher_machine=None, 
               python="python", cleanup=False, **task_args):
    run_train = make_run_train(
        source_train=source_train, source_test=source_test,
        restart=restart, source_model=source_model, save_model=save_model,
        source_pbasis=source_pbasis, source_arg=source_arg, 
        workdir=".", share_folder=share_folder,
        outlog="log.train", dispatcher=dispatcher, resources=resources,
        dpdispatcher_machine=dpdispatcher_machine, 
        python=python, **task_args
    )
    post_train = make_test_train(
        data_paths=["data_train/*","data_test/*"],
        model_file=save_model, output_prefix="test", group_results=True,
        workdir=".", outlog="log.test"
    )
    # concat
    seq = [run_train, post_train]
    # seq = [run_train]
    if cleanup:
        clean_train = make_cleanup(
            ["slurm-*.out", "err", "fin.record", "tag_*finished"],
            workdir=".")
        seq.append(clean_train)
    # make sequence
    return Sequence(
        seq,
        workdir=workdir
    )