import os
import sys
import time
import torch
import argparse
import numpy as np
from pyscf import gto
from model import QCNet
from scf import DeepSCF


def parse_xyz(filename, basis='ccpvtz', verbose=0):
    with open(filename) as fp:
        natoms = int(fp.readline())
        comments = fp.readline()
        xyz_str = "".join(fp.readlines())
    mol = gto.Mole()
    mol.verbose = verbose
    mol.atom = xyz_str
    mol.basis  = basis
    mol.build(0,0,unit="Ang")
    return mol


def solve_mol(mol, model, chkfile=None, verbose=0):
    if verbose:
        tic = time.time()
    cf = DeepSCF(mol, model)
    if chkfile:
        cf.set(chkfile=chkfile)
    ecf = cf.scf()
    if verbose:
        tac = time.time()
        print(f"time of scf: {tac - tic}, converged: {cf.converged}")
    natom = mol.natm
    nao = mol.nao
    nproj = sum(cf.shell_sec)
    meta = np.array([natom, nao, nproj])
    dm = cf.make_rdm1()
    eig = cf.make_eig(dm)
    ehf = cf.energy_tot0(dm)
    return meta, ehf, ecf, dm, eig, cf.converged


def dump_data(dir_name, meta, **data_dict):
    os.makedirs(dir_name, exist_ok = True)
    np.savetxt(os.path.join(dir_name, 'system.raw'), 
               np.reshape(meta, (1,-1)), 
               fmt = '%d', header = 'natom nao nproj')
    for name, value in data_dict.items():
        np.save(os.path.join(dir_name, f'{name}.npy'), value)


def main(xyz_files, model_path, dump_dir=None, group=False, verbose=0):

    model = QCNet.load(model_path)
    if dump_dir is None:
        dump_dir = os.curdir
    if group:
        results = []
        
    for fl in xyz_files:
        mol = parse_xyz(fl, verbose=verbose)
        try:
            result = solve_mol(mol, model, verbose=verbose)
        except Exception as e:
            print(fl, 'failed! error:', e, file=sys.stderr)
            continue
        if not group:
            meta, ehf, ecf, dm, eig, conv = result
            natom, nao, nproj = meta
            sub_dir = os.path.join(dump_dir, os.path.splitext(os.path.basename(fl))[0])
            dump_data(sub_dir, meta,
                e_hf=np.reshape(ehf, (1,1)),
                e_cf=np.reshape(ecf, (1,1)),
                dm_eig=np.reshape(eig, (1, natom, nproj)),
                conv=np.reshape(conv, (1,1)))
        else:
            results.append(result)
            if any(result[0] != results[0][0]):
                print(fl, 'meta does not match! saving previous results only.', file=sys.stderr)
                break
        if verbose:
            print(fl, 'finished')

    if group:
        nframes = len(results)
        meta, ehf, ecf, dm, eig, conv = zip(*results)
        natom, nao, nproj = meta[0]
        dump_data(dump_dir, meta[0],
            e_hf=np.reshape(ehf, (nframes,1)),
            e_cf=np.reshape(ecf, (nframes,1)),
            dm_eig=np.reshape(eig, (nframes, natom, nproj)),
            conv=np.reshape(conv, (nframes,1)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate and save SCF energies and descriptors using given model.")
    parser.add_argument("files", nargs="+", 
                        help="input xyz files")
    parser.add_argument("-m", "--model-path", default='model.pth', 
                        help="path to the trained model, default: model.pth")
    parser.add_argument("-d", "--dump-dir", default='.', 
                        help="dir of dumped files, default: current dir")
    parser.add_argument("-G", "--group", action='store_true',
                        help="group results for all molecules, only works for same system")
    parser.add_argument("-v", "--verbose", default=0, type=int, 
                        help="output calculation information")
    args = parser.parse_args()
    
    main(args.files, args.model_path, args.dump_dir, args.group, args.verbose)