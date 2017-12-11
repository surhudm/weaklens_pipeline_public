import yaml
import numpy as np
import sys
import argparse
from subprocess import call
from weaklens_aroundsource import run_pipeline
from mpi4py import MPI

def getsource_dict(pofz_type, realization, rotation, integral_cut_zl, integral_cut_zdiff, integral_cut_Pth, integral_cut_zmax, include_corrections=True):
    d = {}
    d["type"] = "mockhsc-wide-s16a_v2.0_calib"
    d["filetype"] = "fits"
    d["fullpofz"] = True
    d["pofz_type"] = pofz_type
    d["include_corrections"] = include_corrections
    d['rotation'] = rotation
    d['realization'] = realization
    d['integral_cut_zl'] = integral_cut_zl
    d['integral_cut_zdiff'] = integral_cut_zdiff
    d['integral_cut_Pth'] = integral_cut_Pth
    d['integral_cut_zmax'] = integral_cut_zmax
    return d

#
# Nrun = realization + (rotation*Nreal) + zbin*(Nrot*Nreal) + random*(Nrot*Nreal*Nz)
#
def parse_run_number(Nrun, Nreal, Nrot, Nz):
    realization = (Nrun%Nreal)
    Nrun = (Nrun-realization)/Nreal
    rotation = (Nrun%Nrot)
    Nrun = (Nrun-rotation)/Nrot
    zbin = (Nrun%Nz)
    Nrun = (Nrun-zbin)/Nz
    random = Nrun * 1
    return realization, rotation, zbin, random

def parallel_job(config):
    run_pipeline(config)

if __name__=="__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--Omegam", help="Omegam", type=float)
    parser.add_argument("--dir", help="Directory to output")
    parser.add_argument("--rmax", help="Max radius", type=float)
    parser.add_argument("--rmin", help="Min radius", type=float)
    parser.add_argument("--rbin", help="Radial bin", type=int)
    parser.add_argument("--pofz", help="Photoz type")
    parser.add_argument("--mockdir", help="Mock directory, e.g. Masato-mock-random-full or Masato-mock-random")
    parser.add_argument("--startrun", help="Start Run number", type=int)
    parser.add_argument("--total_realizations", help="Total realizations", type=int)
    parser.add_argument("--total_rotations", help="Total rotations", type=int)
    parser.add_argument("--total_zbins", help="Total zbins", type=int)
    parser.add_argument("--integral_cut_zl", help="Integrated P(z) cut parameter zl", type=float)
    parser.add_argument("--integral_cut_zdiff", help="Integrated P(z) cut parameter zdiff", type=float)
    parser.add_argument("--integral_cut_Pth", help="Integrated P(z) cut parameter Pth", type=float)
    parser.add_argument("--integral_cut_zmax", help="Integrated P(z) cut parameter zmax", type=float)
    
    args = parser.parse_args()
    run_number = rank+args.startrun
    realization, rotation, zbin, random = parse_runnumber(run_number, args.total_realizations, args.total_rotations, args.total_zbins)
    
    dsample = {}
    dsample["type"] = args.mockdir
    dsample["zbin"] = zbin
    dsample['rotation'] = rotation
    dsample['realization'] = realization
    dsample["rannum"] = random
        
    ssample = getsource_dict(args.pofz, realization, rotation, args.integral_cut_zl, args.integral_cut_zdiff, args.integral_cut_Pth, args.integral_cut_zmax)
        
    d = {"lens":dsample, "source":ssample, "outputdir":"%s/%03d_%03d_%02d_%02d/" % (args.dir, realization, rotation, zbin, random), "Omegam": args.Omegam, "rmin":args.rmin, "rmax":args.rmax, "rbin":args.rbin}
    configfile = "%s/ymls/%03d_%03d_%02d_%02d.yml" % (args.dir, realization, rotation, zbin, random)
    call("mkdir -p %s/ymls/" % (args.dir, rotation), shell=1)
    with open(configfile, 'w') as yaml_file:
        yaml.dump(d, yaml_file, default_flow_style=False)
    #run_pipeline(d)
    
    comm.Barrier()
