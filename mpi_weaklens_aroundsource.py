import sys
sys.path.append("install/lib/python2.7/site-packages/")
sys.path.append("install/lib64/python2.7/site-packages/")
import weaklens
import frogress
import pandas
import numpy as np
import sys
from astropy.io import fits as pyfits
import fitsio
import argparse
import yaml
from weaklens_select import lens_select, source_select, get_pofz_array
from scipy.interpolate import interp1d
from mpi4py import MPI
from weaklens_aroundsource import run_pipeline

if __name__ == "__main__":

    comm = MPI.COMM_WORLD

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config", help="Configuration file")
    parser.add_argument("--rmax", help="Max radius, does not override the value in config file", type=float, default=80.0)
    parser.add_argument("--rmin", help="Min radius, does not override the value in config file", type=float, default=0.5)
    parser.add_argument("--rbin", help="Radial bin, does not override the value in config file", type=int, default=25)
    parser.add_argument("--random", help="Radial bin, does not override the value in config file", type=bool, default=False)
    
    args = parser.parse_args()
    
    with open(args.config, 'r') as ymlfile:
        config = yaml.load(ymlfile)
    print config


    if "rmax" not in config:
        config["rmax"] = args.rmax
    if "rmin" not in config:
        config["rmin"] = args.rmin
    if "rbin" not in config:
        config["rbin"] = args.rbin
    if "Omegam" not in config:
        config["Omegam"] = 0.279

    if args.random:
        config["random_rotate_seed"] = comm.rank + 875324

    config["output_dsigma"] = config["output_dsigma"] + "%05d" % comm.rank

    run_pipeline(config)
