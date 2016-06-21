import sys
sys.path.append("install/lib/python2.7/site-packages/")
import weaklens
import frogress
import pandas
import numpy as np
import sys
import pyfits
import fitsio
import argparse
import yaml
from weaklens_select import lens_select, source_select

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--config", help="Configuration file")
parser.add_argument("--rmax", help="Max radius", type=float, default=50.0)
parser.add_argument("--rmin", help="Min radius", type=float, default=0.5)
parser.add_argument("--rbin", help="Radial bin", type=int, default=20)
parser.add_argument("--pairout", help="Output filename with pairs information", default="debug_pairs.dat")

args = parser.parse_args()

with open(args.config, 'r') as ymlfile:
    config = yaml.load(ymlfile)
print config
print args

# Convert arguments to variables
rmax = (args.rmax)
rmin = (args.rmin)
rbin = args.rbin
pairout = args.pairout

wpipe = weaklens.weaklens(rmin, rmax, rbin, pairout, False)

# First get all the lens galaxy data
print config["lens"]
ra, dec, zred, wt = lens_select(config["lens"])

wpipe.allocate_lens_memory(ra.size)
for i in range(ra.size):
    wpipe.process_lens(ra[i], dec[i], zred[i], wt[i])
tree = wpipe.finalize_lenses()


# Ready to pounce on the source data, a million galaxies at a time
itern=0
done = False
sourceargs = config["source"]
sourceargs["iter"] = 0

if sourceargs["fullpofz"]:
    print "Setting up P(z)"
    wpipe.setup_pofz(sourceargs["pofz_zmin"], sourceargs["pofz_zdiff"], sourceargs["pofz_zbins"])
print sourceargs["fullpofz"], type(sourceargs["fullpofz"])

while not done:
    #datagal, sourceargs, Ngal, status, pofz = source_select(sourceargs, chunksize=1000000)
    datagal, sourceargs, Ngal, status, pofz = source_select(sourceargs, chunksize=100000)

    if status!=0 :
        break

    # For every source, query the clusters around it
    #for i in frogress.bar(range(Ngal)):
    #for i in range(Ngal):
    for i in range(10):

        if sourceargs['filetype'] == "fits":
            if sourceargs['include_corrections']:
                ragal, decgal, e1gal, e2gal, wgal, mgal, c1_dp, c2_dp, c1_nb, c2_nb = datagal[i]
                zphotgal = 0.0
            else:
                ragal, decgal, e1gal, e2gal = datagal[i]
                wgal = 1.0
                mgal = 0.0
                c1_dp = 0.0
                c2_dp = 0.0
                c1_nb = 0.0
                c2_nb = 0.0
                zphotgal = 0.0
                print ragal, decgal, e1gal, e2gal

        elif sourceargs['filetype'] == "ascii":
            ragal, decgal, e1gal, e2gal, wgal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal = datagal.ra.values[i], datagal.dec.values[i], datagal.e1.values[i],datagal.e2.values[i],datagal.weight.values[i],datagal.m.values[i],datagal.c1_dp.values[i],datagal.c2_dp.values[i],datagal.c1_nb.values[i],datagal.c2_nb.values[i],datagal.z_b.values[i]

        if sourceargs["fullpofz"]:
            wpipe.process_pofz(np.array(pofz[i], dtype=np.float32))
        else:
            zphotgal = pofz[i]

        #exit(11)
        wpipe.process_source(ragal, decgal, e1gal, e2gal, wgal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal, sourceargs["fullpofz"])
        #break
    #break

wpipe.finalize_results()
