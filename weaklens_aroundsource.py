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
from weaklens_select import lens_select, source_select, get_pofz_array

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--config", help="Configuration file")
parser.add_argument("--rmax", help="Max radius", type=float, default=20.0)
parser.add_argument("--rmin", help="Min radius", type=float, default=0.1)
parser.add_argument("--rbin", help="Radial bin", type=int, default=20)

args = parser.parse_args()

with open(args.config, 'r') as ymlfile:
    config = yaml.load(ymlfile)
print config
print args

# Convert arguments to variables
rmax = (args.rmax)
rmin = (args.rmin)
rbin = args.rbin
pairout = config["outputdir"] + "/pairs.dat"

from subprocess import call
call("mkdir -p %s" % (config["outputdir"]), shell=1)
call("cp %s %s/" % (ymlfile, config["outputdir"]), shell=1)

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
####wpipe.setup_pofz(sourceargs["pofz_zmin"], sourceargs["pofz_zdiff"], sourceargs["pofz_zbins"])
    pofz_zz = get_pofz_array(sourceargs).astype("float32")
    wpipe.setup_pofz_array(pofz_zz)
    ## print pofz_zz.size
    ## pofz_zz_arr = wpipe.doubleArray(pofz_zz.size)
    ## for i in range(pofz_zz.size):
    ##     print i, pofz_zz[i]
    ##     pofz_zz_arr[i] = pofz_zz[i]

#### try:
####     fout = open("Test_%s" % sourceargs["pofz_type"], "w")
#### except:
####     fout = open("Test_debug", "w")

while not done:
    itern = itern + 1
    #datagal, sourceargs, Ngal, status, pofz = source_select(sourceargs, chunksize=1000000)
    datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag = source_select(sourceargs, chunksize=1000000)

    if status!=0 :
        break

    # if itern==2:
    #     break

    # For every source, query the clusters around it
    try:
        sys.stdout.write("\n %s: \n" % (sourceargs['fits_list'][sourceargs['iter']-1]))
    except:
        sys.stdout.write("\n Passed: \n")

    for i in frogress.bar(range(Ngal)):
        if not datagalflag[i] or not pofzflag[i]:
            continue
    #for i in range(Ngal):
    #for i in range(10):

        if sourceargs['filetype'] == "fits":
            if sourceargs['include_corrections']:
                ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, c1_nb, c2_nb = datagal[i]
                zphotgal = 0.0
                e2gal = -e2gal
                c2_dp = -c2_dp
                c2_nb = -c2_nb
            else:
                ragal, decgal, e1gal, e2gal, wgal, rms_egal = datagal[i]
                if sourceargs['type'] != "hsc-wide-josh":
                    #print "\nbefore:", ragal, decgal, e1gal, e2gal, wgal
                    wgal = 1./(0.365**2+wgal**2)
                    e2gal = -e2gal
                    #print "\nafter:", ragal, decgal, e1gal, e2gal, wgal
                e1gal = np.float64(e1gal)
                e2gal = np.float64(e2gal)
                mgal = 0.0*ragal
                c1_dp = 0.0*ragal
                c2_dp = 0.0*ragal
                c1_nb = 0.0*ragal
                c2_nb = 0.0*ragal
                zphotgal = 0.0*ragal
                rms_egal = np.float64(rms_egal)
                #print ragal, decgal, e1gal, e2gal, wgal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal

        elif sourceargs['filetype'] == "ascii" and sourceargs['type'] == "alexie-test":
            ragal, decgal, zphotgal, e1gal, e2gal, wgal = datagal.ra.values[i], datagal.dec.values[i], datagal.z_b.values[i],datagal.e1.values[i],datagal.e2.values[i],datagal.weight.values[i]
            e2gal = -e2gal
            mgal, c1_dp, c2_dp, c1_nb, c2_nb = ragal*0.0, ragal*0.0, ragal*0.0, ragal*0.0, ragal*0.0
            rms_egal = ragal*0.0+1./np.sqrt(2.)
            
        elif sourceargs['filetype'] == "ascii":
            ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal = datagal.ra.values[i], datagal.dec.values[i], datagal.e1.values[i],datagal.e2.values[i],datagal.weight.values[i],datagal.m.values[i],datagal.c1_dp.values[i],datagal.c2_dp.values[i],datagal.c1_nb.values[i],datagal.c2_nb.values[i],datagal.z_b.values[i]

        if sourceargs["fullpofz"]:
            if np.isnan(pofz[i][0]):
                continue
            else:
                pass
            status = wpipe.process_pofz(np.array(pofz[i], dtype=np.float32))
            if status:
                continue
        else:
            zphotgal = np.float64(pofz[i])

        wpipe.process_source(ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal, sourceargs["fullpofz"])
    wpipe.finalize_results()

wpipe.finalize_results()
