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


def run_pipeline(config, verbose=0):
    # Convert arguments to variables
    rmax = config["rmax"]
    rmin = config["rmin"]
    rbin = config["rbin"]

    if "oversample_pofz" not in config:
        config["oversample_pofz"] = 1

    if "zdiff" not in config:
        config["zdiff"] = 0.0

    if "output_dsigma" not in config:
        config["output_dsigma"] = "Dsigma.dat"
    
    pairout = config["outputdir"] + "/" + config["output_dsigma"]

    if "outputpairfile" not in config:
        config["outputpairfile"] = ""
    elif config["outputpairfile"] != "":
        config["outputpairfile"] = config["outputdir"] + "/" + config["outputpairfile"]

    
    from subprocess import call
    call("mkdir -p %s" % (config["outputdir"]), shell=1)
    # call("cp %s %s/" % (ymlfile, config["outputdir"]), shell=1)
    
    wpipe = weaklens.weaklens(rmin, rmax, rbin, pairout, False, config["Omegam"], config["outputpairfile"])
    
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
        pofz_zz = get_pofz_array(sourceargs).astype("float32")

        if config["oversample_pofz"] > 1:
            new_size = config["oversample_pofz"]*(pofz_zz.size-1) + 1
            old_pofz_zz = pofz_zz * 1.0
            pofz_zz = np.linspace(old_pofz_zz[0], old_pofz_zz[-1], new_size)

        wpipe.setup_pofz_array(pofz_zz)
        if verbose:
            print pofz_zz.size
            pofz_zz_arr = wpipe.doubleArray(pofz_zz.size)
            for i in range(pofz_zz.size):
                print i, pofz_zz[i]
                pofz_zz_arr[i] = pofz_zz[i]
    elif "photoz_estimate" not in sourceargs:
        sourceargs["photoz_estimate"] = "_best"

    # Do we want to randomly rotate sources for shape noise estimates
    if "random_rotate_seed" in config:
        wpipe.setup_random_rotate(config["random_rotate_seed"])
    
    while not done:
        itern = itern + 1
        datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag = source_select(sourceargs, chunksize=1000000)
    
        if status!=0 :
            break
        else:
            print "\n", datagalflag.size, pofzflag.size
    
        
        # Uncomment the following two rows for debugging purposes: 
        # if itern==2:
        #     break
    
        # For every source, query the clusters around it
        try:
            sys.stdout.write("\n %s: \n" % (sourceargs['fits_list'][sourceargs['iter']-1]))
        except:
            sys.stdout.write("\n Passed: \n")
    
        #for i in frogress.bar(range(Ngal)):
        for i in range(Ngal):

            if not datagalflag[i] or not pofzflag[i]:
                continue

            if sourceargs['filetype'] == "fits":
                if sourceargs['include_corrections']:
                    if sourceargs['format'] == "des":
                        ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, c1_nb, c2_nb = datagal[i]
                    elif sourceargs["format"] == "hsc":
                        ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, R2gal = datagal[i]
                        c1_nb = 0.0
                        c2_nb = 0.0
                    zphotgal = 0.0
                    e1gal = np.double(e1gal)
                    e2gal = -np.double(e2gal)
                    c2_dp = -c2_dp
                    c2_nb = -c2_nb
                
            elif sourceargs['filetype'] == "ascii":
                ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal, R2gal = datagal.ra.values[i], datagal.dec.values[i], datagal.e1.values[i],datagal.e2.values[i],datagal.weight.values[i],datagal.m.values[i],datagal.c1_dp.values[i],datagal.c2_dp.values[i],datagal.c1_nb.values[i],datagal.c2_nb.values[i],datagal.z_b.values[i], datagal.R2gal.values[i]
    
            if sourceargs["fullpofz"]:
                if np.isnan(pofz[i][0]):
                    continue
                else:
                    pass

                if config["oversample_pofz"] > 1:
                    spl = interp1d(old_pofz_zz, np.array(pofz[i], dtype=np.float32), kind="linear")
                    pofz_array = spl(pofz_zz)
                else:
                    pofz_array = np.array(pofz[i], dtype=np.float32)

                status = wpipe.process_pofz(pofz_array)
                if status:
                    continue
            else:
                zphotgal = np.float64(pofz[i])
    
            if "integral_cut_zl" in sourceargs:
                wpipe.process_source(ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal, R2gal, sourceargs["fullpofz"], sourceargs["integral_cut_zl"], sourceargs["integral_cut_zdiff"], sourceargs["integral_cut_Pth"], sourceargs["integral_cut_zmax"])
            elif "zdiff" in sourceargs and not sourceargs["fullpofz"]:
                wpipe.process_source(ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal, R2gal, sourceargs["fullpofz"], -99.0, sourceargs["zdiff"], 0.0, 99.0)
            else:
                wpipe.process_source(ragal, decgal, e1gal, e2gal, wgal, rms_egal, mgal, c1_dp, c2_dp, c1_nb, c2_nb, zphotgal, R2gal, sourceargs["fullpofz"])
        wpipe.finalize_results()
    
    wpipe.finalize_results(writeok=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--config", help="Configuration file")
    parser.add_argument("--rmax", help="Max radius, does not override the value in config file", type=float, default=80.0)
    parser.add_argument("--rmin", help="Min radius, does not override the value in config file", type=float, default=0.5)
    parser.add_argument("--rbin", help="Radial bin, does not override the value in config file", type=int, default=25)
    
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


    run_pipeline(config)
