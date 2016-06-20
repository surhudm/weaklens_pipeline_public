import numpy as np
import pyfits
import pandas
import sys
import fitsio

def lens_select(lensargs):
    if lensargs['type'] == "lowz_rcsrestrict":
        ra, dec, zred, lam, wt = np.loadtxt("DataStore/rcsrestrict_dr12_v5_lowz_All_catalog.dat", unpack=1)

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]

def source_select(sourceargs, chunksize):

    if sourceargs['type'] == "rcslens" and sourceargs['filetype'] == "fits":
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits'] = fitsio.FITS('DataStore/rcslens_data.fits')
            sourceargs['nrows'] = sourceargs['fits'][1].read_header()['naxis2']
            if sourceargs["fullpofz"]:
                sourceargs['fitspofz'] = fitsio.FITS('DataStore/rcslens_pofz.fits')

        datagal = 0
        status = (itern*chunksize>=sourceargs['nrows'])
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz

        wheremask = fits[1].where("MAG_extcorr_r < %.5f && #row>%d && #row<=%d" % (sourceargs['maglim'], itern*chunksize, (itern+1)*chunksize))
        datagal = sourceargs['fits'][1]['ra','dec','e1','e2','weight','m','c1_dp','c2_dp','c1_nb','c2_nb','z_b'][wheremask]
        if sourceargs["fullpofz"]:
            pofz = sourceargs['fitspofz'][1][wheremask]
            pofz = pofz.view(pofz.dtype[0]).reshape(pofz.shape + (-1,))

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]

    if sourceargs['type'] == "rcslens" and filetype == "ascii":
        itern = sourceargs['iter']
        if itern == 0:
            import pandas
            sourceargs['dfchunks'] = pandas.read_csv("rcslens_ptestimate_zbest.tsv", names=(['ra','dec','e1','e2','weight','m','c1_dp','c2_dp','c1_nb','c2_nb','z_b']), usecols=(range(1, 12)), chunksize=chunksize, delim_whitespace=1, skiprows=1)
        datagal = sourceargs['dfchunks'].next()

        pofz = 0
        if sourceargs["fullpofz"]:
            print "Photo-z pdf not implemented in ascii"
            exit(11)
        if sourceargs["maglim"]:
            print "Maglim not implemented"
            exit(11)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = datagal.ra.size
        status = 0
    return datagal, sourceargs, Ngal, status, pofz
