import numpy as np
import pyfits
import pandas
import sys
import fitsio
import glob

def lens_select(lensargs):
    if lensargs['type'] == "redmapper":
        hdulist = pyfits.open("/home/surhud/DataStore/redmapper/redmapper_dr8_public_v5.10_catalog.fits")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z_lambda"].astype('float64')
        lamda = data["lambda"]
        wt = ra/ra
        
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt

    if lensargs['type'] == "redmapper-random":
        ra, dec, zred, lamda, wt = np.loadtxt("/home/surhud/DataStore/redmapper/redmapper_public_v5.10_randoms_%05d.dat" % (lensargs["rannum"]), unpack=1) 
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d randoms \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]

    if lensargs['type'] == "sdss-vagc-lbg":
        ra, dec, zred, mstel, mag = np.loadtxt("/home/surhud/DataStore/LBG/Isolated_v0.dat", unpack=1)
        mstel = np.log10(mstel)
        wt = ra/ra

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (mstel>lensargs["xmstelmin"]) & (mstel<=lensargs["xmstelmax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]

    if lensargs['type'] == "lowz_rcsrestrict":
        ra, dec, zred, lam, wt = np.loadtxt("DataStore/rcsrestrict_dr12_v5_lowz_All_catalog.dat", unpack=1)

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]

    if lensargs['type'] == "camira":
        ra, dec, zred, lam = np.loadtxt("/home/surhud/DataStore/camira/%s/%s/camira_catalog.dat" % (lensargs['hsc-release'], lensargs['version']), unpack=1, usecols=(0, 1, 2, 3))
        wt = lam/lam

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lam>lensargs["lammin"]) & (lam<=lensargs["lammax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]


def source_select(sourceargs, chunksize):

    if sourceargs['type'] == "rcslens" and sourceargs['filetype'] == "fits":
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits'] = fitsio.FITS('/home/surhud/DataStore/RCSlens/rcslens_data.fits')
            sourceargs['nrows'] = sourceargs['fits'][1].read_header()['naxis2']
            if sourceargs["fullpofz"]:
                sourceargs['fitspofz'] = fitsio.FITS('/home/surhud/DataStore/RCSlens/rcslens_pofz.fits')

        datagal = 0
        status = (itern*chunksize>=sourceargs['nrows'])
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz

        wheremask = sourceargs['fits'][1].where("MAG_extcorr_r < %.5f && #row>%d && #row<=%d" % (sourceargs['maglim'], itern*chunksize, (itern+1)*chunksize))
        datagal = sourceargs['fits'][1]['ra','dec','e1','e2','weight','m','c1_dp','c2_dp','c1_nb','c2_nb'][wheremask]
        if sourceargs["fullpofz"]:
            pofz = sourceargs['fitspofz'][1][wheremask]
            pofz = pofz.view(pofz.dtype[0]).reshape(pofz.shape + (-1,))
        else:
            pofz = sourceargs['fits'][1]['z_b'][wheremask]

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        
        print "pofzflag and datagalflag not implemented yet"
        exit(11)

    if sourceargs['type'] == "rcslens" and sourceargs["filetype"] == "ascii":
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

        
        print "pofzflag and datagalflag not implemented yet"
        exit(11)

    if sourceargs['type'] == "hsc-wide" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:

        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("/home/surhud/DataStore/S15B_v2.1/*_tracts/*_%s.fits" % (sourceargs['pofz_type']))
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0, 0

        #print "=========="
        #print (sourceargs['fits_list'][itern].split("_pz_pdf_%s" % (sourceargs['pofz_type'])))
        #print "=========="

        sourceargs['fits'] = fitsio.FITS(sourceargs['fits_list'][itern].split("_pz_pdf_%s" % (sourceargs['pofz_type']))[0] + ".fits" )
        sourceargs['fitspofz'] = fitsio.FITS(sourceargs['fits_list'][itern])
        ntot = sourceargs['fits'][1].read_header()['naxis2']
        datagalflag = np.zeros(ntot)==np.zeros(ntot)
        pofzflag = np.zeros(ntot)==np.zeros(ntot)

        #wheremask = sourceargs['fits'][1].where()
        datagal = sourceargs['fits'][1]['ira','idec','ishape_hsm_regauss_e1','ishape_hsm_regauss_e2', 'ishape_hsm_regauss_sigma'][:]
        for_gal_flag = sourceargs['fits'][1]['imag_cmodel', 'iflux_cmodel', 'iflux_cmodel_err', 'ishape_hsm_regauss_resolution', 'ishape_hsm_regauss_e1', 'ishape_hsm_regauss_e2', 'ishape_hsm_regauss_sigma'][:]
        
        for i in range(ntot):
            imag_cmodel, iflux_cmodel, iflux_cmodel_err, ishape_hsm_regauss_resolution, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_sigma = for_gal_flag[i]
            datagalflag[i] = (imag_cmodel < 25.0) & (iflux_cmodel/iflux_cmodel_err >=10) & (ishape_hsm_regauss_resolution> 0.3) & (ishape_hsm_regauss_e1**2 + ishape_hsm_regauss_e2**2 < 4) & (ishape_hsm_regauss_sigma<0.4)

        if sourceargs['pofz_type'] == "demp":
            for_pofz_flag = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']), "%s_photoz_mode" % (sourceargs['pofz_type']), "%s_photoz_conf_median" % (sourceargs['pofz_type']) ][:]
        elif sourceargs['pofz_type'] == "mizuki":
            for_pofz_flag = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']), "%s_reduced_chisq" % (sourceargs['pofz_type']), "%s_photoz_variance" % (sourceargs['pofz_type']) ][:]
        elif sourceargs['pofz_type'] == "nnpz":
            for_pofz_flag = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']), "%s_photoz_flag" % (sourceargs['pofz_type']), "%s_is_clean" % (sourceargs['pofz_type']) ][:]
        elif sourceargs['pofz_type'] == "mlz":
            for_pofz_flag = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']), "%s_photoz_stddev_mean" % (sourceargs['pofz_type']), "%s_photoz_conf_mean" % (sourceargs['pofz_type']) ][:]
        else:
            print "Photoz not supported yet!"
            exit(11)

        for i in range(ntot):
            if sourceargs["pofz_type"] == "demp":
                z_median, z_mode, z_conf = for_pofz_flag[i]
                pofzflag[i] = (np.absolute(z_mode - z_median)/(1+z_mode) < 0.8) & (z_conf > 0.09)
            elif sourceargs['pofz_type'] == "mizuki":
                z_median, red_chi2, z_var = for_pofz_flag[i]
                pofzflag[i] = (red_chi2<5) & (z_var<0.45)
            elif sourceargs['pofz_type'] == "nnpz":
                z_median, photoz_flag, is_clean = for_pofz_flag[i]
                pofzflag[i] = (is_clean==1) & (photoz_flag<=1)
            elif sourceargs['pofz_type'] == "mlz":
                z_median, stddev, conf_mean = for_pofz_flag[i]
                pofzflag[i] = (stddev<3) & (conf_mean>0.13)

        pofz = sourceargs['fitspofz'][1]['P(z)'][:]
        pofz = pofz.reshape((pofz.shape[0], -1,))

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        status = 0
        
    if sourceargs['type'] == "hsc-wide" and sourceargs['filetype'] == "fits" and  not sourceargs["fullpofz"]:

        itern = sourceargs['iter']
        if itern == 0:
            #sourceargs['fits_list'] = glob.glob("/home/surhud/DataStore/S15B_v2.1/*_tracts/*_pz.fits")
            sourceargs['fits_list'] = glob.glob("/home/surhud/DataStore/S15B_v2.1/WIDE12H_tracts/*_pz.fits")
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0, 0

        sourceargs['fits'] = fitsio.FITS(sourceargs['fits_list'][itern].split("_pz")[0] + ".fits" )
        sourceargs['fitspofz'] = fitsio.FITS(sourceargs['fits_list'][itern])
        ntot = sourceargs['fits'][1].read_header()['naxis2']
        datagalflag = np.zeros(ntot)==np.zeros(ntot)
        pofzflag = np.zeros(ntot)==np.zeros(ntot)

        #wheremask = sourceargs['fits'][1].where()
        datagal = sourceargs['fits'][1]['ira','idec','ishape_hsm_regauss_e1','ishape_hsm_regauss_e2', 'ishape_hsm_regauss_sigma'][:]
        for_gal_flag = sourceargs['fits'][1]['imag_cmodel', 'iflux_cmodel', 'iflux_cmodel_err', 'ishape_hsm_regauss_resolution', 'ishape_hsm_regauss_e1', 'ishape_hsm_regauss_e2', 'ishape_hsm_regauss_sigma'][:]
        
        for i in range(ntot):
            imag_cmodel, iflux_cmodel, iflux_cmodel_err, ishape_hsm_regauss_resolution, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_sigma = for_gal_flag[i]
            datagalflag[i] = (imag_cmodel < 25.0) & (iflux_cmodel/iflux_cmodel_err >=10) & (ishape_hsm_regauss_resolution> 0.3) & (ishape_hsm_regauss_e1**2 + ishape_hsm_regauss_e2**2 < 4) & (ishape_hsm_regauss_sigma<0.4)

        pofz = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']) ][:]
        if sourceargs['pofz_type'] == "demp":
            for_pofz_flag = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']), "%s_photoz_mode" % (sourceargs['pofz_type']), "%s_photoz_conf_median" % (sourceargs['pofz_type']) ][:]
        elif sourceargs['pofz_type'] == "mizuki":
            for_pofz_flag = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']), "%s_reduced_chisq" % (sourceargs['pofz_type']), "%s_photoz_variance" % (sourceargs['pofz_type']) ][:]
        elif sourceargs['pofz_type'] == "nnpz":
            for_pofz_flag = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']), "%s_photoz_flag" % (sourceargs['pofz_type']), "%s_is_clean" % (sourceargs['pofz_type']) ][:]
        elif sourceargs['pofz_type'] == "mlz":
            for_pofz_flag = sourceargs['fitspofz'][1]["%s_photoz_median" % (sourceargs['pofz_type']), "%s_photoz_stddev_mean" % (sourceargs['pofz_type']), "%s_photoz_conf_mean" % (sourceargs['pofz_type']) ][:]
        else:
            print "Photoz not supported yet!"
            exit(11)

        for i in range(ntot):
            if sourceargs["pofz_type"] == "demp":
                z_median, z_mode, z_conf = for_pofz_flag[i]
                pofzflag[i] = (np.absolute(z_mode - z_median)/(1+z_mode) < 0.8) & (z_conf > 0.09)
            elif sourceargs['pofz_type'] == "mizuki":
                z_median, red_chi2, z_var = for_pofz_flag[i]
                pofzflag[i] = (red_chi2<5) & (z_var<0.45)
            elif sourceargs['pofz_type'] == "nnpz":
                z_median, photoz_flag, is_clean = for_pofz_flag[i]
                pofzflag[i] = (is_clean==1) & (photoz_flag<=1)
            elif sourceargs['pofz_type'] == "mlz":
                z_median, stddev, conf_mean = for_pofz_flag[i]
                pofzflag[i] = (stddev<3) & (conf_mean>0.13)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        status = 0

    if sourceargs['type'] == "hsc-wide-josh" and sourceargs['filetype'] == "fits":
        # Glob up all the tracts
        if sourceargs["fullpofz"]:
            print "P(z) not implemented for Josh's catalog"
            exit(11)
        else:
            pass

        itern = sourceargs['iter']
        if itern == 0:
            #sourceargs['fits'] = fitsio.FITS("/home/surhud/DataStore/RF/hsc_s15b_wl_catalog_v0.fits")
            sourceargs['fits'] = fitsio.FITS("/home/surhud/DataStore/rf_wide12h_hsc_s15b_wl_catalog_v0.fits")
            sourceargs['nrows'] = sourceargs['fits'][1].read_header()['naxis2']

        datagal = 0
        status = (itern*chunksize>=sourceargs['nrows'])
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz

        wheremask = sourceargs['fits'][1].where("#row>%d && #row<=%d" % (itern*chunksize, (itern+1)*chunksize))
        datagal = sourceargs['fits'][1]['ra','dec','e1','e2','weight'][wheremask]
        if sourceargs["fullpofz"]:
            print "Full p(z) not supported yet"
        else:
            pofz = sourceargs['fits'][1]['z_med'][wheremask]

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
            
    return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag
