import numpy as np
import pyfits
import pandas
import sys
import fitsio
import glob

def lens_select(lensargs):
    if lensargs['type'] == "dr12-lowz":
        hdulist = pyfits.open("DataStore/DR12/galaxy_DR12v5_LOWZ_North.fits.gz")
        data = hdulist[1].data
        ra_n = data["ra"]
        dec_n = data["dec"]
        zred_n = data["z"].astype('float64')
        wt_n = data["weight_systot"].astype("float64")
        idx = (zred_n>lensargs["zmin"]) & (zred_n<=lensargs["zmax"])
        ra_n = ra_n[idx]
        dec_n = dec_n[idx]
        zred_n = zred_n[idx]
        wt_n = wt_n[idx]

        hdulist = pyfits.open("DataStore/DR12/galaxy_DR12v5_LOWZ_South.fits.gz")
        data = hdulist[1].data
        ra_s = data["ra"]
        dec_s = data["dec"]
        zred_s = data["z"].astype('float64')
        wt_s = data["weight_systot"].astype("float64")
        idx = (zred_s>lensargs["zmin"]) & (zred_s<=lensargs["zmax"])
        ra_s = ra_s[idx]
        dec_s = dec_s[idx]
        zred_s = zred_s[idx]
        wt_s = wt_s[idx]

        ra = np.append(ra_n, ra_s)
        dec = np.append(dec_n, dec_s)
        zred = np.append(zred_n, zred_s)
        wt = np.append(wt_n, wt_s)
        
        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "dr12-cmass":
        hdulist = pyfits.open("DataStore/DR12/galaxy_DR12v5_CMASS_North.fits.gz")
        data = hdulist[1].data
        ra_n = data["ra"]
        dec_n = data["dec"]
        zred_n = data["z"].astype('float64')
        wt_n = data["weight_systot"].astype("float64")
        idx = (zred_n>lensargs["zmin"]) & (zred_n<=lensargs["zmax"])
        ra_n = ra_n[idx]
        dec_n = dec_n[idx]
        zred_n = zred_n[idx]
        wt_n = wt_n[idx]

        hdulist = pyfits.open("DataStore/DR12/galaxy_DR12v5_CMASS_South.fits.gz")
        data = hdulist[1].data
        ra_s = data["ra"]
        dec_s = data["dec"]
        zred_s = data["z"].astype('float64')
        wt_s = data["weight_systot"].astype("float64")
        idx = (zred_s>lensargs["zmin"]) & (zred_s<=lensargs["zmax"])
        ra_s = ra_s[idx]
        dec_s = dec_s[idx]
        zred_s = zred_s[idx]
        wt_s = wt_s[idx]

        ra = np.append(ra_n, ra_s)
        dec = np.append(dec_n, dec_s)
        zred = np.append(zred_n, zred_s)
        wt = np.append(wt_n, wt_s)
        
        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "alexie-test":
        ra, dec, zred, wt = np.loadtxt("DataStore/alexie-test/bossdr12_GAMA15_zbin1.txt", unpack=1)
        wt = ra/ra
        
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
        return ra, dec, zred, wt

    if lensargs['type'] == "redmapper-smallR":
        hdulist = pyfits.open("DataStore/redmapper-smallR/dr8_run_redmapper_v5.10-0.75_lgt5_prunedcatalog.fit")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z_lambda"].astype('float64')
        lamda = data["lambda_chisq"]
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

    if lensargs['type'] == "redmapper-largeR":
        hdulist = pyfits.open("DataStore/redmapper-largeR/dr8_run_redmapper_v5.10-1.25_lgt5_prunedcatalog.fit")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z_lambda"].astype('float64')
        lamda = data["lambda_chisq"]
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

    if lensargs['type'] == "redmapper":
        hdulist = pyfits.open("DataStore/redmapper/redmapper_dr8_public_v5.10_catalog.fits")
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
        ra, dec, zred, lamda, wt = np.loadtxt("DataStore/redmapper/redmapper_public_v5.10_randoms_%05d.dat" % (lensargs["rannum"]), unpack=1) 
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lamda>lensargs["lammin"]) & (lamda<=lensargs["lammax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        lamda = lamda[idx]
        wt = wt[idx]
        sys.stdout.write("Selecting %d randoms \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]

    if lensargs['type'] == "sdss-vagc-lbg":
        ra, dec, zred, mstel, mag = np.loadtxt("DataStore/LBG/Isolated_v0.dat", unpack=1)
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
        ra, dec, zred, lam = np.loadtxt("DataStore/camira/%s/%s/camira_catalog.dat" % (lensargs['hsc-release'], lensargs['version']), unpack=1, usecols=(0, 1, 2, 3))
        wt = lam/lam

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (lam>lensargs["lammin"]) & (lam<=lensargs["lammax"])
        sys.stdout.write("Selecting %d lenses \n" % (np.sum(idx)))
        return ra[idx], dec[idx], zred[idx], wt[idx]


def source_select(sourceargs, chunksize):

    if sourceargs['type'] == "hsc-wide-s16a_v2.0" and sourceargs['filetype'] == "fits" and not sourceargs["fullpofz"]:

        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("DataStore/S16A_v2.0/*_tracts/*_pz.fits")
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pz")[0] + ".fits")
        data = hdulist[1].data
        ira = data["ira"]
        idec = data["idec"]

        ishape_hsm_regauss_e1 = data["ishape_hsm_regauss_e1"]
        ishape_hsm_regauss_e2 = data["ishape_hsm_regauss_e2"]
        ishape_hsm_regauss_derived_bias_c1 = data["ishape_hsm_regauss_derived_bias_c1"]
        ishape_hsm_regauss_derived_bias_c2 = data["ishape_hsm_regauss_derived_bias_c2"]
        ishape_hsm_regauss_derived_bias_m = data["ishape_hsm_regauss_derived_bias_m"]
        ishape_hsm_regauss_derived_weight = data["ishape_hsm_regauss_derived_weight"] 
        ishape_hsm_regauss_derived_rms_e = data["ishape_hsm_regauss_derived_rms_e"] 
        zeroarr = ishape_hsm_regauss_derived_bias_c1 * 0.0
        datagalflag = data["weak_lensing_flag"]

        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_bias_m, ishape_hsm_regauss_derived_bias_c1, ishape_hsm_regauss_derived_bias_c2, zeroarr, zeroarr]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern])
        data = hdulist[1].data
        pofz = data["%s_photoz_best" % (sourceargs['pofz_type']) ]
        pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<0.1)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = pofz.size
        status = 0

    if sourceargs['type'] == "hsc-wide-s16a_v2.0" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("DataStore/S16A_v2.0/*_tracts/*_%s.fits" % (sourceargs['pofz_type']))
            sourceargs['fits_list_2'] = glob.glob("DataStore/S16A_v2.0/*_tracts/*_pz.fits")
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pz_pdf_%s" % (sourceargs['pofz_type']))[0] + ".fits")
        data = hdulist[1].data
        ira = data["ira"]
        idec = data["idec"]

        ishape_hsm_regauss_e1 = data["ishape_hsm_regauss_e1"]
        ishape_hsm_regauss_e2 = data["ishape_hsm_regauss_e2"]
        ishape_hsm_regauss_derived_bias_c1 = data["ishape_hsm_regauss_derived_bias_c1"]
        ishape_hsm_regauss_derived_bias_c2 = data["ishape_hsm_regauss_derived_bias_c2"]
        ishape_hsm_regauss_derived_bias_m = data["ishape_hsm_regauss_derived_bias_m"]
        ishape_hsm_regauss_derived_weight = data["ishape_hsm_regauss_derived_weight"] 
        ishape_hsm_regauss_derived_rms_e = data["ishape_hsm_regauss_derived_rms_e"] 
        zeroarr = ishape_hsm_regauss_derived_bias_c1 * 0.0
        datagalflag = data["weak_lensing_flag"]

        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_bias_m, ishape_hsm_regauss_derived_bias_c1, ishape_hsm_regauss_derived_bias_c2, zeroarr, zeroarr]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pdf_%s" % sourceargs['pofz_type'])[0] + ".fits")
        data = hdulist[1].data
        if sourceargs['pofz_type'] == "frankenz":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "mizuki":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "mlz":
            pofzflag = (datagalflag == True)
        elif sourceargs['pofz_type'] == "ephor":
            pofzflag = (datagalflag == True)

        sourceargs['fitspofz'] = fitsio.FITS(sourceargs['fits_list'][itern])

        pofz = sourceargs['fitspofz'][1]['P(z)'][:]
        pofz = pofz.reshape((pofz.shape[0], -1,))

        hdulist = pyfits.open(sourceargs['fits_list_2'][itern])
        data = hdulist[1].data
        pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<0.1)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        status = 0

    return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag

    '''
    if sourceargs['type'] == "hsc-wide" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:
        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("DataStore/S15B_v2.1/*_tracts/*_%s.fits" % (sourceargs['pofz_type']))
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pz_pdf_%s" % (sourceargs['pofz_type']))[0] + ".fits")
        data = hdulist[1].data
        imag_cmodel = data["imag_cmodel"]
        iflux_cmodel = data["iflux_cmodel"]
        iflux_cmodel_err = data["iflux_cmodel_err"]
        ishape_hsm_regauss_resolution = data["ishape_hsm_regauss_resolution"]
        ishape_hsm_regauss_e1 = data["ishape_hsm_regauss_e1"]
        ishape_hsm_regauss_e2 = data["ishape_hsm_regauss_e2"]
        ishape_hsm_regauss_sigma = data["ishape_hsm_regauss_sigma"] 
        ira = data["ira"]
        idec = data["idec"]

        datagalflag = (imag_cmodel < 25.0) & (iflux_cmodel/iflux_cmodel_err >=10) & (ishape_hsm_regauss_resolution> 1./3.) & (ishape_hsm_regauss_e1**2 + ishape_hsm_regauss_e2**2 < 4) & (ishape_hsm_regauss_sigma<0.4)
        #datagalflag = (imag_cmodel < 24.5) & (iflux_cmodel/iflux_cmodel_err >=10) & (iflux_cmodel/iflux_cmodel_err <80) & (ishape_hsm_regauss_resolution> 0.3) & (ishape_hsm_regauss_e1**2 + ishape_hsm_regauss_e2**2 < 4) & (ishape_hsm_regauss_sigma<0.4)
        #datagalflag = (imag_cmodel < 27.0)

        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_sigma]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pdf_%s" % sourceargs['pofz_type'])[0] + ".fits")
        data = hdulist[1].data
        if sourceargs['pofz_type'] == "demp":
            z_median = data["%s_photoz_median" % (sourceargs['pofz_type'])]
            z_mode = data["%s_photoz_mode" % (sourceargs['pofz_type'])]
            z_conf = data["%s_photoz_conf_median" % (sourceargs['pofz_type'])]
            pofzflag = (np.absolute(z_mode - z_median)/(1+z_mode) < 0.8) & (z_conf > 0.09)
        elif sourceargs['pofz_type'] == "mizuki":
            red_chi2 = data["%s_reduced_chisq" % (sourceargs['pofz_type'])]
            z_var = data["%s_photoz_variance" % (sourceargs['pofz_type'])]
            pofzflag = (red_chi2<5) & (z_var<0.45)
        elif sourceargs['pofz_type'] == "nnpz":
            photoz_flag = data["%s_photoz_flag" % (sourceargs['pofz_type'])]
            is_clean = data["%s_is_clean" % (sourceargs['pofz_type'])]
            pofzflag = (is_clean==1) & (photoz_flag<=1)
        elif sourceargs['pofz_type'] == "mlz":
            stddev = data["%s_photoz_stddev_mean" % (sourceargs['pofz_type'])]
            conf_mean = data["%s_photoz_conf_mean" % (sourceargs['pofz_type'])]
            pofzflag = (stddev<3) & (conf_mean>0.13)
        elif sourceargs['pofz_type'] == "frrfv0":
            pofzflag = (datagalflag == True)

        sourceargs['fitspofz'] = fitsio.FITS(sourceargs['fits_list'][itern])

        pofz = sourceargs['fitspofz'][1]['P(z)'][:]
        pofz = pofz.reshape((pofz.shape[0], -1,))

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        status = 0

    if sourceargs['type'] == "alexie-test" and sourceargs['filetype'] == "ascii" and  not sourceargs["fullpofz"]:

        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['dfchunks'] = pandas.read_csv("DataStore/alexie-test/GAMA15H_sources_test_original2.cat", names=(['ra','dec','z_b','jk','e1','e2','weight']), chunksize=chunksize, delim_whitespace=1)
        try:
            datagal = sourceargs['dfchunks'].next()
        except:
            status = 1
            return 0.0, sourceargs, 0.0, status, 0, 0, 0

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = datagal.ra.size
        status = 0

        pofz = datagal.z_b.values
        datagalflag = (datagal.ra.values==datagal.ra.values)
        pofzflag = (datagal.ra.values==datagal.ra.values)

    if sourceargs['type'] == "hsc-wide" and sourceargs['filetype'] == "fits" and  not sourceargs["fullpofz"]:

        itern = sourceargs['iter']
        if itern == 0:
            sourceargs['fits_list'] = glob.glob("DataStore/S15B_v2.1/*_tracts/*_pz.fits")
            #sourceargs['fits_list'] = glob.glob("/home/surhud/DataStore/S15B_v2.1/WIDE12H_tracts/*_pz.fits")
            
        if itern == len(sourceargs['fits_list']):
            status = 1

        datagal = 0
        status = not (itern<len(sourceargs['fits_list']))
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0, 0

        # First apply the flags
        hdulist = pyfits.open(sourceargs['fits_list'][itern].split("_pz")[0] + ".fits" )
        data = hdulist[1].data
        imag_cmodel = data["imag_cmodel"]
        iflux_cmodel = data["iflux_cmodel"]
        iflux_cmodel_err = data["iflux_cmodel_err"]
        ishape_hsm_regauss_resolution = data["ishape_hsm_regauss_resolution"]
        ishape_hsm_regauss_e1 = data["ishape_hsm_regauss_e1"]
        ishape_hsm_regauss_e2 = data["ishape_hsm_regauss_e2"]
        ishape_hsm_regauss_sigma = data["ishape_hsm_regauss_sigma"] 
        ira = data["ira"]
        idec = data["idec"]

        datagalflag = (imag_cmodel < 25.0) & (iflux_cmodel/iflux_cmodel_err >=10) & (ishape_hsm_regauss_resolution> 0.3) & (ishape_hsm_regauss_e1**2 + ishape_hsm_regauss_e2**2 < 4) & (ishape_hsm_regauss_sigma<0.4)
        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_sigma]).T

        hdulist = pyfits.open(sourceargs['fits_list'][itern])
        data = hdulist[1].data
        z_median = data["%s_photoz_median" % (sourceargs['pofz_type'])]
        if sourceargs['pofz_type'] == "demp":
            z_mode = data["%s_photoz_mode" % (sourceargs['pofz_type'])]
            z_conf = data["%s_photoz_conf_median" % (sourceargs['pofz_type'])]
            pofzflag = (np.absolute(z_mode - z_median)/(1+z_mode) < 0.8) & (z_conf > 0.09)
        elif sourceargs['pofz_type'] == "mizuki":
            red_chi2 = data["%s_reduced_chisq" % (sourceargs['pofz_type'])]
            z_var = data["%s_photoz_variance" % (sourceargs['pofz_type'])]
            pofzflag = (red_chi2<5) & (z_var<0.45)
        elif sourceargs['pofz_type'] == "nnpz":
            photoz_flag = data["%s_photoz_flag" % (sourceargs['pofz_type'])]
            is_clean = data["%s_is_clean" % (sourceargs['pofz_type'])]
            pofzflag = (is_clean==1) & (photoz_flag<=1)
        elif sourceargs['pofz_type'] == "mlz":
            stddev = data["%s_photoz_stddev_mean" % (sourceargs['pofz_type'])]
            conf_mean = data["%s_photoz_conf_mean" % (sourceargs['pofz_type'])]
            pofzflag = (stddev<3) & (conf_mean>0.13)

        pofz = data["%s_photoz_median" % (sourceargs['pofz_type']) ]

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = pofz.size
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
            sourceargs['fits'] = fitsio.FITS("DataStore/RF/hsc_s15b_wl_catalog_v0.fits")
            #sourceargs['fits'] = fitsio.FITS("/home/surhud/DataStore/rf_wide12h_hsc_s15b_wl_catalog_v0.fits")
            sourceargs['nrows'] = sourceargs['fits'][1].read_header()['naxis2']

        datagal = 0
        status = (itern*chunksize>=sourceargs['nrows'])
        Ngal = 0
        pofz = 0
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0, 0

        wheremask = sourceargs['fits'][1].where("#row>%d && #row<=%d" % (itern*chunksize, (itern+1)*chunksize))
        datagal = sourceargs['fits'][1]['ra','dec','e1','e2','weight'][wheremask]
        if sourceargs["fullpofz"]:
            print "Full p(z) not supported yet"
        else:
            pofz = sourceargs['fits'][1]['z_med'][wheremask]

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = np.shape(datagal)[0]
        forflag = np.arange(Ngal)
        datagalflag = (forflag==forflag)
        pofzflag = (forflag==forflag)
            
    '''

def get_pofz_array(sourceargs):
    hdulist = pyfits.open("DataStore/S16A_v2.0/pz_pdf_bins_%s.fits" % sourceargs['pofz_type'])
    data = hdulist[1].data
    return data["bins"]
