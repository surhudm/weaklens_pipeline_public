import numpy as np
from astropy.io import fits as pyfits
import pandas
import sys
import fitsio
import glob
from sqlitedict import SqliteDict
import gnupg
import os

gpg = gnupg.GPG()

def FITS_to_PGP(message):
    """
    Turns a string stored into an FITS comment back into a proper
    PGP message
    """
    s = "-----BEGIN PGP MESSAGE-----\n\n"
    s += message
    s += "\n-----END PGP MESSAGE-----\n"
    return s

# This is the function which you will use to decrypt the string
def decrypt_string(encrypted_string):
    #string = getpass.getpass(prompt="Enter passphrase to unlock the private key:\n")
    decrypted_data = gpg.decrypt(encrypted_string) #, passphrase=string)

    print('ok: ', decrypted_data.ok)
    print('status: ', decrypted_data.status)
    #print 'stderr: ', decrypted_data.stderr
    #print 'Decrypted string (additive m value for the catalog): ', decrypted_data.data
    return decrypted_data.data

def lens_select(lensargs):
    """
    lens_select(dictionary lensargs)

    Selects lens sample

    :Parameters:
    
    -    lensargs : Dictionary read from the config file which contains various
         cuts on the lens sample

    :Returns:
    
    -    ra : Right ascension in degrees
    -    dec : Declination in degrees
    -    zred : Redshift
    -    weight : Weights for the lenses

    :Adding new samples:
    
    Add new samples by copy pasting one of the lensargs block and then modifying them according to your will
    """
    if lensargs["type"] == "generic":
        hdulist = pyfits.open(lensargs["fname"])
        data = hdulist[1].data
        ra = data["ra"].astype("float64")
        dec = data["dec"].astype("float64")
        z = data["z"].astype("float64")
        wt = data["wt"].astype("float64")
         
        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, z, wt

    # Use generic text file with the columns: 
    # ra (deg), dec (deg), z, lens_wt
    # Any lines beginning with "#" are ignored
    if lensargs["type"] == "generic-text":
        data = pandas.read_csv(lensargs["fname"], delim_whitespace=1, header=None, names=(["ra", "dec", "z", "wt"]), comment="#")
        ra = data["ra"].values.astype("float64")
        dec = data["dec"].values.astype("float64")
        z = data["z"].values.astype("float64")
        wt = data["wt"].values.astype("float64")
         
        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, z, wt
    if lensargs["type"] == "gama":
        hdulist = pyfits.open("./DataStore/gama/G3CFOFGroup.fits")
        data = hdulist[1].data
        ra = data["IterCenRA"].astype("float64") #ra of galaxies in the lense
        dec = data["IterCenDec"].astype("float64") #dec of galaxies in the lense
        z = data["IterCenZ"].astype("float64") #redshift of galaxies in the lense
        vdisp=data["VelDisp"].astype("float64") 
        idx = (z>lensargs["zmin"]) & (z<=lensargs["zmax"]) & (vdisp>lensargs["vdispmin"]) & (vdisp<=lensargs["vdispmax"])
        ra = ra[idx]
        dec = dec[idx]
        z = z[idx]
        vdisp = vdisp[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, z, wt

    if lensargs["type"] == "mgii-absorbers-ran":
        hdulist = pyfits.open("DataStore/INFO_140303_Ref_Mg_II_DR79_match.fits")
        data = hdulist[1].data
        ra = data["ra"].astype("float64")
        dec = data["dec"].astype("float64")
        zabs = data["zabs"].astype("float64")
         
        idx = (zabs>lensargs["zmin"]) & (zabs<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zabs = zabs[idx]

        wt = ra/ra

        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, zabs, wt

    if lensargs["type"] == "mgii-absorbers":
        hdulist = pyfits.open("DataStore/INFO_140303_Mg_II_DR79_match.fits")
        data = hdulist[1].data
        ra = data["ra"].astype("float64")
        dec = data["dec"].astype("float64")
        zabs = data["zabs"].astype("float64")
         
        idx = (zabs>lensargs["zmin"]) & (zabs<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zabs = zabs[idx]

        wt = ra/ra

        sys.stdout.write("Selecting %d samples \n" % (ra.size))

        return ra, dec, zabs, wt

    if lensargs['type'] == "hsc-cmass-subsample-random":
        df = pandas.read_csv("DataStore/HSC_CMASS/CMASS_K+E_v0.99.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values
        wt = (df.weight_noz.values + df.weight_cp.values - 1.0)* df.weight_star.values
        iMag = df.Mag_Wake.values
        idx = (iMag>lensargs["iMagmin"]) & (iMag<=lensargs["iMagmax"]) & (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])

        # First compute the selection function as a function of redshift
        idxall = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        zarr = np.linspace(lensargs["zmin"], lensargs["zmax"], 20)
        selarr = np.zeros(zarr.size-1)
        for i in range(1, zarr.size):
            selarr[i-1] = np.sum(idx & (zred>zarr[i-1]) & (zred<=zarr[i]))*1.0/np.sum(idxall & (zred>zarr[i-1]) & (zred<=zarr[i]))
            print(i, selarr[i-1], "\n=================")
        zarr = zarr[:-1]/2. + zarr[1:]/2.

        from scipy.interpolate import UnivariateSpline
        selspl = UnivariateSpline(zarr, selarr, bbox=(lensargs["zmin"], lensargs["zmax"]), ext=1)

        with SqliteDict("DataStore/HSC_CMASS/randoms/ran%03d.sqlite3" % (lensargs["rannum"])) as data:
            ra = data["ra"]
            dec = data["dec"]
            zred = data["z"].astype('float64')
            wt = data["weight"].astype("float64")

            idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
            ra = ra[idx]
            dec = dec[idx]
            zred = zred[idx]
            wt = wt[idx]

            np.random.seed(lensargs["rannum"])
            rannum = np.random.random(size=ra.size)
            idx = (rannum<selspl(zred))
            ra = ra[idx]
            dec = dec[idx]
            zred = zred[idx]
            wt = wt[idx]

        sys.stdout.write("Selecting %d samples \n" % (ra.size))
        #jackreg = getregions(ra, dec, lensargs["jackregfile"])
        return ra, dec, zred, wt #, jackreg

    if lensargs['type'] == "hsc-lowz-subsample":
        df = pandas.read_csv("DataStore/HSC_LOWZ/LOWZ_K+E_v0.99.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values
        wt = (df.weight_noz.values + df.weight_cp.values - 1.0)* df.weight_star.values
        iMag = df.Mag_Wake.values
        idx = (iMag>lensargs["iMagmin"]) & (iMag<=lensargs["iMagmax"]) & (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])

        if lensargs['selection_write'] == 1:
            idxall = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
            zarr = np.linspace(lensargs["zmin"], lensargs["zmax"], 20)
            selarr = np.zeros(zarr.size-1)
            for i in range(1, zarr.size):
                selarr[i-1] = np.sum(idx & (zred>zarr[i-1]) & (zred<=zarr[i]))*1.0/np.sum(idxall & (zred>zarr[i-1]) & (zred<=zarr[i]))
                print(i, selarr[i-1], "\n=================")
            zarr = zarr[:-1]/2. + zarr[1:]/2.
            dfsel = pandas.DataFrame(zarr, columns=(["zred"]))
            dfsel["sel"] = selarr
            dfsel.to_csv(lensargs["selection"], index=False, sep=" ")

        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        iMag = iMag[idx]

        sys.stdout.write("Selecting %d samples \n" % (ra.size))
        #jackreg = getregions(ra, dec, lensargs["jackregfile"])
        return ra, dec, zred, wt #, jackreg

    if lensargs['type'] == "hsc-cmass-subsample":
        df = pandas.read_csv("DataStore/HSC_CMASS/CMASS_K+E_v0.99.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values
        wt = (df.weight_noz.values + df.weight_cp.values - 1.0)* df.weight_star.values
        iMag = df.Mag_Wake.values
        idx = (iMag>lensargs["iMagmin"]) & (iMag<=lensargs["iMagmax"]) & (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])

        if lensargs['selection_write'] == 1:
            idxall = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
            zarr = np.linspace(lensargs["zmin"], lensargs["zmax"], 20)
            selarr = np.zeros(zarr.size-1)
            for i in range(1, zarr.size):
                selarr[i-1] = np.sum(idx & (zred>zarr[i-1]) & (zred<=zarr[i]))*1.0/np.sum(idxall & (zred>zarr[i-1]) & (zred<=zarr[i]))
                print(i, selarr[i-1], "\n=================")
            zarr = zarr[:-1]/2. + zarr[1:]/2.
            dfsel = pandas.DataFrame(zarr, columns=(["zred"]))
            dfsel["sel"] = selarr
            dfsel.to_csv(lensargs["selection"], index=False, sep=" ")

        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        iMag = iMag[idx]

        sys.stdout.write("Selecting %d samples \n" % (ra.size))
        #jackreg = getregions(ra, dec, lensargs["jackregfile"])
        return ra, dec, zred, wt #, jackreg

    if lensargs['type'] == "SDSS-QSOs-John":
        fname = "DataStore/SDSS_QSOs_zlt1.asc"
        dfreal = pandas.read_csv(fname, delim_whitespace=1, names=(["ra", "dec", "z"]), usecols=([1, 2, 7]))
        ra = dfreal["ra"].values
        dec = dfreal["dec"].values
        zred = dfreal["z"].values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "Masato-mock-forreq-random-full":
        rannum = lensargs['rannum']
        rotation = lensargs['rotation']
        realization = lensargs['realization']
        files = glob.glob("DataStore/Mocks_HOD_full/random_*")
        np.random.seed(rannum)

        for i, fname in enumerate(files):
            tractreg = fname.split("_")[-1]
            realfile = "DataStore/Mocks_HOD_full/r%03d_rotmat%d_%s" % (realization, rotation, tractreg)
            dfreal = pandas.read_csv(realfile, delim_whitespace=1, skiprows=1, names=(["ra", "dec", "z"]), usecols=([1, 2, 3]))
            Nreal = dfreal.ra.values.size

            rowstart = rannum * Nreal
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1+rowstart, names=(["ra", "dec", "z"]), usecols=([0, 1, 2]), nrows=Nreal)
            else:
                dfp = pandas.read_csv(fname, delim_whitespace=1, skiprows=1+rowstart, names=(["ra", "dec", "z"]), usecols=([0, 1, 2]), nrows=Nreal)
                df = df.append(dfp, ignore_index=1)

        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt


    if lensargs['type'] == "Masato-mock-forreq-random":
        rannum = lensargs['rannum']
        rotation = lensargs['rotation']
        realization = lensargs['realization']
        files = glob.glob("DataStore/Mocks_HOD/random_*")
        np.random.seed(rannum)

        for i, fname in enumerate(files):
            tractreg = fname.split("_")[-1]
            realfile = "DataStore/Mocks_HOD/r%03d_rotmat%d_%s" % (realization, rotation, tractreg)
            dfreal = pandas.read_csv(realfile, delim_whitespace=1, skiprows=1, names=(["ra", "dec", "z"]), usecols=([1, 2, 3]))
            Nreal = dfreal.ra.values.size

            rowstart = rannum * Nreal
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1+rowstart, names=(["ra", "dec", "z"]), usecols=([0, 1, 2]), nrows=Nreal)
            else:
                dfp = pandas.read_csv(fname, delim_whitespace=1, skiprows=1+rowstart, names=(["ra", "dec", "z"]), usecols=([0, 1, 2]), nrows=Nreal)
                df = df.append(dfp, ignore_index=1)

        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "Masato-mock-forreq-full":
        realization = lensargs['realization']
        rotation = lensargs['rotation']
        files = glob.glob("DataStore/Mocks_HOD_full/r%03d_rotmat%d_*" % (realization, rotation))
        for i, fname in enumerate(files):
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3]))
            else:
                df = df.append(pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3])), ignore_index=1)
        
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt


    if lensargs['type'] == "Masato-mock-forreq":
        realization = lensargs['realization']
        rotation = lensargs['rotation']
        files = glob.glob("DataStore/Mocks_HOD/r%03d_rotmat%d_*" % (realization, rotation))
        for i, fname in enumerate(files):
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3]))
            else:
                df = df.append(pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3])), ignore_index=1)
        
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "Masato-mock":
        realization = lensargs['realization']
        rotation = lensargs['rotation']
        sample = lensargs['sample']
        zmin = lensargs["zmin"]
        zmax = lensargs["zmax"]
        files = glob.glob("DataStore/hod_gal/%s_%.2f_%.2f/newID/r%03d_rotmat%d_*" % (sample, zmin, zmax, realization, rotation))
        for i, fname in enumerate(files):
            if i==0:
                df = pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3]))
            else:
                df = df.append(pandas.read_csv(fname, delim_whitespace=1, skiprows=1, names=(["id", "ra", "dec", "z"]), usecols=([0, 1, 2, 3])), ignore_index=1)
        
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values

        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = ra/ra

        sys.stdout.write("Selecting %d lenses \n" % (ra.size))
        return ra, dec, zred, wt

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

        if "pofzoutput" in lensargs:
            n, bins = np.histogram(zred, weights=wt, bins=np.linspace(np.min(zred), np.max(zred), 20), normed=1)
            np.savetxt(lensargs["pofzoutput"], np.transpose([bins[1:]/2+bins[:-1]/2, n]))
            exit(11)

        return ra, dec, zred, wt

    if lensargs['type'] == "cfht-cmass-reproduce":
        hdulist = pyfits.open("DataStore/HSC_CMASS/cmass-dr11v1-A-Anderson.dat_wcollind_mstel_xdist.fits")
        data = hdulist[1].data
        ra = data["ra"]
        dec = data["dec"]
        zred = data["z"].astype('float64')
        wt_star = data["weight_star"].astype("float64")
        wt_noz = data["weight_noz"].astype("float64")
        wt_cp = data["weight_cp"].astype("float64")
        wt = wt_star*(wt_noz+wt_cp-1)
        mstar = data["logmass"]
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (mstar>lensargs["mstarmin"]) & (mstar<=lensargs["mstarmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
        
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

        if "pofzoutput" in lensargs:
            n, bins = np.histogram(zred, weights=wt, bins=np.linspace(np.min(zred), np.max(zred), 20), normed=1)
            np.savetxt(lensargs["pofzoutput"], np.transpose([bins[1:]/2+bins[:-1]/2, n]))
            exit(11)

        return ra, dec, zred, wt

    if lensargs['type'] == "voids":

        df = pandas.read_csv("DataStore/Voids/Void_combined_catalog.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.zred.values
        wt = ra/ra
        
        sys.stdout.write("Selecting %d voids \n" % (ra.size))
        #return ra[idx], dec[idx], zred[idx], wt[idx]
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

    if lensargs['type'] == "redmapper_satellite":
        hdulist = pyfits.open("DataStore/redmapper/redmapper_dr8_public_v6.3_members.fits")
        hdulist2 = pyfits.open("DataStore/redmapper/redmapper_dr8_public_v6.3_catalog.fits")
        data = hdulist[1].data
        data2 = hdulist2[1].data
        zlamda=data2["z_lambda"]
        cid=data2["id"]
        ldict={}
        for ii in range(cid.size):
            ldict[cid[ii]]=zlamda[ii]
        mid=data['id']
        ra = data["ra"]
        dec = data["dec"]        
        lamda = data2["lambda"]
        wt = ra/ra
        r = data["r"]
        zred = ra*0.0
        
        for jj in range(mid.size):
            zred[jj]=ldict[mid[jj]].astype('float64')        
            
        idx = (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"]) & (r>lensargs["rmin"]) & (r<=lensargs["rmax"])
        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        #lamda = lamda[idx]
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

        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]

        if "pofzoutput" in lensargs:
            n, bins = np.histogram(zred, weights=wt, bins=np.linspace(np.min(zred), np.max(zred), 20), normed=1)
            np.savetxt(lensargs["pofzoutput"], np.transpose([bins[1:]/2+bins[:-1]/2, n]))
            exit(11)
        return ra, dec, zred, wt

    if lensargs['type'] == "bossdr11-hsc-cosmology":
        if lensargs["sample"] == "LOWZ":
            df = pandas.read_csv("DataStore/BOSS_DR11/LOWZ_K+E_v0.99.dat", delim_whitespace=1)
        elif lensargs["sample"] == "CMASS":
            df = pandas.read_csv("DataStore/BOSS_DR11/CMASS_K+E_v0.99.dat", delim_whitespace=1)
        ra = df.ra.values
        dec = df.dec.values
        zred = df.z.values
        wt = (df.weight_noz.values + df.weight_cp.values - 1.0)* df.weight_star.values
        if 'mag_type' in lensargs:
            if lensargs['mag_type'] == 'Wake':
                iMag = df.Mag_Wake.values
            elif lensargs['mag_type'] == 'Blanton':
                iMag = df.Mag_Blanton.values
        else:
            iMag = df.Mag_Wake.values
        idx = (iMag>lensargs["iMagmin"]) & (iMag<=lensargs["iMagmax"]) & (zred>lensargs["zmin"]) & (zred<=lensargs["zmax"])

        # Selection on spatial distribution. Take galaxies approx delta = 13.5 deg from the edge of each field, which corresponds to ~100 Mpc/h at z=0.15 (the lowest redshift in our samples).
        delta = 13.5
        # HECTOMAP
        idx_s = (236.-delta/np.cos(dec/180.) < ra) & (ra < 249.+delta/np.cos(dec/180.)) & (42.-delta < dec) & (dec < 44.5+delta)
        # GAMA09H
        idx_s |= (128.-delta/np.cos(dec/180.) < ra) & (ra < 142.+delta/np.cos(dec/180.)) & (-2.-delta < dec) & (dec < 5.+delta)
        # WIDE12H
        idx_s |= (175.5-delta/np.cos(dec/180.) < ra) & (ra < 183.+delta/np.cos(dec/180.)) & (-1.5-delta < dec) & (dec < 1.5+delta)
        # GAMA15H
        idx_s |= (210.-delta/np.cos(dec/180.) < ra) & (ra < 226.+delta/np.cos(dec/180.)) & (-1.5-delta < dec) & (dec < 1.5+delta)
        # VVDS
        idx_s |= (330.-delta/np.cos(dec/180.) < ra) & (ra < 342.+delta/np.cos(dec/180.)) & (-1.5-delta < dec) & (dec < 2.+delta)
        # XMM
        idx_s |= (29.-delta/np.cos(dec/180.) < ra) & (ra < 40.+delta/np.cos(dec/180.)) & (-6.5-delta < dec) & (dec < -1.8+delta)

        idx &= idx_s

        ra = ra[idx]
        dec = dec[idx]
        zred = zred[idx]
        wt = wt[idx]
#        print "sum wt:", np.sum(wt), len(ra)
        iMag = iMag[idx]

        sys.stdout.write("Selecting %d samples \n" % (ra.size))
        #jackreg = getregions(ra, dec, lensargs["jackregfile"])
        #np.savetxt("selected_lens_%s.dat" % (lensargs["sample"]), np.array([ra, dec, zred, wt]).transpose())
        return ra, dec, zred, wt #, jackreg

    if lensargs['type'] == "Euclid_Flagship_LightCone":
        # Euclid Flagship mock: lens_bin_<N>.fits in lensargs['dir']
        # Columns: ra_gal, dec_gal, true_redshift_gal
        if 'dir' not in lensargs:
            raise RuntimeError("lens.dir is required for Euclid_Flagship_LightCone")
        lens_bin_ids = lensargs.get('lens_bin_id', [1, 2, 3, 4, 5])
        if isinstance(lens_bin_ids, int):
            lens_bin_ids = [lens_bin_ids]

        ra_list, dec_list, zred_list = [], [], []
        for bin_id in lens_bin_ids:
            filename = os.path.join(lensargs['dir'], 'lens_bin_%d.fits' % bin_id)
            d = pyfits.getdata(filename)
            ra_list.append(d['ra_gal'].astype('float64'))
            dec_list.append(d['dec_gal'].astype('float64'))
            zred_list.append(d['true_redshift_gal'].astype('float64'))

        ra = np.concatenate(ra_list)
        dec = np.concatenate(dec_list)
        zred = np.concatenate(zred_list)

        # Optional redshift cuts (bins may already be pre-filtered)
        if 'zmin' in lensargs:
            sel = zred >= lensargs['zmin']
            ra, dec, zred = ra[sel], dec[sel], zred[sel]
        if 'zmax' in lensargs:
            sel = zred < lensargs['zmax']
            ra, dec, zred = ra[sel], dec[sel], zred[sel]

        wt = np.ones(len(ra), dtype=float)
        sys.stdout.write("Number of lenses (Euclid Flagship LightCone): %s\n" % (ra.size))
        return ra, dec, zred, wt

    if lensargs['type'] == "Euclid_Flagship_LightCone_random":
        # Random points for Euclid_Flagship_LightCone: a single FITS file
        # (lensargs['dir'] is the file path) with the same column names as the
        # lens catalog. The randoms should cover the survey footprint and follow
        # the lens n(z); select the matching lens bin with zmin/zmax.
        if 'dir' not in lensargs:
            raise RuntimeError("lens.dir is required for Euclid_Flagship_LightCone_random")
        d = pyfits.getdata(lensargs['dir'])
        ra = d['ra_gal'].astype('float64')
        dec = d['dec_gal'].astype('float64')
        zred = d['true_redshift_gal'].astype('float64')

        if 'zmin' in lensargs:
            sel = zred >= lensargs['zmin']
            ra, dec, zred = ra[sel], dec[sel], zred[sel]
        if 'zmax' in lensargs:
            sel = zred < lensargs['zmax']
            ra, dec, zred = ra[sel], dec[sel], zred[sel]

        wt = np.ones(len(ra), dtype=float)
        sys.stdout.write("Number of random lenses (Euclid Flagship LightCone): %s\n" % (ra.size))
        return ra, dec, zred, wt

    raise RuntimeError("lens type %s is not defined in lens_select()" % lensargs['type'])


def source_select(sourceargs, chunksize):
    """
    source_select(dictionary lensargs)

    Selects source sample

    :Parameters:
    
    -    sourceargs : Dictionary read from the config file which contains various
         cuts on the source sample
    -    chunksize : The chunksize of sources to be read (only applicable for some ascii readers)

    :Returns:     return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag
    
    -    datagal : array consisting of ra, dec, e1, e2, weight, erms, m, c1, c2, a1, a2
    -    sourceargs : Dictionary read from the config file which contains various
         cuts on the source sample
    -    Ngal : Number of galaxies
    -    status : status of the read
    -    pofz : Record arrays read from the P(z) file for `Ngal` galaxies
    -    datagalflag: Flag array to indicate whether to use the galaxy or not
    -    pofzflag: Flag array to indicate whether to use the photoz of the galaxy or not

    """
    if sourceargs['type'] == "hsc-wide-s19a_v2.0" and sourceargs['filetype'] == "fits" and not sourceargs["fullpofz"]:

        # This variable iterates over the different fields in the catalog
        itern = sourceargs['iter']
        if itern == 0:
            if not sourceargs["blinded"]:
                sys.stderr.write("All source catalogs in S19A_v2.0 are blinded. Aborting\n")
                exit(11)
            else:
                sourceargs["df"] = pandas.read_csv("DataStore/S19A_v2.0/Mapping_from_blinded_catalogs_to_shape_catalogs.dat", delim_whitespace=1)
            if "usefield" in sourceargs:
                restidx = sourceargs["df"].field.values==sourceargs["usefield"]
                sourceargs["df"] = sourceargs["df"][restidx]
                print("Using only field", sourceargs["usefield"])

        datagal = 0
        Ngal = 0
        pofz = 0
        status = not (itern<(sourceargs['df'].field.values.size))
        if status:
            print("Uh oh, problem here")
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        field = sourceargs["df"].field.values[itern]
        tract = sourceargs["df"].tract.values[itern]
        begin = sourceargs["df"].begin.values[itern]
        end = sourceargs["df"].end.values[itern]
        
        # Initialize list of tracts within a given field
        if begin == 0:
            # Define the shape catalog file from Rachel to read from for this field
            if sourceargs["blinded"]:
                # Decrypt the Delta m value, if sourceargs["dm"] not set already
                fname = "DataStore/Blinded_S19A_v2.0/%s_blinded_%d_m.fits" % (field, sourceargs["blindnumber"])
                if "dm" not in sourceargs:
                    try:
                        hdulist = pyfits.open(fname, memmap=True)
                        msg = hdulist[0].header["DM1"]
                        sourceargs['dm'] = float(decrypt_string(FITS_to_PGP(msg)))
                    except:
                        print("GPG decryption failed, check your gpg-agent")
                        exit(12)
                # The commented line is memory intensive
                #sourceargs["fits"] = fitsio.FITS("%s[1][col *, i_hsmshaperegauss_derived_shear_bias_m = i_hsmshaperegauss_derived_shear_bias_m - %f ]" % (fname, sourceargs['dm']))[1]
                # Instead read the file as it is and then subtract off later
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]

        # Read in the photoz file for the current field, tract
        try:
            hdulist = pyfits.open("DataStore/S19A_v2.0/%s_tracts/%d_pz.fits" % (field, tract))
            data = hdulist[1].data
            pofz = data["%s_photoz%s" % (sourceargs['pofz_type'], sourceargs['photoz_estimate']) ]
            pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<sourceargs['photoz_risk_best_cut'])
        except:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        # Read in the shape values for the current field, tract, the row filtering specification indexes rows from 1 - Nrecords
        sourceargs["fits"] = pyfits.open("DataStore/S19A_v2.0/%s_tracts/%d_no_m.fits" % (field, tract))[1].data
        wheremask = sourceargs['fits_mblind'].where("#row>%d && #row<=%d" % (begin, end))

        #datagal = sourceargs['fits'][1]["ira", "idec", "i_hsmshaperegauss_e1", "i_hsmshaperegauss_e2", "i_hsmshaperegauss_derived_shear_bias_c1", "i_hsmshaperegauss_derived_shear_bias_c2", "i_hsmshaperegauss_derived_shear_bias_m", "i_hsmshaperegauss_derived_shape_weight", "i_hsmshaperegauss_derived_rms_e"][wheremask]
        #datagal = sourceargs['fits']["ira", "idec", "i_hsmshaperegauss_e1", "i_hsmshaperegauss_e2", "i_hsmshaperegauss_derived_shape_weight", "i_hsmshaperegauss_derived_rms_e", "i_hsmshaperegauss_derived_shear_bias_m", "i_hsmshaperegauss_derived_shear_bias_c1", "i_hsmshaperegauss_derived_shear_bias_c2"][wheremask]
        ira = sourceargs['fits']["i_ra"]
        idec = sourceargs['fits']["i_dec"]
        ishape_hsm_regauss_e1 = sourceargs['fits']["i_hsmshaperegauss_e1"]
        ishape_hsm_regauss_e2 = sourceargs['fits']["i_hsmshaperegauss_e2"]
        ishape_hsm_regauss_derived_shape_weight = sourceargs['fits']["i_hsmshaperegauss_derived_weight"]
        ishape_hsm_regauss_derived_rms_e = sourceargs['fits']["i_hsmshaperegauss_derived_rms_e"]
        ishape_hsm_regauss_derived_shear_bias_m = sourceargs['fits_mblind']["i_hsmshaperegauss_derived_shear_bias_m"][wheremask] - sourceargs["dm"]
        ishape_hsm_regauss_derived_shear_bias_c1 = sourceargs['fits']["i_hsmshaperegauss_derived_shear_bias_c1"]
        ishape_hsm_regauss_derived_shear_bias_c2 = sourceargs['fits']["i_hsmshaperegauss_derived_shear_bias_c2"]
        ishape_hsm_regauss_resolution = sourceargs['fits']["i_hsmshaperegauss_resolution"]
        ishape_psf_shape11 = sourceargs['fits']["i_sdssshape_psf_shape11"]
        ishape_psf_shape22 = sourceargs['fits']["i_sdssshape_psf_shape22"]
        ishape_psf_shape12 = sourceargs['fits']["i_sdssshape_psf_shape12"]
        ishape_psf_e1 = (ishape_psf_shape11-ishape_psf_shape22)/(ishape_psf_shape11+ishape_psf_shape22)
        ishape_psf_e2 = 2*ishape_psf_shape12/(ishape_psf_shape11+ishape_psf_shape22)
        i_apertureflux_10_mag = sourceargs['fits']["i_apertureflux_10_mag"]
        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_shape_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_shear_bias_m, ishape_hsm_regauss_derived_shear_bias_c1, ishape_hsm_regauss_derived_shear_bias_c2, ishape_hsm_regauss_resolution, i_apertureflux_10_mag, ishape_psf_e1, ishape_psf_e2]).T

        datagalflag = sourceargs['fits']["weak_lensing_flag"]

        sourceargs['iter'] = sourceargs['iter'] + 1

        Ngal = pofzflag.size
        status = 0
        
        return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag

    if sourceargs['type'] == "hsc-wide-s19a_v2.0" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:

        # This variable iterates over the different fields in the catalog
        itern = sourceargs['iter']
        if itern == 0:
            if not sourceargs["blinded"]:
                sys.stderr.write("All source catalogs in S19A_v2.0 are blinded. Aborting\n")
                exit(11)
            else:
                sourceargs["df"] = pandas.read_csv("DataStore/S19A_v2.0/Mapping_from_blinded_catalogs_to_shape_catalogs.dat", delim_whitespace=1)
            if "usefield" in sourceargs:
                restidx = sourceargs["df"].field.values==sourceargs["usefield"]
                sourceargs["df"] = sourceargs["df"][restidx]
                print("Using only field", sourceargs["usefield"])

        datagal = 0
        Ngal = 0
        pofz = 0
        status = not (itern<(sourceargs['df'].field.values.size))
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        field = sourceargs["df"].field.values[itern]
        tract = sourceargs["df"].tract.values[itern]
        begin = sourceargs["df"].begin.values[itern]
        end = sourceargs["df"].end.values[itern]
        
        # Initialize list of tracts within a given field
        if begin == 0:
            # Define the shape catalog file from Rachel to read from for this field
            if sourceargs["blinded"]:
                # Decrypt the Delta m value, if sourceargs["dm"] not set already
                fname = "DataStore/Blinded_S19A_v2.0/%s_blinded_%d_m.fits" % (field, sourceargs["blindnumber"])
                if "dm" not in sourceargs:
                    try:
                        hdulist = pyfits.open(fname, memmap=True)
                        msg = hdulist[0].header["DM1"]
                        sourceargs['dm'] = float(decrypt_string(FITS_to_PGP(msg)))
                    except:
                        print("GPG decryption failed, check your gpg-agent")
                        exit(12)
                # The commented line is memory intensive
                #sourceargs["fits"] = fitsio.FITS("%s[1][col *, i_hsmshaperegauss_derived_shear_bias_m = i_hsmshaperegauss_derived_shear_bias_m - %f ]" % (fname, sourceargs['dm']))[1]
                # Instead read the file as it is and then subtract off later
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]

        # Read in the photoz file for the current field, tract
        try:
        #if 1:
            hdulist = pyfits.open("DataStore/S19A_v2.0/%s_tracts/%d_pz.fits" % (field, tract))
            data = hdulist[1].data
            pofz = data["%s_photoz_best" % (sourceargs['pofz_type']) ]
            pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<sourceargs['photoz_risk_best_cut'])
        except:
        #else:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        # Read in the shape values for the current field, tract
        sourceargs["fits"] = pyfits.open("DataStore/S19A_v2.0/%s_tracts/%d_no_m.fits" % (field, tract))[1].data
        wheremask = sourceargs['fits_mblind'].where("#row>%d && #row<=%d" % (begin, end))

        #datagal = sourceargs['fits']["ira", "idec", "i_hsmshaperegauss_e1", "i_hsmshaperegauss_e2", "i_hsmshaperegauss_derived_shear_bias_c1", "i_hsmshaperegauss_derived_shear_bias_c2", "i_hsmshaperegauss_derived_shape_weight", "i_hsmshaperegauss_derived_rms_e"][wheremask]
        #datagalflag = sourceargs['fits'][1]["weak_lensing_flag"][wheremask]

        ira = sourceargs['fits']["i_ra"]
        idec = sourceargs['fits']["i_dec"]
        ishape_hsm_regauss_e1 = sourceargs['fits']["i_hsmshaperegauss_e1"]
        ishape_hsm_regauss_e2 = sourceargs['fits']["i_hsmshaperegauss_e2"]
        ishape_hsm_regauss_derived_shape_weight = sourceargs['fits']["i_hsmshaperegauss_derived_weight"]
        ishape_hsm_regauss_derived_rms_e = sourceargs['fits']["i_hsmshaperegauss_derived_rms_e"]
        ishape_hsm_regauss_derived_shear_bias_m = sourceargs['fits_mblind']["i_hsmshaperegauss_derived_shear_bias_m"][wheremask] - sourceargs["dm"]
        ishape_hsm_regauss_derived_shear_bias_c1 = sourceargs['fits']["i_hsmshaperegauss_derived_shear_bias_c1"]
        ishape_hsm_regauss_derived_shear_bias_c2 = sourceargs['fits']["i_hsmshaperegauss_derived_shear_bias_c2"]
        ishape_hsm_regauss_resolution = sourceargs['fits']["i_hsmshaperegauss_resolution"]
        ishape_psf_shape11 = sourceargs['fits']["i_sdssshape_psf_shape11"]
        ishape_psf_shape22 = sourceargs['fits']["i_sdssshape_psf_shape22"]
        ishape_psf_shape12 = sourceargs['fits']["i_sdssshape_psf_shape12"]
        ishape_psf_e1 = (ishape_psf_shape11-ishape_psf_shape22)/(ishape_psf_shape11+ishape_psf_shape22)
        ishape_psf_e2 = 2*ishape_psf_shape12/(ishape_psf_shape11+ishape_psf_shape22)
        i_apertureflux_10_mag = sourceargs['fits']["i_apertureflux_10_mag"]
        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_shape_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_shear_bias_m, ishape_hsm_regauss_derived_shear_bias_c1, ishape_hsm_regauss_derived_shear_bias_c2, ishape_hsm_regauss_resolution, i_apertureflux_10_mag, ishape_psf_e1, ishape_psf_e2]).T
        datagalflag = sourceargs['fits']["weak_lensing_flag"]

        try:
            sourceargs['fitspofz'] = fitsio.FITS("DataStore/S19A_v2.0/%s_tracts/%d_pz_pdf_%s.fits" % (field, tract, sourceargs['pofz_type']))[1]

            pofz = sourceargs['fitspofz']['PDF'][:]
            pofz = pofz.reshape((pofz.shape[0], -1,))
        except:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = pofzflag.size
        status = 0

        return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag

    if sourceargs['type'] == "hsc-wide-s16a_v2.0" and sourceargs['filetype'] == "fits" and not sourceargs["fullpofz"]:

        # This variable iterates over the different fields in the catalog
        itern = sourceargs['iter']
        if itern == 0:
            if not sourceargs["blinded"]:
                sourceargs["df"] = pandas.read_csv("DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat.calibrated", delim_whitespace=1)
            else:
                sourceargs["df"] = pandas.read_csv("DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat", delim_whitespace=1)
            if "usefield" in sourceargs:
                restidx = sourceargs["df"].field.values==sourceargs["usefield"]
                sourceargs["df"] = sourceargs["df"][restidx]
                print("Using only field", sourceargs["usefield"])

        datagal = 0
        Ngal = 0
        pofz = 0
        status = not (itern<(sourceargs['df'].field.values.size))
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        field = sourceargs["df"].field.values[itern]
        tract = sourceargs["df"].tract.values[itern]
        begin = sourceargs["df"].begin.values[itern]
        end = sourceargs["df"].end.values[itern]
        
        # Initialize list of tracts within a given field
        if begin == 0:
            # Define the shape catalog file from Rachel to read from for this field
            if sourceargs["blinded"]:
                # Decrypt the Delta m value, if sourceargs["dm"] not set already
                fname = "DataStore/Blinded_S16A_v2.0/%s_blinded_m_%s_%d.fits" % (field, sourceargs["username"], sourceargs["blindnumber"])
                fname_nom = "DataStore/Blinded_S16A_v2.0/%s_no_m.fits" % (field)
                if "dm" not in sourceargs:
                    try:
                        hdulist = pyfits.open(fname, memmap=True)
                        msg = hdulist[0].header["DM1"]
                        sourceargs['dm'] = float(decrypt_string(FITS_to_PGP(msg)))
                    except:
                        print("GPG decryption failed, check your gpg-agent")
                        exit(12)
                # The commented line is memory intensive
                #sourceargs["fits"] = fitsio.FITS("%s[1][col *, ishape_hsm_regauss_derived_shear_bias_m = ishape_hsm_regauss_derived_shear_bias_m - %f ]" % (fname, sourceargs['dm']))[1]
                # Instead read the file as it is and then subtract off later
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["fits"] = fitsio.FITS("%s" % (fname_nom))[1]
            else:
                fname = "DataStore/Calibrated_S16A_v2.0/%s_calibrated.fits" % (field)
                sourceargs["fits"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["dm"] = 0.0

        # Read in the photoz file for the current field, tract
        try:
            hdulist = pyfits.open("DataStore/S16A_v2.0/%s_tracts/%d_pz.fits" % (field, tract))
            data = hdulist[1].data
            pofz = data["%s_photoz%s" % (sourceargs['pofz_type'], sourceargs['photoz_estimate']) ]
            pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<sourceargs['photoz_risk_best_cut'])
        except:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        # Read in the shape values for the current field, tract, the row filtering specification indexes rows from 1 - Nrecords
        wheremask = sourceargs['fits'].where("#row>%d && #row<=%d" % (begin, end))

        #datagal = sourceargs['fits'][1]["ira", "idec", "ishape_hsm_regauss_e1", "ishape_hsm_regauss_e2", "ishape_hsm_regauss_derived_shear_bias_c1", "ishape_hsm_regauss_derived_shear_bias_c2", "ishape_hsm_regauss_derived_shear_bias_m", "ishape_hsm_regauss_derived_shape_weight", "ishape_hsm_regauss_derived_rms_e"][wheremask]
        #datagal = sourceargs['fits']["ira", "idec", "ishape_hsm_regauss_e1", "ishape_hsm_regauss_e2", "ishape_hsm_regauss_derived_shape_weight", "ishape_hsm_regauss_derived_rms_e", "ishape_hsm_regauss_derived_shear_bias_m", "ishape_hsm_regauss_derived_shear_bias_c1", "ishape_hsm_regauss_derived_shear_bias_c2"][wheremask]
        ira = sourceargs['fits']["ira"][wheremask]
        idec = sourceargs['fits']["idec"][wheremask]
        ishape_hsm_regauss_e1 = sourceargs['fits']["ishape_hsm_regauss_e1"][wheremask]
        ishape_hsm_regauss_e2 = sourceargs['fits']["ishape_hsm_regauss_e2"][wheremask]
        ishape_hsm_regauss_derived_shape_weight = sourceargs['fits']["ishape_hsm_regauss_derived_shape_weight"][wheremask]
        ishape_hsm_regauss_derived_rms_e = sourceargs['fits']["ishape_hsm_regauss_derived_rms_e"][wheremask]
        ishape_hsm_regauss_derived_shear_bias_m = sourceargs['fits_mblind']["ishape_hsm_regauss_derived_shear_bias_m"][wheremask] - sourceargs["dm"]
        ishape_hsm_regauss_derived_shear_bias_c1 = sourceargs['fits']["ishape_hsm_regauss_derived_shear_bias_c1"][wheremask]
        ishape_hsm_regauss_derived_shear_bias_c2 = sourceargs['fits']["ishape_hsm_regauss_derived_shear_bias_c2"][wheremask]
        ishape_hsm_regauss_resolution = sourceargs['fits']["ishape_hsm_regauss_resolution"][wheremask]
        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_shape_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_shear_bias_m, ishape_hsm_regauss_derived_shear_bias_c1, ishape_hsm_regauss_derived_shear_bias_c2, ishape_hsm_regauss_resolution]).T
        datagalflag = sourceargs['fits']["weak_lensing_flag"][wheremask]

        sourceargs['iter'] = sourceargs['iter'] + 1

        Ngal = pofzflag.size
        status = 0
        
        return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag

    if sourceargs['type'] == "hsc-wide-s16a_v2.0" and sourceargs['filetype'] == "fits" and  sourceargs["fullpofz"]:

        # This variable iterates over the different fields in the catalog
        itern = sourceargs['iter']
        if itern == 0:
            if not sourceargs["blinded"]:
                sourceargs["df"] = pandas.read_csv("DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat.calibrated", delim_whitespace=1)
            else:
                sourceargs["df"] = pandas.read_csv("DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat", delim_whitespace=1)
            if "usefield" in sourceargs:
                restidx = sourceargs["df"].field.values==sourceargs["usefield"]
                sourceargs["df"] = sourceargs["df"][restidx]
                print("Using only field", sourceargs["usefield"])

        datagal = 0
        Ngal = 0
        pofz = 0
        status = not (itern<(sourceargs['df'].field.values.size))
        if status:
            return datagal, sourceargs, Ngal, status, pofz, 0., 0.

        field = sourceargs["df"].field.values[itern]
        tract = sourceargs["df"].tract.values[itern]
        begin = sourceargs["df"].begin.values[itern]
        end = sourceargs["df"].end.values[itern]
        
        # Initialize list of tracts within a given field
        if begin == 0:
            # Define the shape catalog file from Rachel to read from for this field
            if sourceargs["blinded"]:
                # Decrypt the Delta m value, if sourceargs["dm"] not set already
                fname = "DataStore/Blinded_S16A_v2.0/%s_blinded_m_%s_%d.fits" % (field, sourceargs["username"], sourceargs["blindnumber"])
                fname_nom = "DataStore/Blinded_S16A_v2.0/%s_no_m.fits" % (field)
                if "dm" not in sourceargs:
                    try:
                        hdulist = pyfits.open(fname, memmap=True)
                        msg = hdulist[0].header["DM1"]
                        sourceargs['dm'] = float(decrypt_string(FITS_to_PGP(msg)))
                    except:
                        print("GPG decryption failed, check your gpg-agent")
                        exit(12)
                # The commented line is memory intensive
                #sourceargs["fits"] = fitsio.FITS("%s[1][col *, ishape_hsm_regauss_derived_shear_bias_m = ishape_hsm_regauss_derived_shear_bias_m - %f ]" % (fname, sourceargs['dm']))[1]
                # Instead read the file as it is and then subtract off later
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["fits"] = fitsio.FITS("%s" % (fname_nom))[1]
            else:
                fname = "DataStore/Calibrated_S16A_v2.0/%s_calibrated.fits" % (field)
                sourceargs["fits"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["fits_mblind"] = fitsio.FITS("%s" % (fname))[1]
                sourceargs["dm"] = 0.0

        # Read in the photoz file for the current field, tract
        try:
            hdulist = pyfits.open("DataStore/S16A_v2.0/%s_tracts/%d_pz.fits" % (field, tract))
            data = hdulist[1].data
            pofz = data["%s_photoz_best" % (sourceargs['pofz_type']) ]
            pofzflag = (data["%s_photoz_risk_best" % (sourceargs['pofz_type']) ]<sourceargs['photoz_risk_best_cut'])
        except:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        # Read in the shape values for the current field, tract
        wheremask = sourceargs['fits'].where("#row>%d && #row<=%d" % (begin, end))

        #datagal = sourceargs['fits']["ira", "idec", "ishape_hsm_regauss_e1", "ishape_hsm_regauss_e2", "ishape_hsm_regauss_derived_shear_bias_c1", "ishape_hsm_regauss_derived_shear_bias_c2", "ishape_hsm_regauss_derived_shape_weight", "ishape_hsm_regauss_derived_rms_e"][wheremask]
        #datagalflag = sourceargs['fits'][1]["weak_lensing_flag"][wheremask]

        ira = sourceargs['fits']["ira"][wheremask]
        idec = sourceargs['fits']["idec"][wheremask]
        ishape_hsm_regauss_e1 = sourceargs['fits']["ishape_hsm_regauss_e1"][wheremask]
        ishape_hsm_regauss_e2 = sourceargs['fits']["ishape_hsm_regauss_e2"][wheremask]
        ishape_hsm_regauss_derived_shape_weight = sourceargs['fits']["ishape_hsm_regauss_derived_shape_weight"][wheremask]
        ishape_hsm_regauss_derived_rms_e = sourceargs['fits']["ishape_hsm_regauss_derived_rms_e"][wheremask]
        ishape_hsm_regauss_derived_shear_bias_m = sourceargs['fits_mblind']["ishape_hsm_regauss_derived_shear_bias_m"][wheremask] - sourceargs["dm"]
        ishape_hsm_regauss_derived_shear_bias_c1 = sourceargs['fits']["ishape_hsm_regauss_derived_shear_bias_c1"][wheremask]
        ishape_hsm_regauss_derived_shear_bias_c2 = sourceargs['fits']["ishape_hsm_regauss_derived_shear_bias_c2"][wheremask]
        ishape_hsm_regauss_resolution = sourceargs['fits']["ishape_hsm_regauss_resolution"][wheremask]
        datagal = np.array([ira, idec, ishape_hsm_regauss_e1, ishape_hsm_regauss_e2, ishape_hsm_regauss_derived_shape_weight, ishape_hsm_regauss_derived_rms_e, ishape_hsm_regauss_derived_shear_bias_m, ishape_hsm_regauss_derived_shear_bias_c1, ishape_hsm_regauss_derived_shear_bias_c2, ishape_hsm_regauss_resolution]).T
        datagalflag = sourceargs['fits']["weak_lensing_flag"][wheremask]

        try:
            sourceargs['fitspofz'] = fitsio.FITS("DataStore/S16A_v2.0/%s_tracts/%d_pz_pdf_%s.fits" % (field, tract, sourceargs['pofz_type']))[1]

            pofz = sourceargs['fitspofz']['P(z)'][:]
            pofz = pofz.reshape((pofz.shape[0], -1,))
        except:
            pofz = 0
            pofzflag = np.zeros(end-begin, dtype=bool)

        sourceargs['iter'] = sourceargs['iter'] + 1
        Ngal = pofzflag.size
        status = 0

        return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag

    if sourceargs['type'] == "Euclid_Q1":
        print("Euclid_Q1 not implemented yet")
        exit(11)

    if sourceargs['type'] == "Euclid_RR2":
        print("Euclid_RR2 not implemented yet")
        exit(11)
        
    if sourceargs['type'] == "Euclid_Flagship_LightCone":
        # Euclid Flagship mock: source_bin_<N>.fits in sourceargs['dir']
        # Columns: ra_gal, dec_gal, true_redshift_gal, eps1_gal, eps2_gal, gamma1, gamma2
        # Observed ellipticity (weak-lensing linear approx.): e = eps_intrinsic + gamma
        # Reads in chunks of chunksize rows; use with include_corrections: False
        # and fullpofz: False (true_redshift_gal is used as the source redshift).
        if 'dir' not in sourceargs:
            raise RuntimeError("source.dir is required for Euclid_Flagship_LightCone")
        source_bin_id = sourceargs.get('source_bin_id', 6)
        filename = os.path.join(sourceargs['dir'], 'source_bin_%d.fits' % source_bin_id)

        # Initialize row offset on first call
        if 'source_row_start' not in sourceargs:
            with fitsio.FITS(filename) as f:
                sourceargs['source_nrows'] = f[1].get_nrows()
            sourceargs['source_row_start'] = 0

        row_start = sourceargs['source_row_start']
        nrows_total = sourceargs['source_nrows']
        if row_start >= nrows_total:
            status = 1
            return 0, sourceargs, 0, status, 0, 0, 0

        row_end = min(row_start + chunksize, nrows_total)
        d = fitsio.read(filename, ext=1, rows=range(row_start, row_end))

        e1 = d['eps1_gal'].astype('float64') + d['gamma1'].astype('float64')
        e2 = d['eps2_gal'].astype('float64') + d['gamma2'].astype('float64')
        rms_e = 0.25 * np.ones(len(d['ra_gal']))
        wgal = np.ones(len(d['ra_gal']))  # No weight column in Flagship

        # weaklens_aroundsource flips the sign of e2, so store -e2 here
        datagal = np.array([d['ra_gal'], d['dec_gal'], e1, -e2, wgal, rms_e]).T
        datagalflag = np.ones(len(d['ra_gal']), dtype=bool)
        pofz = d['true_redshift_gal'].astype('float64')
        pofzflag = np.ones(len(d['ra_gal']), dtype=bool)

        Ngal = len(d['ra_gal'])
        status = 0
        sourceargs['source_row_start'] = row_end
        sourceargs['iter'] = sourceargs['iter'] + 1
        sys.stdout.write("Number of sources (Euclid Flagship LightCone bin %d): rows %d-%d of %d\n" % (source_bin_id, row_start, row_end, nrows_total))

        return datagal, sourceargs, Ngal, status, pofz, datagalflag, pofzflag


def get_pofz_array(sourceargs):
    if sourceargs['type'] == "hsc-wide-s16a_v2.0":
        hdulist = pyfits.open("DataStore/S16A_v2.0/pz_pdf_bins_%s.fits" % sourceargs['pofz_type'])
        data = hdulist[1].data
        return data["bins"]
    if sourceargs['type'] == "hsc-wide-s19a_v2.0":
        hdulist = pyfits.open("DataStore/S19A_v2.0/XMM_tracts/8279_pz_pdf_%s.fits" % sourceargs['pofz_type'])
        data = hdulist[2].data
        return data["bins"]
