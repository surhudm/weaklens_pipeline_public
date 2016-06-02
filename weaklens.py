import pandas
from scipy.spatial import cKDTree
import numpy as np
import sys
import pyfits
import argparse
from scipy.interpolate import interp1d
import yaml
from select import lens_select, source_select
import datetime
import time
import frogress

sys.path.append("/home/surhud/aum/install/lib/python2.7/site-packages/")
import cosmology as cc

# Some constants
gee = 4.2994e-9
cee = 299792.458

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--config", help="Configuration file")
parser.add_argument("--rmax", help="Max radius", type=float, default=15.0)
parser.add_argument("--rmin", help="Min radius", type=float, default=0.1)
parser.add_argument("--rbin", help="Radial bin", type=int, default=15)
parser.add_argument("--pairout", help="Output filename with pairs information", default="debug_pairs.dat")

args = parser.parse_args()

with open(args.config, 'r') as ymlfile:
    config = yaml.load(ymlfile)

# Convert arguments to variables
rmax = (args.rmax)
rmin = (args.rmin)
logrmax = np.log10(args.rmax)
logrmin = np.log10(args.rmin)
rbin = args.rbin
logrdiff = (logrmax-logrmin)/rbin
pairout = args.pairout

dcount_sum = np.zeros(rbin)
dsig_sum = np.zeros(rbin)
dsigcross_sum = np.zeros(rbin)
wls_sum = np.zeros(rbin)
wls_ms_sum = np.zeros(rbin)

# Only important parameter here is Omega0 = 0.27
cosmo = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,np.log10(8.0),1.0)

# First get all the lens galaxy data
lensargs = lens_parseargs(config["lens"])
ra, dec, zred, wt = lens_select(lensargs)

ra = ra * np.pi / 180.
dec = dec * np.pi / 180.

c_lra=np.cos(ra)
s_lra=np.sin(ra)
s_ldec=np.sin(dec)
c_ldec=np.cos(dec)

xx = c_ldec*c_lra
yy = c_ldec*s_lra
zz = s_ldec
zzp = c_ldec

# For every zred compute the dlofz in hinv Mpc
dCofzl = zred * 0.0
for i in range(zred.size):
    dCofzl[i] = cosmo.Dcofz(zred[i])

# Ready to pounce on the source data, a million galaxies at a time
itern=0
done = False
sourceargs = source_parseargs(config["source"])

if sourceargs["fullpofz"]:
    pofz_xz = sourceargs["pofz_zmin"] + sourceargs["pofz_zdiff"]*np.arange(sourceargs["pofz_zbins"])
    pofz_dCofzs = zphot_arr * 0.0
    for i in range(zphot_arr.size):
        pofz_dCofzs[i] = cosmo.Dcofz(pofz_xz[i])
        pofz_sigc_inv_fac[i]=1e12*4*M_PI*gee/(cee*cee)*(1+pofz_xz[kk])/(pofz_dCofzs[i])

while (done):
    ragal, decgal, zphotgal, e1gal, e2gal, mgal, wgal, sourceargs, status = getdata(sourceargs, itern, chunksize=1000000)
    if status!=0:
        break

    ragal  = ragal * np.pi/180.0
    decgal = decgal * np.pi/180.0

    dCofzs = zphotgal * 0.0
    for i in range(zphotgal.size):
        dCofzs[i] = cosmo.Dcofz(zphotgal[i])


    # Sigc_inv = 4piG/c^2 dL dLS/dS (1+zl)**2
    # First factor = 4piG/c^2 / dS
    if not sourceargs["fullpofz"]:
        sigc_inv_fac=1e12*4*np.pi*gee/(cee*cee)*(1+zphotgal)/dCofzs

    c_sra=np.cos(ragal)
    s_sra=np.sin(ragal)
    s_sdec=np.sin(decgal)
    c_sdec=np.cos(decgal)

    xxgal = c_sdec*c_sra
    yygal = c_sdec*s_sra
    zzgal = s_sdec
    zzgalp = c_sdec

    tree = cKDTree(zip(xxgal, yygal, zzgal))
    # For every galaxy cluster, query the sources around it
    for i in frogress.bar(range(ra.size)):

        ind = tree.query_ball_point([xx[i], yy[i], zz[i]], rmax/(dCofzl[i]))
        if len(ind)==0:
            # print "skipping", i
            continue
        else:
            # print ind
            pass

        if np.isscalar(ind):
            ind = [ind]
        ind = np.array(ind)
        # print "1:", ind

        # zphot should be larger than the lens redshift plus deltaz
        # idx = (zphotgal[ind]>zred[i]+deltaz)
        # ind = ind[idx]

        # Now compute the actual angles and distances
        cosang = xx[i]*xxgal[ind] + yy[i]*yygal[ind] + zz[i]*zzgal[ind]
        Rp = dCofzl[i] * np.arccos(cosang)
        idx = (Rp>rmin)&(Rp<rmax)

        # Now compute the
        cos_alps_alpl=c_sra[ind][idx]*c_lra[i]+s_sra[ind][idx]*s_lra[i];
        sin_alps_alpl=s_sra[ind][idx]*c_lra[i]-c_sra[ind][idx]*s_lra[i];
        # cosangle=s_ldec[i]*s_sdec[ind][idx]+c_ldec[i]*c_sdec[ind][idx]*cos_alps_alpl;
        cosangle = cosang[idx]
        sinangle=np.sqrt(1-cosangle*cosangle);
        cosphi=c_ldec[i]*sin_alps_alpl/sinangle;
        sinphi=(-s_ldec[i]*c_sdec[ind][idx]+c_ldec[i]*s_sdec[ind][idx]*cos_alps_alpl)/sinangle;
        etan=-e1gal[ind][idx]*(2*cosphi*cosphi-1)+e2gal[ind][idx]*(2*sinphi*cosphi);
        ecross=e1gal[ind][idx]*(2*sinphi*cosphi)+e2gal[ind][idx]*(2*cosphi*cosphi-1);
        # This is CFHT part
        #etan=-e1gal[ind][idx]*(2*cosphi*cosphi-1)-e2gal[ind][idx]*(2*sinphi*cosphi);

        # Sigc_inv = 4piG/c^2 dL dLS/dS (1+zl)**2
        # First factor = 4piG/c^2 / dS
        # Second factor = dL dLS (1+zl)**2
        if sourceargs["fullpofz"]:
            sigc_inv_den = np.zeros(np.sum(idx))
            sigc_inv_num = np.zeros(np.sum(idx))
            for k in range(sigc_inv_den.size):
                sigc_inv_den[k] = np.sum(pofz[ind][idx][k])
                #fac=pz_fac[k]*chi[xx]*(pz_chiz[k]-chi[xx])/(1+pz_xz[k])*(1+z[xx])
                #av_sigc_inv_num+=pz[k]*fac
                pz_idx = (pofz_dCofzs-dCofzl[i] > 0)
                sigc_inv_num[k] = np.sum(pofz[ind][idx][k][pz_idx]*pofz_sigc_inv_fac[pz_idx]*dCofzl[i]*(pofz_dCofzs[pz_idx]-dCofzl[i])/(1.+pofz_xz[pz_idx])*(1.+z[i]))
                sigc_inv[k] = sigc_inv_num/sigc_inv_den
        else:
            sigc_inv = sigc_inv_fac[ind][idx]*dCofzl[i]*(dCofzs[ind][idx]-dCofzl[i])/(1.+zphotgal[ind][idx])*(1+zred[i])

        # Now compute the weights
        wls = wgal[ind][idx]*wt[i]*pow(sigc_inv,2)
        oneplusmcorr = (1+mgal[ind][idx])

        # print "4:", Rp[idx]
        Rp = Rp[idx]

        # Finalize bins and all arrays
        logRp = np.log10(Rp)
        ibin = np.floor((logRp-logrmin)/logrdiff)
        idx = (ibin>=0)&(ibin<rbin)
        ibin = ibin[idx]
        etan = etan[idx]
        ecross = ecross[idx]
        wls = wls[idx]
        sigc_inv = sigc_inv[idx]
        oneplusmcorr = oneplusmcorr[idx]
        for j in range(ibin.size):
            dsig_sum[ibin[j]] += (etan[j]*wls[j]/sigc_inv[j])
            dcount_sum[ibin[j]] += wt[i]
            dsigcross_sum[ibin[j]] += (ecross[j]*wls[j]/sigc_inv[j])
            wls_sum[ibin[j]] += (wls[j])
            wls_ms_sum[ibin[j]] += (wls[j]*oneplusmcorr[j])
            # fout.write("%.5f %.5f\n"%(Rp[j], Theta[j]))


    tstamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
    sys.stderr.write("%s : %d\n"%(tstamp, itern))
    itern = itern + 1


rpmin = 10.0**(logrmin+logrdiff*np.arange(rbin))
rpmax = 10.0**(logrmin+logrdiff*(1+np.arange(rbin)))
np.savetxt(pairout, np.transpose([logrmin+logrdiff*(0.5+np.arange(rbin)), dsig_sum, wls_sum, wls_ms_sum, dsig_sum/wls_ms_sum, wls_sum/np.sum(wt), rpmin/2.+rpmax/2., dsigcross_sum/wls_ms_sum, dcount_sum/np.sum(wt)]), header=" log (r), dSigma_sum, wls_sum, wls_ms_sum, Delta Sigma, Boost factor term, Average r, Delta Sigma cross")
