from astropy.io import fits as pyfits 
from glob import glob
import numpy as np

'''
for pz in ["mlz", "frankenz", "mizuki", "ephor"]:
    hdulist = pyfits.open("pz_pdf_bins_%s.fits" % pz)
    data = hdulist[1].data
    print pz
    print data["bins"]
    print "===="
'''

# Test if the order of the pofz is preserved in the tract files

#for field in ["AEGIS", "GAMA09H", "GAMA15H", "HECTOMAP", "VVDS", "WIDE12H", "XMM"]:

fp = open("DataStore/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat", "w")
fp.write("field tract begin end matched\n")
for field in ["AEGIS", "GAMA09H", "GAMA15H", "HECTOMAP", "VVDS", "WIDE12H", "XMM"]:
    # Read the filename
    #hdulist = pyfits.open("%s.fits" % (field))
    #hdulist = pyfits.open("DataStore/Blinded_S16A_v2.0/hsc-blind/%s_no_m.fits" % (field))
    hdulist = pyfits.open("DataStore/Calibrated_S16A_v2.0/%s_calibrated.fits" % (field))
    data = hdulist[1].data
    object_id = data["object_id"]
    begin = 0
    list_of_tracts = glob("DataStore/S16A_v2.0/%s_tracts/*_pz_pdf_mizuki.fits" % (field))
    list_of_tracts = np.sort([int(tract.split("_pz")[0].split("/")[1]) for tract in list_of_tracts])

    for tract in list_of_tracts:
        hdulist = pyfits.open("%s_tracts/%d_pz_pdf_mizuki.fits" % (field, tract))
        obj_id = hdulist[1].data["object_id"]

        final = begin+obj_id.size    
        if np.any(object_id[begin:final]!=obj_id):
            exit(11)
        else:
            fp.write("%s %d %d %d %d \n" % (field, tract, begin, final, np.all(object_id[begin:final]==obj_id)))
            import sys
            sys.stderr.write("%s %d %d %d %d \n" % (field, tract, begin, final, np.all(object_id[begin:final]==obj_id)))

        begin = final * 1
