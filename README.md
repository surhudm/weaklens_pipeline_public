# weaklens_pipeline
Weak lensing pipeline for HSC

(1) The goal is to measure the HSC weak lensing signal around a set of different
lens catalogs, using a set of shape and photometric redshift catalogs from  HSC.

(2) Let us break this into steps: preprocess the lens and shape catalogs. In the
end as an input to weaklens.py, we need

    (a) For lens catalog (see lens_select in the file weaklens_select.py):
        (i) ra, dec, z, wt

    (b) For source catalog (see source_select in the file weaklens_select.py):
        (i) ra, dec, e1, e2, zmean, zmode, zmedian, zmc, z68up, z68down, z95up, z95down, and possibly p(z)

The code works by selecting lens galaxies, generating a tree of these lens
galaxies. Then for every source, it finds the lenses that are around it in
given radial bins. This allows the code to be very memory efficient as the full
p(z) distribution need not be stored in the memory. There are also less cache
misses.

# Sample configuration file
```
# Matter density
Omegam: 0.31
# Output directory
outputdir: camira/mizuki/zbest_0.70_1.20_10.0_1000.0
# Radial bins
rbin: 20
rmax: 20.0
rmin: 0.1
# Lens configuration (to be used in weaklens_select.py:lens_select)
lens:
  hsc-release: s16a
  lammax: 1000.0
  lammin: 10.0
  type: camira
  version: v1
  zmax: 1.2
  zmin: 0.7
# Source configuration (to be used in weaklens_select.py:source_select)
source:
  # Blinding related
  username: more
  blinded: 1
  blindnumber: 0
  dm: 0.04353
  # Source file related
  type: hsc-wide-s16a_v2.0
  filetype: fits
  format: hsc
  include_corrections: true
  # Photo z related
  pofz_type: mizuki
  fullpofz: false
  photoz_risk_best_cut: 0.5
```

# Calling signature

```
python weaklens_aroundsource.py --config configfile
```
