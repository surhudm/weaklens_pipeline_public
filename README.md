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

# Data directory organization

## For the calibrated/unblinded catalogs

```
mkdir -p DataStore/Calibrated_S16A_v2.0/
cd DataStore/Calibrated_S16A_v2.0/

# Download zip file from https://cmu.app.box.com/s/y3pdl0v8952ddairx59eo2a2f09kfwkh and unzip in this directory

```
The filenames in this directory will be of the type: `fieldname_calibrated.fits`

## For the blinded catalogs

```
mkdir -p DataStore/Blinded_S16A_v2.0/
cd DataStore/Blinded_S16A_v2.0/

# Download zip file from link given by Rachel and unzip in this directory

```
The filenames in this directory will be of the type: `fieldname_no_m.fits` and `fieldname_blinded_m_username_[012].fits`

## The source P(z) and PDF

For convenience Hironao has made the photozs and the photoz pdfs for sources available. Download them from http://master.ipmu.jp/~miyatake/shape_catalog/S16A_v2.0

Dump them as they are in the directory `DataStore/S16A_v2.0`

# Sample configuration file

## For point estimate of photoz
```
# Matter density
Omegam: 0.31
# Output directory
outputdir: camira/mizuki/zbest_0.70_1.20_10.0_1000.0
# Radial bins
rbin: 20
rmax: 20.0
rmin: 0.1
# Whether to output lens source pair information
outputpairfile : "allpairs.dat"
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
  photoz_estimate: _best
  zdiff: 0.0
```
`photoz_estimate` can be `_best` or `_mc` or `_median` or `_mean` that are
supported by the `P(z)` files provided by HSC. If no such option is provided
then the default would be `_best`.

## For using the full P(z)
```
# Matter density
Omegam: 0.31
# Output directory
outputdir: camira/mizuki/zbest_0.70_1.20_10.0_1000.0
# Radial bins
rbin: 20
rmax: 20.0
rmin: 0.1
# Whether to output lens source pair information
outputpairfile : "allpairs.dat"
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
  fullpofz: true
  photoz_risk_best_cut: 0.5
  integral_cut_zl: 0.8
  integral_cut_zdiff: 0.0
  integral_cut_Pth: 0.99
  integral_cut_zmax: 4.0
```

The parameters `integral_cut_zl`, `integral_cut_zdiff`, `integral_cut_Pth`, `integral_cut_zmax` are used to set a cut on the sources to be used.

If `integral_cut_zl=-99.0`, then use all the sources.

If `integral_cut_zl=-49.0`, then use sources if they satisfy
```
\int_{zlens+zdiff}^{zmax} P(z) dz > integral_cut_Pth \,.
```

If `integral_cut_zl>0`,
```
\int_{integral_cut_zl+zdiff}^{zmax} P(z) dz > integral_cut_Pth \,.
```

# Calling signature

```
python weaklens_aroundsource.py --config configfile
```

# FAQ:

1) `python setup.py install` fails complaining :

```
src/gauleg.cpp:9: error: ‘gsl_integration_glfixed_table’ was not declared in this scope
```

Please install `gsl>=v1.14` and then retry. In case you have multiple versions of gsl installed, then try:

```
CPPFLAGS=-I/path/to/gsl/include python setup.py install
```

2) You do not have the mapping from Rachel to Hironao's catalogs. These can be downloaded now from:

```
http://master.ipmu.jp/~miyatake/shape_catalog/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat
http://master.ipmu.jp/~miyatake/shape_catalog/S16A_v2.0/Mapping_from_Rachel_catalogs_to_Hironao_catalogs.dat.calibrated
```

3) See the .circleci files for the configuration which works in order to
compile the code. There are requirements for the python code mentioned in
requirements.txt. They can be installed to a local directory with:

```
pip install --user -r requirements.txt
```
