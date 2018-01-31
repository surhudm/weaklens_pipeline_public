%module weaklens
%include cstring.i
%include cpointer.i
%pointer_class(double,dp)
%apply double& INOUT { double& a };

%include "carrays.i"
%array_class(double, doubleArray);
%feature("autodoc", 1);
%{
    #define SWIG_FILE_WITH_INIT
    #include "weaklens.h"
%}
%include numpy.i
%init %{
import_array();
%}
%apply (float* IN_ARRAY1, int DIM1) {(float* pofz, int zbins)};
%apply (double* IN_ARRAY1, int DIM1) {(double* pofz, int zbins)};
%apply (double* IN_ARRAY1, int DIM1) {(double* pofz_zz, int xNpz)};


%feature("docstring") weaklens::weaklens
"Initializes weaklens object 

:Parameters:

-   xrmin : The minimum of the radial bins
-   xrmax : The maximum of the radial bins
-   xrbins : The total number of the radial bins
-   xoutfile : The output file
-   xverbose : Boolean parameter for verbosity
-   xOmegam : Matter density parameter
-   output_lens_source_pairs : File to output lens source pair information to

:Returns:

-   Weaklens object

    Without any inputs, initializes to the default rmin=0.5, rmax=15.0, rbins=15, verbose=True, outfile=Debug.dat, Omegam=0.279, output_lens_source_pairs=""

:Examples:

>>> import weaklens as wl
>>> a = wl.weaklens(0.5, 15.0, 15, 'Debug.dat', true, 0.279, 'allpairs.dat')
>>> help(a)

"

%feature("docstring") weaklens::allocate_lens_memory
"Allocates memory for lens sample

:Parameters:

-   Nlens : Size of the lens sample (int)

:Examples:

>>> import weaklens as wl
>>> a = wl.weaklens(0.5, 15.0, 15, 'Debug.dat', true, 0.279)
>>> a.allocate_lens_memory(1000)
>>> help(a)

"

%feature("docstring") weaklens::process_lens
"Add a lens to the sample

:Parameters:

-   xra : Ra of the lens (float, in degrees)
-   xdec : Dec of the lens (float, in degrees)
-   xz : Redshift of the lens (float)
-   xwt : Weight of the lens (float)

:Examples:

>>> import weaklens as wl
>>> a = wl.weaklens(0.5, 15.0, 15, 'Debug.dat', true, 0.279)
>>> a.allocate_lens_memory(1000)
>>> a.process_lens(100., 2.0, 0.1, 1.0)
>>> help(a)

"

%feature("docstring") weaklens::finalize_lenses
"Finalize the lens sample, truncates the array if all allocated memory is not
filled and generates a KDTree from the lens points internally for finding
neighbours.

:Parameters:

-   None : 

:Examples:

>>> import weaklens as wl
>>> a = wl.weaklens(0.5, 15.0, 15, 'Debug.dat', true, 0.279)
>>> a.allocate_lens_memory(1000)
>>> for i in range(Ngal):
        a.process_lens(100., 2.0, 0.1, 1.0)
>>> a.finalize_lenses()
>>> help(a)

"

%feature("docstring") weaklens::process_source
"Calculate the weak lensing signal for a source

:Parameters:

-   sra : Source ra (in degrees)
-   sdec : Source declination (in degrees)
-   se1 : Source e1
-   se2 : Source e2
-   swt : Source weight
-   srms_e : Source rms ellipticity
-   smcat : Source multiplicative bias
-   sc1_dp : Source additive bias (c1_dp)
-   sc2_dp : Source additive bias (c2_dp)
-   sc1_nb : Source additive bias (c1_nb)
-   sc2_nb : Source additive bias (c2_nb)
-   szbest : Source best estimate of photoz or the photoz pdf
-   usepdf : Source use pdf or not (false by default)

:Examples:

>>> import weaklens as wl
>>> a = wl.weaklens(0.5, 15.0, 15, 'Debug.dat', true, 0.279)
>>> a.allocate_lens_memory(1000)
>>> for i in range(Ngal):
        a.process_lens(100., 2.0, 0.1, 1.0)
>>> a.finalize_lenses()
>>> a.process_source(100.1, 2.1, 0.21, 0.22, 0.5, 0.365, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, false)
>>> help(a)

"

%feature("docstring") weaklens::finalize_results
"Compile weak lensing signal and dump into output file

:Parameters:

-   None :

:Examples:

>>> import weaklens as wl
>>> a = wl.weaklens(0.5, 15.0, 15, 'Debug.dat', true, 0.279)
>>> a.allocate_lens_memory(1000)
>>> for i in range(Ngal):
        a.process_lens(100., 2.0, 0.1, 1.0)
>>> a.finalize_lenses()
>>> a.process_source(100.1, 2.1, 0.21, 0.22, 0.5, 0.365, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, false)
>>> a.finalize_results()
>>> help(a)

"

%feature("docstring") weaklens::setup_pofz_array
"Set up the redshift array for pofz

:Parameters:

-   zz : Numpy array with redshifts

:Examples:

"

%feature("docstring") weaklens::setup_pofz
"Set up the redshift array for pofz using zmin, zdiff and Nz

:Parameters:

-   zmin : zmin for P(z) array
-   zdiff : Delta z for P(z) array
-   Nz : Number of zbins for P(z) array

:Examples:


"

%feature("docstring") weaklens::process_pofz
"Process pofz distribution for a given source, sets up a spline for sigmacritinverse(zl)

:Parameters:

-   pofz : Numpy array with P(z)

:Examples:


"

%feature("docstring") weaklens::add_pofz
"Add pofz to the pofz distribution

:Parameters:

-   pofz : Numpy array with P(z)

:Examples:


"

%feature("docstring") weaklens::finalize_pofz
"Output pofz to outfile

:Parameters:

-   None :

:Examples:


"

%include "weaklens.h"
