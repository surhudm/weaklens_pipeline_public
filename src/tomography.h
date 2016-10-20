#ifndef TOMOGRAPHY_H
#define TOMOGRAPHY_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <fstream>
#include "cosmology.h"
#include "kdtree2.h"
#include <vector>
#include <ctime>
#include <cstring>
#include <cmath>
#include <boost/multi_array.hpp>
typedef boost::multi_array<float,1> array1dfloat;
#include <boost/timer.hpp>

class tomography;

class tomography
{
    protected:
        cosmology *cosmo;
        bool verbose;
        float rmin;
        float rmax;
        float logrmin;
        float logrmax;
        float logrdiff;
        int rbins;
        float *ra, *dec, *zred, *lam, *wt, *datacen, *c_lra, *s_lra, *c_ldec, *s_ldec, *Dcom;
        float zsrcmin, zsrcmax;
        int Ncen, Nfill;
        float zredmin, zredmax, Dcommin;
        double *sumdgammat_num, *sumdgammacross_num, *sumdcount_num, *sumwls, *sumwls_ms;
        double sumwl, pofzinbin;
        char outfile[1000];
        boost::timer *timer;
        long int Niter;
        bool lenses_alloc, lenses_finalized;
        array1dfloat xper;
        int dim;

        // The spline for inverse sigmacritical
        int Npz;
        double *pz_xz;
        bool pz_alloc;

        double gee;
        double cee;

    public:
        kdtree2::KDTree *tree;
        tomography();
        tomography(float xrmin, float xrmax, int xrbins, float xzsrcmin, float xzsrcmax, char *xoutfile, bool xverbose);
        ~tomography();
        int allocate_lens_memory(int xNcen);
        int process_lens(float xra, float xdec, float xzred, float xwt);
        int finalize_lenses();
        int process_source(float sra, float sdec, float se1, float se2, float swt, float smcat, float sc1_dp, float sc2_dp, float sc1_nb, float sc2_nb, float szbest, bool usepdf);
        int setup_pofz(double pofz_zmin, double pofz_zdiff, int xNpz);
        int process_pofz(double *pofz, int zbins);
        int finalize_results();
        int test_searchrecord();
};

#endif
