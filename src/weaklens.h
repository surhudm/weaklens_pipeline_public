#ifndef WEAKLENS_H
#define WEAKLENS_H

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

class weaklens;

class weaklens
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
        int Ncen, Nfill;
        float zredmin, zredmax, Dcommin;
        double *sumdsig_num, *sumdcross_num, *sumdcount_num, *sumwls, *sumwls_ms;
        double sumwl;
        char outfile[1000];
        boost::timer *timer;
        long int Niter;
        bool lenses_alloc, lenses_finalized;
        array1dfloat xper;
        int dim;

        // The spline for inverse sigmacritical
        bool bool_init_sigcritinv;
        gsl_interp_accel *sigcritinv_acc;
        gsl_spline *sigcritinv_spline;
        int Npz;
        double *pz_xz, *pz_chiz, *pz_fac;
        bool pz_alloc;

        double gee;
        double cee;

    public:
        kdtree2::KDTree *tree;
        weaklens();
        weaklens(float xrmin, float xrmax, int xrbins, char *xoutfile, bool xverbose);
        ~weaklens();
        int allocate_lens_memory(int xNcen);
        int process_lens(float xra, float xdec, float xzred, float xwt);
        int finalize_lenses();
        int process_source(float sra, float sdec, float se1, float se2, float swt, float smcat, float sc1_dp, float sc2_dp, float sc1_nb, float sc2_nb, float szbest, bool usepdf);
        int setup_pofz(double pofz_zmin, double pofz_zdiff, int xNpz);
        int process_pofz(double *pofz, int zbins);
        double sigmacritinverse(double zlens);
        int finalize_results();
        int test_searchrecord();
};

#endif
