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
typedef boost::multi_array<double,1> array1ddouble;
#include <boost/timer.hpp>

class tomography;

class tomography
{
    protected:
        cosmology *cosmo;
        bool verbose;
        double rmin;
        double rmax;
        double logrmin;
        double logrmax;
        double logrdiff;
        int rbins;
        double *ra, *dec, *zred, *lam, *wt, *datacen, *c_lra, *s_lra, *c_ldec, *s_ldec, *Dcom;
        double zsrcmin, zsrcmax;
        int Ncen, Nfill;
        double zredmin, zredmax, Dcommin;
        double *sumdgammat_num, *sumdgammacross_num, *sumdcount_num, *sumwls, *sumwls_ms;
        double sumwl, pofzinbin;
        char outfile[1000];
        boost::timer *timer;
        long int Niter;
        bool lenses_alloc, lenses_finalized;
        array1ddouble xper;
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
        tomography(double xrmin, double xrmax, int xrbins, double xzsrcmin, double xzsrcmax, char *xoutfile, bool xverbose);
        ~tomography();
        int allocate_lens_memory(int xNcen);
        int process_lens(double xra, double xdec, double xzred, double xwt);
        int finalize_lenses();
        int process_source(double sra, double sdec, double se1, double se2, double swt, double smcat, double sc1_dp, double sc2_dp, double sc1_nb, double sc2_nb, double szbest, bool usepdf);
        int setup_pofz(double pofz_zmin, double pofz_zdiff, int xNpz);
        int process_pofz(double *pofz, int zbins);
        int finalize_results();
        int test_searchrecord();
};

#endif
