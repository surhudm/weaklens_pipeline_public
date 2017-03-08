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
typedef boost::multi_array<double,1> array1ddouble;
#include <boost/timer.hpp>

class weaklens;

class weaklens
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
        int Ncen, Nfill;
        double zredmin, zredmax, Dcommin;
        double *sumdsig_num, *sumdcross_num, *sumdcount_num, *sumwls, *sumwls_ms, *sumdsigsq_num, *sumdcrosssq_num, *sumwls_resp;
        double sumwl;
        char outfile[1000];
        boost::timer *timer;
        long int Niter;
        bool lenses_alloc, lenses_finalized;
        array1ddouble xper;
        int dim;

        // The spline for inverse sigmacritical
        bool bool_init_sigcritinv;
        gsl_interp_accel *sigcritinv_acc;
        gsl_spline *sigcritinv_spline;
        int Npz;
        double *pz_xz, *pz_chiz, *pz_fac, *pz_pofz;
        double Nphotgal;
        bool pz_alloc;

        double gee;
        double cee;

    public:
        kdtree2::KDTree *tree;
        weaklens();
        weaklens(double xrmin, double xrmax, int xrbins, char *xoutfile, bool xverbose);
        ~weaklens();
        int allocate_lens_memory(int xNcen);
        int process_lens(double xra, double xdec, double xzred, double xwt);
        int finalize_lenses();
        int process_source(double sra, double sdec, double se1, double se2, double swt, double srms_e, double smcat, double sc1_dp, double sc2_dp, double sc1_nb, double sc2_nb, double szbest, bool usepdf);
        int setup_pofz(double pofz_zmin, double pofz_zdiff, int xNpz);
        int setup_pofz_array(double *pofz_zz, int xNpz);
        int process_pofz(double *pofz, int zbins);
        int add_pofz(double *pofz, int zbins);
        double sigmacritinverse(double zlens);
        int finalize_results();
        int finalize_pofz();
        int test_searchrecord();
};

#endif
