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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class weaklens;

class weaklens
{
    protected:
        cosmology *cosmo;
        bool verbose;
        double R2min;
        double R2max;
        int R2bins;
        double R2diff;
        double rmin;
        double rmax;
        double logrmin;
        double logrmax;
        double logrdiff;
        bool random_rotate;
        int rbins;
        double *ra, *dec, *zred, *lam, *wt, *datacen, *c_lra, *s_lra, *c_ldec, *s_ldec, *Dcom;
        int Ncen, Nfill;
        double zredmin, zredmax, Dcommin;
        double *sumdsig_num, *sumct_num, *sumdcross_num, *sumcx_num, *sumdcount_num, *sumwls, *sumwls_ms, *sumdsigsq_num, *sumdcrosssq_num, *sumwls_resp, *sumpairs, *sumwls_R2firstbin, *sumwls_R2allbins;
        double sumwl;
        char outfile[10000];
        boost::timer *timer;
        long int Niter;
        bool lenses_alloc, lenses_finalized;
        array1ddouble xper;
        int dim;
        FILE *foutput_lens_source_pairs;

        // The spline for inverse sigmacritical
        bool bool_init_sigcritinv;
        gsl_interp_accel *sigcritinv_acc;
        gsl_spline *sigcritinv_spline;
        int Npz;
        double *pz_xz, *pz_chiz, *pz_fac, *pz_pofz, *current_src_pofz;
        double Nphotgal;
        bool pz_alloc;

        double gee;
        double cee;
        double Omegam;

        const gsl_rng_type * rngtype;
        gsl_rng * rng;

    public:
        kdtree2::KDTree *tree;
        weaklens();
        weaklens(double xrmin, double xrmax, int xrbins, char *xoutfile, bool xverbose, double Omegam=0.279, char *output_lens_source_pairs=(char *)"", double xR2min=0.3, double xR2max=1.0, int xR2bins=100);
        ~weaklens();
        int allocate_lens_memory(int xNcen);
        int process_lens(double xra, double xdec, double xzred, double xwt);
        int finalize_lenses();
        int process_source(double sra, double sdec, double se1, double se2, double swt, double srms_e, double smcat, double sc1_dp, double sc2_dp, double sc1_nb, double sc2_nb, double szbest, double sR2, bool usepdf=false, double integral_cut_zl=-99.0, double integral_cut_zdiff=0.0, double integral_cut_Pth=0.0, double integral_cut_zmax=99.0);
        int process_source(double sra, double sdec, double se1, double se2, double swt, double srms_e, double smcat, double sc1_dp, double sc2_dp, double szbest, double sR2, bool usepdf=false, double integral_cut_zl=-99.0, double integral_cut_zdiff=0.0, double integral_cut_Pth=0.0, double integral_cut_zmax=99.0);
        int setup_pofz(double pofz_zmin, double pofz_zdiff, int xNpz);
        int setup_pofz_array(double *pofz_zz, int xNpz);
        int process_pofz(double *pofz, int zbins);
        int add_pofz(double *pofz, int zbins);
        double integrate_pofz(double zmin, double zmax);
        double set_zlmax(double integral_cut_zdiff, double integral_cut_zmax, double integral_cut_Pth);
        double sigmacritinverse(double zlens);
        int finalize_results(bool writeok=false);
        int finalize_pofz();
        int test_searchrecord();
        int setup_random_rotate(int Nseed);
};

#endif
