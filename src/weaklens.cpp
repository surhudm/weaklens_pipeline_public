#include "weaklens.h"

weaklens::weaklens()
{
    rmin=0.5;
    rmax=15.0;
    rbins=15;
    logrmin=log10(rmin);
    logrmax=log10(rmax);
    logrdiff=(logrmax-logrmin)/rbins;
    sprintf(outfile, "Debug.dat");
    verbose=true;

    sumdsig_num = (double *)calloc(rbins,sizeof(double));
    sumdcross_num = (double *)calloc(rbins,sizeof(double));
    sumdcount_num = (double *)calloc(rbins,sizeof(double));
    sumwls = (double *)calloc(rbins,sizeof(double));
    sumwls_ms  =(double *)calloc(rbins,sizeof(double));


    // ========================================================
    // Initialize cosmology
    // ========================================================
    cosmo = new cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0);
    lenses_alloc=false;
    lenses_finalized=false;
    zredmin=1.E30;
    zredmax=-1.E30;
    bool_init_sigcritinv = false;
    pz_alloc=false;
    gee = 4.2994e-9;
    cee = 3.E5;

    if (verbose){
        fprintf(stderr, "Initialized cosmology, and delta sigma sum arrays\n");
        fprintf(stderr, "Options passed rmin:%f rmax:%f rbins:%d outfile:%s \n", rmin, rmax, rbins, outfile);
    }

}

weaklens::weaklens(double xrmin, double xrmax, int xrbins, char *xoutfile, bool xverbose)
{
    rmin=xrmin;
    rmax=xrmax;
    rbins=xrbins;
    logrmin=log10(rmin);
    logrmax=log10(rmax);
    logrdiff=(logrmax-logrmin)/rbins;
    sprintf(outfile, "%s", xoutfile);

    sumdsig_num = (double *)calloc(rbins,sizeof(double));
    sumdcross_num = (double *)calloc(rbins,sizeof(double));
    sumdcount_num = (double *)calloc(rbins,sizeof(double));
    sumwls = (double *)calloc(rbins,sizeof(double));
    sumwls_ms  =(double *)calloc(rbins,sizeof(double));

    verbose=xverbose;

    // ========================================================
    // Initialize cosmology
    // ========================================================
    cosmo = new cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0);
    lenses_alloc=false;
    lenses_finalized=false;
    zredmin=1.E30;
    zredmax=-1.E30;
    bool_init_sigcritinv = false;
    pz_alloc=false;
    gee = 4.2994e-9;
    cee = 3.E5;

    if (verbose){
        fprintf(stderr, "Initialized cosmology, and delta sigma sum arrays\n");
        fprintf(stderr, "Options passed rmin:%f rmax:%f rbins:%d outfile:%s \n", rmin, rmax, rbins, outfile);
    }
}

weaklens::~weaklens(){
    free(sumdsig_num);
    free(sumdcross_num);
    free(sumdcount_num);
    free(sumwls);
    free(sumwls_ms);
    delete cosmo;
    if (verbose){
        fprintf(stderr, "Freed memory related to cosmology and the dsigma arrays\n");
    }

    if (lenses_alloc){
        free(ra      );
        free(dec     );
        free(zred    );
        free(wt      );
        free(datacen );
        free(c_lra   );
        free(s_lra   );
        free(s_ldec  );
        free(c_ldec  );
        free(Dcom    );
        if (verbose){
            fprintf(stderr, "Freed memory related to lenses as lenses were allocated\n");
        }
    }else{
        if (verbose){
            fprintf(stderr, "Did not free memory related to lenses as lenses were not finalized\n");
        }
    }
    if (lenses_finalized){
        delete tree;
        delete timer;
        if (verbose){
            fprintf(stderr, "Freed memory related to tree and timer as lenses were finalized\n");
        }
    }else{
        if (verbose){
            fprintf(stderr, "Did not free memory related to tree and timer as lenses were not finalized\n");
        }
    }
    if (bool_init_sigcritinv){
        gsl_interp_accel_free(sigcritinv_acc);
        gsl_spline_free(sigcritinv_spline);
        bool_init_sigcritinv=false;
    }
    if (pz_alloc){
        delete [] pz_xz;
        delete [] pz_chiz;
        delete [] pz_fac;
        pz_alloc = false;
    }
}

int weaklens::allocate_lens_memory(int xNcen)
{
    Ncen = xNcen;
    ra      =(double *)  calloc(Ncen,sizeof(double));
    dec     =(double *)  calloc(Ncen,sizeof(double));
    zred    =(double *)  calloc(Ncen,sizeof(double));
    wt      =(double *)  calloc(Ncen,sizeof(double));
    datacen =(double *)  calloc(Ncen*3,sizeof(double));
    c_lra   =(double *)  calloc(Ncen,sizeof(double));
    s_lra   =(double *)  calloc(Ncen,sizeof(double));
    s_ldec  =(double *)  calloc(Ncen,sizeof(double));
    c_ldec  =(double *)  calloc(Ncen,sizeof(double));
    Dcom    =(double *)  calloc(Ncen,sizeof(double));
    Nfill=0;
    lenses_alloc=true;
    if (verbose){
        fprintf(stderr, "Initialized lenses Ncen=%d\n", Ncen);
    }
    return 0;
}

int weaklens::process_lens(double xra, double xdec, double xzred, double xwt)
{
    if (Nfill==Ncen){
        fprintf(stderr, "Lens arrays were not allocated correctly, please call allocate_lens_memory once again with the correct number of lenses\n");
        //exit(1);
        return 1;
    }
    
    xra = xra*M_PI/180.;
    xdec = xdec*M_PI/180.;
    ra     [Nfill] = xra; 
    dec    [Nfill] = xdec;
    zred   [Nfill] = xzred;
    wt     [Nfill] = xwt;
    c_lra  [Nfill] = cos(xra);
    s_lra  [Nfill] = sin(xra);
    c_ldec [Nfill] = cos(xdec);
    s_ldec [Nfill] = sin(xdec);
    datacen[3*Nfill] = c_ldec[Nfill]*c_lra[Nfill];
    datacen[3*Nfill+1] = c_ldec[Nfill]*s_lra[Nfill];
    datacen[3*Nfill+2] = s_ldec[Nfill];
    Dcom   [Nfill] = cosmo->Dcofz(double(xzred));
    if (xzred<zredmin) zredmin = xzred;
    if (xzred>zredmax) zredmax = xzred;
    Nfill++;
    if (verbose){
        fprintf(stderr, "Filled %d-th lenses with ra:%f dec:%f\n", Nfill, xra*180./M_PI, xdec*180./M_PI);
    }

    return 0;
}

int weaklens::finalize_lenses()
{
    Dcommin = cosmo->Dcofz(zredmin);
    if (Ncen>Nfill){
        Ncen=Nfill;
        ra      =(double *)  realloc(ra     , Ncen*sizeof(double));
        dec     =(double *)  realloc(dec    , Ncen*sizeof(double));
        zred    =(double *)  realloc(zred   , Ncen*sizeof(double));
        wt      =(double *)  realloc(wt     , Ncen*sizeof(double));
        datacen =(double *)  realloc(datacen, Ncen*3*sizeof(double));
        c_lra   =(double *)  realloc(c_lra  , Ncen*sizeof(double));
        s_lra   =(double *)  realloc(s_lra  , Ncen*sizeof(double));
        s_ldec  =(double *)  realloc(s_ldec , Ncen*sizeof(double));
        c_ldec  =(double *)  realloc(c_ldec , Ncen*sizeof(double));
        Dcom    =(double *)  realloc(Dcom   , Ncen*sizeof(double));
    }

    sumwl = 0.0;
    for(int i=0; i<Ncen; i++){
        sumwl+=wt[i];
    }

    // Ok now use kdtree2
    dim=3;
    xper.resize(boost::extents[dim]);
    for (int j=0; j<dim; j++)
      xper[j] = -1.0;

    // rearrange is true
    timer = new boost::timer;
    tree = new kdtree2::KDTree(datacen,xper,Ncen,dim);
    tree->sort_results = true;
    printf("Tree done %.1lf seconds\n", timer->elapsed());

    Niter=0;

    lenses_finalized=true;

    if (verbose){
        fprintf(stderr, "Finalized %d lenses, now ready to ingest source data\n", Ncen);
    }

    return 0;

}

int weaklens::test_searchrecord(){
}

int weaklens::process_source(double sra, double sdec, double se1, double se2, double swt, double smcat, double sc1_dp, double sc2_dp, double sc1_nb, double sc2_nb, double szbest, bool usepdf=false)
{

    if(std::isnan(sra) || std::isnan(sdec) || std::isnan(se1) || std::isnan(se2) || std::isnan(swt) || std::isnan(smcat) || std::isnan(sc1_dp) || std::isnan(sc2_dp) || std::isnan(sc1_nb) || std::isnan(sc2_nb) || std::isnan(szbest)){
        return 1;
    }
    if(!lenses_finalized){
        fprintf(stderr, "Lenses were not finalized, please call finalize_lenses() before processing sources\n");
        //exit(1);
        return 1;
    }
    Niter++;
    sra *= (M_PI/180.);
    sdec *= (M_PI/180.);
    se1 -= (sc1_dp+sc1_nb);
    se2 -= (sc2_dp+sc2_nb);

    double c_sdec=cos(sdec);
    double s_sdec=sin(sdec);
    double s_sra=sin(sra);
    double c_sra=cos(sra);

    double sx=c_sdec*c_sra;
    double sy=c_sdec*s_sra;
    double sz=s_sdec;
    double dCofzs=cosmo->Dcofz(double(szbest));
    double sigc_inv_fac=1.e12*4.*M_PI*gee/(cee*cee)/dCofzs;

    // Ok now get results
    kdtree2::KDTreeResultVector res;
    std::vector<double> qv(3);
    qv[0]=sx;
    qv[1]=sy;
    qv[2]=sz;
    double dis2 = pow(rmax/Dcommin, 2.);
    tree->r_nearest(qv, dis2, res);
    //0 0.034433 0.260000 435.623138 0.150020
    if (verbose)
        fprintf(stderr, "%d %f %f %f %f\n", res.size(), sqrt(dis2), szbest, Dcommin, zredmin);
    //FILE *dfile =fopen("debug.dat", "w");

    for(int k=0;k<res.size();k++){
        int xx = res.at(k).idx;
        double dis2=sqrt(res.at(k).dis);
        double logrp = log10(dis2*Dcom[xx]);

        if ((!usepdf) && szbest<=zred[xx]) continue;

        if (logrp<=logrmin || logrp>=logrmax) continue;
        //std::cout<<szbest<<" "<<zred[xx]<<std::endl;

        int rpbin = int((logrp-logrmin)/logrdiff);

        double av_sigc_inv = 0.0;
        if (!usepdf){
            av_sigc_inv = sigc_inv_fac*Dcom[xx]*(dCofzs-Dcom[xx])*(1+zred[xx]);
        }else{
            av_sigc_inv = sigmacritinverse(zred[xx]);
        }
        //fprintf(stderr, "zred:%f zbest:%f avsigcinv:%le\n", zred[xx], szbest, av_sigc_inv);


        /// Check thoroughly here
        double wls=wt[xx]*swt*pow(av_sigc_inv,2);
        if (std::isnan(wls)){
            fprintf(stderr, "zred:%f zbest:%f avsigcinv:%le\n", zred[xx], szbest, av_sigc_inv);
            continue;
        }
        double wls_by_av_sigc_inv=wt[xx]*swt*av_sigc_inv;

        ///2->s 1->l
        double cos_alps_alpl=c_sra*c_lra[xx]+s_sra*s_lra[xx];
        double sin_alps_alpl=s_sra*c_lra[xx]-c_sra*s_lra[xx];
        //double cosangle=s_ldec[xx]*s_sdec+c_ldec[xx]*c_sdec*cos_alps_alpl;
        double cosangle=sx*datacen[3*xx] + sy*datacen[3*xx+1] + sz*datacen[3*xx+2];
        double sinangle=sqrt(1-cosangle*cosangle);
        double cosphi=c_ldec[xx]*sin_alps_alpl/sinangle;
        double sinphi=(-s_ldec[xx]*c_sdec+c_ldec[xx]*s_sdec*cos_alps_alpl)/sinangle;
        
        // The following works definitely for Josh's catalog, with e2 = -e2_josh_catalog
        // This checks out with Hironao's code as well, with our sinphi = - sinphi from Hironao.
        double etan=-se1*(2*cosphi*cosphi-1)+se2*(2*sinphi*cosphi);
        double ecross=se2*(2*cosphi*cosphi-1)+se1*(2*sinphi*cosphi);
        
        if (std::isnan(etan) || std::isnan(ecross))
             continue;

        sumdsig_num[rpbin]+=(etan*wls_by_av_sigc_inv);
        sumdcross_num[rpbin]+=(ecross*wls_by_av_sigc_inv);
        sumdcount_num[rpbin]+=(wt[xx]);
        sumwls[rpbin]+=wls;
        sumwls_ms[rpbin]+=(wls*(1+smcat));
        //fprintf(dfile, "%d %.18le %.18le %.18le %.18le %.18le\n", xx, cosangle, sinangle, etan, wls, av_sigc_inv);
        //fprintf(dfile, "%d %.18le %.18le %.18le %.18le %.18le %.18le %.18le %.18le\n", xx, datacen[3*xx], datacen[3*xx+1], datacen[3*xx+2], sx, sy, sz, wls, av_sigc_inv);
        /*
       for (int kk=0;kk<rbins;kk++)
       printf("%le ", sumdsig_num[kk]);
       printf("%le %d", wt[xx], xx);
       printf("\n");
       */

    }
       //exit(11);


    // Now go through the indexing
    //if(Niter%100000==0) printf("%ld done %.1lf seconds\n", Niter, timer->elapsed());

    return 0;
}

int weaklens::finalize_results(){
    fprintf(stderr, "%e %e %e\n", logrmin, logrmax, logrdiff);
    FILE *fout=fopen(outfile, "w");

    for (int i=0;i<rbins;i++)
    {
        double rr = logrmin + (i+0.5)*logrdiff;
        double rrmin = pow(10., logrmin + i*logrdiff);
        double rrmax = pow(10., logrmin + (i+1)*logrdiff);
        fprintf(fout, "%e %le %le %le %le %le %e %le %le %le \n", rr, sumdsig_num[i], sumwls[i], sumwls_ms[i], sumdsig_num[i]/sumwls_ms[i], sumwls[i]/sumwl, rrmin/2+rrmax/2, sumdcross_num[i]/sumwls_ms[i], sumdcount_num[i]/sumwl, sumdcount_num[i]);
    }
    fclose(fout);
    return 0;
}

int weaklens::process_pofz(double *pofz, int zbins)
{

    // Free up spline if already initialized
    if (bool_init_sigcritinv){
        gsl_interp_accel_free(sigcritinv_acc);
        gsl_spline_free(sigcritinv_spline);
        bool_init_sigcritinv=false;
    }

    // This function generates a spline between pofz_zmin and pofz_zmax
#ifndef SWIG
    const int static Nspl=75;
#else
    const int Nspl=75;
#endif
    double xx[Nspl];
    double yy[Nspl];
    for (int i=0; i<Nspl; i++)
    {
        xx[i] = 0.0 + i*(zredmax+0.1)/(Nspl-1.0);

        // av_sigc_inv = av_sigc_inv_num/av_sigc_inv_den
        double av_sigc_inv_num=0.0;
        double av_sigc_inv_den=0.0;
        double chi = cosmo->Dcofz(xx[i]);
        for(int k=0;k<Npz;k++){
            av_sigc_inv_den+=pofz[k];
            if(xx[i]>pz_xz[k]) continue;
            double fac=pz_fac[k]*chi*(pz_chiz[k]-chi)*(1+xx[i]);
            av_sigc_inv_num+=pofz[k]*fac;
        }

        if (av_sigc_inv_den==0.0 || std::isnan(av_sigc_inv_den) || std::isnan(av_sigc_inv_num)){
            yy[i] = 0.0;
	}else{
            yy[i] = av_sigc_inv_num/av_sigc_inv_den;
	}
    }
    sigcritinv_acc = gsl_interp_accel_alloc ();
    sigcritinv_spline = gsl_spline_alloc (gsl_interp_cspline, Nspl);
    bool_init_sigcritinv = true;
    gsl_spline_init (sigcritinv_spline, xx, yy, Nspl);
    
    if(verbose){
    std::cout<<"# "<<"The spline for sigmacritinverse is now initialized"<<std::endl;
    }

    /*
    for(int i=0; i<zbins; i++)
        std::cout<<pofz[i]<<" ";
    std::cout<<std::endl;
    // Now go through the indexing
    */
    //Niter++;
    //if(Niter%100000==0) printf("%ld done %.1lf seconds\n", Niter, timer->elapsed());
    return 0;
}

double weaklens::sigmacritinverse(double zlens)
{
    if (!bool_init_sigcritinv){
        std::cout<<"Spline for sigmacritinverse needs to be initialized before use"<<std::endl;
        exit(11);
    }
    
    double result = 0.0;
    try{
        result = gsl_spline_eval (sigcritinv_spline, zlens, sigcritinv_acc);
    }
    catch (std::exception& e){
        result = 0.0;
    }

    if (result<0.0){
        result = 0.0;
    }

    return result;
}

int weaklens::setup_pofz(double pofz_zmin, double pofz_zdiff, int xNpz)
{
    Npz = xNpz;
    pz_xz = new double [Npz];
    pz_chiz = new double [Npz];
    pz_fac = new double [Npz];
    
    for(int kk=0;kk<Npz;kk++){
        pz_xz[kk]=pofz_zmin+kk*pofz_zdiff;
        pz_chiz[kk]=cosmo->Dcofz(pz_xz[kk]);
        pz_fac[kk]=1.e12*4.*M_PI*gee/(cee*cee)/(pz_chiz[kk]);
    }
    pz_alloc=true;

    return 1;
}
