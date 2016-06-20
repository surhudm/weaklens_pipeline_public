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

    if (verbose){
        fprintf(stderr, "Initialized cosmology, and delta sigma sum arrays\n");
        fprintf(stderr, "Options passed rmin:%f rmax:%f rbins:%d outfile:%s \n", rmin, rmax, rbins, outfile);
    }

}

weaklens::weaklens(float xrmin, float xrmax, int xrbins, char *xoutfile, bool xverbose)
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
}

int weaklens::allocate_lens_memory(int xNcen)
{
    Ncen = xNcen;
    ra      =(float *)  calloc(Ncen,sizeof(float));
    dec     =(float *)  calloc(Ncen,sizeof(float));
    zred    =(float *)  calloc(Ncen,sizeof(float));
    wt      =(float *)  calloc(Ncen,sizeof(float));
    datacen =(float *)  calloc(Ncen*3,sizeof(float));
    c_lra   =(float *)  calloc(Ncen,sizeof(float));
    s_lra   =(float *)  calloc(Ncen,sizeof(float));
    s_ldec  =(float *)  calloc(Ncen,sizeof(float));
    c_ldec  =(float *)  calloc(Ncen,sizeof(float));
    Dcom    =(float *)  calloc(Ncen,sizeof(float));
    Nfill=0;
    lenses_alloc=true;
    if (verbose){
        fprintf(stderr, "Initialized lenses Ncen=%d\n", Ncen);
    }
    return 0;
}

int weaklens::process_lens(float xra, float xdec, float xzred, float xwt)
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
        ra      =(float *)  realloc(ra     , Ncen*sizeof(float));
        dec     =(float *)  realloc(dec    , Ncen*sizeof(float));
        zred    =(float *)  realloc(zred   , Ncen*sizeof(float));
        wt      =(float *)  realloc(wt     , Ncen*sizeof(float));
        datacen =(float *)  realloc(datacen, Ncen*3*sizeof(float));
        c_lra   =(float *)  realloc(c_lra  , Ncen*sizeof(float));
        s_lra   =(float *)  realloc(s_lra  , Ncen*sizeof(float));
        s_ldec  =(float *)  realloc(s_ldec , Ncen*sizeof(float));
        c_ldec  =(float *)  realloc(c_ldec , Ncen*sizeof(float));
        Dcom    =(float *)  realloc(Dcom   , Ncen*sizeof(float));
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

int weaklens::process_source(float sra, float sdec, float se1, float se2, float swt, float smcat, float sc1_dp, float sc2_dp, float sc1_nb, float sc2_nb, float szbest)
{
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
    double gee=4.2994e-9;
    double cee=3.E5;
    double sigc_inv_fac=1.e12*4.*M_PI*gee/(cee*cee)*(1.+szbest)/dCofzs;

    // Ok now get results
    kdtree2::KDTreeResultVector res;
    std::vector<float> qv(3);
    qv[0]=sx;
    qv[1]=sy;
    qv[2]=sz;
    float dis2 = pow(rmax/Dcommin, 2.);
    tree->r_nearest(qv, dis2, res);
    //0 0.034433 0.260000 435.623138 0.150020
    if (verbose)
        fprintf(stderr, "%d %f %f %f %f\n", res.size(), sqrt(dis2), szbest, Dcommin, zredmin);
    //FILE *dfile =fopen("debug.dat", "w");

    for(int k=0;k<res.size();k++){
        int xx = res.at(k).idx;
        float dis2=sqrt(res.at(k).dis);
        float logrp = log10(dis2*Dcom[xx]);

        if (szbest<=zred[xx]) continue;

        if(logrp<=logrmin || logrp>=logrmax) continue;
        //std::cout<<szbest<<" "<<zred[xx]<<std::endl;

        int rpbin=int((logrp-logrmin)/logrdiff);

        double av_sigc_inv= sigc_inv_fac*Dcom[xx]*(dCofzs-Dcom[xx])/(1.+szbest)*(1+zred[xx]);

        /// Check thoroughly here
        double wls=wt[xx]*swt*pow(av_sigc_inv,2);

        ///2->s 1->l
        double cos_alps_alpl=c_sra*c_lra[xx]+s_sra*s_lra[xx];
        double sin_alps_alpl=s_sra*c_lra[xx]-c_sra*s_lra[xx];
        //double cosangle=s_ldec[xx]*s_sdec+c_ldec[xx]*c_sdec*cos_alps_alpl;
        double cosangle=sx*datacen[3*xx] + sy*datacen[3*xx+1] + sz*datacen[3*xx+2];
        double sinangle=sqrt(1-cosangle*cosangle);
        double cosphi=c_ldec[xx]*sin_alps_alpl/sinangle;
        double sinphi=(-s_ldec[xx]*c_sdec+c_ldec[xx]*s_sdec*cos_alps_alpl)/sinangle;

        double etan=-se1*(2*cosphi*cosphi-1)+se2*(2*sinphi*cosphi);
        double ecross=se2*(2*cosphi*cosphi-1)+se1*(2*sinphi*cosphi);

        sumdsig_num[rpbin]+=(etan*wls/av_sigc_inv);
        sumdcross_num[rpbin]+=(ecross*wls/av_sigc_inv);
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
    if(Niter%100000==0) printf("%ld done %.1lf seconds\n", Niter, timer->elapsed());

    return 0;
}

int weaklens::finalize_results(){
    fprintf(stderr, "%e %e %e\n", logrmin, logrmax, logrdiff);
    FILE *fout=fopen(outfile, "w");

    for (int i=0;i<rbins;i++)
    {
        float rr = logrmin + (i+0.5)*logrdiff;
        float rrmin = pow(10., logrmin + i*logrdiff);
        float rrmax = pow(10., logrmin + (i+1)*logrdiff);
        fprintf(fout, "%e %le %le %le %le %le %e %le %le \n", rr, sumdsig_num[i], sumwls[i], sumwls_ms[i], sumdsig_num[i]/sumwls_ms[i], sumwls[i]/sumwl, rrmin/2+rrmax/2, sumdcross_num[i]/sumwls_ms[i], sumdcount_num[i]/sumwl);
    }
    fclose(fout);
    return 0;
}
