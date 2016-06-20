#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include <iostream>
#include <fstream>
#include "kdtree2.h"
#include "weaklens.h"
#include <vector>
#include <ctime>
#include <cstring>
#include <cmath>
typedef boost::multi_array<float,1> array1dfloat;
#include <boost/timer.hpp>

int main(int argc, char** argv)
{
    // ========================================================
    // Read galaxy file
    // ========================================================
    FILE *fin;
    fin=fopen("rcsrestrict_dr12_v5_lowz_All_catalog.dat","r");

    float zmin=0.15;
    float zmax=0.43;
    float rmin=0.1;
    float rmax=15.0;
    int rbins=15;
    char outfile[10000];
    sprintf(outfile, "debug_test.dat");
    weaklens wpipe(rmin, rmax, rbins, outfile, true);

    int Ncen=1000000;
    float *ra, *dec, *zred, *wt;
    ra      =(float *)  calloc(Ncen,sizeof(float));
    dec     =(float *)  calloc(Ncen,sizeof(float));
    zred    =(float *)  calloc(Ncen,sizeof(float));
    wt      =(float *)  calloc(Ncen,sizeof(float));
    int fill=0;
    float xra, xdec, xzred, xlam, xwt;
    while (fscanf(fin, "%e %e %e %e %e", &xra, &xdec, &xzred, &xlam, &xwt) == 5){
        if(xzred<=zmin||xzred>zmax) continue;
        ra     [fill] = xra; 
        dec    [fill] = xdec;
        zred   [fill] = xzred;
        wt     [fill] = xwt;
        fill++;
    }
    Ncen = fill;
    ra      =(float *)  realloc(ra     , Ncen*sizeof(float));
    dec     =(float *)  realloc(dec    , Ncen*sizeof(float));
    zred    =(float *)  realloc(zred   , Ncen*sizeof(float));
    wt      =(float *)  realloc(wt     , Ncen*sizeof(float));

    wpipe.allocate_lens_memory(Ncen);

    for (int i=0; i<Ncen; i++){
        wpipe.process_lens(ra[i], dec[i], zred[i], wt[i]);
    }
    wpipe.finalize_lenses();
    wpipe.test_searchrecord();
}

