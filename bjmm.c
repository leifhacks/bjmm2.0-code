#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include <time.h>

int main(int argc, char *argv[]) {
    // running time counter
    clock_t tStart = clock();
    
    // Compute LUT for H
    int scale=1000000;
    double H1[scale+1];
    for (int i=0; i<=scale; i++) {
        H1[i]=H(i/(2*(double)scale));
        //printf("%d %f %.20f \n",i, i/(2*(double)scale), H1[i]);
    }
    
    // initialization
    int depth = 0;
    int NNflag = 0;
    int HDflag = 0;
    int McEflag = 0;
    double emin[5], emax[5], esteps[5];
    for (int i=0;i<5;i++) {emin[i] = 0.0, emax[i] = 0.0, esteps[i] = 1.0;}
    double pmin, pmax, psteps;
    double lmin, lmax, lsteps;
    double kmin, kmax, ksteps;
    double McElieceD, McElieceR;
    double k=0, w=0;
    double Tmax=0.0, Tmin=1.0, MAXvalues[24], MINvalues[24];
    
    // FDD, BJMM-NN, m=5
    //depth=5;
    //NNflag = 1;
    //kmin=0.419, kmax=0.421, ksteps=0.001;
    //pmin=0.0810, pmax=0.0830, psteps=0.0001;
    //lmin=0.2622, lmax=0.2642, lsteps=0.0001;
    //emin[1]=0.0312, emax[1]=0.0332, esteps[1]=0.0001;
    //emin[2]=0.0146, emax[2]=0.0166, esteps[2]=0.0001;
    //emin[3]=0.0030, emax[3]=0.0050, esteps[3]=0.0001;
    //emin[4]=0.0000, emax[4]=0.0010, esteps[4]=0.0001;
    
    // FDD, BJMM-NN, m=4
    //depth=4;
    //NNflag = 1;
    //kmin=0.423, kmax=0.423, ksteps=0.001;
    //pmin=0.0800, pmax=0.0850, psteps=0.0001;
    //lmin=0.2610, lmax=0.2660, lsteps=0.0001;
    //emin[1]=0.0310, emax[1]=0.0330, esteps[1]=0.0001;
    //emin[2]=0.0140, emax[2]=0.0170, esteps[2]=0.0001;
    //emin[3]=0.0030, emax[3]=0.0050, esteps[3]=0.0001;
    
    // HDD, BJMM-NN, m=4
    //depth=4;
    //NNflag = 1;
    //HDflag=1;
    //kmin=0.475, kmax=0.476, ksteps=0.001;
    //pmin=0.0150, pmax=0.0200, psteps=0.0001;
    //lmin=0.0630, lmax=0.0680, lsteps=0.0001;
    //emin[1]=0.0050, emax[1]=0.0070, esteps[1]=0.0001;
    //emin[2]=0.0005, emax[2]=0.0025, esteps[2]=0.0001;
    //emin[3]=0.0000, emax[3]=0.0010, esteps[3]=0.0001;
    
    // McE, BJMM-NN, m=4
//    depth=4;
//    NNflag = 1;
//    McEflag=1;
//    pmin=0.0060, pmax=0.0160, psteps=0.0001;
//    lmin=0.0350, lmax=0.0450, lsteps=0.0001;
//    emin[1]=0.0010, emax[1]=0.0050, esteps[1]=0.0001;
//    emin[2]=0.0000, emax[2]=0.0020, esteps[2]=0.0001;
//    emin[3]=0.0000, emax[3]=0.0010, esteps[3]=0.0001;
//    McElieceD=0.02, McElieceR=0.775;
    
    // McE, BJMM-NN, m=4
//    depth=4;
//    McEflag=1;
//    pmin=0.0050, pmax=0.0150, psteps=0.0001;
//    lmin=0.0580, lmax=0.0680, lsteps=0.0001;
//    emin[1]=0.0007, emax[1]=0.0047, esteps[1]=0.0001;
//    emin[2]=0.0000, emax[2]=0.0010, esteps[2]=0.0001;
//    emin[3]=0.0000, emax[3]=0.0010, esteps[3]=0.0001;
//    McElieceD=0.02, McElieceR=0.775;
    
    // FDD, BJMM-NN, m=3
    depth=3;
    NNflag = 1;
    kmin=0.420, kmax=0.424, ksteps=0.001;
    pmin=0.0600, pmax=0.0700, psteps=0.0001;
    lmin=0.1800, lmax=0.2000, lsteps=0.0001;
    emin[1]=0.0180, emax[1]=0.0220, esteps[1]=0.0001;
    emin[2]=0.0030, emax[2]=0.0050, esteps[2]=0.0001;
    
    // HDD, BJMM-NN, m=3
    //depth=3;
    //NNflag = 1;
    //HDflag=1;
    //kmin=0.473, kmax=0.475, ksteps=0.001;
    //pmin=0.0150, pmax=0.0200, psteps=0.0001;
    //lmin=0.0640, lmax=0.0690, lsteps=0.0001;
    //emin[1]=0.0030, emax[1]=0.0090, esteps[1]=0.0001;
    //emin[2]=0.0000, emax[2]=0.0050, esteps[2]=0.0001;
    
    // McE, BJMM-NN, m=3
//    depth=3;
//    NNflag = 1;
//    McEflag=1;
//    pmin=0.0050, pmax=0.0150, psteps=0.0001;
//    lmin=0.0370, lmax=0.0480, lsteps=0.0001;
//    emin[1]=0.0010, emax[1]=0.0050, esteps[1]=0.0001;
//    emin[2]=0.0000, emax[2]=0.0020, esteps[2]=0.0001;
//    McElieceD=0.02, McElieceR=0.775;
    
    // McE, BJMM, m=3
//    depth=3;
//    McEflag=1;
//    pmin=0.0000, pmax=0.0180, psteps=0.0001;
//    lmin=0.0420, lmax=0.0620, lsteps=0.001;
//    emin[1]=0.0000, emax[1]=0.0050, esteps[1]=0.0001;
//    emin[2]=0.0000, emax[2]=0.0020, esteps[2]=0.0001;
//    McElieceD=0.02, McElieceR=0.775;

    
    // FDD, BJMM-NN, m=2
//    depth=2;
//    NNflag = 1;
//    kmin=0.410, kmax=0.430, ksteps=0.001;
//    pmin=0.0380, pmax=0.0440, psteps=0.0001;
//    lmin=0.0910, lmax=0.0970, lsteps=0.0001;
//    emin[1]=0.0040, emax[1]=0.0100, esteps[1]=0.0001;
    
    // HDD, BJMM-NN, m=2
    //depth=2;
    //NNflag = 1;
    //HDflag=1;
    //kmin=0.455, kmax=0.461, ksteps=0.0001;
    //pmin=0.0070, pmax=0.0130, psteps=0.0001;
    //lmin=0.0260, lmax=0.0320, lsteps=0.0001;
    //emin[1]=0.0000, emax[1]=0.0050, esteps[1]=0.0001;
    
    // McE, BJMM-NN, m=2
//    depth=2;
//    NNflag = 1;
//    McEflag=1;
//    pmin=0.0030, pmax=0.0130, psteps=0.0001;
//    lmin=0.0100, lmax=0.0300, lsteps=0.0001;
//    emin[1]=0.0000, emax[1]=0.0040, esteps[1]=0.0001;
//    McElieceD=0.02, McElieceR=0.775;
    
    // McE, BJMM, m=2
//    depth=2;
//    McEflag=1;
//    pmin=0.0000, pmax=0.0145, psteps=0.0001;
//    lmin=0.0200, lmax=0.0800, lsteps=0.0001;
//    emin[1]=0.0000, emax[1]=0.0040, esteps[1]=0.0001;
//    McElieceD=0.02, McElieceR=0.775;
    
    
    // McEliece settings
    if (McEflag) {
        kmin=McElieceR;
        kmax=McElieceR;
        ksteps=0.0001;
        w=McElieceD;
    }

    // Optimization for every k
    for (k=kmin; k<=kmax; k = k + ksteps) {
        printf("k = %.4f \n======== \n", k);
        Optimize(k,w,depth,pmin,pmax,psteps,lmin,lmax,lsteps,emin,emax,esteps,NNflag,HDflag,H1,scale,MINvalues,&Tmin);
        if (Tmax < Tmin) {
            Tmax=Tmin;
            for (int i=0; i<24; i++) MAXvalues[i]=MINvalues[i];
        }
    }
    
    // Output
    printf("Tiefe %d ", depth);
    if (HDflag) printf("HD "); else if (McEflag) printf("McE "); else printf("FD ");
    if (NNflag) printf("NN-BJMM \n"); else printf("Plain BJMM \n");
    printf(" k=%.3f w=%f l=%.4f ", MAXvalues[0], MAXvalues[1], MAXvalues[2]);
    for(int i=1; i<depth; i++) printf("R[%d]=%.4f ", i, MAXvalues[i+2]);
    for(int i=0; i<=depth; i++) printf("p[%d]=%.4f ", i, MAXvalues[depth+2+i]);
    for(int i=1; i<=depth; i++) printf("S[%d]=%.4f ", i, MAXvalues[2*depth+3+i]);
    for(int i=1; i<=depth; i++) printf("C[%d]=%.4f ", i, MAXvalues[3*depth+2+i]);
    printf("T=%.5f \n", MAXvalues[4*depth+3]);
    
    // Output time
    printf("\n Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    
    
    // Additional values
    //printf("\n");
    //Optimize(0.4,0,0,pmin,pmax,psteps,lmin,lmax,lsteps,emin,emax,esteps,NNflag,HDflag,H1,scale,MAXvalues);
    
    return 0;

}