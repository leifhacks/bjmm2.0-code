#define scale 1000000

double H(double x) {
	//if (x<0.0 && x > -0.000001) x = 0.0;
    if (x<0.0 || x>1.0) {
        printf("-H Error %f -", x);
        return 0.0;
    }
    if (x==0.0) return 0;
    if (x==1.0) return 0;
    return (-x*log(x)-(1-x)*log(1-x))/log(2);
}

double inverse(double H1[], double y) {
    int left=0, right=scale, mid;
    while(left < right-1) {
        mid= (left+right)/2;
        //printf("%d %d %d \n",left, right, mid);
        if(y < H1[mid]) right=mid-1;
        else left=mid;
    }
    
    //printf("%f \n",fabs(y-H(left/(2*(double)scale))));
    if(fabs(y-H(left/(2*(double)scale))) < fabs(y-H(right/(2*(double)scale))) ) return(left/(2*(double)scale));
    else return (right/(2*(double)scale));
}

double NN(double gamma, double lambda, double H1[]) {
    if (gamma<0.0 || gamma>=1.0 || lambda<0.0) {
        printf("-NN Error %f %f -", gamma, lambda);
        return 2*lambda;
    }
    if (gamma==0.0) {
        return lambda;
    }
    if (gamma < 0.5 && lambda < 1-H(gamma/2) ) 
	   return (1-gamma) * (1 - H( (inverse(H1, (1-lambda) ) - gamma/2 ) / (1 - gamma) ) );
    else return fmin(2*lambda, lambda+H(gamma));
}

double Optimize(double k, double w, int depth, double pmin, double pmax, double psteps, double lmin, double lmax, double lsteps, double emin[5], double emax[5], double esteps[5], int NNflag, int HDflag, double H1[], double MINvalues[24], double *T) {
     
    if (w==0) {
        w = inverse(H1, 1-k);
        if (HDflag) w=w/2;
    }
    int dmin=2, dmax=5, dsteps=1;
    if (depth>0) {
        dmin=depth;
        dmax=depth;
    }
    double S[6], C[6];
    double maximum, Smax;
    double l,invprob;
    double p[6], eps[6], r[6], epsswap;
    *T = 1.0;
 
    for (int d=dmin; d<=dmax; d+=dsteps) {
        for (p[0]=pmin; p[0]<pmax; p[0] = p[0] + psteps) {
            for (l=lmin; l<lmax; l = l + lsteps) {
                invprob = H(w) - H(p[0]/(k+l)) * (k+l) - H((w-p[0])/(1-k-l)) * (1-k-l);
                for (eps[1]=emin[1]; eps[1]<=emax[1]; eps[1] = eps[1] + esteps[1]) {
                    for (eps[2]=emin[2]; eps[2]<=emax[2]; eps[2] = eps[2] + esteps[2]) {
                        for (eps[3]=emin[3]; eps[3]<=emax[3]; eps[3] = eps[3] + esteps[3]) {
                            for (eps[4]=emin[4]; eps[4]<=emax[4]; eps[4] = eps[4] + esteps[4]) {
                                epsswap = eps[d], eps[d] = 0.0;
                                r[d] = 0, r[0] = l;
                                for (int i=0; i<d; i++) {
                                    p[i+1] = p[i]/2 + eps[i+1];
                                }
                                for(int i=1; i<d; i++) {
                                    r[i] = p[i-1] + H(eps[i]/(k+l-p[i-1])) * (k+l-p[i-1]);
                                    S[i] = H(p[i]/(k+l)) * (k+l) - r[i];
                                    if ((i>1) || (NNflag==0)) C[i] = 2*S[i] + r[i] - r[i-1];
                                    else                      C[i] = NN((w-p[i-1])/(1-k-l), S[i]/(1-k-l), H1) * (1-k-l);
                                    //printf("r[%d] = %f, S[%d] = %f, C[%d] = %f \n", i, r[i], i, S[i], i, C[i]);
                                }
                                if (l >= r[1]) {
                                    S[d] = H(2*p[d]/(k+l))*(k+l)/2;
                                    C[d] = 2*S[d] + r[d] - r[d-1];
                                    //printf("S[%d] = %f, C[%d] = %f \n", d, S[d], d, C[d]);
                                    maximum = 0.0;
                                    Smax = 0.0;
                                    for (int i=1; i<=d; i++) {Smax = fmax(Smax, S[i]);}
                                    for (int i=1; i<=d; i++) {maximum = fmax(maximum, C[i]);}
                                    maximum = fmax(maximum,Smax);
                                    if (*T > (maximum + invprob)) {
                                        *T = maximum + invprob;
                                        printf("depth=%d k=%.3f w=%f p=%.4f l=%.4f ", d, k, w, p[0], l);
                                        for(int i=1; i<d; i++) printf("eps[%d]=%.4f ", i, eps[i]);
                                        printf("T=%.5f \n", *T);
                                        MINvalues[0]=k; MINvalues[1]=w; MINvalues[2]=l; // for d=3
                                        for(int i=1; i<d; i++) MINvalues[i+2]=r[i]; // 3,4
                                        for(int i=0; i<=d; i++) MINvalues[d+2+i]=p[i]; // 5,6,7,8
                                        for(int i=1; i<=d; i++) MINvalues[2*d+2+i]=S[i]; // 9,10,11
                                        for(int i=1; i<=d; i++) MINvalues[3*d+2+i]=C[i]; // 12,13,14
                                        MINvalues[4*d+3]= *T;                             // 15
                                    }
                                }
                                eps[d] = epsswap;
                            }
                        }
                    }
                }
            }
        }
    }
 
    return 0.0;
 
}