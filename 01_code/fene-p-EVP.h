/**
# Functions $f_s$ and $f_r$ for the FENE-P model

See [log-conform.h](log-conform.h). */

//#include <iostream>
//#include <fstream>

double L2 = 1.;
double myeps=1e-6;
double maxis;

static void fenep_r (double trA, double tau_pxx, double tau_pxy, double tau_pyy, double tau_qq, double tau00, double * nu, double * eta) {
//   printf("fenep %g %g %g\n",tau_pxx,tau_pxy,tau_pyy);
   double t1=tau_pxx;
   double t2=tau_pxy;
   double t3=tau_pyy;
   double t4=tau_qq;

//   double tauD=sqrt(0.5*(t1*t1+t3*t3+2.0*t2*t2));
   double tauD=sqrt((1.0/6.0)*((t1-t3)*(t1-t3)+(t3-t4)*(t3-t4)+(t4-t1)*(t4-t1))+t2*t2);
   double solid = max(0.0,(tauD-tau00)/(tauD+myeps));
   *eta = solid;
   maxis=0.0;
   
//   for(int i = 1; i<N*N; i++)
//     {
//          if(eta[i] > maxis)
//                maxis = eta[i];
//	       	  };

//   FILE * fp=fopen("myfile", "w");
//   fprintf(fp, "%g", maxis);
//   fflush(fp);


   //*eta = 1.;
   *nu = 1.;
  return;
} 

static void fenep_s (double trA, double tau_pxx, double tau_pxy, double tau_pyy, double tauqq, double tau00, double * nu, double * eta) {
   *eta = 1.;
   *nu = 1.;
  return;
}

event defaults (i = 0) {
  f_s = fenep_s;
  f_r = fenep_r;
}

event init (i = 0) {
#if AXI
  double dim = 3;
#else
  double dim = dimension;
#endif  
  scalar trac = trA;
  foreach()
    trac[] = dim*L2/(dim + L2);
}
