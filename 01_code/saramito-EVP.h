/**
# Functions $f_s$ and $f_r$ for the Saramito EVP model (Axi-symmetric case)

See [log-conform-EVP.h](log-conform-EVP.h). */

double L2 = 1.;
double myeps=1e-6; //Tolerance

static void saramito_r (double trA, double tau_pxx, double tau_pxy, double tau_pyy, double tau_qq, double tau00, double * nu, double * eta) {
   double t1=tau_pxx;
   double t2=tau_pxy;
   double t3=tau_pyy;
   double t4=tau_qq;    // Axi-symmetric case

   double tauD=sqrt((1.0/6.0)*((t1-t3)*(t1-t3)+(t3-t4)*(t3-t4)+(t4-t1)*(t4-t1))+t2*t2);
   double solid = max(0.0,(tauD-tau00)/(tauD+myeps));
   *eta = solid;  // Switch term
   *nu = 1.;
  return;
} 

static void saramito_s (double trA, double tau_pxx, double tau_pxy, double tau_pyy, double tauqq, double tau00, double * nu, double * eta) {
   *eta = 1.;
   *nu = 1.;
  return;
}

event defaults (i = 0) {
  f_s = saramito_s;
  f_r = saramito_r;
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
