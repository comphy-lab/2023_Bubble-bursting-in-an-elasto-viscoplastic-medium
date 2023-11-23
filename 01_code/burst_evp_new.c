#define DT_MAX 0.0005
#define Ldomain 8

# include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "log-conform-EVP.h"
#include "fene-p-EVP.h"
#include "distance.h"
#include "adapt_wavelet_limited.h"

#define tmax 4.5
#define LEVEL 8
#define MAXlevel 11

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-4)                                 // error tolerance in f2 VOF
#define VelErr (1e-2)                               // error tolerances in velocity
#define OmegaErr (1e-3)                             // error tolerances in vorticity

#define tsnap (0.005)

# define B 0.5 // solvent to total viscosity ratio

scalar mupv[], lambdav[], tau0v[];


//u.t[right] = dirichlet(0);
//u.t[left]  = dirichlet(0);

//u.t[top] = dirichlet(0);
//uf.n[bottom] = 0.;

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

double Bond,J,Deb;
char nameOut[80], namepng[80], dumpFile[80];
//face vector D2f[];

int main(int argc, char const *argv[]) {
    
L0 = Ldomain;
origin (-L0/2., 0.);
init_grid (1 << 8);
Bond = 0.001;

J = atof(argv[1]);
Deb = atof(argv[2]);

char comm[80];
sprintf (comm, "mkdir -p intermediate");
system(comm);
sprintf (comm, "mkdir -p 01_pp/png");
system(comm);
sprintf (comm, "mkdir -p 01_pp/pdf");
system(comm);
 
sprintf (dumpFile, "dump");

rho1 = 1., mu1 = 0.01*B;
rho2 = 0.001, mu2 = 0.0002, f.sigma = 1.0;

fprintf(ferr, "J %4.1f De %4.1f \n", J, Deb);

lambda=lambdav;
mup=mupv;
tau0=tau0v;

TOLERANCE = 1e-4;

run();
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= Bond;
}

int refRegion(double x, double y, double z){
  return (y < 1.28 ? MAXlevel+2 : y < 2.56 ? MAXlevel+1 : y < 5.12 ? MAXlevel : MAXlevel-1);
}

event init (t = 0) {
  if (!restore (file = dumpFile)){

    char filename[60];
    sprintf(filename,"Bo%5.4f.dat",Bond);
    FILE * fp = fopen(filename,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", filename);
      return 1;
    }
    coord* InitialShape;
    InitialShape = input_xy(fp);
    fclose (fp);
    scalar d[];
    distance (d, InitialShape);
    while (adapt_wavelet_limited ((scalar *){f, d}, (double[]){1e-8, 1e-8}, refRegion).nf);
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    fractions (phi, f);
  }
}

event properties (i++) {
  foreach() {

    mupv[] = (1. - B)*f[]*mu1/B;
    lambdav[] = Deb*f[];
    tau0v[] = J*f[];
  }
  boundary ({lambdav, mupv, tau0v});
}

event adapt(i++){

  scalar KAPPA[], omega[], D2[];
    face vector D2f[];
  curvature(f, KAPPA);
  adapt_wavelet_limited ((scalar *){f, u.x, u.y, tau_p.x.x, tau_p.x.y, tau_p.y.y, trA, solidreg, KAPPA},
     (double[]){fErr, VelErr, VelErr, VelErr, VelErr, VelErr, fErr, fErr, KErr},
     refRegion);
/*
foreach_face(x) {
  double D11 = 0.5*( (u.y[0,1] - u.y[0,-1] + u.y[-1,1] - u.y[-1,-1])/(2.*Delta) );
  double D22 = (u.y[] + u.y[-1, 0])/(2*max(y, 1e-20));
  double D33 = (u.x[] - u.x[-1,0])/Delta;
  double D13 = 0.5*( (u.y[] - u.y[-1, 0])/Delta + 0.5*( (u.x[0,1] - u.x[0,-1] + u.x[-1,1] - u.x[-1,-1])/(2.*Delta) ) );

  double D2temp = sqrt( sq(D11) + sq(D22) + sq(D33) + 2*sq(D13) );
  D2f.x[] = D2temp;
}

foreach_face(y) {
  double D11 = (u.y[0,0] - u.y[0,-1])/Delta;
  double D22 = (u.y[0,0] + u.y[0,-1])/(2*max(y, 1e-20));
  double D33 = 0.5*( (u.x[1,0] - u.x[-1,0] + u.x[1,-1] - u.x[-1,-1])/(2.*Delta) );
  double D13 = 0.5*( (u.x[0,0] - u.x[0,-1])/Delta + 0.5*( (u.y[1,0] - u.y[-1,0] + u.y[1,-1] - u.y[-1,-1])/(2.*Delta) ) );

  double D2temp = sqrt( sq(D11) + sq(D22) + sq(D33) + 2*sq(D13) );
  D2f.y[] = D2temp;
}

 foreach(){
  D2[] = (D2f.x[]+D2f.y[]+D2f.x[1,0]+D2f.y[0,1])/4.;
  if (D2[] > 0.){
    D2[] = log(D2[])/log(10);
  } else {
    D2[] = -10;
  }
}
boundary({D2,D2f});

//  adapt_wavelet_limited ((scalar *){f, u.x, u.y, KAPPA, D2},
//     (double[]){fErr, VelErr, VelErr, KErr, 1e-3},
//     refRegion);
*/
 }

event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

event end (t = end) {
  fprintf(ferr, "Done: \n");
}

event writedt (i++)
{ static FILE * fp;
  if (i == 0) {
  fp = fopen("timestep.txt","w");
  fprintf(fp, "i: %d dt: %g n: %g INT_MAX: %d\n",i,dt,(tnext-t)/dt,INT_MAX);
  fclose(fp);
              }
  else {
  fp = fopen("timestep.txt","a");
  fprintf(fp,"i: %d dt: %g n: %g INT_MAX: %d\n",i,dt,(tnext-t)/dt,INT_MAX);
  fclose(fp); }
 }
   
event logWriting (i+=100) {
  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += (2*pi*y)*(0.5*(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }
  static FILE * fp;
  if (i == 0) {
    fprintf (ferr, "i dt t ke\n");
    fp = fopen ("log", "w");
    fprintf (fp, "i dt t ke\n");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  else {
    fp = fopen ("log", "a");
    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
  }
  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
  if (ke > 1e3 || ke < 1e-6){
    if (i > 1e2){
      return 1;
    }
  }
}
