/**
	•	@file burst_evp.c
	•	@brief Simulation of two-phase elastoviscoplastic fluid flow using the log-conformation method.

	•	This simulation models the behavior of a elastoviscoplastic fluid using the Navier-Stokes equations coupled with the Saramito model.
	
  The code employs state-of-the-art techniques, including adaptive mesh refinement and the log-conformation method, to efficiently simulate complex fluid dynamics.
*/


#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

// we modify the [log-conform.h](http://basilisk.fr/src/log-conform.h) and [fene-p.h](http://basilisk.fr/src/fene-p.h) files to implement the Saramito model:
#include "log-conform-EVP.h"
#include "saramito-EVP.h"

#include "distance.h"
#include "adapt_wavelet_limited.h"

#define tmax 4.5
#define LEVEL 8
#define MAXlevel 11
#define DT_MAX 0.0005
#define Ldomain 8

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f VOF
#define KErr (1e-4)                                 // error tolerance in curvature
#define VelErr (1e-2)                               // error tolerances in velocity, A
#define OmegaErr (1e-3)                             // error tolerances in vorticity

#define tsnap (0.005)

# define B 0.5 // solvent to total viscosity ratio

scalar mupv[], lambdav[], tau0v[];

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

double Bond, J, Deb;
char nameOut[80], namepng[80], dumpFile[80];

int main(int argc, char const *argv[]) {
    
L0 = Ldomain;
origin (-L0/2., 0.);
init_grid (1 << 8);
Bond = 0.001;

J = atof(argv[1]); // Plastocapillary number
Deb = atof(argv[2]); // Deborah number

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

TOLERANCE = 1e-5;

run();
}

// gravity is added as a body force
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= Bond;
}

// a region-based refinement criterion is used to refine very close to the axis of symmetry. 
int refRegion(double x, double y, double z){
  return (y < 1.28 ? MAXlevel+2 : y < 2.56 ? MAXlevel+1 : y < 5.12 ? MAXlevel : MAXlevel-1);
}

event init (t = 0) {
  if (!restore (file = dumpFile)){
    // read the initial shape from a data file.
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

  scalar KAPPA[], Axx[], Axy[], Ayy[], Aqq[];
  curvature(f, KAPPA);
  foreach()
    {
        Axx[] = f[]*(((1.-B)*0.01)/Deb)*tau_p.x.x[];
        Axy[] = f[]*(((1.-B)*0.01)/Deb)*tau_p.x.y[];
        Ayy[] = f[]*(((1.-B)*0.01)/Deb)*tau_p.y.y[];
        Aqq[] = f[]*(((1.-B)*0.01)/Deb)*tau_qq[];
    }
  adapt_wavelet_limited ((scalar *){f, u.x, u.y, Axx, Axy, Ayy, Aqq, trA, solidreg, KAPPA},
     (double[]){fErr, VelErr, VelErr, VelErr, VelErr, VelErr, VelErr, fErr, fErr, KErr},
     refRegion);
 }

event writingFiles (t = 0; t += tsnap; t <= tmax) {
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

event end (t = end) {
  fprintf(ferr, "Done: \n");
}

// logging on the run data
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
