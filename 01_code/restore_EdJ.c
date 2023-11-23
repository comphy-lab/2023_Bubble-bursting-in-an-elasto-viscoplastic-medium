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
#define MAXlevel 12

// Error tolerancs
#define fErr (1e-3)                                 // error tolerance in f1 VOF
#define KErr (1e-4)                                 // error tolerance in f2 VOF
#define VelErr (1e-2)                               // error tolerances in velocity
#define OmegaErr (1e-3)                             // error tolerances in vorticity

scalar mupv[], lambdav[], tau0v[];

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

double Bond,J,Deb,timew,B,Oh;
char nameOut[80], namepng[80], dumpFile[80], nameIn[80];

int  main(int argc, char const *argv[]) {
  sprintf(nameIn, "%s", argv[1]);
  timew = atof(argv[2]);
  L0 = Ldomain;
  origin (-L0/2., 0.);
  init_grid (1 << 8);
  Bond = 0.001;
  B = atof(argv[3]);
  J = atof(argv[4]);
  Deb = atof(argv[5]);
  
  Oh = 0.01;
  rho1 = 1., mu1 = 0.01*B;
  rho2 = 0.001, mu2 = 0.0002, f.sigma = 1.0;

// fprintf(ferr, "J %4.1f De %4.1f \n", J, Deb);

  lambda=lambdav;
  mup=mupv;
  tau0=tau0v;

  restore(file=nameIn);
 
  scalar D2[];
  face vector D2f[];

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
  D2[] = (D2f.x[]+D2f.y[]+D2f.x[1,0]+D2f.y[0,1])/(4*sqrt(2.));
}
boundary({D2,D2f});

    double EdJ = 0.;
  foreach (reduction(+:EdJ)){
    EdJ += (2*pi*y)*(f[]*(tau0v[]*D2[]))*sq(Delta);}

    static FILE * fp;
    fp = fopen ("EdJ.txt", "a");
    fprintf (fp, "%g %g\n", timew, EdJ);
    fclose (fp);
}
