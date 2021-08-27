#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <fftw.h>
#include <rfftw.h>


#include "main.h"
#include "org.h"
#include "field.h"
#include "minkowski.h"


/* definitions of constants */
#define pi 3.141592653589793238462643

/* declarations */
void Default(parlist*);
void Argument(parlist*,int,char*[]);
void Explain(parlist*,FILE*);

int main(int argc, char *argv[]) 
{
  int i;
  parlist par;			/* parameters */
  float***u;			/* the field */
  float***m;			/* the mask */
  double*t;			/* the threshold values */
  double*b;			/* the bin widths */
  double***v,***vsum,***vsqr;	/* the Minkowski functionals */
  char infname[100];

  Default(&par);
  Argument(&par,argc,argv);

  t=initthreshold(&par,&b);
  m=initmask(&par);
  v=initminkowski(&par);
  vsum=initminkowski(&par);
  vsqr=initminkowski(&par);

  for(i=0;i<par.length;i++) {
    u=initfield(&par);
    createfield(&par,m,u);
    minkowski(&par,m,u,t,b,v);
    addminkowski(&par,v,vsum,vsqr);
    if(i==par.length-1) outputfield(&par,u);
    destroyfield(&par,u);
  }
  
  finishminkowski(&par,vsum,vsqr);
  outputminkowski(&par,t,vsum,vsqr);

  exit(0);
}

/* set default values for global variables */
void Default(parlist*par) {
  int i;

  /* random number generator seeds */
  seeds.i0=20;
  seeds.i1=2;
  seeds.i2=94;

  /* miscellaneous */
  par->length=1;
  par->smooth=0;
  par->inter=1;
  par->contamin=0.5;

  /* data */
  for(i=0;i<3;i++) par->dim[i]=0;
  par->lattice=0;
  par->dataname="";
  par->maskname="";

  /* Minkowski functionals */
  par->lo=0;
  par->hi=0;
  par->bins=0;
  par->width=0;
  par->mininame="";
}

/* evaluate command line into the global parameter list */
void Argument(parlist*par,int argc,char*argv[]) {
  char c;

  while((char)EOF!=(c=getopt(argc,argv,"l:i:s:c:x:y:z:a:d:m:t:T:n:w:o:"))) 
    switch(c) {
    case 'l': par->length   = atoi(optarg); break;
    case 'i': par->inter    = atoi(optarg); break;
    case 's': par->smooth   = atof(optarg); break;
    case 'c': par->contamin = atof(optarg); break;
    case 'x': par->dim[0]   = atoi(optarg); break;
    case 'y': par->dim[1]   = atoi(optarg); break;
    case 'z': par->dim[2]   = atoi(optarg); break;
    case 'a': par->lattice  = atof(optarg); break;
    case 'd': par->dataname =      optarg ; break;
    case 'm': par->maskname =      optarg ; break;
    case 't': par->lo       = atof(optarg); break;
    case 'T': par->hi       = atof(optarg); break;
    case 'n': par->bins     = atoi(optarg); break;
    case 'w': par->width    = atof(optarg); break;
    case 'o': par->mininame =      optarg ; break;
    default : Explain(par,stdout);
    }

  /* calculate parameters from other parameters, if necessary */
  if(!par->lattice) par->lattice=1./(par->inter*par->dim[0]);
  if(!par->width) par->width=(par->hi-par->lo)/par->bins;

  return;
}

/* Explain the meaning of each variable in the global parameter list */
void Explain(parlist*par,FILE*fp) {
  fprintf(fp,"# options\n"
	  "# -l     number of realizations of field       (actual %d)\n"
	  "# -i     interpolation factor                  (actual %d)\n"
	  "# -s     width of smoothing kernel             (actual %g)\n"
	  "# -c     maximum contamination of survey region(actual %g)\n"
	  "# -x,y,z size of grid in pixels                (actual %d,%d,%d)\n"
	  "# -a     lattice constant                      (actual %g)\n"
	  "# -d     input file for data                   (actual %s)\n"
	  "# -m     input file for mask                   (actual %s)\n"
	  "# -t,T   threshold range                       (actual %g,%g)\n"
	  "# -n     number of bins                        (actual %d)\n"
	  "# -w     width of bins                         (actual %g)\n"
	  "# -o     output file for Minkowski functionals (actual %s)\n",
	  par->length,par->inter,par->smooth,par->contamin,
	  par->dim[0],par->dim[1],par->dim[2],par->lattice,
	  par->dataname,par->maskname,
	  par->lo,par->hi,par->bins,par->width,par->mininame);

  if(fp==stdout) exit(99); else return;
}
int sadd(char *as, char *bs, char *cs)
{
  while(*as) *cs++=*as++;
  while(*bs) *cs++=*bs++;
  *cs=0x00;
}
int sdraw(char *as, char *bs, char *cs)
{
  while(*as) {
  while(*bs) {
    *as++;
    *bs++;
  };
  *cs++=*as++;
  }
  *cs=0x00;
}










