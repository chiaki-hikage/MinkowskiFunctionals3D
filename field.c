#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fitsio2.h>
#include <fftw.h>
#include <rfftw.h>
#define pi 3.141592653589793238462643

#include "main.h"
#include "org.h"

float***initmask(parlist*);
float***initfield(parlist*);
void createfield(parlist*,float***,float***);
void outputfield(parlist*,float***);
void destroyfield(parlist*,float***);

void insert_nyquist(int*,float***);
void remove_nyquist(int*,float***);
void westward_ho(int*,float***);
void eastward_ho(int*,float***);
void smoothing(float,int*,float***);
void interpolation(int,int*,float***);
void realspace_smoothing(float,int*,float***);

/* allocate storage space for the random field.  During the course of
   the program, the initial dimensions of the field will increase by
   interpolation and by Fourier transformation.  Hence, the storage
   space is allocated such that it will be able to hold the larger
   fields later on.
*/

float***initfield(parlist*par) {
  int*d=par->dim,in=par->inter;
  int i,j;
  float***u;

  /* allocate storage space for pointers to first level */
  u=(float***)calloc(in*d[0],sizeof(float**));

  /* pointers to second level */
  u[0]=(float**)calloc(in*in*d[0]*d[1],sizeof(float*));
  for(i=1;i<d[0];i++) u[i]=u[i-1]+d[1];

  /* pointers to third level */
  u[0][0]=(float*)malloc(in*in*in*d[0]*d[1]*(d[2]+2)*sizeof(float));
  for(i=1;i<d[0];i++) u[i][0]=u[i-1][0]+d[1]*d[2];
  for(i=0;i<d[0];i++) for(j=1;j<d[1];j++) u[i][j]=u[i][j-1]+d[2];

  return u;
}

float***initmask(parlist*par) {
  int i,j,k,l,linit,flag,*d=par->dim,in=par->inter,n;
  float***m,*mask,c=par->contamin;
  fitsfile*fp;
  int status=0;
  FILE*fp2;
  char inf[200],inf2[200];

  /* attempt to open the mask input file */
  fits_open_file(&fp,par->maskname,READONLY,&status);

  /* the mask is interpolated, adjust dimensions accordingly */
  if(in>1) for(i=0;i<3;i++) d[i]*=in;
  n=d[0]*d[1]*d[2];
  
  /* allocate storage space for the mask */
  m=cube(0,d[0]-1,0,d[1]-1,0,d[2]-1);
  mask=m[0][0];

  if(status==0) {

    /* read mask from file and cleanup */
    fits_read_img(fp,TFLOAT,1,n,NULL,mask,NULL,&status);
    fits_close_file(fp,&status);
    fits_report_error(stderr,status);

    /* count points in mask */
    for(i=linit=0;i<n;i++) if(mask[i]) linit++;
    fprintf(stderr,"%d of %d points (%g%%) inside\n",linit,n,100.*linit/n);

    /* determine distance of points in mask to region outside survey
       boundary. Uses Manhattan metric, which is very fast and
       accurate enough for our purpose of just cutting the edge to a
       certain depth. */

    /* preserve mask information while marking points with unknown distance */
    for(i=0;i<n;i++) if(mask[i]) mask[i]=-1;

    /* loop over relevant distances */
    for(l=flag=1;l<c&&flag;l++)
      for(i=1,flag=0;i<d[0]-1;i++)for(j=1;j<d[1]-1;j++)for(k=1;k<d[2]-1;k++)
        if(m[i][j][k]==l-1) {
          if(m[i+1][j][k]<0) m[i+1][j][k]=l,flag=1;
          if(m[i-1][j][k]<0) m[i-1][j][k]=l,flag=1;
          if(m[i][j+1][k]<0) m[i][j+1][k]=l,flag=1;
          if(m[i][j-1][k]<0) m[i][j-1][k]=l,flag=1;
          if(m[i][j][k+1]<0) m[i][j][k+1]=l,flag=1;
          if(m[i][j][k-1]<0) m[i][j][k-1]=l,flag=1;
        }

    /* count points in the inside region */
    for(i=l=0;i<n;i++) if(mask[i]<0) l++;
    fprintf(stderr,"%d of %d points (%g%%) inside\n",l,n,100.*l/n);

    /* dump the mask for debugging and exit.  Should normally be commented. */
    /* outputfield(par,m); exit(0); */

    /* remove edge from mask */
    for(i=0;i<n;i++) if(mask[i]<0) mask[i]=1; else mask[i]=0;

    /* we are going to need the original dimensions for reading the field */
    if(in>1) for(i=0;i<3;i++) d[i]/=in;

    } else {
    for(i=0;i<n;i++) mask[i]=1;
    /* status=0, fits_clear_errmsg();
    m=(float***)NULL; */
  }

  return m;
}

void createfield(parlist*par,float***m,float***u) {
  int i,n,num,*d=par->dim;
  double mean,sigma;
  float*field=u[0][0];
  float*mask=m[0][0];
  fitsfile*fp;
  int status=0;
  FILE*fp1;
  char inf[200],inf2[200];

  /* attempt to open the data input file */
  fits_open_file(&fp,par->dataname,READONLY,&status);

  
  if(status==0) {

    /* read field from file */
    fits_read_img(fp,TFLOAT,1,d[0]*d[1]*d[2],NULL,field,NULL,&status);
    fits_close_file(fp,&status);
    fits_report_error(stderr,status);

    /* normalize to zero mean and unit variance over the grids where mask=1 
                                                            C.Hikage */
    n=d[0]*d[1]*d[2];
    num=0;
    for(i=0,mean=sigma=0;i<n;i++) if(mask[i]){
      num++,mean+=field[i],sigma+=field[i]*field[i];}
    sumtomean(&mean,&sigma,num); 
    for(i=0;i<n;i++) field[i]=(field[i]-mean)/sigma;
    fprintf(stderr,"Mean:%+12.5e, Sigma:%+12.5e \n",mean,sigma);

    if(par->smooth>0) {
      smoothing(par->smooth/par->lattice,d,u);
      for(i=0,mean=sigma=0;i<n;i++) if(mask[i]){
	num++,mean+=field[i],sigma+=field[i]*field[i];}
      sumtomean(&mean,&sigma,num); 
      fprintf(stderr,"After smoothing\n");
      fprintf(stderr,"Mean:%+12.5e, Sigma:%+12.5e \n",mean,sigma);
    }

  } else {

    /* clear the FITS error from the stack */
    status=0, fits_clear_errmsg();

    /* random field is empty, just add Nyquist plane */
    insert_nyquist(d,u);

    /* generate Gaussian random field in Fourier space */
    for(i=0;i<d[0]*d[1]*(d[2]+2);i++) field[i]=gasdev();

    /* smooth and interpolate, if required, and transform to real space */
    westward_ho(d,u);
    if(par->smooth>0) smoothing(par->smooth/par->lattice,d,u);
    if(par->inter>1) interpolation(par->inter,d,u);
    eastward_ho(d,u);

    /* normalize to zero mean and unit variance */
    n=d[0]*d[1]*d[2];
    for(i=0,mean=sigma=0;i<n;i++) mean+=field[i], sigma+=field[i]*field[i];
    sumtomean(&mean,&sigma,n);
    fprintf(stderr,"Mean:%+12.5e, Sigma:%+12.5e \n",mean,sigma);
    for(i=0;i<n;i++) field[i]=(field[i]-mean)/sigma;

    /* calculate lognormal field and normalize to zero mean
    for(i=0;i<n;i++) field[i]=exp(field[i]);
    for(i=0,mean=sigma=0;i<n;i++) mean+=field[i], sigma+=field[i]*field[i];
    sumtomean(&mean,&sigma,n);
    for(i=0;i<n;i++) field[i]=field[i]/mean-1; */
  }
  return;
}

/* output a random field.  Since this is for debugging purposes only,
   the field is simply written to standard output.  If you use this,
   take care to not write anything else on standard output, and
   redirect it to a file. */

void outputfield(parlist*par,float***u) {
  int*d=par->dim;
  fwrite((void*)(u[0][0]),sizeof(float),d[0]*d[1]*d[2],stdout);
  return;
}

/* Return a processed field to the initial state, so it can be used
   for another run.  Note that this only resets the pointers, not the
   values, so subsequently a new data file has to be read or a new
   random field generated.  */

void destroyfield(parlist*par,float***u) {
  int*d=par->dim,in=par->inter;
  int i,j;

  /* field dimensions */
  for(i=0;i<3;i++) d[i]/=in;

  /* pointers to second level */
  for(i=1;i<d[0];i++) u[i]=u[i-1]+d[1];

  /* pointers to third level */
  for(i=1;i<d[0];i++) u[i][0]=u[i-1][0]+d[1]*d[2];
  for(i=0;i<d[0];i++) for(j=1;j<d[1];j++) u[i][j]=u[i][j-1]+d[2];

}

/* Add and remove a Nyquist plane to a field u of dimensions d,
   respectively.  This is done in place by rearranging values and
   re-setting pointers.
*/

void westward_ho(int*d,float***u) {
  int i,j,k,l=d[0]*d[1]*d[2];
  static rfftwnd_plan plan=0;


  /* add a Nyquist plane to the array */
  insert_nyquist(d,u);

  
  /* plan the Fourier transform */
  if(!plan) plan=rfftw3d_create_plan
              (d[0],d[1],d[2],
	       FFTW_REAL_TO_COMPLEX,
	       FFTW_ESTIMATE|FFTW_IN_PLACE);

  /* compute the fastest Fourier transform of the array */
  rfftwnd_one_real_to_complex(plan,(fftw_real*)u[0][0],(fftw_complex*)NULL);

  /* get the correct normalization */
  for(i=0;i<d[0];i++) for(j=0;j<d[1];j++) for(k=0;k<d[2]+2;k++) u[i][j][k]/=l;

  return;
}

void eastward_ho(int*d,float***u) {
  int i,j,k,l=d[0]*d[1]*d[2];
  static rfftwnd_plan plan=0;

  /* plan the Fourier transform */
  if(!plan) plan=rfftw3d_create_plan
              (d[0],d[1],d[2],
               FFTW_COMPLEX_TO_REAL,
	       FFTW_ESTIMATE|FFTW_IN_PLACE);

  /* compute the reverse fastest Fourier transform of the array */
  rfftwnd_one_complex_to_real(plan,(fftw_complex*)u[0][0],(fftw_real*)NULL);

  /* remove the Nyquist plane */
  remove_nyquist(d,u);

  return;
}

void insert_nyquist(int*d,float***u) {
  int i,j,k;

  /* rearrange the values */
  for(i=d[0]-1;i>=0;i--) for(j=d[1]-1;j>=0;j--) for(k=d[2]-1;k>=0;k--) 
    u[0][0][k+(j+i*d[1])*(d[2]+2)]=u[i][j][k];

  /* calculate new pointers */
  for(i=1;i<d[0];i++) u[i][0]=u[i-1][0]+d[1]*(d[2]+2);
  for(i=0;i<d[0];i++) for(j=1;j<d[1];j++) u[i][j]=u[i][j-1]+(d[2]+2);

  return;
}
void remove_nyquist(int*d,float***u) {
  int i,j,k;

  /* rearrange the values */
  for(i=0;i<d[0];i++) for(j=0;j<d[1];j++) for(k=0;k<d[2];k++)
    u[0][0][k+(j+i*d[1])*d[2]]=u[0][0][k+(j+i*d[1])*(d[2]+2)];

  /* calculate new pointers */
  for(i=1;i<d[0];i++) u[i][0]=u[i-1][0]+d[1]*d[2];
  for(i=0;i<d[0];i++) for(j=1;j<d[1];j++) u[i][j]=u[i][j-1]+d[2];

  return;
}

void smoothing(float s,int*d,float***u) {
  int i,j,k,m;
  double**g;
  double sigma,filter;

  /* calculate smoothing kernels in Fourier space */
  g=(double**)malloc(3*sizeof(double*));

  fprintf(stderr,"smoothing:%+12.5e \n",s/d[0]);
  for(m=0;m<3;m++) {
    g[m]=dvector(0,d[m]);
    sigma=2.*pi*s/d[m];
    for(i=0;i<=d[m]/2;i++) g[m][d[m]-i]=g[m][i]=exp(-.5*pow(i*sigma,2.));
  }
  
  fprintf(stderr,"filter:%+12.5e\n",g[0][d[0]-1]);

  /* apply smoothing kernels in Fourier space */
  for(i=0;i<d[0];i++) for(j=0;j<d[1];j++) for(k=0;k<=d[2]/2;k++) {
    filter=g[0][i]*g[1][j]*g[2][k];
    u[i][j][2*k]*=filter, u[i][j][2*k+1]*=filter;
  }

  /* cleanup */
  for(m=0;m<3;m++) free_dvector(g[m],0,d[m]);
  free((void*)g);
    
}

/* interpolate a three-dimensional field u by moving its Fourier
   transform to a larger array.  Even though the old and new order are
   accessed through different pointer structures, they occupy the same
   area in memory.  Therefore, care must be taken or we will stomp on
   data. */

/* NOTE: I have not been able to get in-place interpolation working
   properly right away.  Therefore, I am using an auxiliary field for
   the time being so I can move on to other parts of the code that
   need more attention.  This introduces a memory overhead of up to
   12.5%, so it should be fixed at some point.  However, Fourier space
   interpolation may not be the way to go anyway, because it can
   introduce considerable ringing for cut fields.  Instead, the fields
   should be generated with higher resolution in the first place.  */

void interpolation(int in,int*d,float***u) {
  int i,j,k;
  float***o;

  /* allocate space for old field and store it */
  o=cube(0,d[0]-1,0,d[1]-1,0,d[2]+1);

  /* store old values */
  for(i=0;i<d[0];i++) for(j=0;j<d[1];j++) for(k=0;k<d[2]+2;k++)
    o[i][j][k]=u[i][j][k];

  /* calculate new pointers */
  for(i=1;i<in*d[0];i++) u[i]=u[i-1]+in*d[1];
  for(i=1;i<in*d[0];i++) u[i][0]=u[i-1][0]+in*d[1]*(in*d[2]+2);
  for(i=0;i<in*d[0];i++) for(j=1;j<in*d[1];j++) u[i][j]=u[i][j-1]+(in*d[2]+2);

  /* empty new field */
  for(i=0;i<in*d[0];i++) for(j=0;j<in*d[1];j++) for(k=0;k<in*d[2]+2;k++)
    u[i][j][k]=0;

  /* store new values in old field */
  for(i=0;i<=d[0]/2;i++) for(j=0;j<=d[1]/2;j++) for(k=0;k<d[2]+2;k++) {
    ;            u[        i][        j][k]=o[     i][     j][k];
    if(i>0     ) u[in*d[0]-i][        j][k]=o[d[0]-i][     j][k];
    if(     j>0) u[        i][in*d[1]-j][k]=o[     i][d[1]-j][k];
    if(i>0&&j>0) u[in*d[0]-i][in*d[1]-j][k]=o[d[0]-i][d[1]-j][k];
  }

  /* forget old field */
  free_cube(o,0,d[0]-1,0,d[1]-1,0,d[2]+1);

  /* calculate new dimensions */
  for(i=0;i<3;i++) d[i]*=in;

  return;

}
// void inplace_interpolation(int in,int*d,float***u) {
//   int i,j,k,offset=d[0]/2;
//   float***o;
//   int n;
// 
//   /* allocate space for pointers only, values are handled in place */
//   o=(float***)calloc(d[0],sizeof(float**));
//   o[0]=(float**)calloc(d[0]*d[1],sizeof(float*));
//   o[0][0]=u[0][0];
//   for(i=1;i<d[0];i++) o[i]=o[i-1]+d[1];
//   for(i=1;i<d[0];i++) o[i][0]=o[i-1][0]+d[1]*(d[2]+2);
//   for(i=0;i<d[0];i++) for(j=1;j<d[1];j++) o[i][j]=o[i][j-1]+(d[2]+2);
// 
//   /* calculate new pointers */
//   for(i=1;i<in*d[0];i++) u[i]=u[i-1]+in*d[1];
//   for(i=1;i<in*d[0];i++) u[i][0]=u[i-1][0]+in*d[1]*(in*d[2]+2);
//   for(i=0;i<in*d[0];i++) for(j=1;j<in*d[1];j++) u[i][j]=u[i][j-1]+(in*d[2]+2);
// 
//   /* empty additional values */
//   for(n=0,i=d[0]*d[1]*(d[2]+2);i<in*d[0]*in*d[1]*(in*d[2]+2);i++,n++)
//     u[0][0][i]=0;
//   fprintf(stderr,"%d new values emptied\n",n);
// 
//   /* move values, first into intermediate positions in empty space ... */
//   for(n=0,i=0;i<d[0];i++) for(j=0;j<d[1];j++) for(k=0;k<d[2]+2;k++,n++)
//     u[i+offset][j][k]=o[i][j][k];
//   fprintf(stderr,"%d values moved\n",n);
// 
//   /* ... and then to final destinations */
//   for(n=0,i=0;i<=d[0]/2;i++) for(j=0;j<=d[1]/2;j++) for(k=0;k<d[2]+2;k++) {
//     ;            u[        i][        j][k]=u[     i+offset][     j][k],n++;
//     if(i>0     ) u[in*d[0]-i][        j][k]=u[d[0]-i+offset][     j][k],n++;
//     if(     j>0) u[        i][in*d[1]-j][k]=u[     i+offset][d[1]-j][k],n++;
//     if(i>0&&j>0) u[in*d[0]-i][in*d[1]-j][k]=u[d[0]-i+offset][d[1]-j][k],n++;
//   }
//   fprintf(stderr,"%d values moved again\n",n);
// 
//   for(n=0,i=0;i<d[0];i++) for(j=0;j<d[1];j++) for(k=0;k<d[2]+2;k++,n++)
//     u[i+offset][j][k]=0;
//   fprintf(stderr,"%d intermediate values emptied\n",n);
// 
//   /* forget old pointers */
//   free((void*)o[0]);
//   free((void*)o);
// 
//   /* calculate new dimensions */
//   for(i=0;i<3;i++) d[i]*=in;
// 
//   return;
// 
// }

#define MAX 8
void realspace_smoothing(float s,int*d,float***u) {
  int i,j,k,m,x,y;
  float*v,*w;
  double*g,norm;

  /* calculate smoothing kernels in real space */
  g=dvector(-MAX,+MAX);
  for(i=-MAX,norm=0;i<=MAX;i++) norm+=g[i]=exp(-.5*pow(i/s,2.));
  for(i=-MAX;i<=MAX;i++) g[i]/=norm;
  fprintf(stderr,"normalization %g (should be %g), border value %g\n",
	  norm,sqrt(2*pi)*s,g[-MAX]);

  /* loop over directions */
  for(m=0;m<3;m++) {
    v=vector(0,d[m]-1); w=vector(0,d[m]-1);
    if(!v||!w) Mecker("not enough space for realspace smoothing");
    
    /* loop over lines along this direction */
    for(j=0;j<d[(m+1)%3];j++) for(k=0;k<d[(m+2)%3];k++) {
      /* retrieve old line from field */
      switch(m) {
      case 0 : for(i=0;i<d[m];i++) v[i]=u[i][j][k]; break;
      case 1 : for(i=0;i<d[m];i++) v[i]=u[k][i][j]; break;
      case 2 : for(i=0;i<d[m];i++) v[i]=u[j][k][i]; break;
      }
      
      /* convolve old line with the kernel to make new line */
      for(i=0;i<d[m];i++) for(w[i]=0,y=-MAX;y<=+MAX;y++) if(g[y]>1e-6) {
	x=i-y; if(x>=d[m]) x-=d[m]; if(x<0) x+=d[m];
	w[i]+=g[y]*v[x];
      }
      
      /* store new line to field */
      switch(m) {
      case 0 : for(i=0;i<d[m];i++) u[i][j][k]=w[i]; break;
      case 1 : for(i=0;i<d[m];i++) u[k][i][j]=w[i]; break;
      case 2 : for(i=0;i<d[m];i++) u[j][k][i]=w[i]; break;
      }
    }
    
    free_vector(v,0,d[m]-1); free_vector(w,0,d[m]-1);
  } /* end of loop over directions */
  
  /* cleanup */
  free_dvector(g,-MAX,+MAX);
}
