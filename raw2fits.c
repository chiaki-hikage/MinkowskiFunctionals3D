#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fitsio2.h>

typedef struct parlist {
  int fortran;
  int dim[3];
  char*output;
} parlist;

void Default(parlist*);
void Argument(parlist*,int,char*[]);

int main(int argc,char*argv[]) {
  parlist par;
  float*u,min,max;
  fitsfile*fp;
  long d[3],n,nread;
  int i,test,status=0;

  /* parse command line */
  Default(&par);
  Argument(&par,argc,argv);

  /* allocate storage space */
  n=par.dim[0]*par.dim[1]*par.dim[2];
  u=(float*)malloc(sizeof(float)*n);
  if(!u) exit(1);

  /* discard the first value of a Fortran array */
  /* if(par.fortran) {
    fread((void*)&test,sizeof(int),1,stdin);
    fprintf(stderr,"read array size %d from Fortran array\n",test);
    if(test!=sizeof(float)*n) 
      fprintf(stderr,
	      "first value of Fortran array is not the expected size %d \n"
	      "make sure field dimensions and byte ordering are correct\n",
              n);
   } */

  /* read raw data from standard input */
  for(i=0;i<n;i++){
    fscanf(stdin,"%f \n",&u[i]);
  }

  nread=fread((void*)u,sizeof(float),n,stdin);
  fprintf(stderr,"read %ld values into array, expected %ld values\n",nread,n); 

  /* determine the range of values */
  fprintf(stderr,
	  "determining range of data ...\n"
	  "check byte ordering if this throws a floating exception\n");
  for(i=1,min=max=u[0];i<n;i++) {
    if(u[i]<min) min=u[i];
    if(u[i]>max) max=u[i];
  }
  fprintf(stderr,"data range from %g to %g\n",min,max);

  /* open fits file */
  fits_create_file(&fp,par.output,&status);

  /* write header - note that FITS uses Fortran ordering for the array */
  d[0]=par.dim[2], d[1]=par.dim[1], d[2]=par.dim[0];
  fits_create_img(fp,FLOAT_IMG,3,d,&status);

  /* write data */
  fits_write_img(fp,TFLOAT,1,n,u,&status);

  /* cleanup */
  fits_close_file(fp,&status);
  fits_report_error(stderr,status);

  exit(0);
}

void Default(parlist*par) {
  par->output="";
  par->dim[0]=par->dim[1]=par->dim[2]=0;
  par->fortran=0;
}

void Argument(parlist*par,int argc,char*argv[]) {
  char c;
  
  while((char)EOF!=(c=getopt(argc,argv,"d:x:y:z:F"))) 
    switch(c) {
    case 'd': par->output  =      optarg ; break;
    case 'x': par->dim[0]  = atoi(optarg); break;
    case 'y': par->dim[1]  = atoi(optarg); break;
    case 'z': par->dim[2]  = atoi(optarg); break;
    case 'F': par->fortran = 1           ; break;
    }
  return;
}
