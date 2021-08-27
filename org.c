#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "main.h"

#define pi 3.141592653589793238462643

typedef struct IMMENSE {unsigned long l,r;} immense;
typedef struct GREAT {unsigned short l,c,r;} great;

/* Declarations */
int jacobi(double**,int,double*,double**);
float newtonroot(float);
void sort(unsigned long, float*);
void indexx(unsigned long, float*, unsigned long*);
double gasdev();
void sumtomean(double*,double*,int);
double errf(double);
double gaussf(double);
double gam(double);
double gaminc(double,double);
double hermite(int,double);
double omega(int);
double qromb(double(*)(double),double,double);
void polint(double*,double*,int,double,double*,double*);
double trapzd(double(*)(double),double,double,int);

int fileopen(FILE**,char*,char*);
void fileclose(FILE**);

void Mecker(char*);

int*ivector(int,int);
unsigned long*lvector(int,int);
float*vector(int,int);
double*dvector(int,int);
char**cmatrix(int,int,int,int);
float**matrix(int,int,int,int);
int**imatrix(int,int,int,int);
double**dmatrix(int,int,int,int);
int***icube(int,int,int,int,int,int);
float***cube(int,int,int,int,int,int);
double***dcube(int,int,int,int,int,int);
void free_ivector(int*,int,int);
void free_vector(float*,int,int);
void free_lvector(unsigned long*,int,int);
void free_dvector(double*,int,int);
void free_cmatrix(char**,int,int,int,int);
void free_imatrix(int**,int,int,int,int);
void free_matrix(float**,int,int,int,int);
void free_dmatrix(double**,int,int,int,int);
void free_icube(int***,int,int,int,int,int,int);
void free_cube(float***,int,int,int,int,int,int);
void free_dcube(double***,int,int,int,int,int,int);

/* Declarations for random number generator */
struct SEEDS { int i0, i1, i2; };
double rand1(void);
struct SEEDS seeds;

/* Computes all eigenvalues and eigenvectors of a real symmetric
   matrix a, which is of size n by n. On output, elements of a above
   the diagonal are destroyed.  d returns the eigenvalues of a.  v is
   a matrix with that same dimensions as a whose columns contain, on
   output, the normalized eigenvectors of a.  returns the number of
   Jacobi rotations which were required. */
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau);
int jacobi(double**a,int n,double*d,double**v) 
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  b=dvector(0,n-1);
  z=dvector(0,n-1);
  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  for (i=1;i<=50;i++) {
    sm=0.0;
    for(ip=0;ip<n;ip++) for(iq=ip+1;iq<n;iq++) sm+=fabs(a[ip][iq]);
    if(sm==0.0) {
      free_dvector(z,0,n-1);
      free_dvector(b,0,n-1);
      return i;
    }
    if(i<4) tresh=0.2*sm/(n*n); else tresh=0.0;
    for (ip=0;ip<n;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
	    && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((float)(fabs(h)+g) == (float)fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=0;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<n;j++) {
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=0;j<n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	}
      }
    }
    for(ip=0;ip<n;ip++) b[ip]+=z[ip], d[ip]=b[ip], z[ip]=0.0;
  }
  Mecker("Too many iterations in routine jacobi");
  return 42;			/* keep gcc happy */
}
#undef ROTATE

/* find the threshold giving a Gaussian v0 using Newton's method */
#define JMAX 20
#define XMAX 4
#define XACC 1e-7
#define FACC 1e-7
float newtonroot(float v0)
{
  int j;
  float df,dx,f,x;
  
  x=0;
  for(j=1;j<=JMAX;j++) {
    /* function and derivative */
    f=.5-.5*errf(x/sqrt(2.)); df=-exp(-.5*x*x)/sqrt(2*pi);
    /* one Newton step */
    dx=(f-v0)/df; x-=dx; if(fabs(f-v0)<FACC) return x;
    /* catch asymptotic values */
    if(x<-XMAX) return -XMAX; if(x>XMAX) return XMAX;
  }
  Mecker("Maximum number of iterations exceeded in rtnewt");
  return 0.0;
}
#undef JMAX
#undef XMAX
#undef XACC

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void sort(unsigned long n, float arr[])
{
  unsigned long i,ir=n,j,k,l=1;
  int jstack=0,*istack;
  float a,temp;
  
  istack=ivector(1,NSTACK);
  for(;;) {
    if(ir-l<M) {
      for(j=l+1;j<=ir;j++) {
	a=arr[j];
	for(i=j-1;i>=l;i--) {
	  if(arr[i]<=a) break;
	  arr[i+1]=arr[i];
	}
	arr[i+1]=a;
      }
      if(jstack==0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    }else {
      k=(l+ir)>>1;
      SWAP(arr[k],arr[l+1])
	if(arr[l]>arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    }
      if(arr[l+1]>arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  }
      if(arr[l]>arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for(;;) {
	do i++;while(arr[i]<a);
	do j--;while(arr[j]>a);
	if(j<i)break;
	SWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      jstack+=2;
      if(jstack>NSTACK) Mecker("NSTACK too small insort.");
      if(ir-i+1>=j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      }else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI

 
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void indexx(unsigned long n, float arr[], unsigned long indx[])
{
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	float a;

	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) Mecker("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software z!0(0. */ 

/* Returns a normally distributed deviate with zero mean and unit
   variance, using a uniform deviates.  */
double gasdev()
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if(iset==0) {
    do {
      v1=2*rand1()-1; v2=2*rand1()-1; rsq=v1*v1+v2*v2;
    } while(rsq>=1||rsq==0);
    fac=sqrt(-2*log(rsq)/rsq); 
    gset=v1*fac; iset=1;
    return v2*fac;
  } else {
    iset=0; return gset;
  }
}

void sumtomean(double*sum,double*err,int length)
{
  *sum/=length; 
  if(length>1) {
    *err=(*err-*sum**sum*length)/(length-1);
    *err=*err<=0?0:sqrt(*err);
  } else {
    *err=0;
  }
}
double errf(double x)
{
  if(x> 10) return  1; if(x<-10) return -1;
  if(fabs(x)<1e-6) return x*2/sqrt(pi);
  return qromb(gaussf,0,x);
}
double gaussf(double x)
{
  return 2/sqrt(pi)*exp(-x*x);
}
double gaminc(double n,double x)
{
  if(n<=0)  return -666;
  if(n==.5) return sqrt(pi)*(1.-errf(sqrt(x)));
  if(n==1)  return exp(-x);
  return (n-1)*gaminc(n-1,x)+pow(x,n-1)*exp(-x);
}
double gam(double x)
{
  if(x<=0)  return -666;
  if(x==.5) return sqrt(pi);
  if(x==1)  return 1;
  return (x-1)*gam(x-1);
}
double hermite(int n,double x)
{
  if(n<0)  return 0;
  if(n==0) return exp(-.5*x*x)/sqrt(2*pi);
  return x*hermite(n-1,x)-(n-1)*hermite(n-2,x);
}
double omega(int k)
{
  switch(k) {
  case 0 : return 1;
  case 1 : return 2;
  case 2 : return pi;
  case 3 : return pi/.75;
  default : return pow(pi,k/2.)/gam(1+k/2.);
  }
}

#define EPS 1.0e-6
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

double qromb(double (*func)(double), double a, double b)
{
  double ss,dss;
  double s[JMAXP],h[JMAXP+1];
  int j;
  
  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=trapzd(func,a,b,j);
    if (j >= K) {
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    h[j+1]=0.25*h[j];
  }
  fprintf(stderr,"Too many steps in routine qromb");
  return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*b;
  
  dif=fabs(x-xa[1]);
  c=dvector(1,n);
  b=dvector(1,n);
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    b[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-b[i];
      if ( (den=ho-hp) == 0.0) fprintf(stderr,"Error in routine polint");
      den=w/den;
      b[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : b[ns--]));
  }
  free_dvector(b,1,n);
  free_dvector(c,1,n);
}

#define FUNC(x) ((*func)(x))

double trapzd(double (*func)(double), double a, double b, int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) {
      sum += FUNC(x);
    }
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}
#undef FUNC

/* open and close files; filenames starting with '-' are interpreted
   as stdin and stdout */
int fileopen(FILE**fp,char*name,char*mode)
{
  switch(name[0]) {
  case 0  : *fp = (FILE*)NULL; return 0;
  case '-': 
    switch(mode[0]) { 
    case 'r': *fp= stdin; return 1;
    case 'w': *fp=stdout; return 1;
    default : return 0;
    }
  default : if((*fp=fopen(name,mode))) return 1; else return 0;
  }
}
void fileclose(FILE**fp)
{
  if(*fp==stdin||*fp==stdout) return;
  fclose(*fp); *fp=(FILE*)NULL; return;
}
  
/* complain about an error and exit */
void Mecker(char*meckertext)
{
  fprintf(stderr,"%s\n",meckertext);
  exit(666);
}

/* various memory allocations */
int*ivector(int nl,int nh)
{
  int*v;
  
  v=(int*)calloc((unsigned) (nh-nl+1),sizeof(int));
  if(!v) return v; else return v-nl;
}
unsigned long*lvector(int nl,int nh)
{
  unsigned long*v;
  
  v=(unsigned long*)calloc((unsigned) (nh-nl+1),sizeof(unsigned long));
  if(!v) return v; else return v-nl;
}
float*vector(int nl,int nh)
{
  float*v;
  
  v=(float*)calloc((unsigned) (nh-nl+1),sizeof(float));
  if (!v) return v; else return v-nl;
}
double*dvector(int nl,int nh)
{
  double*v;
  v=(double*)calloc((size_t)(nh-nl+1),(size_t)sizeof(double));
  if (!v) return v; else return v-nl;
}
int**imatrix(int nrl,int nrh,int ncl,int nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;
  
  /* allocate pointers to rows */
  m=(int**)calloc((nrow),sizeof(int*));
  if(!m) return m; else m-=nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(int*)calloc((nrow*ncol),sizeof(int));
  if(!m[nrl]){ free((void*)m); return (int**)NULL; } else m[nrl]-=ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}
float**matrix(int nrl,int nrh,int ncl,int nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float**m;
  
  /* allocate pointers to rows */
  m=(float**)calloc((nrow),sizeof(float*));
  if(!m) return m; else m-=nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(float*)calloc((nrow*ncol),sizeof(float));
  if(!m[nrl]){ free((void*)m); return (float**)NULL; } else m[nrl]-=ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}
char**cmatrix(int nrl,int nrh,int ncl,int nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  char **m;
  
  /* allocate pointers to rows */
  m=(char**)calloc((nrow),sizeof(char*));
  if(!m) return m; else m-=nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(char*)calloc((nrow*ncol),sizeof(char));
  if(!m[nrl]){ free((void*)m); return (char**)NULL; } else m[nrl]-=ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}
double**dmatrix(int nrl,int nrh,int ncl,int nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;
  
  /* allocate pointers to rows */
  m=(double**)calloc((nrow),sizeof(double*));
  if(!m) return m; else m-=nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double*)calloc((nrow*ncol),sizeof(double));
  if(!m[nrl]){ free((void*)m); return (double**)NULL; } else m[nrl]-=ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}
int***icube(int il,int ih,int jl,int jh, int kl, int kh)
{
  long i,j,ni=ih-il+1,nj=jh-jl+1,nk=kh-kl+1;
  int ***c;
  
  /* allocate pointers to first coordinate */
  c=(int***)calloc((ni),sizeof(int**));
  if(!c) return c;
  c-=il;

  /* allocate second coordinate and set pointers to them */
  c[il]=(int**)calloc((ni*nj),sizeof(int*));
  if(!c[il]){ free((void*)c); return (int***)NULL; }
  c[il]-=jl;
  for(i=il+1;i<=ih;i++) c[i]=c[i-1]+nj;
  
  /* allocate third coordinate and set pointers to them */
  c[il][jl]=(int*)calloc((ni*nj*nk),sizeof(int));
  if(!c[il][jl]){ free((void*)c[il]); free((void*)c); return (int***)NULL; }
  c[il][jl]-=kl;
  for(i=il+1;i<=ih;i++) c[i][jl]=c[i-1][jl]+nj*nk;
  for(i=il;i<=ih;i++) for(j=jl+1;j<=jh;j++) c[i][j]=c[i][j-1]+nk;
  
  return c;
}
float***cube(int il,int ih,int jl,int jh, int kl, int kh)
{
  long i,j,ni=ih-il+1,nj=jh-jl+1,nk=kh-kl+1;
  float ***c;
  
  /* allocate pointers to first coordinate */
  c=(float***)calloc((ni),sizeof(float**));
  if(!c) return c;
  c-=il;

  /* allocate second coordinate and set pointers to them */
  c[il]=(float**)calloc((ni*nj),sizeof(float*));
  if(!c[il]){ free((void*)c); return (float***)NULL; }
  c[il]-=jl;
  for(i=il+1;i<=ih;i++) c[i]=c[i-1]+nj;
  
  /* allocate third coordinate and set pointers to them */
  c[il][jl]=(float*)calloc(ni*nj*nk,sizeof(float));
  if(!c[il][jl]){ free((void*)c[il]); free((void*)c); return (float***)NULL; }
  c[il][jl]-=kl;
  for(i=il+1;i<=ih;i++) c[i][jl]=c[i-1][jl]+nj*nk;
  for(i=il;i<=ih;i++) for(j=jl+1;j<=jh;j++) c[i][j]=c[i][j-1]+nk;
  
  return c;
}
double***dcube(int il,int ih,int jl,int jh, int kl, int kh)
{
  long i,j,ni=ih-il+1,nj=jh-jl+1,nk=kh-kl+1;
  double ***c;
  
  /* allocate pointers to first coordinate */
  c=(double***)calloc((ni),sizeof(double**));
  if(!c) return c;
  c-=il;

  /* allocate second coordinate and set pointers to them */
  c[il]=(double**)calloc((ni*nj),sizeof(double*));
  if(!c[il]){ free((void*)c); return (double***)NULL; }
  c[il]-=jl;
  for(i=il+1;i<=ih;i++) c[i]=c[i-1]+nj;
  
  /* allocate third coordinate and set pointers to them */
  c[il][jl]=(double*)calloc((ni*nj*nk),sizeof(double));
  if(!c[il][jl]){ free((void*)c[il]); free((void*)c); return (double***)NULL; }
  c[il][jl]-=kl;
  for(i=il+1;i<=ih;i++) c[i][jl]=c[i-1][jl]+nj*nk;
  for(i=il;i<=ih;i++) for(j=jl+1;j<=jh;j++) c[i][j]=c[i][j-1]+nk;
  
  return c;
}
void free_ivector(int*v,int nl,int nh)
{
  free((void*) (v+nl));
}
void free_lvector(unsigned long*v,int nl,int nh)
{
  free((void*) (v+nl));
}
void free_vector(float*v,int nl,int nh)
{
  free((void*) (v+nl));
}
void free_dvector(double*v,int nl,int nh)
{
  free((void*) (v+nl));
  return;
}
void free_cmatrix(char**m,int nrl,int nrh,int ncl,int nch)
{
  free((void*) (m[nrl]+ncl));
  free((void*) (m+nrl));
}
void free_imatrix(int**m,int nrl,int nrh,int ncl,int nch)
{
  free((void*) (m[nrl]+ncl));
  free((void*) (m+nrl));
}
void free_matrix(float**m,int nrl,int nrh,int ncl,int nch)
{
  free((void*) (m[nrl]+ncl));
  free((void*) (m+nrl));
}
void free_dmatrix(double**m,int nrl,int nrh,int ncl,int nch)
{
  free((void*) (m[nrl]+ncl));
  free((void*) (m+nrl));
}
void free_icube(int***c,int il,int ih,int jl,int jh, int kl, int kh)
{
  free((void*)(c[il][jl]+kl));
  free((void*)(c[il]+jl));
  free((void*)(c+il));  
}
void free_cube(float***c,int il,int ih,int jl,int jh, int kl, int kh)
{
  free((void*)(c[il][jl]+kl));
  free((void*)(c[il]+jl));
  free((void*)(c+il));  
}
void free_dcube(double***c,int il,int ih,int jl,int jh, int kl, int kh)
{
  free((void*)(c[il][jl]+kl));
  free((void*)(c[il]+jl));
  free((void*)(c+il));  
}

/* --------------------------------------------------------------- */
/*                                                                 */
/*   Pseudo-randomgenerator (A. Starkov) based on D. Lehmer's      */
/*   residual method                                               */
/*                                                                 */
/*   Parameters : none                                             */
/*   Result     : identical distributed pseudo-randomgenerated     */
/*                numbers in [0,1]                                 */
/* --------------------------------------------------------------- */
double rand1(void)
{
  long int  k;
  short     ii0, ii1;
  
  const double        rm = 2.9802322E-08;             /* constants */
  const long int       m = 16384;
  const long int      j0 = 11973;
  const long int      j1 = 17784;  /* j0, j1, j2 have to be of the */
  const long int      j2 = 710;    /* type long (!), otherwise the */
                                   /* multiplication i0*j0, etc.   */
                                   /* will be only int-multiplic.  */

  k   = seeds.i0*j0;
  ii0 = (short) k;                    /* lower 2 Byte of k         */
  k   = 2 * (k >> 16);                /* (higher 2 Byte of k) * 2  */
  if (ii0 < 0)
    {
      ii0 = (short)(ii0 + 32768L);
      k++;
    }

  k   = j1*seeds.i0 + j0*seeds.i1 + k;
  ii1 = (short) k;                    /* lower 2 Byte of k         */
  k   = 2 * (k >> 16);                /* (higher 2 Byte of k) * 2  */
  if (ii1 < 0)
    {
      ii1 = (short)(ii1 + 32768L);
      k++;
    }

  k  = j2*seeds.i0 + j0*seeds.i2 + j1*seeds.i1 + k;
  seeds.i0 = ii0;                      /* store i0, i1, i2           */
  seeds.i1 = ii1;
  seeds.i2 = 1023 & (short)k;          /* 1023 & (lower 2 Byte of k) */
                             /* order of mult. equal, bec. i2 < 1023 */
  k   = 2 * m * seeds.i2 + seeds.i1;   
  return k * rm;
  }
