
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define D 3

#define pi 3.141592653589793238462643
#define infinity 6e66
#define epsilon 1e-10
#define false 0
#define true 1

#include "main.h"
#include "org.h"

/* declarations */
double***initminkowski(parlist*);
double*initthreshold(parlist*,double**);
double*vfrthreshold(parlist*,double**,float***,float***);
void minkowski(parlist*,float***,float***,double*,double*,double***);
void outputminkowski(parlist*,double*,double***,double***);
void addminkowski(parlist*,double***,double***,double***);
void finishminkowski(parlist*,double***,double***);

int threshold_crofton(double,double*,int);
int threshold_lowerkoenderink(double,double*,double*,int);
int threshold_upperkoenderink(double,double*,double*,int);

/* allocate storage space for Minkowski functionals */

double***initminkowski(parlist*par) {
  return dcube(0,2,0,par->bins,0,D);
}

/* allocate storage space for and initialize threshold bins */

double*initthreshold(parlist*par,double**bin) {
  int j,b=par->bins;
  double*thr;

  thr=dvector(0,b);
  for(j=0;j<=b;j++) thr[j]=par->lo+j*(par->hi-par->lo)/b;

  if(bin) {
    *bin=dvector(0,b);
    for(j=0;j<=b;j++) (*bin)[j]=par->width;
  }

  return thr;
}

double*vfrthreshold(parlist*par,double**bin,float***u,float***m) {
  int i,j,b=par->bins,*d=par->dim;
  long  inthre;
  unsigned long *indx,nmesh,nout;
  double*thr,*vfr,dthr,min,max;
  float *field,*mask,fnthre,vfrac,nthre;
  if(m) mask=m[0][0]; else mask=(float*)NULL;

  field=u[0][0];
  nmesh=d[0]*d[1]*d[2];
  min=infinity, max=-infinity;
  for(i=0;i<nmesh;i++) if(!mask||mask[i]==1) {
    if(field[i]<min) min=field[i];
    if(field[i]>max) max=field[i];
  }
  fprintf(stderr,"data range from %g to %g\n",min,max);
  nout=0;
  for(i=0;i<d[0]*d[1]*d[2];i++) if(mask&&mask[i]!=1) {
    field[i]=min;  
    nout++;
  }
  indx=lvector(1,nmesh);
  indexx(nmesh,field,indx);
  fprintf(stderr,"nmesh %d nout %d\n",nmesh,nout);

  thr=dvector(0,b);
  vfr=dvector(0,b);
  *bin=dvector(0,b);
  for(j=0;j<=b;j++) {
    vfr[j]=par->lo+j*(par->hi-par->lo)/b;
    vfrac=.5-.5*errf(vfr[j]/sqrt(2.));
    nthre=(1.-vfrac)*(nmesh-nout)+nout+0.5;
    inthre=(long)nthre;
    fnthre=nthre-(float)inthre;
    if (nthre < 0) nthre=0.;
    thr[j]=field[indx[inthre]]+(field[indx[inthre+1]]-field[indx[inthre]])*fnthre;
    if(j>0) {
      dthr=thr[j]-thr[j-1];
      if (dthr < 1.e-8) dthr=1.e-8;
      (*bin)[j-1]=dthr;
      (*bin)[j]=dthr;
      /* fprintf(stderr,"%g \n",(*bin)[j-1]);*/
    }
    /* fprintf(stderr,"%g %d %g\n",nthre,inthre,field[indx[inthre]]);*/
  }
  for (j=0;j<=b;j++) {
    fprintf(stderr,"%g %g\n",thr[j],(*bin)[j]);
  }
  free_lvector(indx,1,nmesh);
  return thr;
}

void addminkowski(parlist*par,double***v,double***vsum,double***vsqr) {
  int i,j,k;

  for(i=0;i<3;i++) for(j=0;j<=par->bins;j++) for(k=0;k<=D;k++)
    vsum[i][j][k]+=v[i][j][k], vsqr[i][j][k]+=v[i][j][k]*v[i][j][k];
}

void finishminkowski(parlist*par,double***vsum,double***vsqr) {
  int i,j,k;

  for(i=0;i<3;i++) for(j=0;j<=par->bins;j++) for(k=0;k<=D;k++)
    sumtomean(&vsum[i][j][k],&vsqr[i][j][k],par->length);
}


/* mask definitions for subsets of the neighbourhood of a cell */
#define EVERYTHING 0x7ffffff
#define CENTER     0x0002000

/* Calculate the Minkowski functionals of field u and return them in v */

void minkowski(parlist*par,float***m,float***u,double*thr,double*bin,double***v) {
  int b=par->bins,*d=par->dim;
  double a=par->lattice;
  float*field,*mask;
  int i,j,k,x,y,z,xm,ym,zm,xp,yp,zp,idx[27],idy[27],id[27];
  int neighbour,lo,hi,inside,curvature[16],euler[256],o[8],pattern;
  int nvolume,nline,nsquare,ncube,nderivative;
  double min,max,u1,u2,u3,u11,u22,u33,u12,u23,u31,grad,inv[3];
  double mu,sig,tau;
  int line[3][3]={{0,13,22},{0,13,16},{0,13,14}};
  int square[3][5]={{0,13,14,16,17},{0,13,14,22,23},{0,13,16,22,25}};
  int cube[9]={0,13,22,16,25,14,23,17,26};

  /* access the three-dimensional arrays that hold the field and the
     mask as vectors.  This is generally considered a dirty trick, but
     speeds up calculations. */
  field=u[0][0];
  if(m) mask=m[0][0]; else mask=(float*)NULL;

  /* calculate the neighbourhood values used for testing whether a
     subset falls completely within the mask */
  for(i=0;i<3;i++) for(j=1;j<=2;j++) line[i][0]|=1<<line[i][j];
  for(i=0;i<3;i++) for(j=1;j<=4;j++) square[i][0]|=1<<square[i][j];
  for(j=1;j<=8;j++) cube[0]|=1<<cube[j];

  /* calculate the contribution of all possible configurations of a
     square to the curvature */
  for(i=0;i<16;i++) {
    for(j=0;j<4;j++) if(i&(1<<j)) o[j]=1; else o[j]=0;
    curvature[i]=0;
    curvature[i]+=(o[0]+o[1]+o[2]+o[3]);
    curvature[i]-=2*((o[0]|o[1])+(o[0]|o[2])+(o[2]|o[3])+(o[1]|o[3]));
    curvature[i]+=4*(o[0]|o[1]|o[2]|o[3]);
  }

  /* calculate the contribution of all possible configurations of a
     cube to the Euler characteristic */
  for(i=0;i<256;i++) {
    for(j=0;j<8;j++) if(i&(1<<j)) o[j]=1; else o[j]=0;
    euler[i]=0;
    euler[i]-=(o[0]+o[1]+o[2]+o[3]+o[4]+o[5]+o[6]+o[7]);
    euler[i]+=2*((o[0]|o[1])+(o[0]|o[2])+(o[2]|o[3])+(o[1]|o[3])+
		 (o[4]|o[5])+(o[4]|o[6])+(o[6]|o[7])+(o[5]|o[7])+
		 (o[0]|o[4])+(o[1]|o[5])+(o[2]|o[6])+(o[3]|o[7]));
    euler[i]-=4*((o[0]|o[1]|o[2]|o[3])+(o[4]|o[5]|o[6]|o[7])+
		 (o[0]|o[1]|o[4]|o[5])+(o[2]|o[3]|o[6]|o[7])+
		 (o[0]|o[2]|o[4]|o[6])+(o[1]|o[3]|o[5]|o[7]));
    euler[i]+=8*(o[0]|o[1]|o[2]|o[3]|o[4]|o[5]|o[6]|o[7]);
  }

  /* set the outside of the mask to the minimum field value */
  min=infinity, max=-infinity;
  for(i=0;i<d[0]*d[1]*d[2];i++) if(!mask||mask[i]==1) {
    if(field[i]<min) min=field[i];
    if(field[i]>max) max=field[i];
  }
  fprintf(stderr,"data range from %g to %g\n",min,max);
  for(i=0;i<d[0]*d[1]*d[2];i++) if(mask&&mask[i]!=1) field[i]=min;  

  /* set cell counters to zero. */
  nvolume=0;
  nline=nsquare=ncube=0;
  nderivative=0;

  /* empty Minkowski functionals */
  for(i=0;i<3;i++) for(j=0;j<=b;j++) for(k=0;k<=D;k++) v[i][j][k]=0;

  /* triple loop over grid.  Indices of neighbouring cells are
     determined as high up as possible in order to speed up things. */
  for(x=0;x<d[0];x++) {
    if(x==0) xm=d[0]-1; else xm=x-1; if(x==d[0]-1) xp=0; else xp=x+1;
    idx[ 0]=idx[ 1]=idx[ 2]=idx[ 3]=idx[ 4]=idx[ 5]=idx[ 6]=idx[ 7]=idx[ 8]=xm;
    idx[ 9]=idx[10]=idx[11]=idx[12]=idx[13]=idx[14]=idx[15]=idx[16]=idx[17]=x ;
    idx[18]=idx[19]=idx[20]=idx[21]=idx[22]=idx[23]=idx[24]=idx[25]=idx[26]=xp;
    for(i=0;i<27;i++) idx[i]*=d[1];

    for(y=0;y<d[1];y++) {
      if(y==0) ym=d[1]-1; else ym=y-1; if(y==d[1]-1) yp=0; else yp=y+1;
      idy[0]=idy[1]=idy[2]=idy[ 9]=idy[10]=idy[11]=idy[18]=idy[19]=idy[20]=ym;
      idy[3]=idy[4]=idy[5]=idy[12]=idy[13]=idy[14]=idy[21]=idy[22]=idy[23]=y ;
      idy[6]=idy[7]=idy[8]=idy[15]=idy[16]=idy[17]=idy[24]=idy[25]=idy[26]=yp;
      for(i=0;i<27;i++) idy[i]=(idy[i]+idx[i])*d[2];

      for(z=0;z<d[2];z++) if(!m||m[x][y][z]==1) {
	if(z==0) zm=d[2]-1; else zm=z-1; if(z==d[2]-1) zp=0; else zp=z+1;
	id[0]=id[3]=id[6]=id[ 9]=id[12]=id[15]=id[18]=id[21]=id[24]=zm;
	id[1]=id[4]=id[7]=id[10]=id[13]=id[16]=id[19]=id[22]=id[25]=z ;
	id[2]=id[5]=id[8]=id[11]=id[14]=id[17]=id[20]=id[23]=id[26]=zp;
	for(i=0;i<27;i++) id[i]+=idy[i];

	/* determine how the neighbourhood of the cell intersects with
	   the mask.  For each cell within the mask, one bit of an
	   integer value is set.  This makes it very convenient to
	   check whether a pre-defined subset of the neighbourhood
	   falls into the mask - one simply ANDs the mask integer with
	   a known value describing the subset.  */
	if(!mask)
	  neighbour=EVERYTHING;
	else
	  for(i=26,neighbour=0;i>=0;i--)
	    if(!mask||mask[id[i]]==1) neighbour|=1<<i;

	/* determine the range of values and the corresponding
	   threshold indices.  */
	for(i=1,min=infinity,max=-infinity;i<=8;i++)
	  if(!mask||mask[id[cube[i]]]==1) {
	    if(field[id[cube[i]]]<min) min=field[id[cube[i]]];
	    if(field[id[cube[i]]]>max) max=field[id[cube[i]]];
	  }
	lo=threshold_crofton(min,thr,b); if(lo<0) lo=0;
	hi=threshold_crofton(max,thr,b); if(hi>b) hi=b;

	/* calculate contribution to volume.  Note that at first, only
	   the bin at the threshold directly below the value is
	   incremented, the other bins are cumulated at the end of the
	   calculation.  This improves both speed and accuracy. */
	nvolume++;
        j=threshold_crofton(field[id[13]],thr,b);
	if(j>=0) v[0][j][0]++, v[1][j][0]++;

	/* calculate contribution to Minkowski functionals using Crofton's
	   formula */
	for(i=0;i<3;i++) if((neighbour&line[i][0])==line[i][0]) {
	  nline++;
	  for(j=lo;j<=hi;j++) {
	    for(k=1,inside=0;k<=2;k++)
	      if(field[id[line[i][k]]]>thr[j]) inside++;
	    if(inside==1) v[0][j][1]++;
	  }
	}

	for(i=0;i<3;i++) if((neighbour&square[i][0])==square[i][0]) {
	  nsquare++;
	  for(j=lo;j<=hi;j++) {
	    for(k=1,pattern=0;k<=4;k++)
	      if(field[id[square[i][k]]]>thr[j]) pattern|=1<<(k-1);
	    v[0][j][2]+=curvature[pattern];
	  }
	}	

	if((neighbour&cube[0])==cube[0]) {
	  ncube++;
	  for(j=lo;j<=hi;j++) {
	    for(k=1,pattern=0;k<=8;k++)
	      if(field[id[cube[k]]]>thr[j]) pattern|=1<<(k-1);
	    v[0][j][3]+=euler[pattern];
	  }
	}

	/* The following would be much nicer, but work only in theorie ...
	for(i=0;i<3;i++) if((neighbour&triangle[i][0])==triangle[i][0]) {
	  ntriangle++;
	  for(j=lo;j<=hi;j++) {
	    for(k=1,inside=0;k<=3;k++)
	      if(field[id[triangle[i][k]]]>thr[j]) inside++;
	    if(inside==1) v[0][j][2]++;
	    if(inside==2) v[0][j][2]--;
	  }
	}
	gamma=2*atan(sqrt(2))/pi;
	for(i=0;i<2;i++) if((neighbour&pyramid[i][0])==pyramid[i][0]) {
	  npyramid++;
	  for(j=lo;j<=hi;j++) {
	    for(k=1,inside=0;k<=4;k++)
	      if(field[id[pyramid[i][k]]]>thr[j]) inside++;
	    if(inside==1) v[0][j][3]+=2.-3.*gamma;
	    if(inside==2) v[0][j][3]+=2.-4.*gamma;
	    if(inside==3) v[0][j][3]+=2.-3.*gamma;
	  }
	}
	*/

	/* calculate contribution to Minkowski functionals using
	   Koenderink invariants.  This requires calculation of
	   derivatives up to second order and can therefore only be done
	   if the whole neighbourhood lies within the mask.  */
	if(neighbour==EVERYTHING) {
	  nderivative++;

	  /* calculate derivatives */
	  u1=(field[id[22]]-field[id[ 4]])/(2*a);
	  u2=(field[id[16]]-field[id[10]])/(2*a);
	  u3=(field[id[14]]-field[id[12]])/(2*a);
	  u11=(field[id[22]]-2*field[id[13]]+field[id[ 4]])/(a*a);
	  u22=(field[id[16]]-2*field[id[13]]+field[id[10]])/(a*a);
	  u33=(field[id[14]]-2*field[id[13]]+field[id[12]])/(a*a);
	  u12=(field[id[25]]-field[id[19]]-field[id[ 7]]+field[id[1]])/(4*a*a);
	  u23=(field[id[17]]-field[id[15]]-field[id[11]]+field[id[9]])/(4*a*a);
	  u31=(field[id[23]]-field[id[ 5]]-field[id[21]]+field[id[3]])/(4*a*a);
	  
	  /* add contribution to appropriate threshold */
	  lo=threshold_upperkoenderink(field[id[13]],thr,bin,b);
	  hi=threshold_lowerkoenderink(field[id[13]],thr,bin,b);
	  if(lo<=hi&&(grad=u1*u1+u2*u2+u3*u3)>1e-10) {
	    inv[0]=sqrt(grad) / 6;
	    inv[1]=
	      (2*u1*u2*u12-u11*(u2*u2+u3*u3)+
	       2*u2*u3*u23-u22*(u3*u3+u1*u1)+
	       2*u3*u1*u31-u33*(u1*u1+u2*u2)) 
	      / (2*(grad)) / (3*pi);
	    inv[2]=
	      (u1*u1*(u22*u33-u23*u23)+2*u1*u2*(u23*u31-u12*u33)+
	       u2*u2*(u33*u11-u31*u31)+2*u2*u3*(u31*u12-u23*u11)+
	       u3*u3*(u11*u22-u12*u12)+2*u3*u1*(u12*u23-u31*u22))
	      / pow(grad,1.5) / (4*pi);
	    for(j=lo;j<=hi;j++) for(k=0;k<3;k++) v[1][j][k+1]+=inv[k];
	  }	

	  /* also, gather some low-order statistics of the derivatives
	     (see below) */
	  v[2][0][2]+=u1+u2+u3, v[2][0][3]+=u1*u1+u2*u2+u3*u3;
	}

	/* use the space reserved for theoretical expectation values of
	   Minkowski functionals to gather some low-order statistics of
	   the random field.  */
	v[2][0][0]+=field[id[13]], v[2][0][1]+=field[id[13]]*field[id[13]];

      }	/* loop over z */
    } /* loop over y */
  } /* loop over x */

  fprintf(stderr,"%d volume, %d line, %d square, %d cube, %d derivative\n",
	  nvolume,nline,nsquare,ncube,nderivative);

  /* cumulate the volume bins */
  for(i=0;i<2;i++) for(j=b;j>0;j--) v[i][j-1][0]+=v[i][j][0];

  /* apply various normalization factors, such that the Minkowski
     functionals returned correspond to densities per volume */
  for(j=0;j<=b;j++) {
    /* volume */
    for(i=0;i<2;i++) v[i][j][0]/=nvolume;
    /* Crofton */
    v[0][j][1]/=nline  *a    *3. ;
    v[0][j][2]/=nsquare*a*a  *6.;
    v[0][j][3]/=ncube  *a*a*a*8.;
    /* Koenderink */
    for(k=1;k<=D;k++) v[1][j][k]/=nderivative*bin[j];
  }

  /* calculate low-order statistics of the field and its first
     derivatives.  Remember, we used the first few slots for the
     theoretical Minkowski to store sums.  */
  mu=v[2][0][0],sig=v[2][0][1],tau=v[2][0][3];
  mu/=nvolume;
  sig=(sig-nvolume*mu*mu)/(nvolume-1);
  tau/=D*nderivative;

  /* use Tomita's formulae to calculate theoretical expectation values */
  fprintf(stderr,"parameters %g %g %g\n",mu,sig,tau);
  for(j=0;j<=b;j++) {
    v[2][j][0]=.5-.5*errf((thr[j]-mu)/sqrt(2*sig));
    for(k=1;k<=D;k++) 
      v[2][j][k]=gam(1.+(D-k)/2.)*gam(1.+k/2.)/gam(1.+D/2.)*
	pow(tau/sig/2/pi,k/2.)*hermite(k-1,(thr[j]-mu)/sqrt(sig));
  }

  return;
}

/* Output the Minkowski functionals stored in v */

void outputminkowski(parlist*par,double*thr,double***v,double***e) {
  int i,j,k,b=par->bins;
  FILE*fp;

  /* output results to files */
  if(fileopen(&fp,par->mininame,"w")) {
    Explain(par,fp);
    for(j=0;j<=b;j++) {
      fprintf(fp,"%+13.8e",thr[j]);
      for(i=0;i<=2;i++) for(k=0;k<=D;k++) fprintf(fp," %+13.8e",v[i][j][k]);
      for(i=0;i<=2;i++) for(k=0;k<=D;k++) 
	if(par->length>1)
	  fprintf(fp," %+13.8e",e[i][j][k]);
	else
	  fprintf(fp," 0");
      fprintf(fp,"\n");
    }
    fileclose(&fp);
  }
}

/* find and return the index j of the threshold bin that lies directly
   below the value u, that is t_j<u and u<t_j+1.  Necessary for
   plaquette counts for Crofton's formula.  If none of the bins
   satisfies this, return a value below or at the top of the allowed
   range, respectively.  */
int threshold_crofton(double u,double*thr,int b) {
  int i,j,k;

  /* value lies outside or at top of valid threshold range */
  if(u<thr[0]) return -1; if(thr[b]<=u) return b;

  /* bracket the threshold index and search by partitioning */
  for(i=0,j=b;j-i>1;) {
    k=(i+j)/2; if(thr[k]<=u) i=k; else j=k;
  }
  return i;
}
/* Same as above, but use the lower end of the bin instead of its
   center, i.e. find the largest index j that still satisfies
   t_j-b_j/2<=u.  Necessary for the delta function for Koenderink
   invariants.  If none of the bins satisfies this, return a value
   outside the allowed range.  */
int threshold_lowerkoenderink(double u,double*thr,double*bin,int b) {
  int i,j,k;

  /* value lies outside or at top of valid threshold range */
  if(thr[0]-bin[0]/2>u)  return -1;
  if(thr[b]-bin[b]/2<=u) return b;

  /* bracket the threshold index and search by partitioning */
  for(i=0,j=b;j-i>1;) {
    k=(i+j)/2; if(thr[k]-bin[k]/2<=u) i=k; else j=k;
  }
  return i;
}
/* Exact opposite of the above, find the smallest index j that
   satisfies t_j+b_j/2>u.  Also used for the delta function for
   Koenderink invariants. */
int threshold_upperkoenderink(double u,double*thr,double*bin,int b) {
  int i,j,k;

  /* value lies outside or at bottom of valid threshold range */
  if(thr[b]+bin[b]/2<=u) return b+1;
  if(thr[0]+bin[0]/2>u)  return 0;

  /* bracket the threshold index and search by partitioning */
  for(i=0,j=b;j-i>1;) {
    k=(i+j)/2; if(thr[k]+bin[k]/2>u) j=k; else i=k;
  }
  return j;
}

