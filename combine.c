//
// usage: ./combine file1 file2 file3 ... > output
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXCHARS 666
#define MAXENTRY 42

int main(int argc,char*argv[]) {
  int i,j,m,n=argc-1,flag;
  char**line,*next;
  double*sum,*sqr,value;
  FILE**fp;

  /* allocate storage space */
  fp=(FILE**)malloc(sizeof(FILE*)*n);
  line=(char**)malloc(sizeof(char*)*n);
  for(i=0;i<n;i++) line[i]=(char*)malloc(sizeof(char)*MAXCHARS);
  sum=(double*)malloc(sizeof(double)*MAXENTRY);
  sqr=(double*)malloc(sizeof(double)*MAXENTRY);

  /* open all input files */
  for(i=0;i<n;i++) fp[i]=fopen(argv[i+1],"r");

  /* set the flag and loop until the flag is unset by hitting an EOF */
  flag=1; while(flag)

    /* loop over input files */
    for(i=0;i<n;i++) {

      /* read one line and break if EOF is reached */
      if(!fgets(line[i],MAXCHARS,fp[i])) { flag=0; break; }

      if(line[i][0]=='#') {

	/* commented lines are output, provided they are unique */
	if(i==0||strcmp(line[i],line[i-1])) printf("%s",line[i]);

      } else {

	/* uncommented lines are parsed, and the values averaged and output */

	if(i==0) for(j=0;j<MAXENTRY;j++) sum[j]=0;

	for(next=line[i],j=0;(*next!='\n')&&(j<MAXENTRY);j++) {
	  value=strtod(next,&next);
	  sum[j]+=value; sqr[j]+=value*value;
	}
	m=j;

	if(i==n-1) {
	  if(n>1) for(j=1;j<=m/2;j++)
	    sum[j]/=n, sqr[j]=sqrt((sqr[j]-sum[j]*sum[j]*n)/(n-1));

	  printf("%+13.8e",sum[0]/n);
	  for(j=1;j<=m/2;j++) printf(" %+13.8e",sum[j]);
	  for(j=1;j<=m/2;j++)
	    if(n>1) printf(" %+13.8e",sqr[j]); else printf(" 0");
	  printf("\n");
	}
	
      }

    }

  exit(0);
}

