#include <stdio.h>
#include <stdlib.h>

void die(char*);

int main(int nar,char*ar[]) {

  int i;
  char in[4],out[4];

  /* loop till end of file is reached */
  for(i=0;fread((void*)in,sizeof(char),4,stdin)==4;i++) {
    out[3]=in[0], out[2]=in[1], out[1]=in[2], out[0]=in[3];
    if(fwrite((void*)out,sizeof(char),4,stdout)!=4)
      die("Could not write to standard output\n");
  }

  /* fprintf(stderr,"read and reversed %d values\n",i); */ 

  exit(0);
 
}

void die(char*reason) {
  fprintf(stderr,"death caused by: %s\n",reason);
  exit(666);
}
