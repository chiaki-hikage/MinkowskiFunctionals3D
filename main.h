#undef DEBUG

#define NESYO(flag) (flag)?"yes":"no"

/* parameters read from command line */
typedef struct parlist {

  /* miscellaneous */
  int length;			/* number of realizations of field */
  int inter;			/* interpolation factor */
  double smooth;		/* width of the smoothing kernel */
  double contamin;		/* maximum contamination of survey region */

  /* data */
  int dim[3];			/* size of grid in pixels */
  double lattice;		/* lattice constant */
  char*dataname,*maskname;	/* input files for data and mask */

  /* Minkowski functionals */
  double lo,hi;			/* threshold range */
  int bins; double width;	/* number and width of bins */
  char*mininame;		/* output file for Minkowski functionals */

} parlist;

extern void Explain(parlist*,FILE*);
