
/* Declarations */
extern int jacobi(double**,int,double*,double**);
extern void sort(unsigned long, float*);
extern void indexx(unsigned long, float*, unsigned long*);
extern float newtonroot(float);
extern double gasdev();
extern void sumtomean(double*,double*,int);
extern double errf(double);
extern double gaussf(double);
extern double gam(double);
extern double gaminc(double,double);
extern double hermite(int,double);
extern double omega(int);
extern double qromb(double(*)(double),double,double);
extern void polint(double*,double*,int,double,double*,double*);
extern double trapzd(double(*)(double),double,double,int);
extern int fileopen(FILE**,char*,char*);
extern void fileclose(FILE**);
extern void Mecker(char*);
extern int*ivector(int,int);
extern unsigned long*lvector(int,int);
extern float*vector(int,int);
extern double*dvector(int,int);
extern char**cmatrix(int,int,int,int);
extern int**imatrix(int,int,int,int);
extern float**matrix(int,int,int,int);
extern double**dmatrix(int,int,int,int);
extern int***icube(int,int,int,int,int,int);
extern float***cube(int,int,int,int,int,int);
extern double***dcube(int,int,int,int,int,int);
extern void free_ivector(int*,int,int);
extern void free_lvector(unsigned long*,int,int);
extern void free_vector(float*,int,int);
extern void free_dvector(double*,int,int);
extern void free_cmatrix(char**,int,int,int,int);
extern void free_imatrix(int**,int,int,int,int);
extern void free_matrix(float**,int,int,int,int);
extern void free_dmatrix(double**,int,int,int,int);
extern void free_icube(int***,int,int,int,int,int,int);
extern void free_cube(float***,int,int,int,int,int,int);
extern void free_dcube(double***,int,int,int,int,int,int);

/* for random number generator */
struct SEEDS { int i0, i1, i2; };
extern struct SEEDS seeds;
extern double rand1(void);
