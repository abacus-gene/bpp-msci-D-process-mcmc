/* minsub.h */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#define square(a) ((a)*(a))
#define min2(a,b) ((a)<(b)?(a):(b))
#define max2(a,b) ((a)>(b)?(a):(b))
#define swap2(a,b,y) { y=a; a=b; b=y; }

void error2(char* message);
FILE* gfopen(char* filename, char* mode);
int scanfile(FILE* fin, int* nrecords, int* nx, int* HasHeader, char line[], int ifields[]);
int matout(FILE* file, double x[], int n, int m);
int matout2(FILE* fout, double x[], int n, int m, int wid, int deci);

int zero(double x[], int n);
double sum(double x[], int n);
int xtoy(double x[], double y[], int n);
int abyx(double a, double x[], int n);
int identity(double x[], int n);
double distance(double x[], double y[], int n);
double innerp(double x[], double y[], int n);
double norm(double x[], int n);

int H_end (double x0[], double x1[], double f0, double f1, double e1, double e2, int n);

double LineSearch2 (double(*fun)(double x[],int n), double *f, double x0[], 
    double p[], double h, double limit, double e, double space[], int n);

int ming2 (FILE *fout, double *f, double (*fun)(double x[], int n),
    int (*dfun)(double x[], double *f, double dx[], int n),
    double x[], double xb[][2], double space[], double e, int n);

#define MAXNFIELDS 320000
