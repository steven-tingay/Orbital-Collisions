#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <time.h>

#define MAX_PDF 1000
#define X_LOWER -3.0
#define X_UPPER 2.0
#define XV_LOWER -2.0
#define XV_UPPER 7.0

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


typedef struct fragment 
{
	long double areaMass;
	long double mass;
	long double deltaV;
	long double charLen;
	long double csArea;
}fragment;	
typedef struct linkedList
{
	struct fragment** head;
	int numbBins;
	int numbNodes;
	long double binCLen;
}linkedList;

FILE *xfopen (const char *fn, const char *mode);
int badmode (const char *s);
int xfclose (FILE *fp);
int xfexists (char *fn);
char *fnwoext (char *nm, char *fn);

void singleSim(long double, long double,int);

void multiSim(long double,long double,int);

void printFile(linkedList**);

long double initiateAM(linkedList** list,int);

void assignDelV(linkedList**,int);

void addVelVal(long double**,long double,int,int);

long double assignAM(linkedList**,int,int,long double**,long double**);

void deAllocLL(linkedList**);

linkedList** createLL(long double,double);

long double charLenUp(long double,double);

void print(linkedList**);

int createScript(long double**,int,char*,char*);

void createDArray(int,int,long double**);

void deAllocA(long double**,int);

void createEq(long double*);

void getRBParam(long double,long double*);

void getSCParam(long double,long double*);

void addVal(long double**,long double*,int);

void setStanVal(long double**,int,long double, long double);

void normalise(long double**,int);

void createCDF(long double**,long double**, int,long double,long double);

long double sampleCDF(long double**);
