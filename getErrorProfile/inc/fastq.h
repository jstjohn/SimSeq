/*
	Author John St. John
	07/12/2010

	Define a few structures and functions for
	reading and writing fastq format.

	Also defines a few functions for converting
	between various phred and ascii phred
	format.
*/

#ifndef _FASTQ_H
#define _FASTQ_H
#include <ctype.h>
#include "common.h"
#include "linefile.h"
#include <math.h>
#include <string.h>

/* Global Constants */
#define MAX_PHRED 93
#define MIN_PHRED 0
#define PHRED_33 0
#define PHRED_64 1
#define MAX_ID_LENGTH 256
#define MAX_SEQ_LENGTH 1000

/* Data Structures */
struct fastqItem {
	char *id;
	char *seq;
	char *score;
	int len;
};

/* Function Prototypes */
void convPhred33ToPhred64( struct fastqItem * );
void convPhred64ToPhred33( struct fastqItem * );
void phred33ToPhred64( char *, int );
void phred64ToPhred33( char *, int );
int  phred64ToPhred( char );
int  phred33ToPhred( char );
double phredToDouble( int );
double phred33ToDouble( char );
double phred64ToDouble( char );
int doubleToPhred( double );
char phredToPhred33( int );
char phredToPhred64( int );
boolean fastqItemNext(struct lineFile *, struct fastqItem *);
struct fastqItem * allocFastqItem();
void freeFastqItem(struct fastqItem * fq);
void printFastqItem(FILE *fp, struct fastqItem *fq);

#endif

