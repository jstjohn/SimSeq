/*
	Author John St. John
	07/12/2010

	Define a few structures and functions for
	reading and writing fastq format.

	Also defines a few functions for converting
	between various phred and ascii phred
	format.
*/

#include <ctype.h>
#include "common.h"
#include "linefile.h"
#include <math.h>
#include <string.h>
#include "fastq.h"


/* Function Prototypes */
/*
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
*/


struct fastqItem * allocFastqItem()
{
	struct fastqItem *fq = (struct fastqItem *) malloc(sizeof(struct fastqItem ));
	fq->id = (char *)malloc(sizeof(char) * MAX_ID_LENGTH);
	fq->seq = (char *)malloc(sizeof(char) * MAX_SEQ_LENGTH);
	fq->score = (char*)malloc(sizeof(char) * MAX_SEQ_LENGTH);
	return fq;
}

void freeFastqItem(struct fastqItem * fq)
{
	if(fq)
	{
		if(fq->id)
		{
			free(fq->id);
			fq->id = NULL;
		}
		if(fq->seq)
		{
			free(fq->seq);
			fq->seq = NULL;
		}
		if(fq->score)
		{
			free(fq->score);
			fq->score = NULL;
		}
		free(fq);
		fq = NULL;
	}
}


void printFastqItem(FILE *fp, struct fastqItem *fq)
{
	int i;
	fprintf(fp,"@%s\n",fq->id);
	for(i=0; i < fq->len; i++)
	{
		fprintf(fp,"%c",fq->seq[i]);
	}
	fprintf(fp,"\n+\n");
	for(i=0; i< fq->len; i++)
	{
		fprintf(fp,"%c",fq->score[i]);
	}
	fprintf(fp,"\n");
}



boolean fastqItemNext(struct lineFile *lf, struct fastqItem *fq)
{
	int start = 0;
	char *line = NULL;
	boolean gotId = FALSE;
	boolean gotSeq = FALSE;
	boolean gotScore = FALSE;
	//Get the ID line
	gotId = lineFileNextReal(lf,&line); //grab next non-blank non-comment line, should have @
	if (! gotId) return FALSE;
	line = skipLeadingSpaces(line);
	eraseTrailingSpaces(line);
	if (line[0] != '@')
		errAbort("ERROR: %s doesn't seem to be fastq format.  "
				"Expecting '@' start of line %d, got %c.",
				lf->fileName, lf->lineIx, line[0]);
	//copyID(line,fq->id,1,strlen(line));
	//sprintf(fq->id, "%s",line+1);
	strcpy(fq->id,line+1);
	//get the Sequence
	int seqLen = 0;
	while (TRUE)
	{
		gotSeq = lineFileNextReal(lf,&line);
		if (! gotSeq) lineFileUnexpectedEnd(lf);
		line = skipLeadingSpaces(line);
		eraseTrailingSpaces(line);
		if (line[0] == '+') break;
		int partLen = strlen(line);
		start = seqLen;
		seqLen += partLen;
		//strAdd(line,fq->seq,0,partLen,start);
		//sprintf((fq->seq)+start,"%s",line);
		strcpy((fq->seq)+start,line);
		if (line[0] == '@')
			errAbort("ERROR: %s doesn't seem to be fastq format.  "
					"Expecting sequence at line %d, got %c.",
					lf->fileName, lf->lineIx, line[0]);
	}
	fq->len=seqLen;
	//skip the +, just checked in previous loop

	//get the Score, check that len(score) == len(seq)
	seqLen=0;
	while(TRUE)
	{
		//the problem is that @ and + are valid score characters
		gotScore = lineFileNextReal(lf,&line);
		if(! gotScore) break;
		line = skipLeadingSpaces(line);
		eraseTrailingSpaces(line);
		int partLen = strlen(line);
		if (line[0] == '@' && seqLen+partLen > fq->len)
		{
			lineFileReuse(lf); //start with this line next time
			break;
		}
		if (line[0] == '+' && seqLen+partLen > fq->len)
			errAbort("ERROR: %s doesn't seem to be fastq format.  "
					"Expecting sequence at line %d, got %c.",
					lf->fileName, lf->lineIx, line[0]);
		//strAdd(line,fq->score,0,partLen,start);
		//sprintf((fq->score)+start,"%s",line);
		start=seqLen;
		seqLen += partLen;
		strcpy((fq->score)+start,line);
	}
	if (seqLen != fq->len)
		errAbort("ERROR: %s has sequences and score strings that "
				"are not the same length. Problem at line %d.",
				lf->fileName, lf->lineIx);
	return TRUE;
}

void convPhred33ToPhred64( struct fastqItem *fq )
{
	phred33ToPhred64(fq->seq,fq->len);
}//end phred33To64

void convPhred64ToPhred33(struct fastqItem *fq)
{
	phred64ToPhred33(fq->seq,fq->len);
}

void phred33ToPhred64( char * p33, int l )
{
	int i;
	for(i=0;i<l;i++)
	{
		p33[i] = phredToPhred64(phred33ToPhred(p33[i]));
	}
}

void phred64ToPhred33( char * p64, int l)
{
	int i;
	for(i=0;i<l;i++)
	{
		p64[i] = phredToPhred33(phred64ToPhred(p64[i]));
	}
}

char phredToPhred33( int p )
{
	if (p > MAX_PHRED) p=MAX_PHRED;
	else if (p < MIN_PHRED) p=MIN_PHRED;
	return ((char) (p + 33));
}

int phred33ToPhred( char p )
{
	return ((int)p) - 33;
}

int phred64ToPhred( char p )
{
	return ((int)p) - 64;
}

char phredToPhred64( int p )
{
	if (p > MAX_PHRED) p=MAX_PHRED;
	else if (p < MIN_PHRED) p=MIN_PHRED;
	return ((char) (p + 64));
}

int doubleToPhred( double p )
/* formula: -10 log10(p) */
{
	double res = -10.0 * log10(p);
	if (res > MAX_PHRED) res=MAX_PHRED;
	else if (res < MIN_PHRED) res=MIN_PHRED;
	return ((int) (res +0.5));  //guarenteed >= 0
}

double phredToDouble( int p )
/* formula: 10^(-p/10)  */
{
	return 1.0/pow(10,((double)p)/10.0);
}

/* Some Functions that can now be made from  a combination of existing functions */

double phred33ToDouble( char p )
{
	return phredToDouble(phred33ToPhred(p));
}

double phred64ToDouble( char p )
{
	return phredToDouble(phred64ToPhred(p));
}



