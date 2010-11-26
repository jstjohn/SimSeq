/* getErrorProfile - Generate the error profile given a bam alignment reference. */
#include "common.h"
#include <ctype.h>
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "fastq.h" //methods for converting phred scores
#include "dnaLoad.h"
#include "dnaseq.h"

//static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

static struct hash *refHash;
static char * refName;
static char * sam;
static int readLen;
static boolean phred33 = FALSE;
static boolean phred64 = FALSE;
void usage()
/* Explain usage and exit. */
{
  errAbort(
	   "getErrorProfile - Generate the error profile given a sam sequence on stdin\n"
	   "usage:\n"
	   "   samtools view file.bam | getErrorProfile [options] >error_profile.txt\n"
	   "options:\n"
	   "   -ref=FILE_NAME (required)\n"
	   "   -readLen=NUMBER (required)\n"
	   "   -sam=FILE_NAME (optional) headless sam file, otherwise use stdin by default.\n"
	   "   -phred33 if input sam has ascii phred-33 quality scores (required)\n"
	   "   -phred64 if input sam has ascii phred-64 quality scores (other option)\n"
	   );
}

static struct optionSpec options[] = {
  {"ref",OPTION_STRING},
  {"readLen",OPTION_INT},
  {"phred33",OPTION_BOOLEAN},
  {"phred64",OPTION_BOOLEAN},
  {"sam",OPTION_STRING},
  {NULL, 0},
};

char complementSingle(char base)
{
  switch(base)
    {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';
    default: return 'N';
    }
}


int baseIndex(char base)
/*bases indexed alphabetically*/
{
  switch(base)
    {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    case 'a': return 0;
    case 'c': return 1;
    case 'g': return 2;
    case 't': return 3;
    default: return -1;
    }
}

char baseFromIndex(int index)
{
  switch(index)
    {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    default:
      errAbort("Invalid base index %d\n",index);
      return 'N'; //dummy to make compiler happy
    }
}

char baseFromMutationIndex(int from, int to)
{
  switch(from)
    {
    case 0:
      switch (to)
	{
	case 0: return 'C';
	case 1: return 'G';
	case 2: return 'T';
	default: return 'N';
	}
    case 1:
      switch (to)
	{
	case 0: return 'A';
	case 1: return 'G';
	case 2: return 'T';
	default: return 'N';
	}
    case 2:
      switch (to)
	{
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'T';
	default: return 'N';
	}
    case 3:
      switch (to)
	{
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	default: return 'N';
	}
    default: return 'N';
    }

}

int mutationIndex(char from, char to)
{
  to = toupper(to);
  switch(toupper(from))
    {
    case 'A': return baseIndex(to)-1;
    case 'C':
      if (to == 'A') return 0;
      else return baseIndex(to)-1;
    case 'G':
      if (to == 'T')
	return baseIndex(to)-1;
      else return baseIndex(to);
    case 'T':
      return baseIndex(to);
    default: 
      errAbort("Invalid 'from' base %c\n",from);
      return -1;
    }

}

char *strrev(char *s,int n)
{
  int i=0;
  while (i<n/2)
    {
      *(s+n) = *(s+i);       //uses the null character as the temporary storage.
      *(s+i) = *(s + n - i -1);
      *(s+n-i-1) = *(s+n);
      i++;
    }
  *(s+n) = '\0';

  return s;
}


void getErrorProfile()
/* getErrorProfile - Generate the error profile given a bam alignment reference. */
{
  //open stdin and read in sam format
  struct lineFile  *samlf = NULL;
  if (sam == NULL) samlf = lineFileStdin(TRUE);
  else samlf = lineFileOpen(sam,TRUE);
  //struct dnaSeq *dnaLoadAll(char *fileName);
  //struct dnaSeq
  ///* A dna sequence in one-character per base format. */
  //{
  //	struct dnaSeq *next;  /* Next in list. */
  //	char *name;           /* Name of sequence. */
  //	DNA *dna;             /* Sequence base by base. */
  //	int size;             /* Size of sequence. */
  //	Bits* mask;           /* Repeat mask (optional) */
  //};
  //
  //struct hash *dnaSeqHash(struct dnaSeq *seqList);
  struct dnaSeq *seqs = dnaLoadAll(refName);
  refHash = dnaSeqHash(seqs);
  //char *line = NULL;
  char *words[12];
  int i,j,k;
  unsigned long long *** phred = (unsigned long long ***)needMem(sizeof(unsigned long long **)*readLen);//will hold a position specific hist
  unsigned long long *** mutation = (unsigned long long ***)needMem(sizeof(unsigned long long **)*readLen);//will hold a position specific hist
  for(i=0;i<readLen;i++)
    {
      phred[i] = (unsigned long long **)needMem(sizeof(unsigned long long*) * 4); //phred scores range from 0-93
      mutation[i] = (unsigned long long **)needMem(sizeof(unsigned long long*) *4);
      for(j=0;j<4;j++)
	{
	  phred[i][j] = (unsigned long long *)needLargeZeroedMem(sizeof(unsigned long long) * 94); //phred scores range from 0-93
	  mutation[i][j] = (unsigned long long *)needLargeZeroedMem(sizeof(unsigned long long) * 3); //3 bases each can change to
	}
    }
  unsigned int junk = 0;
  //chr2_485_762_0:0:0_0:0:0_1619fe	99	chr1	10022	1	76M	=	10044	98	CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:412	XM:i:0	XO:i:0	XG:i:0	MD:Z:76
  while(lineFileChopNextTab(samlf, words, 12))
    {   	//only interested in first 11 positions, 12th will hold rest?
      DNA *seq = words[9];
      junk++;
      struct dnaSeq *ref = NULL;
      char *score = words[10];
      char *chrom = words[2];
      int leftPos = atoi(words[3])-1; //0 based left most position
      int flag = atoi(words[1]);
      int strand_mask = 16;//10000 in binary
      //int paired_mask = 1;
      char refChar;
      char readChar;
      unsigned pscore = 0;
      unsigned pos = strlen(seq); 
      int forward = flag & strand_mask; //bit mask everything other than 0x0010
      if (forward)
	{
	  for(i=0;i<readLen;i++)
	    {
	      ref = hashFindVal(refHash, chrom);
	      if(!ref) errAbort("Sequence name %s not found in reference\n",chrom);
	      refChar = toupper(ref->dna[leftPos+i]);
	      //else refChar = complementSingle(ref->dna[leftPos+readLen-1-i]);
	      readChar = toupper(seq[i]);
	      if (refChar == 'N' || readChar == 'N') continue; //skip N's
	      //handle phred histogram
	      if (phred33) pscore = phred33ToPhred(score[i]);
	      else if (phred64) pscore = phred64ToPhred(score[i]);
	      phred[i][baseIndex(readChar)][pscore]++; //increment hist
	      //handle mutation index
	      if (refChar != readChar)
		mutation[i][baseIndex(refChar)][mutationIndex(refChar,readChar)]++;
	    } //loop over read length
	}
      else
	{ //sequence and score are complemented and/or reversed
	  for(i=readLen-1; i >= 0; i--)
	    {
	      --pos; //start at last index + 1, decrement at beginning
	      //i stores the 0- desired readlen index
	      //pos stores the true readLen - (true-desired) index, this gets us the end of a read
	      ref = hashFindVal(refHash, chrom);
	      if(!ref) errAbort("Sequence name %s not found in reference\n",chrom);
	      refChar = toupper(ref->dna[leftPos+i]);
	      readChar = toupper(seq[pos]);

	      if (refChar == 'N' || readChar == 'N') continue; //skip N's

	      if (phred33) pscore = phred33ToPhred(score[pos]);
	      else if (phred64) pscore = phred64ToPhred(score[pos]);
	      phred[i][baseIndex(complementSingle(readChar))][pscore]++; //increment hist for sequenced base
	      if (refChar != readChar)
		mutation[i][baseIndex(complementSingle(refChar))][mutationIndex(complementSingle(refChar),complementSingle(readChar))]++; //increment for sequenced error

	    }




	}

    }//end while

  //print out the read file.
  /**
   * The error profile will have to be interesting
   * We need a position specific mutation spectrum
   * and a position specific quality histogram.
   *
   * The First line of the file should specify the number of positions
   * of reads.
   *
   * The mutation spectrum should be on a single line for each position
   * A line will be the integer number of observances of each error type.
   * It should be in order representing a 4x4 [A,C,G,T] matrix read
   * across rows from top left to bottom right. Thus the numbers will
   * be [a->c,a->g,a->t,c->a,c->g,c->t,...]. For simplicity these numbers
   * will be tab delineated. The first number on each line will be the
   * position [0,...,n].
   *
   * The Quality histogram will be harder to parse because there could be
   * quite a few quality values in the histogram, and the number of those
   * values could vary based on base and position.
   * Position\tBase\tNumber(N)\tItem1:Count1,...,ItemN:CountN
   **/
  printf("#Number of reads:\n");
  printf("%d\n\n",readLen);
  printf("#Mutation Spectrum:\n");
  for(i=0;i<readLen;i++)
    {
      printf("%d",i);
      for(j=0;j<4;j++)
	{
	  for(k=0;k<3;k++)
	    printf("\t%llu", mutation[i][j][k]);
	}
      printf("\n");
    }
  printf("\n#Quality Histograms:\n");
  for(i=0;i<readLen;i++)
    {
      for(j=0;j<4;j++)
	{
	  unsigned int nonZero[94];
	  int count = 0;
	  for(k=0;k<94;k++)
	    if(phred[i][j][k] != 0)
	      nonZero[count++] = k;
	  printf("%d\t%c\t%d\t",i,baseFromIndex(j),count);
	  for(k=0;k<count-1;k++) printf("%d:%llu,",nonZero[k],phred[i][j][nonZero[k]]);
	  printf("%d:%llu\n",nonZero[count-1],phred[i][j][nonZero[count-1]]); //last element no comma
	}
    }
  printf("\n\n");
}

int main(int argc, char *argv[])
/* Process command line. */
{
  optionInit(&argc, argv, options);
  //optionMustExist("readLen");
  //optionMustExist("ref");
  readLen = optionInt("readLen",0);//should crash if not set
  refName = optionVal("ref",NULL);
  phred33 = optionExists("phred33");
  phred64 = optionExists("phred64");
  sam = optionVal("sam",NULL);
  if(phred33 && phred64)
    {
      errAbort("Options phred33 and phred64 are mutually exclusive!\n");
    }
  if(!(phred33 || phred64)) usage();
  if(!readLen || !refName) usage();
  getErrorProfile();
  return 0;
}

