/* getErrorProfile - Generate the error profile given a bam alignment reference. */
#include "common.h"
#include <ctype.h>
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "fastq.h" //methods for converting phred scores
#include "dnaLoad.h"
#include "dnaseq.h"
#include <stdbool.h>

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
    case 'N': return 4;
    case 'a': return 0;
    case 'c': return 1;
    case 'g': return 2;
    case 't': return 3;
    case 'n': return 4;
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
    case 4: return 'N';
    default:
      errAbort("Invalid base index %d\n",index);
      return 'N'; //make compiler happy...
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

bool invalid(char c){
  //check that c is in A C G T
  switch(toupper(c)){
  case 'A': return false;
  case 'C': return false;
  case 'G': return false;
  case 'T': return false;
  default: return true;
  }
}

bool hasAny(unsigned long long ** c){
  int i,j;
  for(i=0;i<4;i++){
    for(j=0;j<5;j++){
      if(c[i][j] > 0) return true;
    }
  }
  return false;
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
  //For every position store the off diagonal substitution matrix + the frequency of shifting to N
  unsigned long long **** mutation = (unsigned long long ****)needMem(sizeof(unsigned long long ***)*readLen);//will hold a position specific hist
  for(i=0;i<readLen;i++)
    {
      
      mutation[i] = (unsigned long long ***)needMem(sizeof(unsigned long long**) *61); //quality scores could range from 0-60
      for(j=0;j<=60;j++){
	mutation[i][j] = (unsigned long long **)needMem(sizeof(unsigned long long *)*4); //4 reference bases
	for(k=0;k<4;k++){
	  mutation[i][j][k] = (unsigned long long *)needMem(sizeof(unsigned long long)*5); //|{A,C,G,T,N}| = 5 possible read bases
	}
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
      unsigned int pseudo = 1; //pseudo count to give to all non N characters with at least one phred score in the sequence
      int leftPos = atoi(words[3])-1; //0 based left most position
      unsigned int flag = atoi(words[1]);
      unsigned int strand_mask = 16;//0x10 in hex
      //int paired_mask = 1;
      char refChar;
      char readChar;
      unsigned pscore = 0;
      unsigned pos = strlen(seq); 
      unsigned int forward = flag & strand_mask; //bit mask everything other than 0x0010
      if (forward)
	{
	  for(i=0;i<readLen;i++)
	    {
	      ref = hashFindVal(refHash, chrom);
	      if(!ref) errAbort("Sequence name %s not found in reference\n",chrom);
	      refChar = toupper(ref->dna[leftPos+i]);
	      //else refChar = complementSingle(ref->dna[leftPos+readLen-1-i]);
	      readChar = toupper(seq[i]);
	      if(invalid(refChar)) continue; //skip non-nucleotides in the reference
	      //handle phred histogram
	      if (phred33) pscore = phred33ToPhred(score[i]);
	      else if (phred64) pscore = phred64ToPhred(score[i]);
	      if(mutation[i][pscore][0][0] < pseudo){ //add in pseudocounts for everything other than N
		int tmp1,tmp2;
		for(tmp1=0;tmp1<4;tmp1++){
		  for(tmp2=0;tmp2<4;tmp2++){
		    mutation[i][pscore][tmp1][tmp2]=pseudo;
		  }
		}
	      }
	      mutation[i][pscore][baseIndex(refChar)][baseIndex(readChar)]++;
	    } //loop over read length
	}
      else
	{ //sequence and score are complemented and/or reversed
	  for(i=0; i <readLen ; i++)
	    {
	      --pos; //start at last index + 1, decrement at beginning
	      //i stores the 0- desired readlen index
	      //pos stores the true readLen - (true-desired) index, this gets us the end of a read
	      ref = hashFindVal(refHash, chrom);
	      if(!ref) errAbort("Sequence name %s not found in reference\n",chrom);
	      refChar = toupper(ref->dna[leftPos+i]);
	      readChar = toupper(seq[pos]);

	      if(invalid(refChar)) continue; //skip non-nucleotides in the reference
	      
	      if (phred33) pscore = phred33ToPhred(score[pos]);
	      else if (phred64) pscore = phred64ToPhred(score[pos]);

	      if(mutation[i][pscore][0][0] < pseudo){ //add in pseudocounts for everything other than N, haven't done this already
		int tmp1,tmp2;
		for(tmp1=0;tmp1<4;tmp1++){
		  for(tmp2=0;tmp2<4;tmp2++){
		    mutation[i][pscore][tmp1][tmp2]=pseudo;
		  }
		}
	      }

	      mutation[i][pscore][baseIndex(complementSingle(refChar))][baseIndex(complementSingle(readChar))]++; //increment for sequenced error

	    }




	}

    }//end while

  //print out the error histogram file.
  printf("#Pos\tPhred\tA->A\tA->C\tA->G\tA->T\tA->N\tC->A\tC->C\tC->G\tC->T\tC->N\tG->A\tG->C\tG->G\tG->T\tG->N\tT->A\tT->C\tT->G\tT->T\tT->N\n");
  int p;
  for(i=0;i<readLen;i++){
      for(p=0;p<=60;p++){
	if(hasAny(mutation[i][p])){ //if there are any phred scores to report here...
	  printf("%d\t%d",i,p); //print pos\tphred
	  for(j=0;j<4;j++){
	    for(k=0;k<5;k++)
	      printf("\t%llu", mutation[i][p][j][k]);
	  }
	  printf("\n");
	}
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

