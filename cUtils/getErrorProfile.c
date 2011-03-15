/* getErrorProfile - Generate the error profile given a bam alignment reference. */
#include "common.h"
#include <ctype.h>
#include <math.h>
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "fastq.h" //methods for converting phred scores
#include "dnaLoad.h"
#include "dnaseq.h"
#include <stdbool.h>
#include "sam.h"
#define FIRST_PREV_PHRED (0)

//static char const rcsid[] = "$Id: newProg.c,v 1.30 2010/03/24 21:18:33 hiram Exp $";

static struct hash *refHash;
static char * refName;
static char * sam;
static int readLen;
static boolean phred33 = FALSE;
static boolean phred64 = FALSE;
static int plist[10];
static int plen=8;
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
      "   -pscores=p1,...,pN up to 10 phred scores to model separated by comma with no spaces(default: '0,5,10,15,20,25,30,35,40') \n"
      "   -sam=FILE_NAME (optional) headless sam file, otherwise use stdin by default.\n"
      "   -phred33 if input sam has ascii phred-33 quality scores (required)\n"
      "   -phred64 if input sam has ascii phred-64 quality scores (other option)\n"
  );
}

static struct optionSpec options[] = {
    {"ref",OPTION_STRING},
    {"pscores",OPTION_STRING},
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

//given a pscore, looks at the global list of
//pscores and finds the nearest one.
unsigned int nearestIdxFromPscore(unsigned int pscore){
  int i;
  for(i=0;i<plen-1;i++){
    if (pscore == plist[i]) return i;
    if (pscore < plist[i+1]){ //find out which its closer to i or i+1
      if(abs(pscore - plist[i])<abs(pscore-plist[i+1]))
        return i;
      else
        return i+1;
    }
  }
  return plen-1;
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

void intrev(int *s, int n){
  int i=0;
  while (i<n/2){
    *(s+n) = *(s+i);
    *(s+i) = *(s + n - i -1);
    *(s+n-i-1) = *(s+n);
    i++;
  }
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

void getErrorProfile(char *reflist)
/* getErrorProfile - Generate the error profile given a bam alignment reference. */
{

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
  unsigned long long insertions = 0; //for now just count these and write to stderr as stats
  unsigned long long deletions = 0;
  struct dnaSeq *seqs = dnaLoadAll(refName);
  refHash = dnaSeqHash(seqs);
  int i,j,k,l;
  //For every position store the complete substitution matrix A->A...T->N
  unsigned long long **** mutation = (unsigned long long ****)needMem(sizeof(unsigned long long ***)*readLen);//will hold a position specific hist
  for(i=0;i<readLen;i++)
  {
    mutation[i] = (unsigned long long ****)needMem(sizeof(unsigned long long***) *plen); //previous quality scores
    for(j=0;j<plen;j++){
      mutation[i][j] = (unsigned long long ***)needMem(sizeof(unsigned long long **)*plen); //current quality score
      for(k=0;k<plen;k++){
        mutation[i][j][k] = (unsigned long long **)needMem(sizeof(unsigned long long *)*4); //4 reference bases
        for(k=0;l<4;l++){
          mutation[i][j][k][l] = (unsigned long long *)needMem(sizeof(unsigned long long)*5); //|{A,C,G,T,N}| = 5 possible read bases
        }
      }
    }
  }
  unsigned int junk = 0;
  //chr2_485_762_0:0:0_0:0:0_1619fe	99	chr1	10022	1	76M	=	10044	98	CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT	~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:412	XM:i:0	XO:i:0	XG:i:0	MD:Z:76



  /**
   * Start working on BAM alignment
   */
  //open stdin and read in sam format
  samfile_t *fp;
  if (sam == NULL){
    fp = samopen("-","rb",reflist);
  }else{
    fp = samopen(sam,"rb",reflist);
  }
  assert(fp);
  if(verboseOut)
    fprintf(stderr,"Processing bam alignment...\n");
  bam1_t *b;
  bam1_core_t *c;

  unsigned int *cigar;

  b = bam_init1();
  c = &b->core;
  DNA seq[10000];
  DNA refChar,readChar;
  uint8_t *tmpqual;
  uint8_t *tmpseq;
  uint8_t qual[1000];
  unsigned int ref_pos;
  unsigned int read_pos;
  unsigned int op;
  int ret;
  char *chrom;
  unsigned int prev_score_idx = 0;
  unsigned int score_idx = 0;

  while ((ret = samread(fp, b)) >= 0){

    chrom = fp->header->target_name[c->tid];
    ref = hashFindVal(refHash, chrom);
    if(!ref) errAbort("Sequence name %s not found in reference, but this read maps: %s\n",chrom,bam1_qname(b));




    //deal with the seq and qual string, reverse complement if needed
    tmpqual = bam1_qual(b);
    tmpseq = bam1_seq(b);
    for(i=0;i<c->l_qseq;i++){
      qual[i] = tmpqual[i];
      seq[i] = bam_nt16_rev_table[bam1_seqi(tmpseq, i)];
    }
    if(c->flag & BAM_FREVERSE){//rc read and rev qual
      reverseComplement(seq,c->l_qseq);
      intrev(qual,c->l_qseq);
    }
    cigar=bam1_cigar(b); //go through the cigar string until the alignment is fully explained
    k=0;
    op = cigar[0] & BAM_CIGAR_MASK;
    ref_pos = c->pos; //0based
    read_pos = 0; //start at the beginning of the read
    prev_score_idx = 0;
    score_idx = 0; //at first position it should be 0
    for(l = cigar[0] >> BAM_CIGAR_SHIFT; l>0; --l){
      // for each base with a 1-1 mapping:
      // 1. check if correct, incorrect, or unmapped
      // 2. increment stats
      if(l==0){
        k++;
        if(k>=c->n_cigar) break; //done
        op = cigar[k] & BAM_CIGAR_MASK;
        l = cigar[k] >> BAM_CIGAR_SHIFT;
      }
      refChar = toupper(ref->dna[ref_pos]);
      readChar = toupper(seq[read_pos]);
      score_idx = nearestIdxFromPscore(qual[read_pos]);

      /***
       * parse this position of the alignment
       */
      //let op==BAM_SOFTCLIP be a mapping position.
      if( op == BAM_CINS ){ //insertion to the reference
        //TODO: handle insertions
        read_pos++; //inc read not ref
      }else if(op == BAM_CDEL){ //deletion from the reference
        //TODO: handle deletions
        ref_pos++; //inc ref not read
      }else if(op == BAM_CREF_SKIP){ //skip stuff in the ref
        //TODO: handle... whatever this is (basically a long deletion?).
        ref_pos++; //inc ref not read
      }else{// 1-1 mapping (for the purposes of this, a clipped read is fair game)
        if(!invalid(refChar)){
          if(read_pos == 0){ //no prev phred score!
            mutation[FIRST_PREV_PHRED][score_idx][baseIndex(refChar)][baseIndex(readChar)]++;
          }else{ //standard situation
            mutation[prev_score_idx][score_idx][baseIndex(refChar)][baseIndex(readChar)]++;
          }
        }
        read_pos++;
        ref_pos++;
      }


      prev_score_idx = score_idx;
    } //end for loop over this read
  }//end while loop over reads
  bam_destroy1(b);

  samclose(fp);


  //print out the error histogram file.
  printf("#Pos\tPrev_Phred\tCurr_Phred\tA->A\tA->C\tA->G\tA->T\tA->N\tC->A\tC->C\tC->G\tC->T\tC->N\tG->A\tG->C\tG->G\tG->T\tG->N\tT->A\tT->C\tT->G\tT->T\tT->N\n");
  int pp,pc;
  for(i=0;i<readLen;i++){
    for(pp=0;pp<plen;pp++){
      if(i==0 && pp != FIRST_PREV_PHRED)
        continue; //nothing to print for this value of pp at this position
      for(pc=0;pc<plen;pc++){
        printf("%d\t%d\t%d",i,plist[pp],plist[pc]); //print pos\tprev_pred\tcurr_phred
        for(j=0;j<4;j++){
          for(k=0;k<5;k++)
            printf("\t%llu", mutation[i][pp][pc][j][k]);
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
  char *pscoreString;
  readLen = optionInt("readLen",0);//should crash if not set
  refName = optionVal("ref",NULL);
  phred33 = optionExists("phred33");
  phred64 = optionExists("phred64");
  sam = optionVal("sam",NULL);
  pscoreString = optionVal("pscores","0,5,10,15,20,25,30,35,40");
  if(phred33 && phred64)
  {
    errAbort("Options phred33 and phred64 are mutually exclusive!\n");
  }
  if(!(phred33 || phred64)) usage();
  if(!readLen || !refName) usage();
  char *pscores[10];
  chopString(pscoreString, ",", pscores, 10);
  int i;
  for(i=0;i<10;i++){
    if(strlen(pscores[i])==0){
      plen = i;
      break;
    }else{
      plist[i] = atoi(pscores[i]);
    }
  }

  getErrorProfile();
  return 0;
}

