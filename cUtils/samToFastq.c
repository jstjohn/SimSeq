/**
 * samToFastq
 *  Author: John St. John
 *  Date: 2/20/2010
 *  
 * 
 *
 *
 */




#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include "common.h"
#include "options.h"
#include "dnautil.h"
#include "fastq.h"
#include "sam.h"

bool verboseOut = false;



void usage()
/* Explain usage and exit. */
{
  errAbort(
      "samToFastq: Converts a **Name Sorted (read pairs must be grouped at least)** sam or bam file into fastq files. The fastq reads are reverse complemented as specified in the alignment.\n"
      "usage:\n"
      "\tsamToFastq [required options] [options] alignment1.bam(sorted by name, sam is ok too) [alignment2.bam, ..., alignmentN.bam] \n"
      "\n**required** options:\n"
      "\t-single=FILE\tEither writes only the unpaired reads if read1 and read2 is also supplied, or all reads in interleaved format (if bam is sorted!).\n"
      "\t===and/or===\n"
      "\t-read1=FILE\tFile to hold the first reads of a pair.\n"
      "\t ==and==\n"
      "\t-read2=FILE\tFile to hold the second reads of a pair. \n"
      "\noptions:\n"
      "\t-refList=FILE\t Required if using a sam file without a header, use the .fai file made by samtools index for this, or make a file of tab seperated 'ref\tlen' on new lines.\n"
      "\t-noUnmated\tonly print mated reads\n"
      "\t-phred64\t output ascii phred+64 strings rather than ascii phred+33 strings \n"
      "\t-h, -help\tdisplay this message and exit\n"
      "\t-verbose\tOutput verbose debug messages to stderr.\n\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"read1",OPTION_STRING},
    {"read2",OPTION_STRING},
    {"single",OPTION_STRING},
    {"refList",OPTION_STRING},
    {"phred64",OPTION_BOOLEAN},
    {"noUnmated",OPTION_BOOLEAN},
    {"verbose",OPTION_BOOLEAN},
    {"help",OPTION_BOOLEAN},
    {"h",OPTION_BOOLEAN},
    {NULL, 0}
}; //end options()


inline void fillFqItem(struct fastqItem *fq, bam1_t *b, bool phred64){
  bam1_core_t *c;
  c = &b->core;
  if(c->flag & BAM_FREAD1){
    sprintf(fq->id,"%s/1",bam1_qname(b));
  }else if(c->flag & BAM_FREAD2){
    sprintf(fq->id,"%s/2",bam1_qname(b));
  }else{
    strcpy(fq->id,bam1_qname(b));
  }
  fq->len = c->l_qseq;
  int i;
  for(i=0;i<c->l_qseq;i++)
    fq->seq[i] = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)];
  if(phred64)
    for(i=0;i<c->l_qseq;i++)
        fq->score[i] = phredToPhred64( (bam1_qual(b))[i] );
  else
    for(i=0;i<c->l_qseq;i++)
      fq->score[i] = phredToPhred33( (bam1_qual(b))[i] );
  fq->seq[c->l_qseq] = '\0';

  if(c->flag & BAM_FREVERSE){
    reverseComplementFastqItem(fq);
  }
}


void processBamFile(samfile_t *fp, FILE *single, FILE *read1, FILE *read2, bool noUnmated, bool phred64)
/*Iterate through all bam reads, a la samtools flagstat, and calculate our stats */
{
  if(verboseOut)
    fprintf(stderr,"Processing bam alignment...\n");
  bam1_t *b;
  bam1_t *b_prev;
  bam1_core_t *c;
  int ret;
  bool print_single_file = false;
  bool print_pair_file = false;
  struct fastqItem *fq = allocFastqItem();
  struct fastqItem *fq2 = allocFastqItem();
  if(read1 != NULL && read2 != NULL){
    print_pair_file = true;
  }
  if(single != NULL && print_pair_file && noUnmated){
    print_single_file = false; //don't need to print singles if print paired and no unmated
  }else if(single != NULL){
    print_single_file = true;
  }
  b = bam_init1();
  b_prev = bam_init1();
  c = &b->core;
  bool second_wait = false; //waiting for a second matching read to show up
  while ((ret = samread(fp, b)) >= 0)
  {
    //are we printing single reads, and pairs?
    if(c->flag & BAM_FSECONDARY) continue; //skip non-primary alignments;
    if(!second_wait){
      bam_copy1(b_prev,b);
      second_wait = true;
    }else{
      //check to see if read ids match
      if(strcmp(bam1_qname(b),bam1_qname(b_prev))==0){
        //print out b and b_prev
        fillFqItem(fq,b,phred64);
        fillFqItem(fq2,b_prev,phred64);
        if(print_pair_file){
          if(c->flag & BAM_FREAD1){
            printFastqItem(read1,fq);
            printFastqItem(read2,fq2);
          }else{ //reads go in other files
            printFastqItem(read1,fq2);
            printFastqItem(read2,fq);
          }
        }else{//print both to single file
          if(c->flag & BAM_FREAD1){
            printFastqItem(single,fq);
            printFastqItem(single,fq2);
          }else{ //reads go in other files
            printFastqItem(single,fq2);
            printFastqItem(single,fq);
          }
        }
        //set second_wait to false
        second_wait = false;
      }else{ //reads don't match, print the single
        if(print_single_file && !noUnmated){
          fillFqItem(fq,b_prev,phred64);
          printFastqItem(single,fq);
        }
        //copy b into b_prev, and second_wait = true
        second_wait = true;
        bam_copy1(b_prev,b);
      }

    }//end second_wait == true
  }//end while loop over reads
  if(second_wait){//hit last read, and it couldn't be a match to the prev
    if(print_single_file && !noUnmated){
      fillFqItem(fq,b_prev,phred64);
      printFastqItem(single,fq);
    }
  }
  bam_destroy1(b);
  freeFastqItem(fq);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *single = NULL;
  char *read1 = NULL;
  char *read2 = NULL;
  char *reflist = NULL;
  bool noUnmated = false;
  bool help = false;
  bool phred64 = false;
  FILE *fs = NULL;
  FILE *f1 = NULL;
  FILE *f2 = NULL;
  optionInit(&argc, argv, options);
  verboseOut = optionExists("verbose");
  noUnmated = optionExists("noUnmated");
  help = (optionExists("h") || optionExists("help"));
  if(help) usage();
  single = optionVal("single",NULL);
  read1 = optionVal("read1",NULL);
  read2 = optionVal("read2",NULL);
  reflist = optionVal("refList",NULL);
  phred64 = optionExists("phred64");
  if(single==NULL && (read1==NULL || read2==NULL)){
    //must supply single and/or (read1 and read2)
    fprintf(stderr,"must supply single and/or (read1 and read2)\n");
    usage();
  }else if((read1 == NULL && read2 != NULL) || (read2==NULL && read1 != NULL)){
    //must supply both read1 and read2 if either are supplied
    fprintf(stderr,"must supply both read1 and read2 if either are supplied\n");
    usage();
  }
  if(single != NULL){
    fs = fopen(single,"w");
  }
  if(read1 != NULL){
    f1 = fopen(read1,"w");
    f2 = fopen(read2,"w");
  }
  int i;
  for(i=1;i<argc;i++){
    samfile_t *fp;
    if(strcasecmp(argv[i]+(strlen(argv[i])-4),"sam")){
      fp = samopen(argv[i],"r",reflist);
    }else{
      fp = samopen(argv[i],"rb",reflist);
    }
    assert(fp);
    //process the bam file, write output
    processBamFile(fp,fs,f1,f2,noUnmated,phred64);
    samclose(fp);
  }
  if(single != NULL){
      fclose(fs);
    }
    if(read1 != NULL){
      fclose(f1);
      fclose(f2);
   }

  return 0;
} //end main()


