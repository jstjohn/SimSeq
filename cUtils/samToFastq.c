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
bam_header_t *header;



void usage()
/* Explain usage and exit. */
{
  errAbort(
      "samToFastq: Converts a sam or bam file into fastq files. The fastq reads are reverse complemented as specified in the alignment.\n"
      "usage:\n"
      "\tsamToFastq [required options] [options] alignment1.bam(sorted by name, sam is ok too) [alignment2.bam, ..., alignmentN.bam] \n"
      "\n**required** options:\n"
      "\t-single=FILE\tEither writes only the unpaired reads if read1 and read2 is also supplied, or all reads in interleaved format (if bam is sorted!).\n"
      "\t===and/or===\n"
      "\t-read1=FILE\tFile to hold the first reads of a pair.\n"
      "\t ==and==\n"
      "\t-read2=FILE\tFile to hold the second reads of a pair. \n"
      "\noptions:\n"
      "\t-noUnmated\tonly print mated reads\n"
      "\t-h, -help\tdisplay this message and exit\n"
      "\t-verbose\tOutput verbose debug messages to stderr.\n\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"read1",OPTION_STRING},
    {"read2",OPTION_STRING},
    {"single",OPTION_STRING},
    {"noUnmated",OPTION_BOOLEAN},
    {"verbose",OPTION_BOOLEAN},
    {"help",OPTION_BOOLEAN},
    {"h",OPTION_BOOLEAN},
    {NULL, 0}
}; //end options()


inline void fillFqItem(struct fastqItem *fq, bam1_t *b){
  bam1_core_t *c;
  c = &b->core;
  if(c->flag & BAM_FREAD1){
    sprintf(fq->id,"%s/1",bam1_qname(b));
  }else if(c->flag & BAM_FREAD2){
    sprintf(fq->id,"%s/2",bam1_qname(b));
  }else{
    strcpy(fq->id,bam1_qname(b));
  }

  strcpy(fq->score,bam1_qual(b));
  fq->len = c->l_qseq;
  int i;
  for(i=0;i<c->l_qseq;i++){
    fq->seq[i] = bam1_seqi(bam1_seq(b), i);
  }
  fq->seq[c->l_qseq] = '\0';

  if(c->flag & BAM_FREVERSE){
    reverseComplementFastqItem(fq);
  }
}


void processBamFile(bamFile fp, FILE *single, FILE *read1, FILE *read2, bool noUnmated)
/*Iterate through all bam reads, a la samtools flagstat, and calculate our stats */
{
  if(verboseOut)
    fprintf(stderr,"Processing bam alignment...\n");
  bam1_t *b;
  bam1_core_t *c;
  int ret;
  bool print_single_file = false;
  bool print_pair_file = false;
  struct fastqItem *fq = allocFastqItem();
  if(read1 != NULL && read2 != NULL){
    print_pair_file = true;
  }
  if(single != NULL && print_pair_file && noUnmated){
    print_single_file = false; //don't need to print singles if print paired and no unmated
  }else if(single != NULL){
    print_single_file = true;
  }
  b = bam_init1();
  c = &b->core;
  while ((ret = bam_read1(fp, b)) >= 0)
  {
    //are we printing single reads, and pairs?
    if(c->flag & BAM_FSECONDARY) continue; //skip non-primary alignments;
    fillFqItem(fq,b);
    if(print_single_file && print_pair_file){
      if(c->flag & BAM_FPAIRED){
        if(c->flag & BAM_FREAD1)
          printFastqItem(read1,fq);
        else
          printFastqItem(read2,fq);
      }else{//single read
        printFastqItem(single,fq);
      }
    }else if(print_pair_file){
      //just print paired reads to the pair files
      if(c->flag & BAM_FPAIRED){
        if(c->flag & BAM_FREAD1)
          printFastqItem(read1,fq);
        else
          printFastqItem(read2,fq);
      }
    }else{
      //print everything to singles
      if(noUnmated){
        if(c->flag & BAM_FPAIRED)
          printFastqItem(single,fq);
      }else{
        printFastqItem(single,fq);
      }
    }
  }//end while loop over reads
  bam_destroy1(b);
  freeFastqItem(fq);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *single = NULL;
  char *read1 = NULL;
  char *read2 = NULL;
  bool noUnmated = false;
  bool help = false;
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
    bamFile fp = bam_open(argv[i],"r");
    assert(fp);
    header = bam_header_read(fp);
    bam_init_header_hash(header);
    //process the bam file, write output
    processBamFile(fp,fs,f1,f2,noUnmated);
    bam_header_destroy(header);
    bam_close(fp);
  }

  return 0;
} //end main()


