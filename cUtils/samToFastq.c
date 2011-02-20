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
      "\t-verbose\tOutput verbose debug messages to stderr.\n\n"
  );
}//end usage()


static struct optionSpec options[] = {
    /* Structure holding command line options */
    {"read1",OPTION_STRING},
    {"read2",OPTION_STRING},
    {"single",OPTION_STRING},
    {"verbose",OPTION_BOOLEAN},
    {NULL, 0}
}; //end options()



void processBamFile(bamFile fp, FILE *single, FILE *read1, FILE *read2)
/*Iterate through all bam reads, a la samtools flagstat, and calculate our stats */
{
  if(verboseOut)
    fprintf(stderr,"Processing bam alignment...\n");
  bam1_t *b;
  bam1_core_t *c;
  int ret;
  bool print_singles = false;
  b = bam_init1();
  c = &b->core;
  while ((ret = bam_read1(fp, b)) >= 0)
  {
    //are we printing single reads?
    if()
  }//end while loop over reads
  bam_destroy1(b);
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *single = NULL;
  char *read1 = NULL;
  char *read2 = NULL;
  FILE *fs = NULL;
  FILE *f1 = NULL;
  FILE *f2 = NULL;
  optionInit(&argc, argv, options);
  verboseOut = optionExists("verbose");
  single = optionVal("single",NULL);
  read1 = optionVal("read1",NULL);
  read2 = optionVal("read2",NULL);
  if(single==NULL && (read1==NULL || read2==NULL)){
    //must supply single and/or (read1 and read2)
    errAbort("must supply single and/or (read1 and read2)\n");
  }else if((read1 == NULL && read2 != NULL) || (read2==NULL && read1 != NULL)){
    //must supply both read1 and read2 if either are supplied
    errAbort("must supply both read1 and read2 if either are supplied\n");
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
    processBamFile(fp,fs,f1,f2);
    bam_header_destroy(header);
    bam_close(fp);
  }

  return 0;
} //end main()


