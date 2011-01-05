/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package simseq;

import java.io.File;
import java.util.Random;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;


/**
 *
 * @author john
 */
public class Main{
    private static String lastUpdate = "1.5.2011";
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        try {
            // TODO code application logic here
            String def_read_len = "100";
            String def_mate_sheer_mean = "500";
            String def_mate_sheer_stdev = "50";
            String def_library_ins_mean = "200";
            String def_library_ins_stdev = "20";
            String def_duplicate_p = "0.0";
            String def_mate_non_biotin_pulldown_p = "0.3";
            String def_read_prefix = "SimSeq_";
            String def_read_number = "1000000";

            Options options = new Options();
            //REQUIRED ARGS
            options.addOption("o", "out", true, "Filename for output sam file (REQUIRED)");
            options.addOption("r", "reference",true, "Reference genome sequence in uncompressed fasta format (REQUIRED)");


            //OPTIONAL ARGS
            options.addOption("d", "dip", true, "If diploid data desired, path to diploid file. (format: chrom pos(0 based) altChar");
            options.addOption("e", "error", true, "If simulated read error desired, path to read error file.");
            options.addOption(null, "error2", true, "If you desire a seperate error distribution to be applied to the second read, then provide a path to that error profile with this option");
            options.addOption("p", "read_prefix", true, "Prefix for simulated reads. Default: "+def_read_prefix);
            options.addOption("1", "read1_length", true, "Integer length of first read. Default: "+def_read_len);
            options.addOption("2", "read2_length", true, "Integer length of second read. Default: "+def_read_len);
            options.addOption("n", "read_number",true, "Integer number of reads you would like to sample. Default: "+def_read_number);
            options.addOption("u", "duplicate_probability", true, "probability of generating a duplicate. Default: "+def_duplicate_p);
            options.addOption("l", "insert_size",true,"mean library insert size for either mate-paired or paired-end. Default: "+def_library_ins_mean);
            options.addOption("s", "insert_stdev",true,"mean library insert stdev for either mate-paired or paired-end. Default: "+def_library_ins_stdev);
            options.addOption("m", "mate_pair",false,"Perform mate-pair rather than paired end run.");
            options.addOption(null,"mate_frag",true, "If using a mate-pair library, what is your desired loop fragmentation size? Default: "+def_mate_sheer_mean);
            options.addOption(null,"mate_frag_stdev",true,"If using a mate-pair library, what is your desired loop fragmentation size standard deviation? Default: "+def_mate_sheer_stdev);
            options.addOption(null,"mate_pulldown_error_p",true, "If using a mate-pair library, what is the probability that a read does not include the biotin marker? Default: "+def_mate_non_biotin_pulldown_p);
            options.addOption(null,"phred64",false,"Output phred+64 quality string rather than phred+33");
            options.addOption("h", "help",false,"Print this usage message.");
            options.addOption(null,"debug",false,"Write debug info to stderr.");
            CommandLineParser parser = new PosixParser();
            HelpFormatter formatter = new HelpFormatter();
            CommandLine cmd = parser.parse(options, args);
            
            if(cmd.hasOption('h')){
                formatter.printHelp( "java -jar -Xmx2048m SimSeq.jar [required options] [options]\nLast Updated: "+lastUpdate, options );
                System.exit(1);
            }
            if(!cmd.hasOption('o') || !cmd.hasOption('r')){
                throw new ParseException("Missing required argument, try -h or --help for help.");
            }
            File out = new File(cmd.getOptionValue('o'));
            File ref  = new File(cmd.getOptionValue('r'));
            int read1_length = Integer.parseInt(cmd.getOptionValue('1', def_read_len));
            int read2_length = Integer.parseInt(cmd.getOptionValue('2', def_read_len));
            int mate_sheer_mean = Integer.parseInt(cmd.getOptionValue("mate_frag", def_mate_sheer_mean));
            int mate_sheer_stdev = Integer.parseInt(cmd.getOptionValue("mate_frag_stdev", def_mate_sheer_stdev));
            int library_ins_mean = Integer.parseInt(cmd.getOptionValue('l', def_library_ins_mean));
            int library_ins_stdev = Integer.parseInt(cmd.getOptionValue('s', def_library_ins_stdev));
            int read_number = Integer.parseInt(cmd.getOptionValue('n',def_read_number));
            double duplicate_p = Double.valueOf(cmd.getOptionValue('u', def_duplicate_p));
            double mate_non_biotin_pulldown_p = Double.valueOf(cmd.getOptionValue("mate_pulldown_error_p",def_mate_non_biotin_pulldown_p));
            String read_prefix = cmd.getOptionValue('p', def_read_prefix);
            boolean debug = cmd.hasOption("debug");
            AddDiploid dadd = null;
            AddError eadd = null;
            AddError eadd2 = null;
            boolean phred33 = true;
            if(cmd.hasOption("phred64"))phred33=false;
            if(cmd.hasOption('d')){
                System.err.println("WARNING: diploid sites addition is not currently implemented");
                dadd = new AddDiploid(cmd.getOptionValue('d'),debug);
            }
            if(cmd.hasOption("error") && !cmd.hasOption("error2")){
                eadd = new AddError(cmd.getOptionValue('e'),Math.max(read1_length, read2_length),phred33,debug);
            }
            else if(cmd.hasOption("error2"))
            {
                if(!cmd.hasOption("error"))
                    throw new ParseException("You must supply an error file with '-e' or '--error' if you also want to supply '--error2'");
                eadd = new AddError(cmd.getOptionValue("error"),read1_length,phred33,debug);
                eadd2 = new AddError(cmd.getOptionValue("error2"),read2_length,phred33,debug);
            }
            SamWriter swrite = new SamWriter(out);
            //get seq length
            FastaReader far = new FastaReader(ref);
            int totalLen = 0;
            for(FastaRecord rec :far){
                totalLen += rec.getSeq().length();
            }
            far = new FastaReader(ref);
            StringBuilder qual1sb = new StringBuilder();
            StringBuilder qual2sb = new StringBuilder();
            StringBuilder dummy1sb = new StringBuilder();
            StringBuilder dummy2sb = new StringBuilder();

            char perfQual = QualityUtil.getPhred33ScoreFromPhredScore(93);
            for(int i=0;i<read1_length;i++){
                qual1sb.append(perfQual);
                dummy1sb.append("A");
            }
            for(int i=0;i<read2_length;i++){
                qual2sb.append(perfQual);
                dummy2sb.append("A");
            }
            String qual1 = qual1sb.toString();
            String qual2 = qual2sb.toString();
            String seq1tmp = dummy1sb.toString();
            String seq2tmp = dummy2sb.toString();

            SamRecord sr1 = new SamRecord(Integer.toString(read1_length)+"M", read_prefix, seq1tmp, qual1, true);
            SamRecord sr2 = new SamRecord(Integer.toString(read2_length)+"M", read_prefix, seq2tmp, qual2, false);

            for(FastaRecord rec :far){
                sr1.setRef(rec.getID());
                sr2.setRef(rec.getID());
                SeqSampler sampler = new SeqSampler(rec.getSeq());
                //find out how many reads to sample from this sequence
                int len = rec.getSeq().length();
                long num = Math.round(((double)len/(double)totalLen)*(double)read_number);
                Random r = new Random();
                String tmpseq1 = null,
                                tmpseq2 = null,
                                tmpqual1 = null,
                                tmpqual2 = null;
                boolean duplicate = false;
                if(cmd.hasOption('m')){//mate pair run?
                    for(long i = 0; i < num; i++){
                        /*
                         *     public void MPSample(SamRecord sr1,
                         *         SamRecord sr2, String qual1, String qual2,
                         *         int mate_ins, int mate_stdev, int read_ins, int read_stdev,
                         *         int read1_len, int read2_len, double p_bad_pulldown, Random r){
                         */
                        sampler.MPSample(sr1, sr2,
                                qual1, qual2,
                                library_ins_mean,
                                library_ins_stdev,
                                mate_sheer_mean,
                                mate_sheer_stdev,
                                read1_length,
                                read2_length,
                                mate_non_biotin_pulldown_p, r);
                        if(dadd != null){
                            dadd.AddDiploidRead(sr1);
                            dadd.AddDiploidRead(sr2);
                        }
                        tmpseq1 = tmpseq2 = tmpqual1 = tmpqual2 = null;
                        duplicate = false;
                        if(r.nextFloat()< duplicate_p){
                            duplicate = true;
                            tmpseq1 = sr1.seqLine.toString();
                            tmpqual1 = sr1.qualLine.toString();
                            tmpseq2 = sr2.seqLine.toString();
                            tmpqual2 = sr2.qualLine.toString();
                        }
                        if(eadd != null && eadd2 == null){
                                eadd.AddErrorRead(sr1);
                                eadd.AddErrorRead(sr2);
                        }else if(eadd2 != null){
                                eadd.AddErrorRead(sr1);
                                eadd2.AddErrorRead(sr2);
                        }

                        swrite.write(sr1);
                        swrite.write(sr2);

                        while(duplicate){
                            //mark duplicate
                            i++;
                            if(i>=num){
                                duplicate = false;
                                sr1.duplicate = false;
                                sr2.duplicate = false;
                                break;
                            }//double check we aren't duplicating over the number of reads
                            sr1.seqIndex++; //increment the counters
                            sr2.seqIndex++; //increment the counters
                            sr1.duplicate = duplicate;
                            sr2.duplicate = duplicate;
                            //restore original reads and qual prior to error
                            sr1.seqLine = new StringBuilderDNA(tmpseq1);
                            sr1.qualLine = new StringBuilder(tmpqual1);
                            sr2.seqLine = new StringBuilderDNA(tmpseq2);
                            sr2.qualLine = new StringBuilder(tmpqual2);
                            if(eadd != null && eadd2 == null){
                                eadd.AddErrorRead(sr1);
                                eadd.AddErrorRead(sr2);
                            }else if(eadd2 != null){
                                eadd.AddErrorRead(sr1);
                                eadd2.AddErrorRead(sr2);
                            }
                            swrite.write(sr1);
                            swrite.write(sr2);
                            if(r.nextFloat()< duplicate_p){
                                duplicate = true;
                            }else{//exit loop
                                duplicate = false;
                                sr1.duplicate = false;
                                sr2.duplicate = false;
                                break;
                            }
                        }



                    }

                }else{ //regular run?
                    for(long i = 0; i < num; i++){
                        /*
                         *     public void PESample(SamRecord sr1,
                         *          SamRecord sr2, String qual1, String qual2,
                         *          int mean_ins, int stdev,
                         *          int read1_len, int read2_len, Random r){
                         */
                        sampler.PESample(sr1, sr2,
                                qual1, qual2,
                                library_ins_mean,
                                library_ins_stdev,
                                read1_length,
                                read2_length, r);
                        if(dadd != null){
                            dadd.AddDiploidRead(sr1);
                            dadd.AddDiploidRead(sr2);
                        }
                        tmpseq1 = tmpseq2 = tmpqual1 = tmpqual2 = null;
                        duplicate = false;
                        if(r.nextFloat()< duplicate_p){
                            duplicate = true;
                            tmpseq1 = sr1.seqLine.toString();
                            tmpqual1 = sr1.qualLine.toString();
                            tmpseq2 = sr2.seqLine.toString();
                            tmpqual2 = sr2.qualLine.toString();
                        }
                        if(eadd != null && eadd2 == null){
                                eadd.AddErrorRead(sr1);
                                eadd.AddErrorRead(sr2);
                        }else if(eadd2 != null){
                                eadd.AddErrorRead(sr1);
                                eadd2.AddErrorRead(sr2);
                        }

                        swrite.write(sr1);
                        swrite.write(sr2);

                        while(duplicate){
                            //mark duplicate
                            i++;
                            if(i>=num){
                                duplicate = false;
                                sr1.duplicate = false;
                                sr2.duplicate = false;
                                break;
                            }//double check we aren't duplicating over the number of reads
                            sr1.seqIndex++; //increment the counters
                            sr2.seqIndex++; //increment the counters
                            
                            sr1.duplicate = duplicate;
                            sr2.duplicate = duplicate;
                            //restore original reads and qual prior to error
                            sr1.seqLine = new StringBuilderDNA(tmpseq1);
                            sr1.qualLine = new StringBuilder(tmpqual1);
                            sr2.seqLine = new StringBuilderDNA(tmpseq2);
                            sr2.qualLine = new StringBuilder(tmpqual2);
                            if(eadd != null && eadd2 == null){
                                eadd.AddErrorRead(sr1);
                                eadd.AddErrorRead(sr2);
                            }else if(eadd2 != null){
                                eadd.AddErrorRead(sr1);
                                eadd2.AddErrorRead(sr2);
                            }
                            swrite.write(sr1);
                            swrite.write(sr2);
                            if(r.nextFloat()< duplicate_p){
                                duplicate = true;
                            }else{//exit loop
                                duplicate = false;
                                sr1.duplicate = false;
                                sr2.duplicate = false;
                                break;
                            }
                        }

                    }
                }
            }

            swrite.close();
            





            
            //

        } catch (Exception ex) {
            throw new RuntimeException(ex);
            //System.err.println("Program failure. Reason: "+ex.getStackTrace());
        }
    }

}



