/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package simseq;
import java.util.Random;

/**
 *
 * @author jstjohn
 */
public class SeqSampler{
    private String seq;
    private int len;
    public SeqSampler(String sequence){
        this.seq = sequence;
        this.len = sequence.length();
    }
    public void PESample(SamRecord sr1,
            SamRecord sr2, String qual1, String qual2,
            int mean_ins, int stdev,
            int read1_len, int read2_len, Random r){
        int p,l;
        if(mean_ins > this.len){
            throw new RuntimeException("Fasta file has shorter sequence than the requested mean library length");
        }
        //increment the seq index
        sr1.seqIndex++;
        sr2.seqIndex++;

        do{
            p = r.nextInt(this.len);
            l = ((int)(r.nextGaussian()*((double)stdev))) + mean_ins;
        }while((p+l)>this.len || l < read1_len || l < read2_len);
        sr1.qualLine.replace(0, qual1.length(), qual1);
        sr2.qualLine.replace(0, qual2.length(), qual2);

        if(r.nextBoolean()){ //sample from forward strand
            sr1.seqLine.replace(seq.substring(p, p+read1_len));
            sr2.seqLine.replace(seq.substring(p+l-read2_len,p+l));
            //strand settings
            sr1.mate_reverse_strand = true;
            sr1.query_reverse_strand = false;
            sr1.pos = p+1;
            sr1.mpos = p+l-read2_len+1;
            sr1.isize = l;
            sr2.isize = -l;
        }else{//sample from reverse strand
            sr1.seqLine.replace(seq.substring(p+l-read1_len,p+l));
            sr2.seqLine.replace(seq.substring(p,p+read2_len));
            //strand settings
            sr1.mate_reverse_strand = false;
            sr1.query_reverse_strand = true;
            sr1.pos = p+l-read1_len+1;
            sr1.mpos = p+1;
            sr1.isize = -l;
            sr2.isize = l;
        }
        //things that can be inferred in 2 from 1
        sr2.mate_reverse_strand  = sr1.query_reverse_strand;
        sr2.query_reverse_strand = sr1.mate_reverse_strand;
        sr2.pos = sr1.mpos;
        sr2.mpos = sr1.pos;
    }


    
    public void MPSample(SamRecord sr1,
            SamRecord sr2, String qual1, String qual2,
            int mate_ins, int mate_stdev, int read_ins, int read_stdev,
            int read1_len, int read2_len, double p_bad_pulldown, Random r){
        int i,b,l,p;
        if(mate_ins > this.len){
            throw new RuntimeException("Fasta file has shorter sequence than the requested mean library length");
        }
        //increment the seq index
        sr1.seqIndex++;
        sr2.seqIndex++;

        //grab our mate insert, our sheared mate loop, and position
        do{
            p = r.nextInt(this.len);
            l = ((int)(r.nextGaussian()*((double)mate_stdev))) + mate_ins;
            i = ((int)(r.nextGaussian()*((double)read_stdev))) + read_ins;
        }while( i < read1_len  ||
                i < read2_len  ||
                (p+l)>this.len ||
                l < read1_len  ||
                l < read2_len  ||
                i > p);

        
        b = r.nextInt(i); //the location of biotin
        //initialize sr1 and sr2
        sr1.qualLine.replace(0, qual1.length(), qual1);
        sr2.qualLine.replace(0, qual2.length(), qual2);
        sr1.chimeric = sr2.chimeric = false;
        sr1.mate_unmapped = sr2.mate_unmapped = false;
        sr1.query_unmapped = sr2.query_unmapped = false;
        sr1.proper_pair = sr2.proper_pair = true;
        
        boolean rev = r.nextBoolean();
        if(rev){ //set up flags for mp library
            sr2.query_reverse_strand =
                sr1.query_reverse_strand =
                sr1.mate_reverse_strand =
                sr2.query_reverse_strand =
                sr2.mate_reverse_strand = true;
        }else{
            sr2.query_reverse_strand = 
                sr1.query_reverse_strand =
                sr1.mate_reverse_strand =
                sr2.query_reverse_strand = 
                sr2.mate_reverse_strand = false;
        }
        if(b==0 || b == i-1 || p_bad_pulldown > r.nextDouble()){
            if(rev){//reverse strand sample
                sr2.seqLine.replace(seq.substring(p, p+read2_len));
                sr1.seqLine.replace(seq.substring(p+i-read1_len,p+i));
                sr1.pos=p+i-read1_len+1;
                sr1.mpos=p+read2_len;
                sr1.isize = -i;
            }else{
                sr1.seqLine.replace(seq.substring(p, p+read1_len));
                sr2.seqLine.replace(seq.substring(p+i-read2_len,p+i));
                sr1.mpos=p+i-read1_len+1;
                sr1.pos=p+read2_len;
                sr1.isize = i;
            }

            //infer sr2 things from sr1
            sr2.isize = -1*sr1.isize;
            sr2.pos=sr1.mpos;
            sr2.mpos=sr1.pos;
        }else{
            int rstart;
            if(rev){ //sample from reverse
                rstart = Math.abs(Math.min(0,i-b-read2_len));
                sr2.seqLine.replace(rstart,read2_len,seq.substring(p,p+read2_len-rstart));
                if(rstart != 0){
                    sr2.seqLine.replace(0,rstart,seq.substring(p+l-rstart,p+l));
                    sr2.chimeric = true;
                    sr2.query_unmapped = true;
                    sr1.mate_unmapped = true;
                    sr2.proper_pair=sr1.proper_pair=false;
                    sr1.isize=sr2.isize=0;
                }
                rstart = Math.abs(Math.min(0,b-read1_len));
                sr1.seqLine.replace(0,rstart,seq.substring(p+l-b,p+l-b+rstart));
                if(rstart != 0){
                    sr1.seqLine.replace(rstart, read1_len, seq.substring(p,p+read1_len-rstart));
                    sr1.chimeric = true;
                    sr1.query_unmapped = true;
                    sr2.mate_unmapped = true;
                    sr2.proper_pair=sr1.proper_pair=false;
                    sr1.isize=sr2.isize=0;
                }
                if(!sr1.chimeric && !sr2.chimeric){
                    sr2.isize=l-read1_len;
                    sr1.isize = - sr2.isize;
                }
                
            }else{//sample from forward
                rstart = Math.abs(Math.min(0,i-b-read1_len));
                sr1.seqLine.replace(rstart,read1_len,seq.substring(p,p+read1_len-rstart));
                if(rstart != 0){
                    sr1.seqLine.replace(0,rstart,seq.substring(p+l-rstart,p+l));
                    sr1.chimeric = true;
                    sr1.query_unmapped = true;
                    sr2.mate_unmapped = true;
                    sr1.proper_pair=sr2.proper_pair=false;
                    sr1.isize=sr2.isize=0;
                }
                rstart = Math.abs(Math.min(0,b-read2_len));
                sr2.seqLine.replace(0,rstart,seq.substring(p+l-b,p+l-b+rstart));
                if(rstart != 0){
                    sr2.seqLine.replace(rstart, read1_len, seq.substring(p,p+read2_len-rstart));
                    sr2.chimeric = true;
                    sr2.query_unmapped = true;
                    sr1.mate_unmapped = true;
                    sr1.proper_pair=sr2.proper_pair=false;
                    sr2.isize=sr1.isize=0;
                }
                if(!sr1.chimeric && !sr2.chimeric){
                    sr1.isize=l-read1_len;
                    sr2.isize = - sr1.isize;
                }
            }
        }
    }
}
