
package simseq;

/**
 * Represents a fastq record, fairly literally, i.e. without any conversion.
 */
public class SamRecord {
    //Sam Flags
    public boolean paired;
    public boolean proper_pair;
    public boolean query_unmapped;
    public boolean mate_unmapped;
    public boolean query_reverse_strand;
    public boolean mate_reverse_strand;
    public boolean first;
    public boolean second;
    public boolean duplicate;
    public int pos;//position
    public int mpos;//mate position
    public int isize;
    public long seqIndex;
    public String cigar;
    public String seqHeaderPrefix;
    public String refname;
    public StringBuilderDNA seqLine;
    public StringBuilder qualLine;
    public boolean chimeric;
    
    public SamRecord(String cigar_string, String header_prefix, String dummyseq, 
            String dummyqual, boolean first_read){
        this.seqIndex=0;
        this.cigar=cigar_string;
        this.seqHeaderPrefix=header_prefix;
        this.first = first_read;
        this.second = ! this.first;
        this.paired = true;
        this.seqLine = new StringBuilderDNA(dummyseq);
        this.qualLine = new StringBuilder(dummyqual);
    }
    public void setRef(String ref_name){
        refname = ref_name;
    }
    public void initMP(){
    }
    public void initPE(){
        this.proper_pair = true; //no chance of chimeric reads
        this.query_unmapped = false;
        this.mate_unmapped = false;
    }
    public int getFlag(){
        int                     flag  = 0;
        if(paired)              flag |= 0x001;
        if(proper_pair)         flag |= 0x002;
        if(query_unmapped)      flag |= 0x004;
        if(mate_unmapped)       flag |= 0x008;
        if(query_reverse_strand)flag |= 0x0010;
        if(mate_reverse_strand) flag |= 0x0020;
        if(first)               flag |= 0x0040;
        if(second)              flag |= 0x0080;
        if(duplicate)           flag |= 0x0400;
        return flag;
    }
}
