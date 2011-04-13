/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package simseq;


import java.io.BufferedWriter;
import java.io.IOException;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;

public class SamWriter {
    private final File file;
    private final PrintWriter writer;
    private final boolean informative_id;
    public SamWriter(final File file, boolean inf_id) throws IOException {
        this.file = file;
        this.informative_id = inf_id;
        this.writer = new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
       // this.writer = new PrintWriter(IoUtil.openFileForBufferedWriting(file));
    }
    public void write(final SamRecord rec) throws IOException {
        
        //<QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> \
        //[<TAG>:<VTYPE>:<VALUE> [...]]
        
        //print the informative id if desired
        if(this.informative_id){
            writer.print(rec.refname);
            writer.print("_");
            long start = Math.min(rec.pos, rec.mpos);
            writer.print(start);
            writer.print("_");
            writer.print(start+Math.abs(rec.isize));
            if(rec.first && rec.query_reverse_strand || rec.second && rec.mate_reverse_strand){
                writer.print("_R_");
            }else{// from forward strand
                writer.print("_F_");
            }
        }else{
            writer.print(rec.seqHeaderPrefix);//print prefix
        }
        writer.print(Long.toHexString(rec.seqIndex)); //qname
        writer.print("\t");
        writer.print(Integer.toString(rec.getFlag()));//flag
        writer.print("\t");
        writer.print(rec.refname);
        writer.print("\t");
        writer.print(Integer.toString(rec.pos));//pos
        writer.print("\t");
        writer.print("255"); //mapq=255 means no map quality information
        writer.print("\t");
        writer.print(rec.cigar); //cigar
        writer.print("\t");
        writer.print("="); //name of mate's reference(MRNM) = for same
        writer.print("\t");
        writer.print(Integer.toString(rec.mpos));//mpos
        writer.print("\t");
        writer.print(Integer.toString(rec.isize));//isize
        writer.print("\t");
        writer.print(rec.seqLine.toString());//seq
        writer.print("\t");
        writer.print(rec.qualLine.toString());
        writer.print("\t");
        if(rec.chimeric) writer.print("YC:Z:"+rec.c_cigar);
        else writer.print("YC:Z:"+rec.cigar);
        writer.print("\n");


        if (writer.checkError()) {
            throw new IOException("Error in writing file "+file);
        }
    }

    public void close() {
        writer.close();
    }
}
