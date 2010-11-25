/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package simseq;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 *
 * @author jstjohn
 */
public class FastaReader implements Iterator<FastaRecord>, Iterable<FastaRecord>, Closeable {
    final private File fastaFile;
    final private BufferedReader reader;
    private FastaRecord nextRecord;
    private String nextID;
    private int line=1;

    public FastaReader(final File file) throws FileNotFoundException, IOException {
        fastaFile = file;
        reader = new BufferedReader(new InputStreamReader(new FileInputStream(fastaFile)));
        String tline;
        while(null != (tline = reader.readLine()) && !tline.startsWith(">")){
            line++;
        }
        if(tline != null){
            String tmp = tline.substring(1);
            nextID=tmp.split("\\s")[0];
        }
        else nextID = null;
        nextRecord = readNextRecord();
    }

    private FastaRecord readNextRecord(){
        try {

            // read until the first id
            String tline;
            //fast forward to first id
            if(nextID==null){
               return null;
            }

            final String seqHeader = nextID;
            if (seqHeader == null) return null;
            if (StringUtil.isBlank(seqHeader)) {
                throw new IOException(error("Missing sequence id"));
            }

            // Read sequence line
            StringBuilder seqBuilder = new StringBuilder();
            while((tline = reader.readLine()) != null && !tline.startsWith(">")){
                line++;
                if(!StringUtil.isBlankOrComment(tline)){
                    seqBuilder.append(tline);
                }
            }
            if(tline==null){
                nextID=null;
            }else if(!tline.startsWith(">")){ //it should start with > or be null
                throw new RuntimeException("Unknown error at line: "+Integer.toString(line));
            }else{//good id line
            String tmp = tline.substring(1);
            nextID=tmp.split("\\s")[0];
        }
            final String seqLine = seqBuilder.toString();
            checkLine(seqLine,"dna sequence");

            final FastaRecord frec = new FastaRecord(seqHeader, seqLine);
            return frec;

        } catch (IOException e) {
            throw new RuntimeException(String.format("Error reading '%s'", fastaFile.getAbsolutePath()),e);
        }
    }

    public boolean hasNext() { return nextRecord != null; }

    public FastaRecord next(){
        if (!hasNext()) {
            throw new NoSuchElementException("next() called when !hasNext()");
        }
        final FastaRecord rec = nextRecord;
        nextRecord = readNextRecord();
        return rec;
    }

    public void remove() { throw new UnsupportedOperationException("Unsupported operation"); }

    public Iterator<FastaRecord> iterator() { return this; }

    public int getLineNumber() { return line ; }
    public File getFile() { return fastaFile ; }

    public void close() {
        try {
            reader.close();
        } catch (IOException e) {
            throw new RuntimeException("IO problem in file "+fastaFile.getAbsolutePath(),e);
        }
    }

    private void checkLine(final String line, final String kind) {
        if (line == null) {
            throw new RuntimeException(error("File is too short - missing "+kind+" line"));
        }
        if (StringUtil.isBlank(line)) {
            throw new RuntimeException(error("Missing "+kind));
        }
    }

    private String error(final String msg) {
        return msg + " at line "+line+" in "+fastaFile.getAbsolutePath();
    }
}
