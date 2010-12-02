/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package simseq;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Random;

/**
 *
 * @author john
 */
public class AddError {

    private int[][][][] mp;
    private Random rand = new Random();
    private boolean phred33;

    public AddError(String errorFname, int rlen, boolean phred_33, boolean debug) throws IOException {
        /*
         * Constructor, simply read in the error file and get ready to
         * add error to any reads that come by.
         */
        phred33 = phred_33; //for now this is not an option
        parseErrorFile(errorFname, rlen);
    }


    private int[] parseLine(String line){
        //turns a tab spaced line into an int array
        String parts[] = line.split("\\s");
        int i;
        int[] res = new int[parts.length];
        for(i=0;i<parts.length;i++){
            res[i] = Integer.parseInt(parts[i]);
        }
        return res;
    }
    private void parseErrorFile(String fname, int rlen) throws FileNotFoundException, IOException {
        BufferedReader fh = new BufferedReader(new FileReader(fname));
        String s;

        //now that we have the length lets allocate our storage

        mp = new int[rlen][61][4][6]; //61 phred scores, 4 ref bases, 5 subst bases + vert_cumsum

        int pos = 0;
        //read the mutation spectrum
        //#Pos	Phred	A->A	A->C	A->G	A->T	A->N	C->A    \
        //  C->C	C->G	C->T	C->N	G->A	G->C	G->G	\
        //  G->T        G->N	T->A	T->C	T->G	T->T	T->N
        while (((s = fh.readLine()) != null) && (pos < rlen)) {
            if (!StringUtil.isBlankOrComment(s)) {
                int parts[] = parseLine(s);
                pos = parts[0];
                if (pos >= rlen) break;
                int phred = parts[1];
                if (phred > 60) throw new IOException("Phred score: "+Integer.toString(phred)+" in error profile is greater than 60!");
                int i = 2;
                for(int base = 0; base < 4; base++){
                    mp[pos][phred][base][0] = parts[i];
                    for (int j = 1; j < 5; j++) {

                        //store the cumulative frequency
                       mp[pos][phred][base][j] = mp[pos][phred][base][j - 1] + parts[j + i];
                    }
                    i+= 5; //increment to the next set of 5 X->A,...,X->N
                }
            }
        }

        //calculate the cumulative count of each occurance
        //for all positions, calculate the cumulative count of observations
        //under each phred score
        for(int i = 0; i < rlen; i++){
            for(int base = 0; base < 4; base++){
                mp[i][0][base][5] = mp[i][0][base][4];
                for(int j = 1;j <= 60; j++){
                    mp[i][j][base][5] = mp[i][j][base][4]+mp[i][j-1][base][5];
                }
            }
        }
    }




    public void AddErrorRead(SamRecord rec) {
        //reverse it from reference direction to sequencing direction
        if(rec.query_reverse_strand){
            rec.seqLine.reverseComplement();
            rec.qualLine.reverse();
        }

        //Now add error to the sequence returning the modified sequence
        for (int i = 0; i < rec.seqLine.length(); i++) {
            int base = indexBase(rec.seqLine.charAt(i));
            if(base!=indexBase('N')){
                //choose our phred score
                int r;
                r = rand.nextInt(mp[i][60][base][5]);//total for this position
                int phred = 60;
                for(int p = 0; p <= 60; p++){
                    if(r < mp[i][p][base][5]){
                        phred = p;
                        break;
                    }
                }

                //now choose our base given this phred score, ref base, and pos
                r = rand.nextInt(mp[i][phred][base][4]);
                int ind = 5;
                for(int c = 0; c < 5; c++){
                    if(r < mp[i][phred][base][c]){
                        ind = c;
                        break;
                    }
                }
                if(phred33) rec.qualLine.setCharAt(i,QualityUtil.getPhred33ScoreFromPhredScore(phred));
                else rec.qualLine.setCharAt(i, QualityUtil.getPhred64ScoreFromPhredScore(phred));
                rec.seqLine.sub(i, baseIndex(ind));
            } else { //we have an N
                if(phred33)rec.qualLine.setCharAt(i, QualityUtil.getPhred33ScoreFromPhredScore(0));
                else rec.qualLine.setCharAt(i, QualityUtil.getPhred64ScoreFromPhredScore(0));
            }
        }

        //reverse it back to reference direction
        if(rec.query_reverse_strand){
            rec.seqLine.reverseComplement();
            rec.qualLine.reverse();
        }

        //return rec;
    }

    private char baseIndex(int i){
        switch(i){
            case 0 : return 'A';
            case 1 : return 'C';
            case 2 : return 'G';
            case 3 : return 'T';
            default: return 'N';
        }
    }
    
    private int indexBase(char c){
        switch(Character.toUpperCase(c)){
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default : return 4;
        }
    }

}
