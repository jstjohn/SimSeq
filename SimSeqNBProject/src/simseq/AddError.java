/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package simseq;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

/**
 *
 * @author john
 */
public class AddError {

    private int[][][] mpA;
    private int[][][] mpC;
    private int[][][] mpG;
    private int[][][] mpT;
    private Random rand = new Random();
    //private boolean phred33;

    public AddError(String errorFname, int rlen, boolean debug) throws IOException {
        /*
         * Constructor, simply read in the error file and get ready to
         * add error to any reads that come by.
         */
        //phred33 = phred_33; //for now this is not an option
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

        mpA = new int[rlen][62][7]; //61 phred scores + sum and 5 ints + sum + vert_cumsum
        mpC = new int[rlen][62][7];
        mpG = new int[rlen][62][7];
        mpT = new int[rlen][62][7];

        int pos = 0;
        //read the mutation spectrum
        while (((s = fh.readLine()) != null) && (pos < rlen)) {
            if (!StringUtil.isBlankOrComment(s)) {
                int parts[] = parseLine(s);
                pos = parts[0];
                if (pos >= rlen) break;
                int phred = parts[1];
                if (phred > 60) throw new IOException("Phred score: "+Integer.toString(phred)+" in error profile is greater than 60!");
                int i = 2;
           
                mpA[pos][phred][0] = parts[2];
                for (int j = 1; j < 5; j++) {
                    mpA[pos][phred][j] = mpA[pos][phred][j - 1] + parts[j + 2];
                }
                mpA[pos][phred][5] = sum(parts, i, i+5);//total
                mpA[pos][61][6] += mpA[pos][phred][5];//increment total phred for this pos

                
                mpC[pos][phred][0] = parts[7];
                for (int j = 1; j < 5; j++) {
                    mpC[pos][phred][j] = mpC[pos][phred][j - 1] + parts[j + 7];
                }
                mpC[pos][phred][5] = sum(parts, 7, 12);//total
                mpC[pos][61][6] += mpC[pos][phred][5];//increment total phred for this pos

                mpG[pos][phred][0] = parts[12];
                for (int j = 1; j < 5; j++) {
                    mpG[pos][phred][j] = mpG[pos][phred][j - 1] + parts[j + 12];
                }
                mpG[pos][phred][5] = sum(parts, 12, 17);//total
                mpG[pos][61][6] += mpG[pos][phred][5];//increment total phred for this pos

                mpT[pos][phred][0] = parts[17];
                for (int j = 1; j < 5; j++) {
                    mpT[pos][phred][j] = mpT[pos][phred][j - 1] + parts[j + 17];
                }
                mpT[pos][phred][5] = sum(parts, 17, 22);//total
                mpT[pos][61][6] += mpT[pos][phred][5];//increment total phred for this pos
            }
        }//done with mutation spectrum

        //calculate the cumulative count of each occurance
        for(int i = 0; i < rlen; i++){
            mpA[i][0][6] = mpA[i][0][5];
            mpC[i][0][6] = mpC[i][0][5];
            mpG[i][0][6] = mpG[i][0][5];
            mpT[i][0][6] = mpT[i][0][5];
            for(int j = 1;j <= 60; j++){
                mpA[i][j][6] = mpA[i][j][5]+mpA[i][j-1][6];
                mpC[i][j][6] = mpC[i][j][5]+mpC[i][j-1][6];
                mpT[i][j][6] = mpT[i][j][5]+mpT[i][j-1][6];
                mpG[i][j][6] = mpG[i][j][5]+mpG[i][j-1][6];
            }
        }
    }

    private int sum(int[] a, int from, int to) {
        int total = 0;
        for (int i = from; i < to; i++) {
            total += a[i];
        }
        return total;
    }



    public void AddErrorRead(SamRecord rec) {

        StringBuilderDNA seq = rec.seqLine; //pointer
        StringBuilder score = rec.qualLine; //pointer


        //reverse it from reference direction to sequencing direction
        if(rec.query_reverse_strand){
            seq.reverseComplement();
            score.reverse();
        }

        //Now add error to the sequence returning the modified sequence
        for (int i = 0; i < seq.length(); i++) {
            int[][] mpp = null;
            switch (Character.toUpperCase(seq.charAt(i))) {
                case 'A':
                    mpp = mpA[i];
                    break;
                case 'C':
                    mpp = mpC[i];
                    break;
                case 'G':
                    mpp = mpG[i];
                    break;
                case 'T':
                    mpp = mpT[i];
                    break;
                default:
                    mpp = null;
                    break;
            }//end switch
            if(mpp != null){
                //choose our phred score
                int r;
                r = rand.nextInt(mpp[61][6]);//total for this position
                int phred = 60;
                for(int p = 0; p <= 60; p++){
                    if(r < mpp[p][6]){
                        phred = p;
                        break;
                    }
                }

                //now choose our base given this phred score, ref base, and pos
                r = rand.nextInt(mpp[phred][5]);
                int ind = 5;
                for(int c = 0; c < 5; c++){
                    if(r < mpp[phred][c]){
                        ind = c;
                        break;
                    }
                }
                //if(phred33) //for now this is not an option
                score.setCharAt(i,QualityUtil.getPhred33ScoreFromPhredScore(phred));
                //else score.setCharAt(i, QualityUtil.getPhred64ScoreFromPhredScore(phred));
                seq.sub(i, baseIndex(ind));
            } else { //we have an N
                score.setCharAt(i, QualityUtil.getPhred33ScoreFromPhredScore(0));
            }
        }

        //reverse it back to reference direction
        if(rec.query_reverse_strand){
            seq.reverseComplement();
            score.reverse();
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
