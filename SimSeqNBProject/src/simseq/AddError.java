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

    private float[][] mpA;
    private float[][] mpC;
    private float[][] mpG;
    private float[][] mpT;
    private ArrayList<float[]> qhA;
    private ArrayList<float[]> qhC;
    private ArrayList<float[]> qhG;
    private ArrayList<float[]> qhT;
    private Random rand = new Random();

    public AddError(String errorFname, boolean debug) throws IOException {
        /*
         * Constructor, simply read in the error file and get ready to
         * add error to any reads that come by.
         */
        parseErrorFile(errorFname);
    }

    private void parseErrorFile(String fname) throws FileNotFoundException, IOException {
        BufferedReader fh = new BufferedReader(new FileReader(fname));
        String s;
        int rlen = 0;
        //Get the read length
        while ((s = fh.readLine()) != null) {
            if (!StringUtil.isBlankOrComment(s)) {
                rlen = Integer.parseInt(s);
                break;
            }
        }
        if (rlen <= 0) {
            throw new IOException("Invalid read length: " + Integer.toString(rlen));
        }

        //now that we have the length lets allocate our storage

        mpA = new float[rlen][3];
        mpC = new float[rlen][3];
        mpG = new float[rlen][3];
        mpT = new float[rlen][3];
        qhA = new ArrayList<float[]>(rlen);
        qhC = new ArrayList<float[]>(rlen);
        qhG = new ArrayList<float[]>(rlen);
        qhT = new ArrayList<float[]>(rlen);


        //read the mutation spectrum
        int i = 0;
        while (((s = fh.readLine()) != null) && (i < rlen)) {
            if (!StringUtil.isBlankOrComment(s)) {
                String parts[] = s.split("\\s");
                float total = sumWords(parts, 1, 4);
                mpA[i][0] = Float.parseFloat(parts[1]) / total;
                for (int j = 1; j < 3; j++) {
                    mpA[i][j] = mpA[i][j - 1] + Float.parseFloat(parts[j + 1]) / total;
                }
                total = sumWords(parts, 4, 7);
                mpC[i][0] = Float.parseFloat(parts[4]) / total;
                for (int j = 1; j < 3; j++) {
                    mpC[i][j] = mpC[i][j - 1] + Float.parseFloat(parts[j + 4]) / total;
                }
                total = sumWords(parts, 7, 10);
                mpG[i][0] = Float.parseFloat(parts[7]) / total;
                for (int j = 1; j < 3; j++) {
                    mpG[i][j] = mpG[i][j - 1] + Float.parseFloat(parts[j + 7]) / total;
                }
                total = sumWords(parts, 10, 13);
                mpT[i][0] = Float.parseFloat(parts[10]) / total;
                for (int j = 1; j < 3; j++) {
                    mpT[i][j] = mpT[i][j - 1] + Float.parseFloat(parts[j + 10]) / total;
                }
                i++;
            }
        }//done with mutation spectrum

        //read and allocate the mutation histogram
        i = 0;
        while (((s = fh.readLine()) != null) && (i < (rlen * 4))) {
            if (!StringUtil.isBlankOrComment(s)) {

                String[] parts = s.split("\\s");
                char base = parts[1].charAt(0);
                String[] hist = parts[3].split(",");
                int[] phred = new int[hist.length];
                int[] counts = new int[hist.length];
                for (int j = 0; j < hist.length; j++) {
                    String[] tmp = hist[j].split(":");
                    phred[j] = Integer.parseInt(tmp[0]);
                    counts[j] = Integer.parseInt(tmp[1]);
                }
                int maxPhred = phred[phred.length - 1];
                float[] cumProb = new float[maxPhred + 1];
                int total = sum(counts);
                for (int j = 0; j < phred.length; j++) {
                    cumProb[phred[j]] = ((float) counts[j]) / total;
                }
                for (int j = 1; j < cumProb.length; j++) {
                    cumProb[j] += cumProb[j - 1];
                }
                switch (Character.toUpperCase(base)) {
                    case 'A':
                        qhA.add(cumProb);
                        break;
                    case 'C':
                        qhC.add(cumProb);
                        break;
                    case 'G':
                        qhG.add(cumProb);
                        break;
                    case 'T':
                        qhT.add(cumProb);
                        break;
                    default:
                        throw new IOException("Bad base in error distribution ");
                }



            }
        }
    }

    private int sum(int[] a) {
        int total = 0;
        for (int i = 0; i < a.length; i++) {
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
            float[] qhp = null;
            float[] mpp = null;
            switch (Character.toUpperCase(seq.charAt(i))) {
                case 'A':
                    qhp = qhA.get(i);
                    mpp = mpA[i];
                    break;
                case 'C':
                    qhp = qhC.get(i);
                    mpp = mpC[i];
                    break;
                case 'G':
                    qhp = qhG.get(i);
                    mpp = mpG[i];
                    break;
                case 'T':
                    qhp = qhT.get(i);
                    mpp = mpT[i];
                    break;
                default:
                    qhp = null;
                    mpp = null;
                    break;
            }//end switch

            if (qhp != null && mpp != null) {
                float rf = rand.nextFloat();
                for (int k = 0; k < qhp.length; k++) {
                    if (rf <= qhp[k]) {
                        score.setCharAt(i, QualityUtil.getPhred33ScoreFromPhredScore(k));
                        break;
                    }
                }// end get score from distribution
                rf = rand.nextFloat();
                if (rf <= QualityUtil.getErrorProbabilityFromAscii33PhredScore(score.charAt(i))) {
                    rf = rand.nextFloat();
                    for (int k = 0; k < 3; k++) {
                        if (rf <= mpp[k]) {
                            seq.sub(i, mutateChar(seq.charAt(i), k));
                            break;
                        }
                    }
                }
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

    private int sumWords(String[] s, int from, int to) {
        int total = 0;
        for (int i = from; i < to; i++) {
            total += Integer.parseInt(s[i]);
        }
        return total;
    }

    private char mutateChar(char ori, int i) {
        switch (Character.toUpperCase(ori)) {
            case 'A':
                switch (i) {
                    case 0:
                        return 'C';
                    case 1:
                        return 'G';
                    default:
                        return 'T';
                }
            case 'C':
                switch (i) {
                    case 0:
                        return 'A';
                    case 1:
                        return 'G';
                    default:
                        return 'T';
                }
            case 'G':
                switch (i) {
                    case 0:
                        return 'A';
                    case 1:
                        return 'C';
                    default:
                        return 'T';
                }
            default:
                switch (i) {
                    case 0:
                        return 'A';
                    case 1:
                        return 'C';
                    default:
                        return 'G';
                }
        }
    }
}
