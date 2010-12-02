package simseq;

import java.util.HashMap;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author john
 */
public final class QualityUtil {
    private static final float[] errorProbabilityByPhredScore;
    private static final HashMap<Character,Integer>asciiPhred33ToPhredScore;
    private static final HashMap<Character,Integer>asciiPhred64ToPhredScore;
    private static final char[] phredScoreToAsciiPhred33;
    private static final char[] phredScoreToAsciiPhred64;
    static {
        errorProbabilityByPhredScore = new float[101];
        asciiPhred33ToPhredScore = new HashMap<Character,Integer>(101);
        phredScoreToAsciiPhred33 = new char[101];
        asciiPhred64ToPhredScore = new HashMap<Character,Integer>(101);
        phredScoreToAsciiPhred64 = new char[101];
        for (int i=0; i<errorProbabilityByPhredScore.length; ++i) {
            errorProbabilityByPhredScore[i] = (float) (1d / Math.pow(10d, i / 10d));
            asciiPhred33ToPhredScore.put(((char)(i+33)), i);
            phredScoreToAsciiPhred33[i] = ((char)(i+33));
            asciiPhred64ToPhredScore.put(((char)(i+64)), i);
            phredScoreToAsciiPhred64[i] = ((char)(i+64));
        }


    }

    /** Given a phred score between 0 and 100 returns the probability of error. */
    public static float getErrorProbabilityFromPhredScore(final int i) {
        return errorProbabilityByPhredScore[i];
    }

    /** Gets the phred score for any given probability of error. */
    public static int getPhredScoreFromErrorProbability(final float probability) {
        return (int) Math.round(-10 * Math.log10(probability));
    }

    public static float getErrorProbabilityFromAscii33PhredScore(final char c){
        return getErrorProbabilityFromPhredScore(asciiPhred33ToPhredScore.get(c));
    }

    public static char getPhred33ScoreFromPhredScore(final int i){
        return phredScoreToAsciiPhred33[i];
    }

    public static char getPhred33ScoreFromErrorProbability(final float probability){
        return phredScoreToAsciiPhred33[getPhredScoreFromErrorProbability(probability)];
    }

    public static float getErrorProbabilityFromAscii64PhredScore(final char c){
        return getErrorProbabilityFromPhredScore(asciiPhred64ToPhredScore.get(c));
    }

    public static char getPhred64ScoreFromPhredScore(final int i){
        return phredScoreToAsciiPhred64[i];
    }

    public static char getPhred64ScoreFromErrorProbability(final float probability){
        return phredScoreToAsciiPhred64[getPhredScoreFromErrorProbability(probability)];
    }
}
