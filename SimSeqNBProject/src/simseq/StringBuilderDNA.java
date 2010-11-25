/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package simseq;

/**
 *
 * @author john
 */
public class StringBuilderDNA{
    private StringBuilder psb;
    public StringBuilderDNA(StringBuilder sb){
        this.psb = sb;
    }
    public StringBuilderDNA(String s){
        this.psb = new StringBuilder(s);
    }
    @Override
    public String toString(){
        return this.psb.toString();
    }
    public void sub(int i,char c){
        this.psb.setCharAt(i, c);
    }
    public void reverseComplement(){
        this.psb = this.psb.reverse();
        for(int i=0;i<psb.length();i++){
            this.psb.setCharAt(i, this.complement(this.psb.charAt(i)));
        }
    }
    public void reverse(){
        this.psb = this.psb.reverse();
    }
    public void replace(String nstr){
        this.psb.replace(0, nstr.length(), nstr);
    }
    public void replace(int from, int to, String nstr){
        this.psb.replace(from, to, nstr);
    }
    private char complement(char c){
        switch(Character.toUpperCase(c)){
            case 'A': return 'T';
            case 'T': return 'A';
            case 'G': return 'C';
            case 'C': return 'G';
            default: return 'N';
        }
    }
    public char charAt(int i){
        return this.psb.charAt(i);
    }
    public int length(){
        return this.psb.length();
    }

}
