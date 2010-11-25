/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package simseq;

/**
 *
 * @author jstjohn
 */
public class FastaRecord {
    private final String seq;
    private final String id;
    public FastaRecord(final String idline,final String sequence){
        this.seq=sequence;
        this.id=idline;
    }
    public String getID(){
        return this.id;
    }
    public String getSeq(){
        return this.seq;
    }
}
