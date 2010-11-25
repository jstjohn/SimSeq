/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package simseq;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 *
 * @author john
 */
public class AddDiploid {
    /*Constructor*/

    private HashMap<String, HashMap<Integer, Character>> hmp;

    public AddDiploid(String dipFname, boolean debug) throws IOException {
        readDiploid(dipFname);
    }
    public void AddDiploidRead(SamRecord rec){
        
    }
    private void readDiploid(String dipfName) throws IOException {
        hmp = new HashMap<String, HashMap<Integer, Character>>(22);
        BufferedReader fh = new BufferedReader(new FileReader(dipfName));
        String s;
        while ((s = fh.readLine()) != null) {
            String[] f = s.split("\\s");
            //chr1	840	T
            if (f.length == 3) {
                if (!hmp.containsKey(f[0])) {
                    hmp.put(f[0], new HashMap<Integer,Character>());
                }
                hmp.get(f[0]).put(Integer.parseInt(f[1]), f[2].charAt(0));
            }
        }
    }
}
