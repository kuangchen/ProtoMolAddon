/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package jpmolconfig;

import java.awt.Toolkit;
import javax.swing.JFrame;
import javax.swing.JPanel;

/**
 *
 * @author Chris
 */
public class JPmolConfig {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {        

        // TODO code application logic here
        
        JFrame frame = new JFrame("JProtomol Config generator");
        
        //ClassLoader cl = frame.getClass().getClassLoader();
        //frame.setIconImage(Toolkit.getDefaultToolkit().getImage(cl.getResource("icon.gif")));

        JPanel panel = new configJPanel();
        //frame.setSize(256,256);
        frame.setDefaultCloseOperation(frame.EXIT_ON_CLOSE);
        frame.getContentPane().add(panel);
        frame.pack();
        frame.setVisible(true);//show();
    }

}
