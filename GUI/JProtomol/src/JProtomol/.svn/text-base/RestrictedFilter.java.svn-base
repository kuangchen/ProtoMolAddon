/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package JProtomol;
import java.io.*;
import javax.swing.filechooser.FileFilter;
/**
 *
 * @author chris
 */
public class RestrictedFilter extends FileFilter {
     public boolean accept(File file)
     {
         if(file.isDirectory())
             return false;
         return true;
     }

     public String getDescription()
     {
         return "this folder only";
     }

}
