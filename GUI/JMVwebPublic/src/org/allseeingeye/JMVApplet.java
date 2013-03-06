package org.allseeingeye;

import java.applet.*;
import java.awt.*;
import java.io.*;
import java.net.*; 

import javax.media.opengl.*;
import com.sun.opengl.util.*;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

public class JMVApplet extends Applet {
  private Animator animator;

  public void init() {
    //JMolecularViewer jmv;
    String rport, hostip, rgbs, myip;
    float [] rgb = {0.5f, 0.5f, 0.5f};
    
    //inline parameters
    rport = getParameter("randport");
    if(rport == null) rport = "0xCE11";
    hostip = getParameter("hostip");
    if(hostip == null) hostip = "localhost";
    myip = getParameter("myip");
    if(myip == null) myip = "localhost";
    rgbs = getParameter("red");
    if(rgbs != null) rgb[0] = Float.parseFloat(rgbs);
    rgbs = getParameter("green");
    if(rgbs != null) rgb[1] = Float.parseFloat(rgbs);
    rgbs = getParameter("blue");
    if(rgbs != null) rgb[2] = Float.parseFloat(rgbs);
    //
    setLayout(new BorderLayout());
    GLCanvas canvas = new GLCanvas();
    jmv = new JMolecularViewer( true );
    //jmv.port = Integer.parseInt(rport);
    canvas.addGLEventListener(jmv);//new JMolecularViewer());
    canvas.setSize(getSize());
    add(canvas, BorderLayout.CENTER);
    animator = new FPSAnimator(canvas, 60);
    jmv.port = Integer.parseInt(rport);
    jmv.host = hostip;
    for(int i=0;i<3;i++) jmv.back_rgb[i] = rgb[i];
    jmv.myip = myip;
    jmv.applet = true;
  }
  
  JMolecularViewer jmv;
  
  public void change_port(String inpPort) {
     int tPort = 0;
     if(inpPort!=null){
         try{
            tPort = Integer.valueOf(inpPort);
         }catch(Exception efl){
         }
     }
     if(tPort <= 0 || tPort > 60000){
        jmv.port++;
        jmv.mol.port++;
     }else{
        jmv.port = tPort;
        jmv.mol.port = tPort;
     }
     //close port/reset flags
     jmv.mol.closePort();
     jmv.mol.clearFlags();
     jmv.jopengl.data_available = false;
     //reset view
     jmv.mview.resetView();
     //
     //if(inpPort!=null) System.out.println("Inline data = "+inpPort);
  }
  
  public String get_port() {
      return String.valueOf(jmv.port);
  }

  public void start() {
    animator.start();
  }

  public void stop() {
    animator.stop();

  }
  
 
}
