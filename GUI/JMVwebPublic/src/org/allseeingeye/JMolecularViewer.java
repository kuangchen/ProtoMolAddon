package org.allseeingeye;

import java.awt.*;
import java.awt.event.*;

import com.sun.opengl.util.*;
import java.net.*;
import java.io.*;
import java.util.*;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;

import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLCanvas;
import javax.media.opengl.GLEventListener;
import com.sun.opengl.util.j2d.*;
import java.awt.geom.*;

import javax.swing.JOptionPane;
import javax.swing.JFileChooser;

import javax.media.opengl.GLCapabilities;
/**
 * JMolecularViewer.java 
 * author: Chris Sweet, James Sweet.
 *
 */

public class JMolecularViewer implements GLEventListener, MouseListener, MouseMotionListener, MouseWheelListener 
{
        public static Frame frame=null;
        
	public static void main(String[] args)
	{
		frame = new Frame("LCLS Molecular Viewer v3.0 - University of Notre Dame");
                frame.setIconImage(Toolkit.getDefaultToolkit().getImage("icon.gif"));
                
                //anti- alias
                GLCapabilities capabilities = new GLCapabilities();
                capabilities.setSampleBuffers(true);
                capabilities.setNumSamples(4);
                GLCanvas canvas = new GLCanvas(capabilities);

		//GLCanvas canvas = new GLCanvas();
		
		canvas.addGLEventListener(new JMolecularViewer( false ));
		frame.add(canvas);
		frame.setSize(640, 480);
                frame.setMinimumSize(new Dimension (320, 240));
		final Animator animator = new Animator(canvas);
		frame.addWindowListener(new WindowAdapter()
		{
			public void windowClosing(WindowEvent e)
			{
				// Run this on another thread than the AWT event queue to
				// make sure the call to Animator.stop() completes before
				// exiting
				new Thread(new Runnable()
				{
					public void run()
					{
						animator.stop();
						System.exit(0);
					}
				}).start();
			}
		});
		// Center frame
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
		animator.start();
	}

    public JMolecularViewer( boolean isApplet ) {

    }


        public boolean applet = false;
        public Molecule mol = null;
        public String host = "localhost";
        public String myip = "localhost";
        public int port = 0xCE11;
        private int host_wait = 0;
        //OpenGL
        public JOpenGL jopengl = null;
        //Mouse variables
        private int prevMouseX, prevMouseY;
        private boolean mouseRButtonDown = false;
        private boolean mouseButtonDown = false;
        //view control
        public QuartView mview = null;
        static float rangle=0.0f;
        //background color
        public float [] back_rgb = { 0.2305f, 0.3477f, 0.5938f};
        //display mode
        public int display_mode = 0;
        private boolean was_rotating;
        //save image
        private boolean captureImage = false;
        private File capFile = null;
        //timing diagnostics
        private long dispCount = 0, startTime, renderTime;
               
        public void init(GLAutoDrawable drawable)
	{
                //create OpenGL
                jopengl = new JOpenGL(drawable, back_rgb, myip, applet);
                // Create and start the comms thread
                mol = new Molecule(host,port); //new molecule class 
                Thread thread = new commsThread();
                thread.start();
                //
                //listners and view
                mview = new QuartView(640,480);
                drawable.addMouseListener(this);
                drawable.addMouseMotionListener(this);
                drawable.addMouseWheelListener(this);
                //
                
	}
        
	public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height)
	{
            //reshape
            jopengl.reshape(drawable, x, y, width, height, mview);
                
	}
	
        
	public void display(GLAutoDrawable drawable)
	{
            
            //find start time
            startTime = System.currentTimeMillis ();
            //draw molecule
            jopengl.display(drawable, mol, mview, rangle, 
                                mouseButtonDown && !mouseRButtonDown, display_mode);
            //draw zoom
            jopengl.zoom_buttons(drawable, mol, mview, host, port, host_wait, display_mode);
            //Capture?
            if(captureImage){
                captureImage = false;
                //save
		if( capFile != null ) {
                    try {
                        Screenshot.writeToFile(capFile,drawable.getWidth(), drawable.getHeight());
                        //System.out.println("File saved");
                    } catch (Exception ex) {//IO
                        System.out.println("Error saving file: " + ex);
                    }
		}
            }
        }
        
	public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean deviceChanged)
	{}
        
//****Communication thread code************************************************/

        // This class extends Thread
        class commsThread extends Thread {
        // This method is called when the thread runs
            public void run() {
                while(true){
                    //get data
                    try{
                        socketCom();
                    }catch (IOException e){
                        //System.out.println("Failed");  
                        mol.data_valid = false;
                        mol.metadata_valid = false;
                        mol.molRad = 0.0f;
                        mol.port_open = false; //closed if failed DONT NEED
                        host_wait++;            //screen graphics for waiting
                        if(host_wait > 4) host_wait = 0;
                        //close port?
                        mol.closePort();
                        //mol.selected_atom = 0;  //remove selection
                    }
                    //then wait so no comms over-run
                    try{
                        Thread.sleep(50);
                    }catch (InterruptedException ie){
                    }
                }
            }
        } 
        
        public void socketCom() throws IOException {
            
                if(!mol.port_open) mol.openPort();
                if(!mol.metadata_valid) mol.getMetadata();
                if(mol.metadata_valid){ 
                    //mol.openPort();
                    mol.getData(); 
                    //System.out.println("Read data "+String.valueOf(mol.data_valid));
                    if(!mol.data_valid) System.out.println("Read data error");
                    //
                }else{
                    System.out.println("Metadata read error");
                    mol.selected_atom = 0;
                }
                //

        }

      // Methods required for the implementation of MouseListener
      public void mouseEntered(MouseEvent e) {}
      public void mouseExited(MouseEvent e) {}

      public void mousePressed(MouseEvent e) {
        //
        buttons2D(e.getX(), e.getY());   //check for 2D button press
        if(jopengl.zflag[0]) mview.scaleMult(1.15f);
        if(jopengl.zflag[1]) mview.scaleMult(0.87f);
        if(jopengl.zflag[2]){
            display_mode++;
            if(display_mode > 3) display_mode = 0;
            if(display_mode == 1 && mol.a_carb_count == 0) display_mode = 3;
        }
        //capture screen
        if(!applet && jopengl.zflag[3] && !captureImage){
            //get file name
            JFileChooser chooser = new JFileChooser();
            chooser.setDialogTitle("Save PNG File");
            int result = chooser.showSaveDialog(frame);
            if(result == JFileChooser.APPROVE_OPTION){
                capFile = chooser.getSelectedFile();
                captureImage = true;
            }else{
                capFile = null;
                if(result != JFileChooser.CANCEL_OPTION)
                    JOptionPane.showMessageDialog(frame, "Error saving file.");
            }
        }
        //change host/port
        if(!applet && jopengl.zflag[4]) setHostPort();
        //non 2D
        if(!(jopengl.zflag[0] || jopengl.zflag[1] || jopengl.zflag[2] || jopengl.zflag[3] || jopengl.zflag[4])){
            mview.set_start_v(e.getX(), e.getY());
            prevMouseX = e.getX();
            prevMouseY = e.getY();
            mouseButtonDown = true;
            if ((e.getModifiers() & e.BUTTON3_MASK) != 0) {
              mouseRButtonDown = true;
            }else{
                //test if was rotating/disable select if so
                if(Math.abs(mview.angle_inc)>1.0f) was_rotating = true;
                else was_rotating = false;
                //
                mview.last_angle  = mview.current_rot_a;
                mview.angle_inc = 0.0f;
            }
            //save current mouse data
            mview.scaleSave();
        }
        //rotation
        rangle = 0.0f;
     }

      public void mouseReleased(MouseEvent e) {
          mouseButtonDown = false;
          boolean mouse_stationary = (Math.abs(prevMouseX - e.getX()) < 5 && 
                        Math.abs(prevMouseY - e.getY()) < 5);
          if ((e.getModifiers() & e.BUTTON3_MASK) != 0) {
              //about box
              if(mouse_stationary && !applet){
                  JOptionPane.showMessageDialog(frame,"LCLS Java Molecular Viewer\nVersion 3.0 July 2008\nAuthors: Chris Sweet, James Sweet.\n"+
                          "http://www.nd.edu/~lcls/\nhttp://protomol.sourceforge.net",
                          "About JMV",JOptionPane.PLAIN_MESSAGE);
              }
              mouseRButtonDown = false;
          }else{
              //link to website?
              if(mouse_stationary && !was_rotating){
                    jopengl.selMouseX = e.getX();jopengl.selMouseY = e.getY();
                    jopengl.select = true;
                    //getSelection(glp, e.getX(), e.getY());
                    //try{
                    //    Desktop dsk = Desktop.getDesktop();
                    //    dsk.browse(new URI("http://www.nd.edu/~lcls/index.html"));
                    //}catch(Exception ex){
                    //}
              }
          }
          //jopengl.zflag[0] = jopengl.zflag[1] = false;
      }

      public void mouseClicked(MouseEvent e) {
      }

      // Methods required for the implementation of MouseMotionListener
      public void mouseDragged(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        Dimension size = e.getComponent().getSize();

        float thetaY = -360.0f * ( (float)(prevMouseX-x)/(float)size.width);
        float thetaX = -360.0f * ( (float)(prevMouseY-y)/(float)size.height);

        if(mouseRButtonDown)      
            mview.scaleSum(-thetaX/50.0f - thetaY/50.0f);
        else
            mview.set_current_v(e.getX(), e.getY());
      }

      public void mouseMoved(MouseEvent e) {
            buttons2D(e.getX(), e.getY());
      }
      
      public void mouseWheelMoved(MouseWheelEvent e) {
            mview.scaleMult(1.0f + e.getScrollAmount() * e.getWheelRotation()/20.0f);
      }  

      private void buttons2D(float x, float y){

            int numButtons = 3;
            if(!applet) numButtons = 5;
            for(int i=0;i<numButtons;i++) jopengl.zflag[i] = false;
            if(x > jopengl.zrect[0][0] && x < jopengl.zrect[0][1]){
                for(int i=0;i<numButtons;i++){
                    if(y < jopengl.zrect[i][2] && y > jopengl.zrect[i][3])
                        jopengl.zflag[i] = true;
                }
            }
      }
      
      //Host/port Utilities
      private void setHostPort(){
            Object[] possibilities = {"new connection", null, null, null, null, null};
            int noEntries = 0;
            try{
                BufferedReader in = new BufferedReader(new FileReader("hosts.info"));
                String fst;
                while((fst = in.readLine()) != null) possibilities[++noEntries] = fst;
                in.close();
            }catch(IOException er){
            }
            //
            String s = (String)JOptionPane.showInputDialog(
                    frame,
                    "Select host/port",
                    "Select new host and port for simulation",
                    JOptionPane.PLAIN_MESSAGE,
                    null,
                    possibilities,
                    "new connection");
            //get host/port
            if(s != null && s.length() > 0 && s.compareTo("new connection")==0){
                s = (String)JOptionPane.showInputDialog(
                    frame,
                    "Set new host:port, e.g. myUrl:12345",
                    "Set host and port for simulation",
                    JOptionPane.QUESTION_MESSAGE,
                    null,
                    null,
                    host+":"+port);
            }
            //Parse string
            if(s != null && s.length() > 0){
                int index = s.lastIndexOf(':');
                if(index > 0){
                    host = s.substring(0,index);
                    String s2 = s.substring(index + 1);
                    try{
                        port = Integer.valueOf(s2);
                    }catch(Exception efl){
                        port = 0xCE11;
                    }
                    mol.port = port;
                    mol.hostname = host;
                    mol.closePort();
                    mol.clearFlags();
                    jopengl.data_available = false;
                    //reset view
                    mview.resetView();
                    //save new connection
                    boolean hpExists = false;   //exists?
                    for(int i=1;i<6;i++){
                        if(possibilities[i] != null && s.compareTo((String)possibilities[i])==0){
                            hpExists = true;
                            break;
                        }
                    }
                    if(!hpExists){
                        for(int i=5;i>1;i--) possibilities[i] = possibilities[i-1]; //ripple down
                        possibilities[1] = s;
                        try{
                            FileWriter fileWriter = new FileWriter("hosts.info",false); //false for not append
                            BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
                            for(int i=1;i<6;i++){
                                if(possibilities[i] != null){
                                    bufferedWriter.write((String)possibilities[i]);
                                    bufferedWriter.newLine();
                                }else{
                                    break;
                                }
                            }
                            bufferedWriter.close();
                        }catch(IOException ew){
                        }
                    }
                    //
                }
            }
      }
      
}

