package org.allseeingeye;

import java.awt.*;
import java.awt.event.*;

import com.sun.opengl.util.*;
import javax.media.opengl.glu.GLU;
import java.net.*;
import java.io.*;
import java.util.*;

import javax.media.opengl.GL;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.glu.GLU;
import javax.media.opengl.glu.GLUquadric;
import com.sun.opengl.util.j2d.*;
import java.awt.geom.*;
import javax.swing.*;
//
import com.sun.opengl.util.texture.*;
//
import java.nio.*;

/**
 *
 * @authors Chris Sweet & James Sweet
 */
public class JOpenGL {
        //open gl defines
        float front_shininess[] = { 60.0f };
        float front_specular[] = { 0.7f, 0.7f, 0.7f, 1.0f };
        float ambient0[] = { 0.4f, 0.4f, 0.4f, 1.0f };
        float diffuse0[] = { 1.0f, 1.0f, 1.0f, 1.0f };
        float ambient1[] = { 0.2f, 0.2f, 0.2f, 1.0f };
        float diffuse1[] = { 0.5f, 0.5f, 0.5f, 1.0f };
        float position0[] = { 1.0f, 1.0f, 1.0f, 0.0f };
        float position1[] = { -1.0f, -1.0f, 1.0f, 0.0f };
        float lmodel_ambient[] = { 0.5f, 0.5f, 0.5f, 1.0f };
        float lmodel_twoside[] = { GL.GL_TRUE };
        //axes etc.
        private int axes1;
        private int [] spheres;
        private int [] cylinders;
        private int tube_sphere;
        private int arrow;
        private TextRenderer renderer, rendererFlat;
        //zoom buttons, for 640x480
        int [][] zrect = {{554, 620, 22, 2},{554, 620, 44, 24},{554, 620, 66, 46},{554, 620, 88, 68},{554, 620, 110, 90}};
        Boolean [] zflag = {false, false, false, false, false};
        //background
        private Texture texture = null;
        private Texture [] buttonTex = {null, null, null, null, null, null, null, null, null, null};
        //Select
        public boolean select = false;
        public int selMouseX, selMouseY;
        int [][] boundingRect;
        private float [] a_rad = { 1.7f, 1.2f, 1.55f, 1.52f, 1.85f, 1.76f, 2.00f };
        //
        public Boolean data_available = false;
        //
        public float rear_clip = 100;
        //
        boolean isApplet=false; 
        int numButtons = 3;
        
        public JOpenGL(GLAutoDrawable drawable, float [] back_rgb, String myip, boolean applet){
            
                float [][] gl_col = { { 0.5f, 0.1f, 0.1f, 1.0f },
                                        { 0.1f, 0.5f, 0.1f, 1.0f },
                                        { 0.1f, 0.1f, 0.5f, 1.0f }};
                float [][] a_color = {
                    { 0.20f, 0.20f, 0.20f, 1.00f }, /* dark grey */
                    { 0.50f, 0.50f, 0.50f, 1.00f }, /* grey */
                    { 0.10f, 0.10f, 0.80f, 1.00f }, /* blue */
                    { 0.80f, 0.15f, 0.15f, 1.00f }, /* red */
                    { 0.00f, 0.80f, 0.80f, 1.00f },  /* purple - unknown */
                    { 0.35f, 0.04f, 0.34f, 1.00f },   /* dark purp */
                    { 1.00f, 0.87f, 0.00f, 1.00f },   /* gold */
                    { 1.00f, 0.00f, 0.00f, 1.00f }   /* bright red */
                }; 
                
                //applet?
                isApplet = applet;
                if(!isApplet) numButtons = 5;
                //text
                //renderer = new TextRenderer(new Font("SansSerif", Font.BOLD, 25));
                rendererFlat = new TextRenderer(new Font("SansSerif", Font.BOLD, 16));
		// Use debug pipeline
		// drawable.setGL(new DebugGL(drawable.getGL()));
		
		GL gl = drawable.getGL();
                GLU glu = new GLU();
                //
                drawable.setAutoSwapBufferMode(true);
		System.err.println("INIT GL IS: " + gl.getClass().getName());
		
		// Enable VSync
		gl.setSwapInterval(2);

                gl.glColor4f( 1.0f, 1.0f, 1.0f, 1.0f );
                gl.glClearColor( back_rgb[0], back_rgb[1], back_rgb[2], 1.0f );
                gl.glClearDepth( 1.0f );

                gl.glDepthFunc( GL.GL_LESS );
                gl.glEnable( GL.GL_DEPTH_TEST );
                //gl.glEnable( GL.GL_LINE_SMOOTH );
                //gl.glEnable( GL.GL_POINT_SMOOTH );
                //gl.glEnable( GL.GL_POLYGON_SMOOTH );
                //
                gl.glBlendFunc( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA );
                gl.glEnable( GL.GL_BLEND );

                /* lights */
                gl.glLightfv( GL.GL_LIGHT0, GL.GL_AMBIENT, ambient0, 0 );
                gl.glLightfv( GL.GL_LIGHT0, GL.GL_DIFFUSE, diffuse0, 0 );
                gl.glLightfv( GL.GL_LIGHT0, GL.GL_POSITION, position0, 0 );
                gl.glLightfv( GL.GL_LIGHT1, GL.GL_AMBIENT, ambient1, 0 );
                gl.glLightfv( GL.GL_LIGHT1, GL.GL_DIFFUSE, diffuse1, 0 );
                gl.glLightfv( GL.GL_LIGHT1, GL.GL_POSITION, position1, 0 );
                gl.glLightModelfv( GL.GL_LIGHT_MODEL_AMBIENT, lmodel_ambient, 0 );
                gl.glLightModelfv( GL.GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside, 0 );
                gl.glEnable( GL.GL_LIGHTING );
                gl.glEnable( GL.GL_LIGHT0 );
                gl.glEnable( GL.GL_LIGHT1 );

                gl.glEnable( GL.GL_NORMALIZE );
                gl.glEnable( GL.GL_CULL_FACE );
                gl.glShadeModel( GL.GL_SMOOTH );

                gl.glMaterialfv( GL.GL_FRONT_AND_BACK, GL.GL_SHININESS, front_shininess, 0 );
                gl.glMaterialfv( GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, front_specular, 0 );
                gl.glDrawBuffer( GL.GL_BACK );
                //create axes object
                axes1 = gl.glGenLists(1);
                gl.glNewList(axes1, GL.GL_COMPILE);
                //cylinders
                GLUquadric quadObj = glu.gluNewQuadric();
                for(int i=0;i<3;i++){
                    gl.glMaterialfv( GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, gl_col[i], 0);   
                    glu.gluSphere( quadObj, 0.08f, 10, 10 );
                    glu.gluCylinder( quadObj, 0.08f, 0.08f, 0.5f, 10, 10 );
                    gl.glTranslatef( 0.0f, 0.0f, 0.5f ); 
                    glu.gluCylinder( quadObj, 0.2f, 0.0f, 0.5f, 10, 10 );
                    gl.glTranslatef( 0.0f, 0.0f, -0.5f );
                    //
                    if(i==0) gl.glRotatef(-90f, 1.0f, 0.0f, 0.0f);  
                    else gl.glRotatef(90f, 0.0f, 1.0f, 0.0f);  
                }
                gl.glEndList();
                //spheres
                spheres = new int[7];
                for(int i=0;i<7;i++){
                    spheres[i] = gl.glGenLists(1);
                    gl.glNewList(spheres[i], GL.GL_COMPILE);
                    gl.glMaterialfv( GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, a_color[i], 0 );
                    glu.gluSphere( quadObj, 0.25f * a_rad[i], 16, 16 );
                    gl.glEndList();
                }
                //cylinders
                cylinders = new int[5];
                float length = 1.2f, bond_radius = 0.1f;
                for(int i=0;i<5;i++){
                    cylinders[i] = gl.glGenLists(1);
                    gl.glNewList(cylinders[i], GL.GL_COMPILE);
                    gl.glMaterialfv( GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, a_color[i], 0 );
                    glu.gluCylinder( quadObj, bond_radius, bond_radius, 0.5 * length, 8, 1 );
                    gl.glEndList();
                }
                //tube sphere
                tube_sphere = gl.glGenLists(1);
                gl.glNewList(tube_sphere, GL.GL_COMPILE);
                gl.glMaterialfv( GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, a_color[5], 0 );
                glu.gluSphere( quadObj, 0.44f, 16, 16 );
                gl.glEndList();
                //Arrow
                arrow = gl.glGenLists(1);
                gl.glNewList(arrow, GL.GL_COMPILE);
                gl.glMaterialfv( GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, a_color[7], 0 );
                    glu.gluPartialDisk(quadObj,0.47f,0.53f,20,10,0f,270f);
                    glu.gluPartialDisk(quadObj,0.42f,0.58f,20,10,250f,10f);
                    glu.gluPartialDisk(quadObj,0.44f,0.56f,20,10,250f,20f);
                    glu.gluPartialDisk(quadObj,0.46f,0.54f,20,10,250f,30f);
                    glu.gluPartialDisk(quadObj,0.48f,0.52f,20,10,250f,40f);
                    gl.glRotatef( 180f, 0.0f, 1.0f, 0.0f );	//draw other side!
                    glu.gluPartialDisk(quadObj,0.47f,0.53f,20,10,0f,-270f);
                    glu.gluPartialDisk(quadObj,0.42f,0.58f,20,10,-250f,-10f);
                    glu.gluPartialDisk(quadObj,0.44f,0.56f,20,10,-250f,-20f);
                    glu.gluPartialDisk(quadObj,0.46f,0.54f,20,10,-250f,-30f);
                    glu.gluPartialDisk(quadObj,0.48f,0.52f,20,10,-250f,-40f);
                gl.glEndList();
                //Background texture?
                if (texture == null) {
                    try {
                        //System.err.println("Loading background.");
                        if(isApplet){
                            TextureData texData = TextureIO.newTextureData(new URL(myip+"/bground.jpg"),false,"JPG");
                            texture = TextureIO.newTexture(texData);
                        }else{
                            File file = new File("bground.jpg");
                            texture = TextureIO.newTexture(file, true);
                        }
                        //System.err.println("Texture estimated memory size = " + texture.getEstimatedMemorySize());
                    } catch (IOException e) {
                        //e.printStackTrace();
                        //System.err.println("No background file.");
                    }
                }
                //button textures
                ClassLoader cl = this.getClass().getClassLoader();
                String [] bFile = {"buttonZoomPOff","buttonZoomPOn","buttonZoomMOff","buttonZoomMOn","buttonEyeOff","buttonEyeOn","buttonCamOff","buttonCamOn","buttonConOff","buttonConOn"};
                for(int i=0;i<numButtons*2;i++){
                    try {
                        buttonTex[i] = TextureIO.newTexture(cl.getResourceAsStream(bFile[i] + ".jpg"), false, "JPG");//
                    } catch (Exception e) {
                        System.err.println("No texture "+bFile[i] + ".jpg. "+ e.getMessage());
                    }
                }
                //
                boundingRect = new int[6][4];
        }
        
	public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height, QuartView mview)
	{
		GL gl = drawable.getGL();
		GLU glu = new GLU();

		if (height <= 0) // avoid a divide by zero error!
				height = 1;
		final float h = (float) width / (float) height;
		gl.glViewport(0, 0, width, height);
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glLoadIdentity();
		glu.gluPerspective(45.0f, h, 0.1, rear_clip);
		gl.glMatrixMode(GL.GL_MODELVIEW);
		gl.glLoadIdentity();
                //reset view sphere
                mview.sphereRad(width, height);
                //System.out.println("width "+ String.valueOf(width));
                
	}

	public void display(GLAutoDrawable drawable, Molecule mol, 
                                QuartView mview, float rangle, boolean lMouseDown, int display_mode)
	{
                int width, height;
                
                GL gl = drawable.getGL();
		GLU glu = new GLU();
                width = drawable.getWidth();
                height = drawable.getHeight();
                //Select?
                if(select){
                    select = false;
                    if(mol.data_valid || data_available)
                        getSelection(gl, selMouseX, selMouseY, width, height, mol, mview, display_mode);
                }
                // Clear the drawing area
                gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
                // Reset the current matrix to the "identity"
                gl.glLoadIdentity();
                //background?
                if (texture != null) {
                    //2D view
                    gl.glMatrixMode (GL.GL_PROJECTION);
                    gl.glLoadIdentity ();
                    gl.glOrtho (0., 1, 0., 1, -1., 1.);
                    gl.glMatrixMode (GL.GL_MODELVIEW);
                    gl.glLoadIdentity ();
                    gl.glDisable(GL.GL_DEPTH_TEST); 
                    //
                    texture.enable();
                    texture.bind();
                    gl.glTexEnvi(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_REPLACE);
                    TextureCoords coords = texture.getImageTexCoords();
                    //draw textured rectangle
                    gl.glBegin(GL.GL_QUADS);
                    gl.glTexCoord2f(coords.left(), coords.bottom());
                    gl.glVertex3f(0, 0, 0);
                    gl.glTexCoord2f(coords.right(), coords.bottom());
                    gl.glVertex3f(1, 0, 0);
                    gl.glTexCoord2f(coords.right(), coords.top());
                    gl.glVertex3f(1, 1, 0);
                    gl.glTexCoord2f(coords.left(), coords.top());
                    gl.glVertex3f(0, 1, 0);
                    gl.glEnd();
                    texture.disable();
                    //3D view
                    gl.glMatrixMode(GL.GL_PROJECTION);
                    gl.glLoadIdentity();
                    glu.gluPerspective(45.0f, (float)width/(float)height, 0.1, rear_clip);
                    gl.glMatrixMode(GL.GL_MODELVIEW);
                    gl.glLoadIdentity();
                    gl.glEnable(GL.GL_DEPTH_TEST);
                }
                //
                if(mol.data_valid || data_available){
                    data_available = true;
                    //Find radius of molecule
                    mol.molRadius();
                    //set rear clipping plane
                    float rClip = mol.molRad * 3.0f + mol.molRad * mview.scale_xyz;
                    if(rClip != rear_clip){
                        rear_clip = rClip;
                        if (height <= 0) height = 1;// avoid a divide by zero error!
                        final float h = (float) width / (float) height;
                        glu.gluPerspective(45.0f, h, 0.1, rear_clip);
                        //System.out.println("Changing view");
                    }
                    //
                    gl.glLoadIdentity();
                    gl.glPushMatrix();
                    //scale for full view of molecule
                    gl.glTranslatef( 0.0f, 0.0f, -mol.molRad * 3.0f  );
                    gl.glPushAttrib( GL.GL_ALL_ATTRIB_BITS );
                    //axes
                    //axes(gl, mol.molRad, width, height);
                    gl.glPushMatrix();
                    //rotations for molecule
                    if(lMouseDown){
                        mview.angle_inc = mview.current_rot_a - mview.last_angle;
                        float temp_ang_inc = Math.abs(mview.angle_inc);
                        if(temp_ang_inc > 50.0f) mview.angle_inc /= (temp_ang_inc / 50.0f);
                        mview.last_angle  = mview.current_rot_a;
                    }else{
                        if(Math.abs(mview.angle_inc) > 0.01f){
                           mview.angle_inc *= 0.99f; //damping
                           mview.set_current_rot();
                        }
                    }
                    mview.rot_a += rangle;  //rotate view?
                    if(mview.rot_a > 360.0f) mview.rot_a -= 360.0f;
                    //set rotation
                    gl.glRotatef(mview.rot_a, mview.rot_v[0], mview.rot_v[1], mview.rot_v[2]);
                    // scale
                    gl.glScalef(mview.scale_xyz,mview.scale_xyz,mview.scale_xyz);
                    //standard CPK or Space fill
                    if(display_mode == 0 || display_mode == 2 || display_mode == 3){
                        for ( int i = 0 ; i < mol.atom_count ; i++ ){
                            gl.glPushMatrix();
                            gl.glTranslatef( mol.coord[i][0], mol.coord[i][1], mol.coord[i][2] );
                            //scale if spacefill
                            if(display_mode == 3) gl.glScalef(4.0f,4.0f,4.0f);
                            if(i == mol.selected_atom - 1){ //selected atom?
                                float selScale = a_rad[mol.atom_color[i]] * 0.5375f; //7.5% larger
                                gl.glPushMatrix();
                                gl.glScalef(selScale,selScale,selScale);
                                gl.glCallList(spheres[6]);
                                gl.glPopMatrix();
                                //\Phi, \Psi if not spacefill
                                if(display_mode != 3 && mol.selected_atom_n != 0){
                                    phiPsiArrows(i, gl,mol,true);
                                }
                            }
                            else gl.glCallList(spheres[mol.atom_color[i]]);
                            gl.glPopMatrix();
                        }
                    }
                    //Bonds
                    if(display_mode == 0 || display_mode == 2){
                        for ( int atomp = 0 ; atomp < mol.bond_count ; atomp++ ){
                            float length, x, y, z;
  
                            x = mol.coord[mol.b[atomp]][0] - mol.coord[mol.a[atomp]][0];
                            y = mol.coord[mol.b[atomp]][1] - mol.coord[mol.a[atomp]][1];
                            z = mol.coord[mol.b[atomp]][2] - mol.coord[mol.a[atomp]][2];
                            length = (float) Math.sqrt( x * x + y * y + z * z );
                            gl.glPushMatrix();
                            gl.glTranslatef( mol.coord[mol.a[atomp]][0], mol.coord[mol.a[atomp]][1], mol.coord[mol.a[atomp]][2] );
                            // rotate by acos( z / |xyz| ) around xyz crossproduct (0, 0, 1]) 
                            gl.glRotatef( (float) Math.acos( z/length ) * 57.2957f, -y, x, 0.0f ); 
                            gl.glTranslatef( 0f, 0f, 0.5f * (length-1.2f));//length );
                            gl.glCallList(cylinders[mol.atom_color[mol.a[atomp]]]);
                            gl.glTranslatef( 0f, 0f, 0.5f * 1.2f);//length );
                            gl.glCallList(cylinders[mol.atom_color[mol.b[atomp]]]);
                            gl.glPopMatrix();
                        }
                    }
                    if(display_mode == 1 || display_mode == 2){
                    //}else{
                        //alpha carbon tube
                        GLUquadric quadObj = glu.gluNewQuadric();
                        float gry[] = { 0.35f, 0.04f, 0.34f, 1.0f };
                        gl.glMaterialfv( GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE, gry, 0 );
                        for(int ac=0;ac<mol.a_carb_count-1;ac++){
                            float length, x, y, z;
                            
                            if(ac==0){
                                gl.glPushMatrix();
                                 gl.glTranslatef( mol.coord[mol.alpha_carbon[ac]][0], mol.coord[mol.alpha_carbon[ac]][1], mol.coord[mol.alpha_carbon[ac]][2] );
                                 gl.glCallList(tube_sphere);
                                gl.glPopMatrix();
                            }
                            x = mol.coord[mol.alpha_carbon[ac]][0] - mol.coord[mol.alpha_carbon[ac+1]][0];
                            y = mol.coord[mol.alpha_carbon[ac]][1] - mol.coord[mol.alpha_carbon[ac+1]][1];
                            z = mol.coord[mol.alpha_carbon[ac]][2] - mol.coord[mol.alpha_carbon[ac+1]][2];
                            length = (float) Math.sqrt( x * x + y * y + z * z );
                            gl.glPushMatrix();
                             gl.glTranslatef( mol.coord[mol.alpha_carbon[ac+1]][0], mol.coord[mol.alpha_carbon[ac+1]][1], mol.coord[mol.alpha_carbon[ac+1]][2] );
                             gl.glRotatef( (float) Math.acos( z/length ) * 57.2957f, -y, x, 0.0f ); 
                             glu.gluCylinder( quadObj, 0.44f, 0.44f, length, 8, 1 );
                             gl.glCallList(tube_sphere);
                            gl.glPopMatrix();
                        }
                    }
                    gl.glPopMatrix();  
                    //axes
                    gl.glClear(GL.GL_DEPTH_BUFFER_BIT);
                    axes(gl, mol.molRad, width, height, mview);
                    //
                    gl.glPopAttrib();
                    gl.glPopMatrix();  
                    //
                    gl.glFlush();
                    //notifyAll();
                    
                }
                //
	}
        
        private void phiPsiArrows(int i, GL gl, Molecule mol, boolean compObj ){
            float length, x, y, z;
            x = mol.coord[mol.selected_atom_n - 1][0] - mol.coord[i][0];
            y = mol.coord[mol.selected_atom_n - 1][1] - mol.coord[i][1];
            z = mol.coord[mol.selected_atom_n - 1][2] - mol.coord[i][2];
            length = (float) Math.sqrt( x * x + y * y + z * z );
            gl.glPushMatrix();
            gl.glRotatef( (float) Math.acos( z/length ) * 57.2957f, -y, x, 0.0f ); 
            gl.glTranslatef( 0f, 0f, 0.5f * length);
            if(compObj) gl.glCallList(arrow);
            else{
                GLU glu = new GLU();
                GLUquadric quadObj = glu.gluNewQuadric();
                gl.glPassThrough((float)(mol.atom_count+4));
                glu.gluDisk(quadObj,0.42f,0.58f,5,5);
                gl.glRotatef( 180f, 0.0f, 1.0f, 0.0f );	//draw other side!
                glu.gluDisk(quadObj,0.42f,0.58f,5,5);
            }
            gl.glPopMatrix();
            x = mol.coord[mol.selected_atom_c - 1][0] - mol.coord[i][0];
            y = mol.coord[mol.selected_atom_c - 1][1] - mol.coord[i][1];
            z = mol.coord[mol.selected_atom_c - 1][2] - mol.coord[i][2];
            length = (float) Math.sqrt( x * x + y * y + z * z );
            gl.glRotatef( (float) Math.acos( z/length ) * 57.2957f, -y, x, 0.0f ); 
            gl.glTranslatef( 0f, 0f, 0.5f * length);
            gl.glRotatef(180f,0f,1f,0f);
            gl.glRotatef(120f,0f,0f,1f);
            if(compObj) gl.glCallList(arrow);
            else{
                GLU glu = new GLU();
                GLUquadric quadObj = glu.gluNewQuadric();
                gl.glPassThrough((float)(mol.atom_count+5));
                glu.gluDisk(quadObj,0.42f,0.58f,5,5);
                gl.glRotatef( 180f, 0.0f, 1.0f, 0.0f );	//draw other side!
                glu.gluDisk(quadObj,0.42f,0.58f,5,5);
            }
        }

        public void zoom_buttons(GLAutoDrawable drawable, Molecule mol, QuartView mview, String host, 
                                    int port, int host_wait, int display_mode){
                String [] hwait = {".","..","...","....","....."};
                int width, height;

                GL gl = drawable.getGL();
		GLU glu = new GLU();
                width = drawable.getWidth();
                height = drawable.getHeight();
                
                //find outline boxes etc
                findSelection(gl, width, height, mol, mview, display_mode);
                //2D view
                gl.glMatrixMode (GL.GL_PROJECTION);
                gl.glLoadIdentity ();
                gl.glOrtho (0., (double)width, 0., (double)height, -1., 1.);
                gl.glMatrixMode (GL.GL_MODELVIEW);
                gl.glLoadIdentity ();
                //
                gl.glDisable(GL.GL_DEPTH_TEST); 
                gl.glDisable(GL.GL_LIGHTING);
                //
                for(int i=0;i<numButtons;i++){
                    zrect[i][0] = width - 69; zrect[i][1] = width - 3;
                    int indX = i*2;
                    if(zflag[i]) indX++;
                    if(buttonTex[indX] != null){
                        buttonTex[indX].enable();
                        buttonTex[indX].bind();
                        gl.glTexEnvi(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_REPLACE);
                        TextureCoords coords = buttonTex[i*2].getImageTexCoords();
                        //draw textured rectangle
                        gl.glBegin(GL.GL_QUADS);
                        gl.glTexCoord2f(coords.left(), coords.bottom()); gl.glVertex2f(zrect[i][0], height - zrect[i][2]);
                        gl.glTexCoord2f(coords.right(), coords.bottom()); gl.glVertex2f(zrect[i][1], height - zrect[i][2]); 
                        gl.glTexCoord2f(coords.right(), coords.top()); gl.glVertex2f(zrect[i][1], height - zrect[i][3]);
                        gl.glTexCoord2f(coords.left(), coords.top()); gl.glVertex2f(zrect[i][0], height - zrect[i][3]);
                        gl.glEnd();
                        //
                        buttonTex[indX].disable();
                    }
                }
                //Progress bar
                gl.glColor4f(0.9f, 0.9f, 0.9f, 1.0f); 
                if(mol.frames > 0 && mol.data_valid && mol.iterations_done == 0){
                    for(int i=0;i<10;i++){
                        if((mol.frames_done*10)/mol.frames > i){
                            gl.glBegin(GL.GL_QUADS);
                            gl.glVertex2f(10+i*11, height-30); gl.glVertex2f(10+i*11, height-40); 
                            gl.glVertex2f(20+i*11, height-40); gl.glVertex2f(20+i*11, height-30);
                            gl.glEnd();
                        }else{
                            gl.glBegin(GL.GL_LINE_LOOP);
                            gl.glVertex2f(11+i*11, height-31); gl.glVertex2f(11+i*11, height-39); 
                            gl.glVertex2f(19+i*11, height-39); gl.glVertex2f(19+i*11, height-31);
                            gl.glEnd();
                        }
                    }
                }
                //outline box
                if(mol.data_valid || data_available){
                    int textWidth = (int) rendererFlat.getBounds("Simulation [ "+String.valueOf(mol.name, 0 , mol.namelen)+" ]").getWidth();
                    gl.glColor4f(0.9f, 0.9f, 0.9f, 1.0f);
                    gl.glBegin(GL.GL_LINE_LOOP);
                        gl.glVertex2i(4, height-4);
                        gl.glVertex2i(4, height-44);
                        gl.glVertex2i(14+textWidth, height-44);
                        gl.glVertex2i(14+textWidth, height-4);
                    gl.glEnd();
                }
                //Box around selected item
                if(mol.selected_atom > 0 && display_mode != 1){
                    gl.glColor4f(1.0f, 0.87f, 0.0f, 1.0f);
                    gl.glBegin(GL.GL_LINE_LOOP);
			gl.glVertex2i(boundingRect[0][2], boundingRect[0][0]);
			gl.glVertex2i(boundingRect[0][2], boundingRect[0][1]);
			gl.glVertex2i(boundingRect[0][3], boundingRect[0][1]);
			gl.glVertex2i(boundingRect[0][3], boundingRect[0][0]);
                    gl.glEnd();
                }
                //
                gl.glEnable(GL.GL_LIGHTING);
                //
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glLoadIdentity();
		glu.gluPerspective(45.0f, (float)width/(float)height, 0.1, rear_clip);
		gl.glMatrixMode(GL.GL_MODELVIEW);
		gl.glLoadIdentity();
                gl.glEnable(GL.GL_DEPTH_TEST);
                //text
                rendererFlat.setColor(0.7f, 0.7f, 0.7f, 1.0f);
                rendererFlat.beginRendering(width, height);
                if(mol.data_valid || data_available){
                                //rendererFlat.draw("Connected to ["+host+":"+
                                //        String.valueOf(port)+"] ", 
                                //            10, height-20);
                    rendererFlat.draw("Simulation [ "+String.valueOf(mol.name, 0 , mol.namelen)+" ]", 
                                            10, height-20);
                    //X,Y,Z
                    //int offs = (int) rendererFlat.getBounds("X").getWidth()/2;
                    rendererFlat.setColor(0.2f, 0.2f, 0.6f, 1.0f);
                    rendererFlat.draw("X",boundingRect[1][3]+2, boundingRect[1][0]+2);  
                    rendererFlat.setColor(0.2f, 0.6f, 0.2f, 1.0f);
                    rendererFlat.draw("Y",boundingRect[2][3]+2, boundingRect[2][0]+2);  
                    rendererFlat.setColor(0.6f, 0.2f, 0.2f, 1.0f);
                    rendererFlat.draw("Z",boundingRect[3][3]+2, boundingRect[3][0]+2);
                    //Temperature/energy
                    //if(data_available){
                    //    rendererFlat.draw("Temperature : "+String.valueOf((int)mol.temperature)+"\u00B0K", 
                    //                        10, height-40);
                    //}
                }else{
                                rendererFlat.draw("Waiting for ["+host+":"+
                                        String.valueOf(port)+"] "+hwait[host_wait],
                                            10, height-20);
                }
                //current mode
                if(mol.iterations_done != 0){
                    rendererFlat.setColor(0.7f, 0.7f, 0.7f, 1.0f);
                    rendererFlat.draw("Mode number: "+
                                        String.valueOf(mol.iterations_done)+".", 
                                            10, height-40);
                }
                //selected atom
                if(mol.selected_atom > 0 && display_mode != 1){
                    rendererFlat.setColor(1.0f, 0.87f, 0.0f, 1.0f);
                    int l;
                    for(l=0;l<4;l++) 
                        if(mol.type[mol.selected_atom-1][l] == 0) break;
                    rendererFlat.draw("Atom    : "+
                                        String.valueOf(mol.type[mol.selected_atom-1], 0, l)+" ("+String.valueOf(mol.selected_atom)+")", 
                                            boundingRect[0][3]+2, (boundingRect[0][0] + boundingRect[0][1])/2+2);  
                    rendererFlat.draw("Charge: "+
                                        String.valueOf(mol.charge[mol.selected_atom-1]), 
                                            boundingRect[0][3]+2, (boundingRect[0][0] + boundingRect[0][1])/2-14);  
                    if(mol.selected_atom_n > 0 && display_mode < 3){
                        rendererFlat.setColor(1.0f, 0.0f, 0.0f, 1.0f);
                        //Calculate dihedral
                        if(mol.selected_atom_c_n > 0 && mol.selected_atom_n_c > 0){
                            //System.out.println("Dihs "+String.valueOf(mol.selected_atom_n_c)+":"+String.valueOf(mol.selected_atom_n)+":"
                            //        +String.valueOf(mol.selected_atom)+":"
                            //        +String.valueOf(mol.selected_atom_c)+":"
                            //        +String.valueOf(mol.selected_atom_c_n)+".");
                            int Phi = calcDihedral(mol.selected_atom_n_c-1,mol.selected_atom_n-1,mol.selected_atom-1,mol.selected_atom_c-1,mol);
                            int Psi = calcDihedral(mol.selected_atom_n-1,mol.selected_atom-1,mol.selected_atom_c-1,mol.selected_atom_c_n-1,mol);
                            //
                            int textWidth = (int) rendererFlat.getBounds("\u03A6:"+String.valueOf(Phi)+"\u00B0").getWidth();
                            rendererFlat.draw("\u03A6:"+String.valueOf(Phi)+"\u00B0",(boundingRect[4][2]+boundingRect[4][3])/2-textWidth/2, boundingRect[4][1]+2);  
                            textWidth = (int) rendererFlat.getBounds("\u03A8"+String.valueOf(Psi)+"\u00B0").getWidth();
                            rendererFlat.draw("\u03A8:"+String.valueOf(Psi)+"\u00B0",(boundingRect[5][2]+boundingRect[5][3])/2-textWidth/2, boundingRect[5][1]+2);  
                        }else{
                            int textWidth = (int) rendererFlat.getBounds("\u03A6").getWidth();
                            rendererFlat.draw("\u03A6",(boundingRect[4][2]+boundingRect[4][3])/2-textWidth/2, boundingRect[4][1]+2);  
                            rendererFlat.draw("\u03A8",(boundingRect[5][2]+boundingRect[5][3])/2-textWidth/2, boundingRect[5][1]+2);  
                        }
                    }
                }
                //Buttons
                //String[] Btext = {"Zoom+","Zoom-","Capt."};
                //for(int i=0;i<3;i++){
                //    if(zflag[i]) rendererFlat.setColor(0.75f, 0.0f, 0.0f, 0.8f);
                //    else rendererFlat.setColor(0.6f, 0.6f, 0.6f, 0.8f);
                //    rendererFlat.draw(Btext[i], width - 60, height-zrect[i][2]+5); 
                //}
                rendererFlat.endRendering();

        }
        
        private int calcDihedral(int a1, int a2, int a3, int a4, Molecule mol){
            float [] a = {0.0f, 0.0f, 0.0f}, b = {0.0f, 0.0f, 0.0f}, c = {0.0f, 0.0f, 0.0f};
            float [] r12 = {0.0f, 0.0f, 0.0f}, r23 = {0.0f, 0.0f, 0.0f}, r34 = {0.0f, 0.0f, 0.0f};
            //find differences
            for(int i=0;i<3;i++){
                r12[i] = mol.coord[a2][i] - mol.coord[a1][i]; 
                r23[i] = mol.coord[a3][i] - mol.coord[a2][i]; 
                r34[i] = mol.coord[a4][i] - mol.coord[a3][i];
            }
            //find cross products
            a = cross3D(r12,r23);
            b = cross3D(r23,r34);
            c = cross3D(r23,a);
            //normalize
            float norma = norm3D(a);
            float normb = norm3D(b);
            float normc = norm3D(c);
            if(norma == 0.0f || normb == 0.0f || normc == 0.0f) return(0);
            //dot product
            float cosPhi = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
            float sinPhi = c[0]*b[0]+c[1]*b[1]+c[2]*b[2];
            return (int)Math.ceil(Math.atan2(sinPhi,cosPhi)*57.2957);
        }
        
        private float [] cross3D(float [] a, float [] b){
            float [] c = {0.0f, 0.0f, 0.0f};
            //find cross products
            c[0] = a[1]*b[2]-b[1]*a[2];
            c[1] = a[2]*b[0]-b[2]*a[0];
            c[2] = a[0]*b[1]-b[0]*a[1];
            return c;
        }
        
        private float norm3D(float [] a){
            float norm = (float)Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
            if(norm > 0){
                for(int i=0;i<3;i++) a[i] /= norm;
            }
            return norm;
        }
        
	public void axes(GL gl, float radius, int width, int height, QuartView mview)
	{
            float width_ratio, height_ratio;
                      
            if(width > height){
                width_ratio = (float)width / (float)height;
                height_ratio = 1.0f;
            }else{
                width_ratio = 1.0f;
                height_ratio = (float)height / (float)width;
            }
            gl.glPushMatrix();
            gl.glTranslatef( -radius*width_ratio, -radius*height_ratio, 0.0f  );  
            gl.glScalef(radius/5f,radius/5f,radius/5f);
            gl.glRotatef(mview.rot_a, mview.rot_v[0], mview.rot_v[1], mview.rot_v[2]);
            gl.glPushMatrix();
            //
            gl.glCallList(axes1);
            //
            gl.glPopMatrix();  
            gl.glPopMatrix(); 
        }
        
     private void getSelection(GL gl, int xPos, int yPos, int width, int height, Molecule mol, QuartView mview, int display_mode)
     {
        GLU glu = new GLU();
        //
        gl.glPushAttrib(GL.GL_ALL_ATTRIB_BITS );//GL.GL_CURRENT_BIT | GL.GL_LIGHTING_BIT);
        // Space for selection buffer
        //final int BSIZE = mol.atom_count * 4;
        ByteBuffer byteBuffer = ByteBuffer.allocateDirect(mol.atom_count * 4 * 4); //4 ints per atom/4 bytes per int
        byteBuffer.order(ByteOrder.nativeOrder());
        ByteBuffer byteBuffer2 = ByteBuffer.allocateDirect(4*4);
        byteBuffer2.order(ByteOrder.nativeOrder());
        //  
        IntBuffer selectBuff = byteBuffer.asIntBuffer();//
        //IntBuffer.allocate(BSIZE);
        IntBuffer viewport = byteBuffer2.asIntBuffer();//IntBuffer.allocate(4);
	// Hit counter and viewport storage
	int hits;
	// Setup selection buffer
	gl.glSelectBuffer(selectBuff.capacity(), selectBuff);
	// Get the viewport
	gl.glGetIntegerv(GL.GL_VIEWPORT, viewport);
	// Switch to projection and save the matrix
	gl.glMatrixMode(GL.GL_PROJECTION);
	gl.glPushMatrix();
	// Change render mode
	gl.glRenderMode(GL.GL_SELECT);
	// Establish new clipping volume to be unit cube around
	// mouse cursor point (xPos, yPos) and extending one pixels
	// in the vertical and horizontal direction
	gl.glLoadIdentity();
        int curSize = 2;
        //if spacefill then large cursor
        if(display_mode == 3) curSize = 20;
	glu.gluPickMatrix(xPos, viewport.get(3) - yPos + viewport.get(1), curSize,curSize, viewport);
	// Apply perspective matrix 
        glu.gluPerspective(45.0f, (float)width/(float)height, 0.1, rear_clip);
	// Draw the scene
        // Clear the drawing area
        gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
        gl.glPushMatrix();
        //scale for full view of molecule
        gl.glTranslatef( 0.0f, 0.0f, -mol.molRad * 3.0f  );
        //set rotation
        gl.glRotatef(mview.rot_a, mview.rot_v[0], mview.rot_v[1], mview.rot_v[2]);
        // scale
        gl.glScalef(mview.scale_xyz,mview.scale_xyz,mview.scale_xyz);
        //
        // Initialize the names stack
	gl.glInitNames();
	gl.glPushName(0);
        //
        for ( int i = 0 ; i < mol.atom_count ; i++ ){
            gl.glPushMatrix();
            gl.glTranslatef( mol.coord[i][0], mol.coord[i][1], mol.coord[i][2] );
            gl.glLoadName(i+1);
            gl.glCallList(spheres[mol.atom_color[i]]);
            gl.glPopMatrix();
        }
        gl.glPopMatrix();  
        //
        //gl.glFlush();
	// Collect the hits
	hits = gl.glRenderMode(GL.GL_RENDER);
	// Restore the projection matrix
	gl.glMatrixMode(GL.GL_PROJECTION);
	gl.glPopMatrix();
	// Go back to modelview for normal rendering
	gl.glMatrixMode(GL.GL_MODELVIEW);
	// If a hit occurred, display the info.
	if(hits > 0){
            //find front atom
            int minz = Integer.MAX_VALUE;//10000;
            int selAt=0, numNames, k=0;
            for(int j=0;j<hits;j++){
                numNames = selectBuff.get(k);   //should allways be one here
                if(selectBuff.get(k+1) < minz){
                    minz = selectBuff.get(k+1);
                    selAt = selectBuff.get(k+3);
                }
                k += numNames + 3;
            }
            //select or de-select
            if(mol.selected_atom == selAt){//selectBuff.get(3))
		mol.selected_atom = 0;
            }else{
		mol.selected_atom = selAt;
                //if alpha carbon find N and C atoms for \Phi, \Psi
                mol.selected_atom_n = mol.selected_atom_c = 0;
                if(mol.type[mol.selected_atom-1][0] == 'C' && mol.type[mol.selected_atom-1][1] == 'A'){    //alpha carbon?
                    //find nearest N and C for Phi/Psi angle position
                    for(int j=0;j<mol.bond_count;j++){
                        int tempAb=-1;
                        if(mol.a[j] == mol.selected_atom-1){
                            tempAb = mol.b[j];
                        }else{
                            if(mol.b[j] == mol.selected_atom-1) tempAb = mol.a[j];
                        }
                        if(tempAb>=0){
                            if(mol.type[tempAb][0]=='N') mol.selected_atom_n = tempAb+1;
                            if(mol.type[tempAb][0]=='C' && mol.type[tempAb][1]==0) mol.selected_atom_c = tempAb+1;
                        }
                    }
                    if(mol.selected_atom_n == 0 || mol.selected_atom_c == 0){
                        mol.selected_atom_n = mol.selected_atom_c = 0;
                    }else{
                        //find additional N-C and C-N bonded atoms for dihedral angle
                        int tempAb1, tempAb2;
                        mol.selected_atom_n_c = mol.selected_atom_c_n = 0;
                        for(int j=0;j<mol.bond_count;j++){
                            tempAb1 = tempAb2 = 0;
                            if(mol.a[j] == mol.selected_atom_n-1) tempAb1 = mol.b[j];
                            if(mol.b[j] == mol.selected_atom_n-1) tempAb1 = mol.a[j];
                            if(mol.a[j] == mol.selected_atom_c-1) tempAb2 = mol.b[j];
                            if(mol.b[j] == mol.selected_atom_c-1) tempAb2 = mol.a[j];
                            if(tempAb1 > 0 && mol.type[tempAb1][0]=='C' && tempAb1 != mol.selected_atom-1)
                                mol.selected_atom_n_c = tempAb1+1;
                            if(tempAb2 > 0 && mol.type[tempAb2][0]=='N' && tempAb2 != mol.selected_atom-1)
                                mol.selected_atom_c_n = tempAb2+1;
                        }
                    }
                }
                //
            }
        }else{
            mol.selected_atom = 0;
        }
        gl.glPopAttrib();    //get attributes back 
     }
        
     private void findSelection(GL gl, int width, int height, Molecule mol, QuartView mview, int display_mode){
        int size, i;//, count;
        float width_ratio, height_ratio;
        
        gl.glPushAttrib(GL.GL_ALL_ATTRIB_BITS );//GL.GL_CURRENT_BIT | GL.GL_LIGHTING_BIT);
        // Initial minimum and maximum values
        for(i=0;i<6;i++){
            boundingRect[i][3] = boundingRect[i][1] = -999999;
            boundingRect[i][2] = boundingRect[i][0] =  999999;
        }
        GLU glu = new GLU();
        GLUquadric quadObj = glu.gluNewQuadric();
        glu.gluQuadricNormals(quadObj, GLU.GLU_SMOOTH);
	glu.gluDeleteQuadric(quadObj);

        // Space for the feedback buffer
        ByteBuffer byteBuffer = ByteBuffer.allocateDirect(256 * 4 * 6); //<100 for sphere resolution 4, 4 spheres
        byteBuffer.order(ByteOrder.nativeOrder());
        //  
        FloatBuffer feedBackBuff = byteBuffer.asFloatBuffer();
	// Set the feedback buffer
	gl.glFeedbackBuffer(feedBackBuff.capacity(),GL.GL_2D, feedBackBuff);
        // Redraw the scene
	// Enter feedback mode
	gl.glRenderMode(GL.GL_FEEDBACK);
        //
        gl.glPushMatrix();
        //*********do draw**************
        glu.gluPerspective(45.0f, (float)width/(float)height, 0.1, rear_clip);
        gl.glLoadIdentity();
        if(mol.selected_atom > 0){
            //scale for full view of molecule
            gl.glPushMatrix();
            //Draw selected atom
            gl.glTranslatef( 0.0f, 0.0f, -mol.molRad * 3.0f  );
            //set rotation
            gl.glRotatef(mview.rot_a, mview.rot_v[0], mview.rot_v[1], mview.rot_v[2]);
            // scale
            gl.glScalef(mview.scale_xyz,mview.scale_xyz,mview.scale_xyz);
            //
            gl.glPushMatrix();
                gl.glTranslatef( mol.coord[mol.selected_atom-1][0], mol.coord[mol.selected_atom-1][1], mol.coord[mol.selected_atom-1][2] );
                gl.glPassThrough((float)mol.selected_atom);
                float curSiz = a_rad[mol.atom_color[mol.selected_atom-1]] * 0.5375f / 1.8f; //10% larger
                if(display_mode == 3) curSiz *= 4.0f;
                glu.gluSphere( quadObj, curSiz, 4, 4 );
                //arrows
                if(mol.selected_atom_n > 0 && display_mode < 3)
                    phiPsiArrows(mol.selected_atom - 1, gl,mol,false);
                //
            gl.glPopMatrix();
            //
            gl.glPopMatrix();
        }
        //draw X,Y,Z
        if(width > height){
                width_ratio = (float)width / (float)height;
                height_ratio = 1.0f;
        }else{
                width_ratio = 1.0f;
                height_ratio = (float)height / (float)width;
        }
        gl.glTranslatef( -mol.molRad*width_ratio, -mol.molRad*height_ratio, -mol.molRad * 3.0f );  
        gl.glScalef(mol.molRad/5f,mol.molRad/5f,mol.molRad/5f);
        gl.glRotatef(mview.rot_a, mview.rot_v[0], mview.rot_v[1], mview.rot_v[2]);
        gl.glPushMatrix();
        gl.glTranslatef( 0.0f, 0.0f, 1.0f );    
        gl.glPassThrough((float)(mol.atom_count+3));
        glu.gluSphere( quadObj, 0.05, 4, 4 );
        gl.glPopMatrix();
        //
        gl.glRotatef(-90f, 1.0f, 0.0f, 0.0f); 
        gl.glPushMatrix();
        gl.glTranslatef( 0.0f, 0.0f, 1.0f );    
        gl.glPassThrough((float)(mol.atom_count+2));
        glu.gluSphere( quadObj, 0.05, 4, 4 );
        gl.glPopMatrix();
        //
        gl.glRotatef(90f, 0.0f, 1.0f, 0.0f);
        gl.glPushMatrix();
        gl.glTranslatef( 0.0f, 0.0f, 1.0f );    
        gl.glPassThrough((float)(mol.atom_count+1));
        glu.gluSphere( quadObj, 0.05, 4, 4 );
        gl.glPopMatrix();
        //*********end draw**************
        gl.glPopMatrix();
        //
	// Leave feedback mode
	size = gl.glRenderMode(GL.GL_RENDER);
        //System.out.println("Size: "+String.valueOf(size));
	// Parse the feedback buffer and get the
	// min and max X and Y window coordinates
	i = 0;
	while(i < size){
            // Search for appropriate token
            if(feedBackBuff.get(i++) == GL.GL_PASS_THROUGH_TOKEN){
                if(feedBackBuff.get(i) >= (mol.atom_count+1) && feedBackBuff.get(i) <= (mol.atom_count+5)){
                    i++;
                    findpolygons((int)feedBackBuff.get(i-1)-mol.atom_count, i, size, feedBackBuff);
                }    
                if(feedBackBuff.get(i) == mol.selected_atom){     //selected molecule
                    i++;
                    findpolygons(0, i, size, feedBackBuff);
                }
            }
        }
        gl.glPopAttrib();    //get attributes back 
     }

     private void findpolygons(int rInd, int i, int size, FloatBuffer feedBackBuff){
        int count;
                    
        // Loop until next token is reached
        while(i < size && feedBackBuff.get(i) != GL.GL_PASS_THROUGH_TOKEN){
            // Just get the polygons
            if(feedBackBuff.get(i) == GL.GL_POLYGON_TOKEN){
                // Get all the values for this polygon
                count = (int)feedBackBuff.get(++i); // How many vertices
                i++;
                for(int j=0; j<count; j++){	// Loop for each vertex
                    // Min and Max X
                    if((int)feedBackBuff.get(i) > boundingRect[rInd][3])
                        boundingRect[rInd][3] = (int)feedBackBuff.get(i);
                    if((int)feedBackBuff.get(i) < boundingRect[rInd][2])
                        boundingRect[rInd][2] = (int)feedBackBuff.get(i);
                    i++;
                    // Min and Max Y
                    if((int)feedBackBuff.get(i) > boundingRect[rInd][1])
                        boundingRect[rInd][1] = (int)feedBackBuff.get(i);
                    if((int)feedBackBuff.get(i) < boundingRect[rInd][0])
                        boundingRect[rInd][0] = (int)feedBackBuff.get(i);
                    i++;
                }
            }else{
                i++;	// Get next index and keep looking
            }
        }
     }

}
