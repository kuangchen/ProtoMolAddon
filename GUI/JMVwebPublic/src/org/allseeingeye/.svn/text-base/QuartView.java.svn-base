package org.allseeingeye;

/* QuartView contains the molecular viewer 'view' control.
 * 
 * Uses quarternions for the product of successive 
 * rotations about axes.
 * @authors Chris Sweet & James Sweet. Spring 2008
 */

public class QuartView {
    
          //rotation angle and vector
          public float [] rot_v = {0.57735f,0.57735f,0.57735f};
          public float rot_a = 0.0f;
          //current values
          public float [] current_rot_v = {0.57735f,0.57735f,0.57735f};
          public float current_rot_a = 0.0f;
          //old values
          private float [] old_rot_v = {0.57735f,0.57735f,0.57735f};
          private float old_rot_a = 0.0f;
          private float sphere_rad, halfWidth, halfHeight;
          private float [] start_v, current_v;
          //convertion
          float rad2deg = 180.0f / (float)Math.PI;
          float deg2rad = (float)Math.PI / 180.0f;
          //rotation product variables
          float [] start_q, current_q, new_q;
          //rate of change
          public float last_angle = 0.0f, angle_inc = 0.0f;
          //scale
          public float scale_xyz = 1.0f, old_scale_xyz = 0.0f;
          
          //initialize
          public QuartView(int w, int h){
              halfWidth = (float)w / 2.0f; halfHeight = (float)h / 2.0f;
              sphere_rad = Math.min(halfWidth, halfHeight);
              if(sphere_rad <= 0.0f) sphere_rad = 1.0f;
              //vectors
              start_v = new float[3];
              current_v = new float[3];
              //quartonians
              start_q = new float[4];
              current_q = new float[4];
              new_q = new float[4];
          }
          
          public void set_current_v(float x, float y){
              float cosRAng, sinRAng, normm1;
              
              current_v = find_z(x,y); //get normalized current vector
              //cross product for rotation axis
              current_rot_v = cross_v(start_v,current_v);
              //find rotation angle (resolve full 2PI radians)
              cosRAng = dot_v(start_v,current_v);
              sinRAng = norm_v(current_rot_v );
              current_rot_a = (float)Math.atan2(sinRAng,cosRAng) * rad2deg;
              //unit rotation vector
              if(Math.abs(sinRAng) < 1e-8) return;
              normm1 = 1.0f / sinRAng;
              current_rot_v = scale_v(current_rot_v,normm1);
              //find product of old and new rotation
              start_q = rotation2Quartonian( old_rot_a, old_rot_v);
              current_q = rotation2Quartonian( current_rot_a, current_rot_v);
              new_q = quartonianProduct(start_q, current_q);
              //put into current rotation
              rot_a = quartonian2Rotation(rot_a, rot_v, new_q);  
              
          }

          public void set_current_rot(){
              
              //find product of old and new rotation
              start_q = rotation2Quartonian( rot_a, rot_v);
              current_q = rotation2Quartonian( angle_inc/2.0f, current_rot_v);
              new_q = quartonianProduct(start_q, current_q);
              //put into current rotation
              rot_a = quartonian2Rotation(rot_a, rot_v, new_q);  
          }
          
          //set starting vector 
          public void set_start_v(float x, float y){
              //get vector on sphere
              start_v = find_z(x,y);
              //save old values
              copy_v(old_rot_v,rot_v);
              old_rot_a = rot_a;
          }
          
          //convert rotation 'theta' about vetctor 'v' to quarternion
          private float [] rotation2Quartonian(float theta, float [] v){
              float halfAngle, sinHalfAng;
              float [] op_q={0.0f,0.0f,0.0f,0.0f};
              
              halfAngle = theta * deg2rad * 0.5f;
              sinHalfAng = (float)Math.sin(halfAngle);
              op_q[0] = (float)Math.cos(halfAngle); //cosine of angle
              op_q[1] = sinHalfAng * v[0];          //quarternion must have norm 1
              op_q[2] = sinHalfAng * v[1];
              op_q[3] = sinHalfAng * v[2];
              return op_q;
          }
          
          //convert quarternion to rotation 'theta' about vetctor 'v'
          private float quartonian2Rotation(float theta, float [] v, float [] q){
              float v_norm, theta_t;
              
              theta_t = (float)Math.acos(q[0]) * rad2deg * 2.0f;
              v[0] = q[1]; v[1] = q[2]; v[2] = q[3];
              v_norm = norm_v(v);
              if(Math.abs(v_norm) < 1e-8) return 0.0f;
              v[0] /= v_norm; v[1] /= v_norm; v[2] /= v_norm; 
              return theta_t;
          }
          
          private float [] quartonianProduct(float [] q1, float [] q2){
              float [] q = {0.0f,0.0f,0.0f,0.0f};
              q[0] = q2[0]*q1[0] - q2[1]*q1[1] - q2[2]*q1[2] - q2[3]*q1[3];
              q[1] = q2[2]*q1[3] - q2[3]*q1[2] + q2[0]*q1[1] + q2[1]*q1[0];
              q[2] = q2[3]*q1[1] - q2[1]*q1[3] + q2[0]*q1[2] + q2[2]*q1[0];
              q[3] = q2[1]*q1[2] - q2[2]*q1[1] + q2[0]*q1[3] + q2[3]*q1[0];
              return q;
          }
          
          //find z coordinate on sphere, return unit vector
          private float [] find_z(float x, float y){
              float x2py2, fact;
              float [] v = {0.0f,0.0f,0.0f};
              
              v[0] = (x - halfWidth) / halfWidth;//sphere_rad ; 
              v[1] = -(y - halfHeight) / halfHeight;//sphere_rad;
              x2py2 = v[0]*v[0]+v[1]*v[1];
              if(x2py2 > 1.0){
                  v[2] = 0.0f;
                  fact = (float)Math.sqrt(1.0/x2py2);
                  v[0] *= fact;
                  v[1] *= fact;
              }else{
                  v[2] = (float)Math.sqrt((double)(1.0 - x2py2));
              }
              return v;
          }
          
          //set shere diamiter if view changes
          public void sphereRad(int w, int h){
              halfWidth = (float)w / 2.0f; halfHeight = (float)h / 2.0f;
              sphere_rad = Math.min(halfWidth, halfHeight);
              if(sphere_rad <= 0.0f) sphere_rad = 1.0f;
          }
          
          private void copy_v(float [] a, float [] b){
              for(int i=0;i<3;i++) a[i] = b[i];
          }
          
          private float dot_v(float [] a, float [] b){
              return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
          }
          
          private float norm_v(float [] a){
              return (float)Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
          }
          
          private float [] cross_v(float [] a, float [] b){
              float [] op_v={0.0f,0.0f,0.0f};
              op_v[0] = a[1]*b[2]-a[2]*b[1];
              op_v[1] = a[2]*b[0]-a[0]*b[2];
              op_v[2] = a[0]*b[1] - a[1]*b[0];
              return op_v;
          }
          
          private float [] scale_v(float [] a, float b){
              float [] op_v={0.0f,0.0f,0.0f};
              op_v[0] = a[0]*b; op_v[1] = a[1]*b; op_v[2] = a[2]*b;
              return op_v;
          }
          
          //scale factor utilities
          public void scaleMult(float sf){
              scale_xyz *= sf;
              if(scale_xyz < 0.05f) scale_xyz = 0.05f;
          }
          
          public void scaleSum(float sm){
              scale_xyz = sm + old_scale_xyz;
              if(scale_xyz < 0.05f) scale_xyz = 0.05f;
          }
          
          public void scaleSave(){
              old_scale_xyz = scale_xyz;
          }
          
          //reset view
          public void resetView(){ 
            rot_a = 0.0f; 
            rot_v[0] = 1.0f;
            rot_v[1] = 1.0f; 
            rot_v[2] = 1.0f;
            scale_xyz = 1.0f;
            old_scale_xyz = 0.0f;
            angle_inc = 0.0f;
          }

}
