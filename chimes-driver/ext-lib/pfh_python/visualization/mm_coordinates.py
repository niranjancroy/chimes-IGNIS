import numpy as np
import gc
from mm_utilities import *

 
#
# given a set of positions, project them to the coordinates they 'would be at'
#   with respect to a given camera (for ray-tracing purposes)
#
def coordinates_project_to_camera(x, y, z, \
        camera_pos=[0.,0.,0.], camera_dir=[0.,0.,1.], \
        screen_distance=1.0 ):
    camera_pos=np.array(camera_pos,dtype='d');
    camera_dir=np.array(camera_dir,dtype='d');
    ## recenter at camera location
    x-=camera_pos[0]; y-=camera_pos[1]; z-=camera_pos[2];
    ## line-of-sight distance along camera viewing angle
    znew = x*camera_dir[0] + y*camera_dir[1] + z*camera_dir[2]; 
    ## get perpendicular axes to determine position in 'x-y' space:
    cperp_x,cperp_y = return_perp_vectors(camera_dir);
    xnew = x*cperp_x[0] + y*cperp_x[1] + z*cperp_x[2];
    ynew = x*cperp_y[0] + y*cperp_y[1] + z*cperp_y[2];

    ## project to 'screen' at fixed distance along z-axis :: 
    znew[znew==0.]=-1.0e-10; ## just to prevent nans
    ## need correction factors for other quantities (like hsml, for example)
    r_corr = screen_distance/znew;
    xnew *= r_corr; ## put positions 'at the screen'
    ynew *= r_corr; ## put positions 'at the screen'
    znew = -znew; ## so values in the direction of camera are *negative* (for later sorting)

    gc.collect()
    return xnew, ynew, znew, r_corr


def coordinates_rotate(x_all, y_all, z_all, theta, phi, coordinates_cylindrical=0):
    ## set viewing angle for image
    ## project into plane defined by vectors perpendicular to the vector of this angle
    x = np.cos(phi)*x_all + np.sin(phi)*y_all + 0.*z_all
    y = -np.cos(theta)*np.sin(phi)*x_all + np.cos(theta)*np.cos(phi)*y_all + np.sin(theta)*z_all
    z =  np.sin(theta)*np.sin(phi)*x_all - np.sin(theta)*np.cos(phi)*y_all + np.cos(theta)*z_all
    if (coordinates_cylindrical==1):
        x=x/abs(x)*np.sqrt(x*x+z*z); y=y; z=z
    gc.collect()
    return x,y,z
  

## simple function to wrap x/y/z coordinates near the box edges, to give 
##   the 'correct' x-x0 in periodic coordinates
def periodic_dx( x, x0, box ):
    b0 = box / 2.
    dx = x - x0
    too_large = (dx > b0)
    dx[too_large] -= box
    too_small = (dx < -b0)
    dx[too_small] += box
    gc.collect()
    return dx
    
def get_closest_periodic_value( x_expected, x_final, x_initial, limiter=1.0 ):
    x_remainder = np.mod(x_expected, 2.*np.pi)
    dx = periodic_dx( x_final, x_remainder, 2.*np.pi ) ## returns wrapped dx from x_remainder
    dx_full = x_expected + dx - x_initial
    
    return dx_full + x_initial
    
    # now put a limiter so we don't allow some crazy number of 'laps':
    sign_dx = 0.*dx_full + 1.
    sign_dx[dx_full<0.]=-1.
    dx_full_abs = (sign_dx*dx_full) 
    dx_full_cycles = dx_full_abs / (2.*np.pi) # number of 'laps'
    dx_full_remainder = np.mod(dx_full_cycles, 1.0)
    
    dx_corrected = np.copy(dx_full)
    lminusone = limiter - 1.
    hi = (dx_full_cycles > limiter) & (dx_full_remainder > lminusone)
    dx_corrected[hi] = 2.*np.pi*sign_dx[hi]*dx_full_remainder[hi]
    hi = (dx_full_cycles > limiter) & (dx_full_remainder < lminusone)
    dx_corrected[hi] = 2.*np.pi*sign_dx[hi]*(1.+dx_full_remainder[hi])
    
    gc.collect()
    return dx_corrected + x_initial


def xyz_to_cylindrical(x,y,z,vx,vy,vz):
    for q in [x,y,z,vx,vy,vz]: q=np.array(q,dtype='d')+1.0e-10

    ## get the angular momentum vector (to use for interpolation)
    r = np.sqrt( x*x + y*y + z*z )
    pos = np.zeros([x.size,3],dtype='d'); 
    pos[:,0]=x/r; pos[:,1]=y/r; pos[:,2]=z/r;
    v = np.sqrt( vx*vx + vy*vy + vz*vz )
    vel = np.zeros([vx.size,3],dtype='d'); 
    vel[:,0]=vx/v; vel[:,1]=vy/v; vel[:,2]=vz/v;
    
    j_hat = cross( pos, vel )
    jj = np.sqrt( j_hat[:,0]*j_hat[:,0] + j_hat[:,1]*j_hat[:,1] + j_hat[:,2]*j_hat[:,2] )
    for j in [0,1,2]: j_hat[:,j] /= jj
    
    s = np.argsort(r)
    r_s = r[s]; 
    j_s = np.zeros((r.size,3))
    j_x = np.cumsum(j_hat[s,0] * jj[s])
    j_y = np.cumsum(j_hat[s,1] * jj[s])
    j_z = np.cumsum(j_hat[s,2] * jj[s])
    j_0 = np.sqrt(j_x*j_x + j_y*j_y + j_z*j_z)
    j_x /= j_0; j_y /= j_0; j_z /= j_0;
    j_hat[s,0] = j_x
    j_hat[s,1] = j_y
    j_hat[s,2] = j_z

    ## define the absolute 'phi' relative to an arbitrary (but fixed) 90-degree rotation of the angular momentum vector j
    x_jhat = 0.*j_hat
    x_jhat[:,0]=0.*j_hat[:,0]; x_jhat[:,1]=-j_hat[:,2]; x_jhat[:,2]=j_hat[:,1]
    jj = np.sqrt( x_jhat[:,0]*x_jhat[:,0] + x_jhat[:,1]*x_jhat[:,1] + x_jhat[:,2]*x_jhat[:,2] )
    for j in [0,1,2]: x_jhat[:,j] /= jj
    ## generate y-vector by cross-product a x b = c (gaurantees right-hand rule)
    y_jhat = cross( j_hat, x_jhat )
    jj = np.sqrt( y_jhat[:,0]*y_jhat[:,0] + y_jhat[:,1]*y_jhat[:,1] + y_jhat[:,2]*y_jhat[:,2] )
    for j in [0,1,2]: y_jhat[:,j] /= jj
    
    z_vec = j_hat
    x_vec = x_jhat
    y_vec = y_jhat
    
    xn = x*x_vec[:,0] + y*x_vec[:,1] + z*x_vec[:,2]
    yn = x*y_vec[:,0] + y*y_vec[:,1] + z*y_vec[:,2]
    zn = x*z_vec[:,0] + y*z_vec[:,1] + z*z_vec[:,2]
    vxn = vx*x_vec[:,0] + vy*x_vec[:,1] + vz*x_vec[:,2]
    vyn = vx*y_vec[:,0] + vy*y_vec[:,1] + vz*y_vec[:,2]
    vzn = vx*z_vec[:,0] + vy*z_vec[:,1] + vz*z_vec[:,2]
    
    R = np.sqrt(xn*xn + yn*yn)
    V_R = (vxn * xn + vyn * yn) / R
    V_Phi = (-vxn * yn + vyn * xn) / R
    Phi = np.arctan2( yn , xn )
    ## reset to a 0-2pi system instead of numpy's -pi,pi system
    lo = (Phi < 0.); Phi[lo] += 2.*np.pi

    gc.collect()
    return R, V_R, zn, vzn, Phi, V_Phi, j_hat


    
def xyz_to_polar(x,y,z,vx,vy,vz):
    for q in [x,y,z,vx,vy,vz]: q=np.array(q,dtype='d')+1.0e-10
    x += 1.0e-5; vy += 1.0e-5;

    ## get the angular momentum vector (to use for interpolation)
    r = np.sqrt( x*x + y*y + z*z )
    pos = np.zeros([x.size,3],dtype='d'); 
    pos[:,0]=x/r; pos[:,1]=y/r; pos[:,2]=z/r;
    v = np.sqrt( vx*vx + vy*vy + vz*vz )
    vel = np.zeros([vx.size,3],dtype='d'); 
    vel[:,0]=vx/v; vel[:,1]=vy/v; vel[:,2]=vz/v;
    
    j_hat = cross( pos, vel )
    jj = np.sqrt( j_hat[:,0]*j_hat[:,0] + j_hat[:,1]*j_hat[:,1] + j_hat[:,2]*j_hat[:,2] )
    for j in [0,1,2]: j_hat[:,j] /= jj
    
    ## get spherical polar coordinates
    v_r = vx*pos[:,0] + vy*pos[:,1] + vz*pos[:,2]

    ## now get the vector perpendicular to both j and r
    phi_hat = cross( j_hat, pos )
    jj = np.sqrt( phi_hat[:,0]*phi_hat[:,0] + phi_hat[:,1]*phi_hat[:,1] + phi_hat[:,2]*phi_hat[:,2] )
    for j in [0,1,2]: phi_hat[:,j] /= jj

    v_phi = vx*phi_hat[:,0] + vy*phi_hat[:,1] + vz*phi_hat[:,2]
    v_phi /= r

    ## define the absolute 'phi' relative to an arbitrary (but fixed) 90-degree rotation 
    ##   of the angular momentum vector j
    x_jhat = 0.*j_hat
    x_jhat[:,0]=0.*j_hat[:,0]; x_jhat[:,1]=-j_hat[:,2]; x_jhat[:,2]=j_hat[:,1]
    jj = np.sqrt( x_jhat[:,0]*x_jhat[:,0] + x_jhat[:,1]*x_jhat[:,1] + x_jhat[:,2]*x_jhat[:,2] )
    for j in [0,1,2]: x_jhat[:,j] /= jj
    
    ## generate y-vector by cross-product a x b = c (gaurantees right-hand rule)
    y_jhat = cross( j_hat, x_jhat )
    jj = np.sqrt( y_jhat[:,0]*y_jhat[:,0] + y_jhat[:,1]*y_jhat[:,1] + y_jhat[:,2]*y_jhat[:,2] )
    for j in [0,1,2]: y_jhat[:,j] /= jj

    ## now project r onto this, to obtain the components in this plane
    x_for_phi = x_jhat[:,0]*pos[:,0] + x_jhat[:,1]*pos[:,1] + x_jhat[:,2]*pos[:,2]
    y_for_phi = y_jhat[:,0]*pos[:,0] + y_jhat[:,1]*pos[:,1] + y_jhat[:,2]*pos[:,2]
    phi = np.arctan2( y_for_phi , x_for_phi )
    ## reset to a 0-2pi system instead of numpy's -pi,pi system
    lo = (phi < 0.); phi[lo] += 2.*np.pi

    gc.collect()
    return r, v_r, phi, v_phi, j_hat
