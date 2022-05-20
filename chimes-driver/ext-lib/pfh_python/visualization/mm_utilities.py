import numpy as np
import gc


## procedure which will return, for a given input vector A_in, 
##    the perpendicular unit vectors B_out and C_out 
##    which form perpendicular axes to A
def return_perp_vectors(a, LOUD=0):
    eps = 1.0e-10
    a = np.array(a,dtype='f');
    a /= np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    for i in range(len(a)):
        if (a[i]==0.): a[i]=eps;
        if (a[i]>=1.): a[i]=1.-eps;
        if (a[i]<=-1.): a[i]=-1.+eps;
    a /= np.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    ax=a[0]; ay=a[1]; az=a[2];

    ## use a fixed rotation of the a-vector by 90 degrees:
    ## (this anchors the solution so it changes *continously* as a changes)
    t0=np.double(math.pi/2.e0);
    bx=0.*ax; by=np.cos(t0)*ay-np.sin(t0)*az; bz=np.sin(t0)*ay+np.cos(t0)*az;
    ## c-sign is degenerate even for 'well-chosen' a and b: gaurantee right-hand 
    ##   rule is obeyed by defining c as the cross product: a x b = c
    cx=(ay*bz-az*by); cy=-(ax*bz-az*bx); cz=(ax*by-ay*bx); 
    B_out=np.zeros(3); C_out=np.zeros(3);
    B_out[:]=[bx,by,bz]; C_out[:]=[cx,cy,cz];
    
    if (LOUD==1):
        print(a )
        print(B_out)
        print(C_out)
        print('a_tot=',ax*ax+ay*ay+az*az)
        print('b_tot=',bx*bx+by*by+bz*bz)
        print('c_tot=',cx*cx+cy*cy+cz*cz)
        print('adotb=',ax*bx+ay*by+az*bz)
        print('adotc=',ax*cx+ay*cy+az*cz)
        print('bdotc=',bx*cx+by*cy+bz*cz)
    return B_out, C_out


def frame_ext(snum,four_char=1):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	gc.collect()
	return ext;



def ok_scan(input,xmax=1.0e10,pos=0):
    if (pos==1):
        gc.collect()
        return (np.isnan(input)==False) & (np.abs(input)<=xmax) & (input > 0.);
    else:
        gc.collect()
        return (np.isnan(input)==False) & (np.abs(input)<=xmax);


def idmatch(id_i, id_f): 
    ## match items in set (i) to those in set (f)
    ##  --- this is particularly speedy when the ids are all unique (case here) --
    index_dict_i = dict((k,i) for i,k in enumerate(id_i));
    index_dict_f = dict((k,i) for i,k in enumerate(id_f));
    inter = set(id_i).intersection(set(id_f));
    indices_i = np.array([ index_dict_i[x] for x in inter ]);
    indices_f = np.array([ index_dict_f[x] for x in inter ]);
    gc.collect()
    return indices_i, indices_f;


def compile_matched_ids( id_i, id_f ):
    sub_id_i, sub_id_f = idmatch(id_i,id_f)
    
    nomatch_from_i_in_f = (id_i > -1.0e40) ## should return all true
    if (sub_id_i.size > 0):
        nomatch_from_i_in_f[sub_id_i] = False ## is matched
    
    nomatch_from_f_in_i = (id_f > -1.0e40) 
    if (sub_id_f.size > 0):
        nomatch_from_f_in_i[sub_id_f] = False
    
    gc.collect()
    return sub_id_i, sub_id_f, nomatch_from_i_in_f, nomatch_from_f_in_i

    
def cross(x,y):
    #return np.cross(x,y,axis=1)
    c=0.*x
    c[:,0] = x[:,1]*y[:,2] - x[:,2]*y[:,1]
    c[:,1] =-x[:,0]*y[:,2] + x[:,2]*y[:,0]
    c[:,2] = x[:,0]*y[:,1] - x[:,1]*y[:,0]
    gc.collect()
    return c


def cosmological_time(a,h=0.71,Omega_M=0.27):
    ## exact solution for a flat universe
    x=Omega_M/(1.-Omega_M) / (a*a*a);
    t=(2./(3.*np.sqrt(1.-Omega_M))) * np.log( np.sqrt(x) / (-1. + np.sqrt(1.+x)) );
    t *= 13.777 * (0.71/h); ## in Gyr
    gc.collect()
    return t;
    
    
def fcor(x):
    gc.collect()
    return np.array(x,dtype='f',ndmin=1)

def vfloat(x):
    gc.collect()
    return x.ctypes.data_as(ctypes.POINTER(ctypes.c_float));

