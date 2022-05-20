#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int los_column_allstars(
    int N_s, // number of input particles/positions
    float* xs, float* ys, float* zs, // positions of sources
    int N_gas, // number of input particles/positions
    float* xg, float* yg, float* zg, // positions of sources
    float* wt, // total weight for 'extinction' part of calculation
    float* hsml, // smoothing lengths for each
    float* rg, // distance of gas from coordinate center [should be sorted on this!]
    float* OUT ) // output vectors with final weights
{
  long i, j, N_KERNEL_TABLE = 1000;
  double hkernel_over_hsml_to_use=1.4, dpi=3.1415926535897932384626433832795, r2_n=0, *Kernel = calloc(N_KERNEL_TABLE+1, sizeof(double));
  for(i=0;i<N_KERNEL_TABLE;i++)  // build a kernel lookup table to save on the calculation below // 
  {
    Kernel[i] = (135./(14.*dpi)) * exp(-(135./14.)*r2_n); // quintic spline we're using now (approximate gaussian for projected, integrate kernel)
    r2_n += hkernel_over_hsml_to_use*hkernel_over_hsml_to_use/(double)N_KERNEL_TABLE; // radius at this point in the table (to save sqrt operation we're interpolating in r^2)
  }
  Kernel[N_KERNEL_TABLE]=0;
  for(i=0;i<N_s;i++) {OUT[i]=0;} // zero before below
  for(i=0;i<N_gas;i++) {hsml[i] = 1./(hkernel_over_hsml_to_use*hsml[i]*hkernel_over_hsml_to_use*hsml[i]);} // this is the number we'll actually need
  for(i=0;i<N_s;i++)
  {
    if(i%1000) {continue;} // skip for speedup
    if(!(i%10000)) {printf(".%ld.",i);fflush(stdout);}
    OUT[i]=0;
    double x=xs[i],y=ys[i],z=zs[i],r=sqrt(x*x+y*y+z*z),xh=x/r,yh=y/r,zh=z/r;
    j = 0; 
    while((j < N_gas) && (rg[j] <= r))
    {
        double rproj = xg[j]*xh + yg[j]*yh + zg[j]*zh;
        if(rproj > 0)
        {
            double dx=rproj*xh-xg[j], dy=rproj*yh-yg[j], dz=rproj*zh-zg[j], dr2=(dx*dx+dy*dy+dz*dz)*hsml[j];
            if(dr2 < 1)
            {
                double q = dr2 * (double)N_KERNEL_TABLE;
                long k = (long)q;
                /*
                printf("xyz_s=%g/%g/%g xyz_g=%g/%g/%g rs=%g rg=%g dr=%g h=%g dX=%g q=%g k=%d Kernel=%g/%g/%g wt=%g \n",
                    x,y,z,xg[j],yg[j],zg[j],r,rproj,sqrt(dr2/hsml[j]),sqrt(1./(hkernel_over_hsml_to_use*hsml[j])),
                    dr2,q,k,Kernel[k],Kernel[k+1],(Kernel[k] + (Kernel[k+1]-Kernel[k])*(q-k)),wt[j]); fflush(stdout);
                */
                OUT[i] += wt[j] * (Kernel[k] + (Kernel[k+1]-Kernel[k])*(q-k)); // ok that's the weighted result
            }
        }
        j++;
    }
  }  
  return 1;
} // closes main program 


