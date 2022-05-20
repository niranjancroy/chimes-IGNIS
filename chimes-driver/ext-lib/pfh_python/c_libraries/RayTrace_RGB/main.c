#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>           // for a better way of getting the time and time remaining
#include <omp.h>


//#define HLSCALE

/* 
    this routine takes a set of points and does an approximate ray-trace through them 
    (not an exact ray-trace, assuming each particle can be pre-integrated along the 
     z-axis, but this provides a tremendous speedup). 
    
    the input arrays MUST BE SORTED along the z axis, with increasing index number being 
    closer to the 'camera' 

    along with the basic information for the image grid (dimensions and size), 
    you send in the particle positions (in x,y), smoothing lengths (hsml), 
      and 'weights' (mass, and 'luminosity' in each of three 'bands')
      
*/


/* routine to allocate a 2D float array */
float** allocate2D(int nrows, int ncols) {
  int i;
  float **dat2;
  /*  allocate array of pointers  */
//  dat2 = malloc( nrows*sizeof(float*));
  dat2 = calloc( nrows, sizeof(float*));

  if(dat2==NULL) {
    printf("\nError allocating memory\n");
    exit(1);
  }
  /*  allocate each row  */
  for(i = 0; i < nrows; i++) {
//    dat2[i] = malloc( ncols*sizeof(float));
    dat2[i] = calloc( ncols, sizeof(float));
  }
  if(dat2[i-1]==NULL) {
    printf("\nError allocating memory\n");
    exit(1);
  }
  return dat2;
}
/* to free dat2 do:
  for(i=0;i<rows;i++)
    free(dat2[i]);
  free(dat2);
*/



/* extremely fast approximation function for the exponential, 
    useful here since fractional accuracy errors are smaller than the kernel sources anyways */
inline double fast_exp(double y) {
    double d;
    *((int*)(&d) + 0) = 0;
    *((int*)(&d) + 1) = (int)(1512775 * y + 1072632447);
    return d;
}

// changed to better fit python wrapper, not IDL //
int raytrace_rgb(
    int N_xy, // number of input particles/positions
    float* x, float* y, // positions (assumed already sorted in z)
    float* hsml, // smoothing lengths for each
    float* Mass, // total weight for 'extinction' part of calculation
    float* wt1, float* wt2, float* wt3, // weights for 'luminosities'
    float KAPPA1, float KAPPA2, float KAPPA3, // opacities for each channel
    float Xmin, float Xmax, float Ymin, float Ymax, // boundaries of output grid
    int Xpixels, int Ypixels, // dimensions of grid
    float hmax_scaling, // scaling of hamx.  was defaulting to 100; making it variable now
    float* OUT0, float* OUT1, float* OUT2, float* OUT3 ) // output vectors with final weights
{
  // print out the input parameters // 
  printf("-- starting ray-trace with %d particles total...\n",N_xy);
  printf("  Xmin=%f...Xmax=%f...Ymin=%f...Ymax=%f\n",Xmin,Xmax,Ymin,Ymax);
  printf("  Xpixels=%d...Ypixels=%d\n",Xpixels,Ypixels);
  printf("  Kappa_1=%f...Kappa_2=%f...Kappa_3=%f...\n",KAPPA1,KAPPA2,KAPPA3);

  double dx, dy, dx_i, dy_i, dx_n, dy_n, i_x_flt, i_y_flt, d_ij, h, hmax, hmin;
  double h2, x2_n, y2_n, r2_n, h2_i, wk, wt_sum, hkernel_over_hsml_to_use, *Kernel; 
  double dpi=3.1415926535897932384626433832795;
  long n,i,j,k,imin,imax,jmin,jmax,N_KERNEL_TABLE,print_modulus,npixels;

  time_t t1, t2;
  double elapsedTime;
  float pct_done, total_est;

  npixels = Xpixels*Ypixels;

  dx = (Xmax - Xmin)/((double)Xpixels);
  dy = (Ymax - Ymin)/((double)Ypixels);
  dx_i = 1./dx; dy_i = 1./dy;

  // hmax = 100.*sqrt(dx*dx+dy*dy); // set this purely to prevent going over too many cells //
  // hmax = 25.*sqrt(dx*dx+dy*dy);
  hmax = hmax_scaling * sqrt(dx*dx + dy*dy);  // cap the smoothing length.  if HLSCALE is defined, then the lum also goes to zero at hmax
  hmin = 0.5*sqrt(dx*dx+dy*dy); // ensures at least one cell 'sees' the particle // 
  
  // pre-define cell positions so we save a step in the loop //
  double x_i[Xpixels]; for(i=0;i<Xpixels;i++) x_i[i]=Xmin+dx*((double)i+0.5);
  double y_j[Ypixels]; for(i=0;i<Ypixels;i++) y_j[i]=Ymin+dy*((double)i+0.5);
  
  // build a kernel lookup table to save on the calculation below // 
  hkernel_over_hsml_to_use = 1.4; 
  //hkernel_over_hsml_to_use = 2.0; 
  // default = 1 (integrate out to kernel), but low-density regions are represented 
  //   more accurately if this is larger (~2); code expense increases though!
  N_KERNEL_TABLE = 1000;
  Kernel = calloc(N_KERNEL_TABLE+1, sizeof(double)); 
  dx_n=(hkernel_over_hsml_to_use*hkernel_over_hsml_to_use)/((double)N_KERNEL_TABLE); r2_n=0.;
  double kernel_spacing_inv = 1./dx_n;
  for(n=0;n<N_KERNEL_TABLE;n++)
  {
   h = sqrt(r2_n); // radius (to save sqrt operation we're interpolating in r^2 //
   // approximate gaussian for projected, integrate kernel: //
        wk = (135./(14.*dpi)) * exp(-(135./14.)*h*h); // quintic spline we're using now   
        //wk = (135./(14.*dpi)) * exp(-(135./14.)*h*h/(1.+h)); 
        // this has more extended tails, designed to reduce artifacts in low-density regions
        //  (where we would really want to use a proper volume-render)
        // wk = 1.91 * exp(-5.56*h*h); // cubic spline 
   // cubic spline kernel
        //h2=(1.-h); if(h<=0.5) {wk=(1.-6.*h*h*h2);} else {wk=2.*h2*h2*h2;} wk*=40./(7.0*dpi);
   Kernel[n] = wk;
   r2_n += dx_n; // radius at this point in the table
  }
  Kernel[N_KERNEL_TABLE]=0;
  
  // zero out the output vectors before the main sum
  for(n=0;n<npixels;n++)
  {
    OUT0[n]=0.0; OUT1[n]=0.0; OUT2[n]=0.0; OUT3[n]=0.0;
  }

#ifdef HLSCALE
    double hc = hmax / 2.0;  // probably makes the most sense to tie this to hmax, which will be where the lum goes to zero
    int lum_power = 2.0;  // higher keeps it closer to 1 to a higher level.  less than 2 and you don't get a smooth turnover
    printf("Will downweight luminosities for h > %g = hmax/%g", hc, hmax/hc);
    double lf_shift = pow(hmax - hc, lum_power);
    double lf_norm = 1.0/(lf_shift*lf_shift);
#else
    double lum_fact = 1.0;
#endif

  // loop over particles //
  // from talking with phil, the way to do this is in parallel is to have each processor operate on a chunk of the particles
  // which effectively means that each processor is doing a slab in space (i.e. at relatively constant z)
  // then you make sure that the slabs are in order, and multiply together the maps

  // furthermore, I'm going to want to use openMP for that so that I don't have to send copies of the particle data to each processor
  // best way would just be to have OUT0, 1, 2, 3 be arrays of arrays, where each thread (in order) handles one of those arrays

  // that's a task for TODAY!
  // while it's great to parallelize the movie making by having each processor do a different chunk of snapshots,
  // there's a ton of wasted procs because of the memory constraints.
  // would be awesome to speed up the raytrace, so let's try to implement the above.
  // basically that means that I want to split the particles up evenly among the threads
  // do the raytrace for each thread
  // then combine them all together at the end

  // ok so here's the deal (I think).
  // each thread should have a luminosity array and an attenuation array.
  // THAT LUMINOSITY ARRAY INCLUDES SELF ATTENUATION <-- very important
  // so then you can apply the atten to the (self atten) lum of the thread before (i.e. further away)
  // then add on your own (self atten) lum and move on to the next proc

  int nthreads, n_per_thread, threads_done;
  double index_norm, split_power;

  // will assign each thread 1/(thread_idx +1)**split_power fraction of each particles
  // for nthreads = 8 and split_power = 0.5, last thread gets 50% as many particles; for split_power = 1, last thread gets 12% as many
  // for nthreads = 16, these numbers get halved and so on
  // some tuning may be required to maximize effeciency:
  // the smaller hmax is, the smaller this number should be (to 0, of course, which should do linear splitting)
  // not making this an option for now though.
  split_power = 0.7;

  // quick parallel section to properly get the number of threads:
  #pragma omp parallel
  {
    nthreads = omp_get_num_threads();
  }
  if(nthreads==1) printf("!! -- WARNING:  called OpenMP routine, but omp_get_num_threads yielded %i !! \n",nthreads);
  // calculate the normalization for my splitting:
  // python:  index_norm = 1./np.sum(1./np.arange(1,nthreads+1)**p)
  index_norm = 0; // for now, this is actually the reciprocal though; then I'll take 1/it in a minute
  for(i=1;i<nthreads+1;i++) index_norm += (1.0 / pow(i, split_power));
  index_norm = 1.0/index_norm;

//  n_per_thread = (int) ceil(N_xy / nthreads);  // round up here and then trim down the number of particles
//  print_modulus = 1+(n_per_thread/ 200);   // update progress bar every half % (200 updates, should be relatively cheap)

  printf("-- allocating %i MB for the %i images\n", sizeof(float) * nthreads * npixels * 4 / (1024*1024), nthreads );
  float **lum_array1, **lum_array2, **lum_array3, **atten_array;
  lum_array1 = allocate2D(nthreads, npixels);
  lum_array2 = allocate2D(nthreads, npixels);
  lum_array3 = allocate2D(nthreads, npixels);
  atten_array = allocate2D(nthreads, npixels);

  // zero the arrays:
  for(i=0;i<nthreads;i++) {
    for(j=0;j<npixels;j++)
      {
        lum_array1[i][j] = 0.0; lum_array2[i][j] = 0.0; lum_array3[i][j] = 0.0; atten_array[i][j] = 0.0;
      }
    }

  printf("--splitting the particles by assigning each thread (1/n)^%.2f the total particle share\n", split_power);
  threads_done = 0;
  // all the varaibles used in the actual calculation have to be private (including all loop indices and everything)
  // note that I only have to include variables that were already declared -- any that I declare within the pragma are automatically private
  #pragma omp parallel private(n, i, j, imin, imax, jmin, jmax, i_x_flt, i_y_flt, h, h2_i, h2, d_ij, dx_n, x2_n, dy_n, y2_n, r2_n, k, wk, pct_done, elapsedTime, total_est, t2)
   {
    // define variables within the parallel because that makes them private
    int tid, start_idx, end_idx, thread_idx;
    int col, barWidth, seconds_left;
    double sum_temp = 0;
    barWidth = 60;
#ifdef HLSCALE
    double lum_fact;  // lum_fact is how much we modify the luminosity by
#endif

    // split the particles up by their z.  remember increasing index number being closer to the 'camera'
    tid = omp_get_thread_num();

    /*
    // do linear splitting (each thread gets an equal number of particles)
    start_idx = tid*n_per_thread;
    end_idx = (tid+1) * n_per_thread;
    */

    // because later (higher index) particles take longer (go over more pixels), front-load the work:
    if(tid==0)
        {start_idx = 0;}
    else
      {
        // python:  start_idx = np.ceil(Np * index_norm * np.sum(1./np.arange(1,ii)**p))
        for(i=1;i<tid+1;i++)    sum_temp += (1.0 / pow(i, split_power));
        start_idx = (int) ceil(N_xy * index_norm * sum_temp);
      }
    // python:  end_idx = start_idx + np.ceil(Np*index_norm*(1.0/ii)**p)
    end_idx = start_idx + (int) ceil(N_xy * index_norm * 1.0/pow(tid+1, split_power));
    if(end_idx>N_xy)    end_idx = N_xy;

    print_modulus = 1 + (end_idx-start_idx)/200;

    // print number of particles in first and last thread
    if((tid==0)||(tid==nthreads-1))      printf("...thread %i has %i particles", tid, end_idx - start_idx);

    // alternatively, print for all threads
    // printf("thread %i checking in for %i -- %i (n = %i)\n", tid, start_idx, end_idx, end_idx - start_idx);
    
    // ok this barrier is just to make the printing nice, but it should be very short
    #pragma omp barrier
    if(tid==nthreads-1)  {
        time(&t1);     // get the starting time with time.h
        printf("\n");  // clear the print from the previous 
      }
    // now do the actual loop, but only bother with my set of particles and save to my image
    for(n=start_idx;n<end_idx;n++) {
        i_x_flt = (x[n] - Xmin) * dx_i;
        i_y_flt = (y[n] - Ymin) * dy_i;
        h = hsml[n];
        h2 = h*h;

        if(h<hmin) h=hmin; // assume 'intrinsic' h is smeared by some fraction of pixel

        h2_i = 1./h2; // here we need the 'real' h (not the expanded search) //
        h *= hkernel_over_hsml_to_use; // make search area larger for kernel //

        if(h > hmax) { // this uses the expanded search radius //
            h=hmax;
            h2 = h*h;
          }

#ifdef HLSCALE
        if(h<=hc) { 
            lum_fact = 1.0; 
          }
        else {
            lum_fact =  lf_norm * pow(pow(h-hc,lum_power) - lf_shift, 2);
          }
#endif

        d_ij=h*dx_i;
        imin=(long)(i_x_flt-d_ij);
        imax=(long)(i_x_flt+d_ij)+1;

        if(imin<0) imin=0;
        if(imax>Xpixels-1) imax=Xpixels-1;

        d_ij=h*dy_i;
        jmin=(long)(i_y_flt-d_ij);
        jmax=(long)(i_y_flt+d_ij)+1;
        if(jmin<0) jmin=0;
        if(jmax>Ypixels-1) jmax=Ypixels-1;

        // now loop over pixels that this particle (once smeared by the kernel) falls on and do the sum (including self-atten)
        for(i=imin;i<imax;i++)
        {
         dx_n = x[n]-x_i[i];
         if (fabs(dx_n) < h)
         {
         x2_n = dx_n*dx_n;
         for(j=jmin;j<jmax;j++)
         {
           dy_n = y[n]-y_j[j];
           y2_n = dy_n*dy_n;
           r2_n = x2_n + y2_n;
           if (r2_n < h2)
           {
            // use the kernel lookup table compiled above for the weighting //
            r2_n *= h2_i*kernel_spacing_inv; // ok now have the separation in units of hsml, then kernel_table //
            k = (long)r2_n;
            wk = lum_fact * h2_i * (Kernel[k] + (Kernel[k+1]-Kernel[k])*(r2_n-k)); // ok that's the weighted result

            k = j + Ypixels*i; // j runs 0-Ypixels-1, so this provides the necessary indexing //
            //wt_sum += wk;

            // first 'extinct' the background
            if(Mass[n]>0.)
            {
                d_ij = Mass[n]*wk;
                lum_array1[tid][k] *= exp(-KAPPA1 * d_ij);
                lum_array2[tid][k] *= exp(-KAPPA2 * d_ij);
                lum_array3[tid][k] *= exp(-KAPPA3 * d_ij);
                atten_array[tid][k] += d_ij;
                // here the surface density extinct the background sources,
                //   with effective 'opacities' KAPPA1/2/3 in each channel
            }
            // now 'contribute' the particles own luminosity
            if(wt1[n] != 0.) lum_array1[tid][k] += wt1[n]*wk; // adds 'surface brightness' of wt1 to total in channel1
            if(wt2[n] != 0.) lum_array2[tid][k] += wt2[n]*wk;
            if(wt3[n] != 0.) lum_array3[tid][k] += wt3[n]*wk;
           } // if (r2_n < h2)
          } // for(j=jmin;j<jmax;j++)
          } // if (x2_n < h2)
        } // for(i=imin;i<imax;i++) -- end of loop over pixels

        // update the progress bar after finishing some fraction of the particles
        if( (!(n%print_modulus)||(n==end_idx-1)) && (tid==nthreads-1)) { // only the last thread prints progress bar
          time(&t2);
          pct_done = (n-start_idx)*100./(end_idx-start_idx);
          elapsedTime = difftime(t2, t1);    // difftime(a,b) returns a-b in seconds
          // total_est = 100./pct_done * elapsedTime;  // scaling linearly
          total_est = elapsedTime * pow(100./pct_done, 2); // scaling quadratically; probably more accurate
          seconds_left = (int) ceil(total_est - elapsedTime);
          printf("\rFinal thread: [");   // reset to the start of the line and reprint the progress bar each time

          // print either a space or a progress bar character, depending on our progress
          for(col=0; col < barWidth; col++)  {
              if(col<=(pct_done/100)*barWidth)  printf("=");
              else printf(" ");
            }
          printf("] %.0f%%, ", pct_done);

          // print either how much time is estimated remaining or how much time we took to get this far
          if(n==end_idx-1) {
              int elapsedTime_int = (int) ceil(elapsedTime);
              printf("%i:%02i tot,  ", elapsedTime_int/60, elapsedTime_int % 60);
              printf("%i/%i threads done   ", threads_done+1, nthreads);
            }
          else {
              printf("%i:%02i left, ", seconds_left/60, seconds_left % 60);
              printf("%i/%i threads done   ", threads_done, nthreads);
            }
          fflush(stdout);
          } // end of progress bar update check

      } // for(n=start_idx;n<end_idx;n++) -- end of particle loop
      threads_done++;  // increment once a thread gets here because it means it's finished it's particle loop

      // make sure all the threads have finished making their image
      #pragma omp barrier
      
      // ok, now I loop over the processors and combine all the images pixel by pixel

      // each pixel is completely independent of each other pixel, so this is a trivial 
      // for loop to parallelize:  just have to evenly split up the pixel indices to each 
      if(tid==0)  printf("\nNow layering images...");

      // need an indix for my pixel arrays to loop over in parallel
      int pidx_min, pidx_max, pidx, pixels_per_thread;
  
      // do linear splitting (each thread gets an equal number of pixels)
      // makes sense because each pixel should take the same amount of time
      pixels_per_thread = npixels / nthreads;
      pidx_min = tid * pixels_per_thread;
      pidx_max = (tid+1) * pixels_per_thread;
      if((tid==nthreads-1)&(pidx_max<npixels)) pidx_max = npixels;  //handle any rounding errors


      // each thread loops over the pixels assigned to it
      for(pidx=pidx_min; pidx<pidx_max; pidx++) {
        // each thread is looping over it's pixels, then looping through each of the arrays
        // because this is openmp, each thread has access to the arrays that the other threads
        // previously calculated from their share of particles
        for(thread_idx=0; thread_idx<nthreads; thread_idx++) {
          // for each pixel, do the following:
          // 1.  attenuate the existing images (i.e. the images behind it)
          // 2.  add on the (self-attenuated) luminosity from the current thread

          // attenuating images behind this thread:
          // (this does nothing if first thread, but that's fine -- 0 * anything = 0)
          OUT1[pidx] *= exp(-KAPPA1 * atten_array[thread_idx][pidx]);
          OUT2[pidx] *= exp(-KAPPA2 * atten_array[thread_idx][pidx]);
          OUT3[pidx] *= exp(-KAPPA3 * atten_array[thread_idx][pidx]);

          // then add on the (self attenuated) luminosity of this pixel from this thread
          OUT1[pidx] += lum_array1[thread_idx][pidx];
          OUT2[pidx] += lum_array2[thread_idx][pidx];
          OUT3[pidx] += lum_array3[thread_idx][pidx];

          // and for completeness, sum up the total attentuation along each pixel too
          // though I don't think that gets used anywhere
          OUT0[pidx] += atten_array[thread_idx][pidx];
        } // end loop over images (threads)
      } // end combination loop over pixels

      // wait until all threads are done with their pixels to free the images
      // #pragma omp barrier
    } // pragma omp parallel
    
    // now I'm totally done with the individual images and can free them:
    free(atten_array);
    free(lum_array1);
    free(lum_array2);
    free(lum_array3);

  // now print how long the entire calculation took
  time(&t2);
  elapsedTime = difftime(t2, t1) / 60.;    //difftime returns first minus second
  printf("done\nray-trace + stack took a total of %g minutes using %i threads\n", elapsedTime, nthreads);
  fflush(stdout);

  return 1;
} // closes main program 


