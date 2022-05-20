#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>                // for gettimeofday()
#include "proto.h"
#include "ngbtree3d.h"

void allocate_3d(int Nsize);
void set_particle_pointer(int Nsize, float* x, float* y, float* z);
void free_memory_3d(int Nsize);

struct particle_3d 
{
  float Pos[3];
} **P3d;


// revised call for python calling //
//int stellarhsml(int argc,void *argv[])
int stellarhsml(int N_in, float* x, float* y, float* z, int DesNgb, float Hmax, float* H_OUT)
{
  float h_guess, h2, xyz[3], dummy[3], h_guess_0;
  int i, ngbfound;
  allocate_3d(N_in);
  float *r2list;
  int *ngblist;

  struct timeval t1, t2;
  double elapsedTime;
  float pct_done, total_est;

  printf("Calculating stellar smoothing lengths for %d particles with Hmax = %g and %d neighbors\n",N_in, Hmax, DesNgb);
//  printf("N=%d\n",N_in); printf("Hmax=%g\n",Hmax); printf("DesNgb=%d\n",DesNgb);

  ngb3d_treeallocate(N_in, 10*N_in);
  set_particle_pointer(N_in, x,y,z);
  ngb3d_treebuild((float **)&P3d[1], N_in, 0, dummy, dummy);
  h_guess = Hmax/150.0e0; h_guess_0=h_guess;

  int print_modulus = (int) ceil(N_in / 20);   // print 20 times in the loop

  gettimeofday(&t1, NULL); // start the timer
  for(i=0;i<N_in;i++)
  {
	  xyz[0]=P3d[i+1]->Pos[0]+1.0e-10;
	  xyz[1]=P3d[i+1]->Pos[1]+1.0e-10;
	  xyz[2]=P3d[i+1]->Pos[2]+1.0e-10;
	  h2=ngb3d_treefind( xyz, DesNgb ,1.04*h_guess, &ngblist, &r2list, Hmax, &ngbfound); 

    if(!(i%print_modulus))
    {
      pct_done = i*100./N_in;
      printf("pct=%g i=%d hmax=%g h_guess=%g h=%g xyz=%g|%g|%g ngb=%d \n",
        pct_done,i,Hmax,h_guess,sqrt(h2),xyz[0],xyz[1],xyz[2],ngbfound); 

      if (i>0)
        {
          // get time elapsed
          gettimeofday(&t2, NULL);
          // compute the elapsed time in minutes
          elapsedTime = (t2.tv_sec - t1.tv_sec) / 60.0;
          // ok, so if it's taken elapsed time to do x% of the calculation, then scaling linearly:
          total_est = 100./pct_done * elapsedTime;
          printf("estimate %g minutes left\n", total_est - elapsedTime);
        }
      fflush(stdout);
    } // end print statement

      H_OUT[i] = sqrt(h2);
      h_guess = H_OUT[i]; // use this value for next guess, should speed things up // 
      //if (h_guess>10.*h_guess_0) h_guess=2.*h_guess_0;
    } 

  ngb3d_treefree();
  free_memory_3d(N_in);
  printf("done\n");
  return 0;
}



void set_particle_pointer(int Nsize, float* x, float* y, float* z)
{
  int i;
  float *pos;
  for(i=1;i<=Nsize;i++)
    {
      P3d[i]->Pos[0] = x[i-1];
      P3d[i]->Pos[1] = y[i-1];
      P3d[i]->Pos[2] = z[i-1];
    }
}


void allocate_3d(int Nsize)
{
  printf("allocating memory...\n");
  int i;
  if(Nsize>0)
    {
      if(!(P3d=malloc(Nsize*sizeof(struct particle_3d *))))
	{
	  printf("failed to allocate memory. (A)\n");
	  exit(0);
	}
      P3d--;   /* start with offset 1 */
      if(!(P3d[1]=malloc(Nsize*sizeof(struct particle_3d))))
	{
	  printf("failed to allocate memory. (B)\n");
	  exit(0);
	}
      for(i=2;i<=Nsize;i++)   /* initiliaze pointer table */
	P3d[i]=P3d[i-1]+1;
    }
  printf("allocating memory...done\n");
}


void free_memory_3d(int Nsize)
{
  if(Nsize>0)
    {
      free(P3d[1]);
      P3d++;
      free(P3d);
    }
}
