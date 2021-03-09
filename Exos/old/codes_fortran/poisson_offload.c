#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <sys/time.h>

#include "parameters.h"

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

static int write_to_bmp(int n, float **tmptab, int iter,
                        float minval, float maxval);

double second();




//#pragma offload target(mic : target_id) \
//                             in(all_Vals : length(MAXSZ)) \
//                             inout(numEs) out(E_vals : length(MAXSZ/2))

int main() {
  int     i, j, k;
  float **u;
  float **uo;
  float **f;
  float   h;
  float   l2;
  double start, end;

  h = 1. / (N+1);

  // Allocation of arrays of pointers
  u  = (float**) malloc((N+1)*sizeof(float*));
  uo = (float**) malloc((N+1)*sizeof(float*));
  f  = (float**) malloc((N+1)*sizeof(float*));

  // allocation of the rows
  for(i=0; i<N+1 ;i++) {
    u [i] = (float*) malloc((N+1)*sizeof(float));
    uo[i] = (float*) malloc((N+1)*sizeof(float));
    f [i] = (float*) malloc((N+1)*sizeof(float));
  }


  // initialization of u0 and f
  for(i = 0; i < N+1; i++) {
    for(j = 0; j < N+1; j++) {
      uo[i][j] = 0;
      f [i][j] = -2.*100. * M_PI * M_PI * sin(10.*M_PI*i*h) * sin(10.*M_PI*j*h);
    }
  }

  k=0;
  start = second();
  do {
    l2 = 0.;

    for(i = 1; i < N; i++) {
      for(j = 1; j < N ;j++) {
        // computation of the new step
        u[i][j] = 0.25 * ( uo[i-1][j] + uo[i+1][j] + uo[i][j-1] + uo[i][j+1] - f[i][j]*h*h);

        // L2 norm
        l2 += (uo[i][j] - u[i][j])*(uo[i][j] - u[i][j]);
      }
    }

    // copy new grid in old grid
    for(i = 0; i < N+1; i++){
      for(j = 0; j < N+1; j++){
        uo[i][j] = u[i][j];
      }
    }

    // outputs
    printf("l2=%.5f\n", sqrt(l2));

    //write_to_bmp(N+1, u, k, -1., 1.);
    k++;
  } while(l2 > eps);
  end=second();

  printf("T=%.5f s for %d steps (%.5f s/step)\n", (end - start), k, (end - start)/k);


  // deallocation of the rows
  for(i=0; i<N+1 ;i++) {
    free(u [i]);
    free(uo[i]);
    free(f [i]);
  }

  // deallocate the pointers
  free(u);
  free(uo);
  free(f);

  return 0;
}

double second()
{
        struct timeval tp;
        struct timezone tzp;
        int i;

        i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}



int write_to_bmp(int n, float **tmptab, int iter,
                 float minval, float maxval) {
  unsigned char *img = NULL;
  int i,j;
  int x,y;
  int w,h;
  int padding, filesize, datasize;
  float v,r,g,b;
  char filename[50];
  FILE *f;
  unsigned char bmpfileheader[14] = {'B','M',0,0,0,0,0,0,0,0,54,0,0,0};
  unsigned char bmpinfoheader[40] = {40, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      1, 0,24, 0};

  w = n; // width
  h = n; // height

  padding = (4 - (w*3) % 4) % 4;

  datasize = (3*w + padding)*h;
  filesize = 54 + datasize;

  img = (unsigned char *)calloc(datasize, 1);

  for(i = 0; i < w; i++){
    for(j = 0; j < h; j++){
      x = i;
      y = (h -1) - j;

      v = ((tmptab[j][i] - minval)/(maxval - minval));
      r = v * 255; // Red channel
      g = v * 255; // Green channel
      b = v * 255; // Red channel

      r = min(r, 255);
      g = min(g, 255);
      b = min(b, 255);

      img[(x + y*w)*3 + y*padding + 2] = (unsigned char)(r);
      img[(x + y*w)*3 + y*padding + 1] = (unsigned char)(g);
      img[(x + y*w)*3 + y*padding + 0] = (unsigned char)(b);
    }
  }

  bmpfileheader[ 2] = (unsigned char)(filesize      );
  bmpfileheader[ 3] = (unsigned char)(filesize >>  8);
  bmpfileheader[ 4] = (unsigned char)(filesize >> 16);
  bmpfileheader[ 5] = (unsigned char)(filesize >> 24);

  bmpinfoheader[ 4] = (unsigned char)(       w      );
  bmpinfoheader[ 5] = (unsigned char)(       w >>  8);
  bmpinfoheader[ 6] = (unsigned char)(       w >> 16);
  bmpinfoheader[ 7] = (unsigned char)(       w >> 24);
  bmpinfoheader[ 8] = (unsigned char)(       h      );
  bmpinfoheader[ 9] = (unsigned char)(       h >>  8);
  bmpinfoheader[10] = (unsigned char)(       h >> 16);
  bmpinfoheader[11] = (unsigned char)(       h >> 24);
  bmpinfoheader[20] = (unsigned char)(datasize      );
  bmpinfoheader[21] = (unsigned char)(datasize >>  8);
  bmpinfoheader[22] = (unsigned char)(datasize >> 16);
  bmpinfoheader[23] = (unsigned char)(datasize >> 24);


  sprintf(filename, "out_%04d.bmp", iter);
  f = fopen(filename, "w");

  fwrite(bmpfileheader, 1, 14,f);
  fwrite(bmpinfoheader, 1, 40,f);

  fwrite(img, 1, datasize, f);

  fclose(f);

  free(img);
  return 0;

}
