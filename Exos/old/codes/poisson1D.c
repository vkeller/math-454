#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define N 16
#define eps 0.005

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define id(width,row,col) (width*row+col)

int write_to_file_binary1D(int n, float *tmptab, int iter, float minval, float maxval) ;


int main() {
  int     i, j, k;
  float *u;
  float *uo;
  float *f;
  float   h;
  float   l2;

  h = 1. / (N+1);

  // Allocation of arrays of pointers
  u  = (float*) malloc((N+1)*(N+1)*sizeof(float*));
  uo = (float*) malloc((N+1)*(N+1)*sizeof(float*));
  f  = (float*) malloc((N+1)*(N+1)*sizeof(float*));

  // initialization of u0 and f
  for(i = 0; i < N+1; i++) {
    for(j = 0; j < N+1; j++) {
      uo[id(N+1,i,j)] = 0;
      f [id(N+1,i,j)] = -2.*100. * M_PI * M_PI * sin(10.*M_PI*i*h) * sin(10.*M_PI*j*h);
printf("%d ",id(N,i,j));
    }
printf("\n");
  }

  k=0;
  do {
    l2 = 0.;

    for(i = 1; i < N; i++) {
      for(j = 1; j < N ;j++) {
        // computation of the new step
        u[id(N+1,i,j)] = 0.25 * ( uo[id(N+1,i-1,j)] + uo[id(N+1,i+1,j)] + uo[id(N+1,i,j-1)] + uo[id(N+1,i,j+1)] - f[id(N+1,i,j)]*h*h);

        // L2 norm
        l2 += (uo[id(N+1,i,j)] - u[id(N+1,i,j)])*(uo[id(N+1,i,j)] - u[id(N+1,i,j)]);
      }
    }

    // copy new grid in old grid
    for(i = 0; i < N+1; i++){
      for(j = 0; j < N+1; j++){
        uo[id(N+1,i,j)] = u[id(N+1,i,j)];
      }
    }

    // outputs
//    printf("[iteration %d] l2=%.5f\n", k, sqrt(l2));
    printf("[iteration %d] l2=%.5f\n", k, l2);

//    write_to_file_binary1D(N+1, u, k, -1., 1.);
    k++;
  } while(l2 > eps);

  printf("Size of N = %d, eps = %f \n",N,eps);

  // deallocation of the rows
    free(u);
    free(uo);
    free(f);

  return 0;
}


int write_to_file_binary1D(int n, float *tmptab, int iter,
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

//  printf("Print matrix here \n");
//  print_matrix(tmptab,n);


  for(i = 0; i < w; i++){
    for(j = 0; j < h; j++){
      x = i;
      y = (h -1) - j;

      v = ((tmptab[id(n,i,j)] - minval)/(maxval - minval));
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

