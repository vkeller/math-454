#include "io_binary.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

static void write_bmp(int w, int h, float **u, unsigned int iter, float minval,
                      float maxval) {
  int i, j;
  int padding, filesize;

  float v, r, g, b;

  FILE *fh;

  char filename[256] = {'0'};
  char prefix[] = "outputs/out_";
  char suffix[] = ".bmp";

  sprintf(filename, "%s%05d%s", prefix, iter, suffix);

  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  int i_start = (prank == 0 ? 0 : 1);
  int i_end = (prank == psize - 1 ? h : h - 1);

  h = i_end - i_start;
  
  int row_size = 3 * w;
  
  // if the file width (3*w) is not a multiple of 4 adds enough bytes to make it
  // a multiple of 4
  padding = (4 - (row_size) % 4) % 4;
  row_size += padding;

  int img_size = row_size * h;
  unsigned char *img = (unsigned char *) malloc(img_size);
 
  for (i = i_start; i < i_end; i++) {
    for (j = 0; j < w; j++) {
      v = ((u[i][j] - minval) / (maxval - minval));

      r = v * 255; // Red channel
      g = v * 255; // Green channel
      b = v * 255; // Red channel

      r = min(r, 255);
      g = min(g, 255);
      b = min(b, 255);

      img[row_size * (i-i_start) + 3 * j + 2] = (unsigned char)(r);
      img[row_size * (i-i_start) + 3 * j + 1] = (unsigned char)(g);
      img[row_size * (i-i_start) + 3 * j + 0] = (unsigned char)(b);
    }
  }

  int total_h = 0;
  MPI_Allreduce(&h, &total_h, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (total_h / psize != h) {
    printf("The dumper cannot handle this case %d / %d != %d\n", total_h, psize, h);
    exit(1);
  }
  
  if (prank != 0) {
    MPI_Gather(img, img_size, MPI_CHAR,
               NULL, 0, MPI_CHAR, 0,
               MPI_COMM_WORLD);

    //free(img);
    return;
  }

  int total_img_size = row_size*total_h;
  filesize = 54 + total_img_size;  
  unsigned char *total_img = (unsigned char *) malloc (total_img_size);
    
  MPI_Gather(img, img_size, MPI_CHAR,
             total_img, img_size, MPI_CHAR,
             0, MPI_COMM_WORLD);
  
  unsigned char bmpfileheader[14] = {'B', 'M', 0, 0,  0, 0, 0,
                                     0,   0,   0, 54, 0, 0, 0};
  unsigned char bmpinfoheader[40] = {40, 0, 0, 0, 0, 0, 0,  0,
                                     0,  0, 0, 0, 1, 0, 24, 0};

  bmpfileheader[2] = (unsigned char)(filesize);
  bmpfileheader[3] = (unsigned char)(filesize >> 8);
  bmpfileheader[4] = (unsigned char)(filesize >> 16);
  bmpfileheader[5] = (unsigned char)(filesize >> 24);

  bmpinfoheader[4] = (unsigned char)(w);
  bmpinfoheader[5] = (unsigned char)(w >> 8);
  bmpinfoheader[6] = (unsigned char)(w >> 16);
  bmpinfoheader[7] = (unsigned char)(w >> 24);
  bmpinfoheader[8] = (unsigned char)(total_h);
  bmpinfoheader[9] = (unsigned char)(total_h >> 8);
  bmpinfoheader[10] = (unsigned char)(total_h >> 16);
  bmpinfoheader[11] = (unsigned char)(total_h >> 24);
  bmpinfoheader[20] = (unsigned char)((filesize - 54));
  bmpinfoheader[21] = (unsigned char)((filesize - 54) >> 8);
  bmpinfoheader[22] = (unsigned char)((filesize - 54) >> 16);
  bmpinfoheader[23] = (unsigned char)((filesize - 54) >> 24);

  if ((fh = fopen(filename, "wb")) == NULL) {
    printf("Could not open file %s.", filename);
    exit(1);
  }

  fwrite(bmpfileheader, 1, 14, fh);
  fwrite(bmpinfoheader, 1, 40, fh);

  fwrite(total_img, 1, total_h * row_size, fh);

  fclose(fh);

  free(img);
  free(total_img);
}

void write_to_file(int m, int n, float **x, unsigned int iter, float minval,
                   float maxval) {
  write_bmp(n, m, x, iter, minval, maxval);
}
