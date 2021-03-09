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

  char filename[256] = {'0'};
  char prefix[] = "outputs/out_";
  char suffix[] = ".bmp";

  sprintf(filename, "%s%05d%s", prefix, iter, suffix);

  int prank, psize;
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  int start = 0;

  // removing the ghosts from the height
  if (psize > 1) {
    h = (prank == 0 || prank == psize - 1 ? h - 1 : h - 2);
    start = (prank == 0 ? 0 : 1);
  }

  // Gathering the size of every processors, this could be done as in the
  // constructor of the Simulation instead
  int * size_per_proc = (int *)malloc(psize*sizeof(int));
  MPI_Allgather(&h, 1, MPI_INT, size_per_proc, 1, MPI_INT,
                MPI_COMM_WORLD);

  // determining the local offset
  int offset_h = 0;
  for (int i = 0; i < prank; ++i) {
    offset_h += size_per_proc[i];
  }

  int total_h = offset_h;
  for (int i = prank; i < psize; ++i) {
    total_h += size_per_proc[i];
  }

  free(size_per_proc);
  
  int row_size = 3 * w;

  // if the file width (3*w) is not a multiple of 4 adds enough bytes to make it
  // a multiple of 4
  padding = (4 - (row_size) % 4) % 4;
  row_size += padding;

  int img_size = row_size * h;
  unsigned char *img = (unsigned char *)malloc(img_size);

  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++) {
      v = ((u[i+start][j] - minval) / (maxval - minval));

      r = v * 255; // Red channel
      g = v * 255; // Green channel
      b = v * 255; // Red channel

      r = min(r, 255);
      g = min(g, 255);
      b = min(b, 255);

      img[row_size * i + 3 * j + 2] = (unsigned char)(r);
      img[row_size * i + 3 * j + 1] = (unsigned char)(g);
      img[row_size * i + 3 * j + 0] = (unsigned char)(b);
    }
  }

  int total_img_size = row_size * total_h;
  filesize = 54 + total_img_size;

  MPI_File fh;
  MPI_Status status;

  // opening a file in write and create mode
  MPI_File_open(MPI_COMM_WORLD, filename,
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  // defining the size of the file
  MPI_File_set_size(fh, filesize);

  if (prank == 0) {
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
    MPI_File_write_at(fh, 0, bmpfileheader, 14, MPI_CHAR, &status);
    MPI_File_write_at(fh, 14, bmpinfoheader, 40, MPI_CHAR, &status);
  }

  int offset = 54 + row_size * offset_h;

  // We also could write that data with a write_at, the view is just to show
  // different possibilities
  MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

  MPI_File_write(fh, img, img_size, MPI_CHAR, &status);
  MPI_File_close(&fh);

  free(img);
}

void write_to_file(int m, int n, float **x, unsigned int iter, float minval,
                   float maxval) {
  write_bmp(n, m, x, iter, minval, maxval);
}
