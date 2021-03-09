#include "io_binary.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

static void write_bmp(int w, int h, float ** u, unsigned int iter, float minval,
                     float maxval) {
  int i, j;
  int padding, filesize;

  float v, r, g, b;

  FILE * fh;

  char filename[8 + 4 + 1] = {'0'};
  char prefix[] = "out_";
  char suffix[] = ".bmp";

  sprintf(filename, "%s%05d%s", prefix, iter, suffix);

  if ((fh = fopen(filename, "wb")) == NULL) {
    printf("Could not open file %s.", filename);
    exit(1);
  }

  int row_size = 3 * w;
  // if the file width (3*w) is not a multiple of 4 adds enough bytes to make it
  // a multiple of 4
  padding = (4 - (row_size) % 4) % 4;
  row_size += padding;

  filesize = 54 + (row_size)*h;

  unsigned char * img;
  img = (unsigned char *)calloc((row_size)*h, 1);

  for (i = 0; i < h; i++) {
    for (j = 0; j < w; j++) {
      v = ((u[h - 1 - i][j] - minval) / (maxval - minval));

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
  bmpinfoheader[8] = (unsigned char)(h);
  bmpinfoheader[9] = (unsigned char)(h >> 8);
  bmpinfoheader[10] = (unsigned char)(h >> 16);
  bmpinfoheader[11] = (unsigned char)(h >> 24);
  bmpinfoheader[20] = (unsigned char)((filesize - 54));
  bmpinfoheader[21] = (unsigned char)((filesize - 54) >> 8);
  bmpinfoheader[22] = (unsigned char)((filesize - 54) >> 16);
  bmpinfoheader[23] = (unsigned char)((filesize - 54) >> 24);

  fwrite(bmpfileheader,1, 14, fh);
  fwrite(bmpinfoheader,1, 40, fh);

  fwrite(img, 1, h * row_size, fh);

  fclose(fh);

  free(img);
}

void write_to_file(int n, float ** x, unsigned int iter, float minval,
                  float maxval) {
  return write_bmp(n, n, x, iter, minval, maxval);
}
