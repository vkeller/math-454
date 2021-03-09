#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))

/*
From : http://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
*/

void print_matrix(float **matrix, int s);

int write_to_file_binary1D(int n, float **tmptab, int iter,
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

void print_matrix(float **matrix, int s)
{
    int p, q;
    for (p = 0; p < s; ++p)
    {
        for (q = 0; q < s; ++q)
            printf("%f ", matrix[p][q]);
        printf("\n");
    }
}

