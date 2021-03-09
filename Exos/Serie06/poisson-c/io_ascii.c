#include "io_ascii.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* -------------------------------------------------------------------------- */
int write_to_file(int n, float ** x, unsigned int iter, float minval,
                  float maxval) {
  FILE * fp;
  char filename[256];
  char prefix[] = "outputs/out_";
  char suffix[] = ".pgm";
  
  int i, j, lines;
  int v;

  sprintf(filename, "%s%05d%s", prefix, iter, suffix);
  
  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "Can't open output file %s!\n", filename);
    exit(1);
  }

  fprintf(fp, "%s\n%s\n", "P2", "# CREATOR: poisson program");
  fprintf(fp, "%d %d\n", n, n);
  fprintf(fp, "%d\n", 255);

  lines = 4;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      v = (int)roundf(255 * (x[i][j] - minval) / (maxval - minval));
      v = v <= 255 ? v : 255;
      fprintf(fp, "%d \n", (int)(v));
      lines++;
    }
  }

  fclose(fp);
  printf("%d lines written to file %s\n", lines, filename);
  return 0;
}

/* -------------------------------------------------------------------------- */
int colormap(int n) {
  FILE * fp;
  int i, j;
  double v;
  fp = fopen("colormap.pgm", "w");

  fprintf(fp, "%s\n%s\n", "P2", "# CREATOR: poisson program");
  fprintf(fp, "%d %d\n", 20, n);
  fprintf(fp, "%d\n", 255);

  for (j = 0; j < n; j++) {
    v = (int)roundf(255.0 * j / (n - 1));
    for (i = 0; i < 20; i++) {
      fprintf(fp, "%d \n", (int)(v));
    }
  }

  fclose(fp);
  return 0;
}
