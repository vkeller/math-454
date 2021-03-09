#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <mpi.h>

#define N 255

#define NORTH 0
#define SOUTH 1

#define eps 0.005

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#define ID(i, j) ((i) * n + (j))

int write_to_file_mpiio(int n,
                        int n_loc,
                        int offset,
                        float *tmptab,
                        int iter,
                        float minval, float maxval,
                        int prank);

int main(int argc, char *argv[]) {
  int     i, j, k;
  float  *u;
  float  *uo;
  float  *f;
  float   h  = 0.;
  float   l2 = 0.;

  MPI_Status status;
  double start, end;

  int n, n_loc, offset;
  int prank, psize;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);
  MPI_Comm_size(MPI_COMM_WORLD, &psize);

  n = N + 1;
  h = 1. / n;

  // Divide n in psize and distribute the excess to the (n % psize) proc
  n_loc = (n / psize) + (prank < n % psize ? 1 : 0);

  { // Computing the offset of where in the global array the local array is
    // located. (this is needed to initialize f) It could be computed locally
    // without communication but this way you see a Allgather
    int i;
    int buf[psize];
    MPI_Allgather(&n_loc, 1, MPI_INT, buf, 1, MPI_INT, MPI_COMM_WORLD);

    offset = 0;
    for(i = 0; i < prank; ++i) {
      offset += buf[i];
    }
  }

  // add 2 for north and south ghost
  n_loc += 2;

  // Allocation of arrays of pointers
  u  = (float *) malloc(n_loc * n * sizeof(float));
  uo = (float *) malloc(n_loc * n * sizeof(float));
  f  = (float *) malloc(n_loc * n * sizeof(float));

  // initialization of u0 and f
  for(i = 1; i < n_loc - 1; i++) {
    for(j = 0; j < n; j++) {
      u [ID(i, j)] = 0;
      uo[ID(i, j)] = 0;
      f [ID(i, j)] = -2.*100. * M_PI * M_PI * sin(10.*M_PI*((i-1) + offset)*h) * sin(10.*M_PI*j*h);
    }
  }

  k=0;

  start = MPI_Wtime();
  do {
    int i_start = 1;
    int i_end   = n_loc - 1;

    l2 = 0.;

    // First synchronize uo
    if(prank > 0) { // send recv with top proc
      MPI_Sendrecv(&uo[ID(1, 0)], n, MPI_FLOAT, prank - 1, NORTH,
                   &uo[ID(0, 0)], n, MPI_FLOAT, prank - 1, SOUTH, MPI_COMM_WORLD, &status);
    } else {
      ++i_start;
    }

    if(prank < psize - 1) { // send recv with bottom proc
      MPI_Sendrecv(&uo[ID(n_loc - 2, 0)], n, MPI_FLOAT, prank + 1, SOUTH,
                   &uo[ID(n_loc - 1, 0)], n, MPI_FLOAT, prank + 1, NORTH, MPI_COMM_WORLD, &status);
    } else {
      --i_end;
    }

    for(i = i_start; i < i_end; i++) {
      for(j = 1; j < n - 1; j++) {
        // computation of the new step
        u[ID(i, j)] = 0.25 * ( uo[ID(i-1, j)] + uo[ID(i+1, j)] +
                               uo[ID(i, j-1)] + uo[ID(i, j+1)] -
                               f[ID(i, j)]*h*h);

        // L2 norm
        l2 += (uo[ID(i, j)] - u[ID(i, j)])*(uo[ID(i, j)] - u[ID(i, j)]);
      }
    }

    // copy new grid in old grid appart from ghost that are communicated
    for(i = 1; i < n_loc - 1; i++) {
      for(j = 0; j < n; j++){
        uo[ID(i, j)] = u[ID(i, j)];
      }
    }


    // reduce the 12 norm to every proc
    MPI_Allreduce(MPI_IN_PLACE, &l2, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    // outputs
    if(prank == 0) {
      printf("l2=%.5f (k=%d)\n", sqrt(l2), k);
    }

    write_to_file_mpiio(n, n_loc - 2, offset,
                        u + n,
                        k,
                        -1, 1,
                        prank);

    k++;
  } while(l2 > eps);
  end = MPI_Wtime();

  if(prank == 0) {
    double t = end - start;
    printf("T=%.5f s (%.5f s/step)\n", t, t/k);
  }

  // deallocate the pointers
  free(u);
  free(uo);
  free(f);

  MPI_Finalize();

  return 0;
}


int write_to_file_mpiio(int n,
                        int n_loc,
                        int offset,
                        float *tmptab,
                        int iter,
                        float minval, float maxval,
                        int prank) {
  unsigned char *img = NULL;

  int i,j;
  int x,y;
  int w, h, h_loc;
  int padding, filesize;

  float v,r,g,b;

  char filename[8+4+1] = { '0' };
  char prefix[] = "out_";
  char suffix[] = ".bmp";
  char number[5];

  MPI_File fh;
  MPI_Status status;

  strcpy(filename,prefix);
  sprintf(number, "%d", iter);
  while(iter<1000){
    strcat(filename,"0");
    if(iter == 0) iter=1;
    iter *= 10;
  }
  strcat(filename,number);
  strcat(filename,suffix);


  w = n; // width
  h = n; // height
  h_loc = n_loc;

  padding = (4-(w*3)%4)%4;

  filesize = 54 + (3*w + padding)*h;
  img = (unsigned char *) calloc((3 * w + padding) * h_loc, 1);

  for(i = 0; i < w; i++){
    for(j = 0; j < h_loc; j++){
      x = i;
      y = (h_loc - 1) - j;

      v = ((tmptab[ID(j, i)] - minval)/(maxval - minval));
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

  unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
  unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
  unsigned char bmppad[3] = {0,0,0};

  bmpfileheader[ 2] = (unsigned char)(filesize    );
  bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
  bmpfileheader[ 4] = (unsigned char)(filesize>>16);
  bmpfileheader[ 5] = (unsigned char)(filesize>>24);

  bmpinfoheader[ 4] = (unsigned char)(       w    );
  bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
  bmpinfoheader[ 6] = (unsigned char)(       w>>16);
  bmpinfoheader[ 7] = (unsigned char)(       w>>24);
  bmpinfoheader[ 8] = (unsigned char)(       h    );
  bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
  bmpinfoheader[10] = (unsigned char)(       h>>16);
  bmpinfoheader[11] = (unsigned char)(       h>>24);
  bmpinfoheader[20] = (unsigned char)( (filesize-54)    );
  bmpinfoheader[21] = (unsigned char)( (filesize-54)>> 8);
  bmpinfoheader[22] = (unsigned char)( (filesize-54)>>16);
  bmpinfoheader[23] = (unsigned char)( (filesize-54)>>24);


  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fh);
  MPI_File_set_size(fh, filesize);

  if(prank == 0) {
    MPI_File_write_at(fh,  0, bmpfileheader, 14, MPI_CHAR, &status);
    MPI_File_write_at(fh, 14, bmpinfoheader, 40, MPI_CHAR, &status);
  }

  offset = (h - offset - h_loc) * (3 * w + padding);
  MPI_File_set_view(fh, 54 + offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

  MPI_File_write(fh, img, (3*w + padding) * h_loc, MPI_CHAR, &status);
  MPI_File_close(&fh);
}
