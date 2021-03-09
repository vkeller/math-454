#ifndef IO_ASCII_H
#define IO_ASCII_H

/// Write a file a PGM file named out_<iter>.pgm containing them data stored in
/// x as a 2-dimensional array. Min and max val are used to re-scale the data
int write_to_file(int n, float ** x, unsigned int iter, float minval,
                  float maxval);

#endif
