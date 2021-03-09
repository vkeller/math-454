#ifndef IO_BINARY_H
#define IO_BINARY_H

/// Write a file a BMP file named out_<iter>.bmp containing the data stored in
/// x as a 2-dimensional array. Min and max val are used to re-scale the data
void write_to_file(int n, float ** x, unsigned int iter, float minval,
                  float maxval);

#endif
