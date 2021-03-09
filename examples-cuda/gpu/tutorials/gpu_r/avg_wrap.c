#include <Rmath.h>
#include <R.h>
#include <Rdefines.h>

void avg_cuda(double *heights, 
	      int nx,  // number of points in x directions
	      int ny,  // number of points in y directions
	      int r, // radius of circular neighborhood
	      double *output);  //output array

// =================================================
//
//              lmr_gpu_wrap
//
// This is a main routine called from R 
//
// =================================================
SEXP avg_wrap (SEXP heights_ptr, 
		   SEXP nx_in, 
		   SEXP ny_in, 
		   SEXP	r_in,   
		   SEXP output_ptr){ 


  // input/output arrays
  double *heights;
  double *output;

  // pass SEXP parameters to local variables
  int nx       = INTEGER_VALUE(nx_in);
  int ny       = INTEGER_VALUE(ny_in);
  double r       = NUMERIC_VALUE(r_in);

  // obtain pointers for the input vectors
  heights = NUMERIC_POINTER(heights_ptr);
  output  = NUMERIC_POINTER(output_ptr);
  
  // call cuda-run C routine
  avg_cuda( heights, nx, ny, r, output );

  // back to R
  return(R_NilValue);
}
