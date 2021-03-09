#ifndef READER_H_
#define READER_H_

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "nbody_bruteforce.h"

#include "parameters.h"

particle_t * read_test_case(const char * restrict fn);
particle_t getMinMax (particle_t * array, int nbr_particles) ;
int get_nbr_particles(const char * restrict fn);

#endif /*READER_H_*/
