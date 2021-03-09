PHPC - N-BODY PROJECT

HOWTO COMPILE AND RUN
=====================

Requirements : 

- a recent compiler (like gcc or intel)

compile on SCITAS clusters :

```
$ module load gcc
$ git clone ssh://git@c4science.ch/source/nbody-phpc.git
$ make
```

run on a SCITAS cluster interactively (reservation name is indicative) :

```
$ salloc -N 1 -n 1 --reservation=phpc-reservation
$ srun ./nbody-code examples/very-small.txt
```
You should see this output (timing is indicative) :

```
$ srun ./nbody-code examples/very-small.txt
====================================================
N-Body 3D simulation code for MATH-454 course EPFL  
Parameters for the Barnes-Hut algorithm:

Gravitational constant : 0.500000
Theta                  : 1.000000
Time step              : 1.000000
Space multiplicator    : 10.000000
Number of iterations   : 4

These parameters can be modified in "parameters.h"

(c) 2020, Vincent Keller (Vincent.Keller@epfl.ch)
====================================================

Read data from file 
Reading file ... OK
Number of particles : 10
BRUTE FORCE simulation starting 
ITERATION 0 
ITERATION 1 
ITERATION 2 
ITERATION 3 
N-Body brute force for 4 particles : 0.000023 [s] 
BARNES-HUT simulation starting 
Creation of the tree ...OK 
Construction of the tree from file ...OK 
Init forces ... OK 
ITERATION 1 
ITERATION 2 
ITERATION 3 
ITERATION 4 
It remains 10 particles in space 
N-Body barnes-hut for 10 particles : 0.000152 [s] 
Simulation finished 
```

The given example is a 10 bodies example. The example is in [Gadget 2](https://wwwmpa.mpa-garching.mpg.de/gadget) format translated into a text file using [this converter](https://github.com/martinsparre/Gadget2Conversion).


