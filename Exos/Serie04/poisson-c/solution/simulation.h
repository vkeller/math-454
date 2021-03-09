#ifndef SIMULATION_H_
#define SIMULATION_H_

void allocate_grids(int n, float *** u, float *** uo, float *** f);
void deallocate_grids(float ***uo, float ***u, float ***f);
void swap_grids(float *** uo, float *** u);
void initialize_grids(int n, float ** uo, float ** u,
                      float ** f, float h);
float compute_step(int n, float ** uo, float ** u, float ** f, float h);

float simulate(int n, float **uo, float **u, float **f, float h, float epsilon, int * k);

#endif // SIMULATION_H_
