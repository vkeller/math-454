#include "nbody_bruteforce.h"

/*
Min and max functions
*/

double max(double x, double y) 
{
	return ((x) > (y) ? (x) : (y));
} 

double min(double x, double y) 
{
	return ((x) < (y) ? (x) : (y));
} 


/*
Implementation of a simple N-Body code in brute force.
The parallelization target is CUDA
Input format is 
*/
void nbodybruteforce (particle_t * array, int nbr_particles, int nbr_iterations) {

	int i,n;
	double step = 1.;
	for (n = 0 ; n  < nbr_iterations ; n++){
		printf("ITERATION %d \n",n);
		for (i = 0 ; i  < nbr_particles ; i++){
			compute_brute_force(&array[i], array, nbr_particles,step);
		}
	}
}

/*
Compute force (brute force method) of particle p2 on particle p1
Update particle p1
*/

void compute_brute_force(particle_t * p1, particle_t * array, int nbr_particles, double step) {
	double x_sep, y_sep, z_sep, dist_sq, grav_base;
	double F_x=0.;
	double F_y=0.;
	double F_z=0.;
	double a_x=0.;
	double a_y=0.;
	double a_z=0.;
	particle_t tmp;

	for (int i = 0 ; i  < nbr_particles ; i++){
		tmp = array[i];
		if (i!=p1->id){
			x_sep = p1->x - tmp.x;
			y_sep = p1->y - tmp.y;
			z_sep = p1->z - tmp.z;
			dist_sq = max((x_sep*x_sep) + (y_sep*y_sep) + (z_sep*z_sep), 0.01);
			grav_base = GRAV_CONSTANT*(p1->m)*(tmp.m)/dist_sq / sqrt(dist_sq);
			F_x += grav_base*x_sep;
			F_y += grav_base*y_sep;
			F_z += grav_base*z_sep;
		}
	}

// F = m a
// a = F/m
// V = a step
// pos = V * step
	a_x = F_x/p1->m;
	a_y = F_y/p1->m;
	a_z = F_z/p1->m;
	p1->vx += a_x * step;
	p1->vy += a_y * step;
	p1->vz += a_z * step;
	
	p1->x += p1->vx * step;
	p1->y += p1->vy * step;
	p1->z += p1->vz * step;

}




