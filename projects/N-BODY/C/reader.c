#include "reader.h"


/*
Reads an XYZ file containing the particles and their masses.
Use : https://github.com/martinsparre/Gadget2Conversion to read an IC's file example from Gadget 2

File format for this function is :

Line 1 		: N = number of particles
Line 2 -> N	: x : y : z : vx : vy : vz : m : ID : V

where :
	(x,y,z)		: position of particle ID
	(vx,vy,vz)	: velocity of particle ID
	m		: mass of particle ID
	ID		: unique ID number
	V		: potential (not used in this simulation)
*/


particle_t * read_test_case(const char * restrict fn)
{
	particle_t * mat;

	int nbr_particles = 0;
	float x,y,z,vx,vy,vz,m,V;
	int id;
	int i;

	FILE *f;

	
	if ((f = fopen(fn, "r")) == NULL) 
	{
		printf("Could not open file");
		exit(1);
	}
	rewind(f);
	fscanf(f, "%d", &nbr_particles);
	printf("Reading file ... ");
	mat = malloc(nbr_particles*sizeof(particle_t));

/*
File format : 
Line 2 -> N	: x : y : z : vx : vy : vz : m : ID : V
*/

	for (i = 0; i<nbr_particles; i++){
		fscanf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f", &x, &y, &z, &vx, &vy, &vz, &m, &id, &V);
		mat[i].x = x;
		mat[i].y = y;
		mat[i].z = z;
		mat[i].vx = vx;
		mat[i].vy = vy;
		mat[i].vz = vz;
		if (m==0.) m = 1.;
		mat[i].m = m;
		mat[i].id = id;
		mat[i].V = V;
		mat[i].node = NULL;
//		printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f", &x, &y, &z, &vx, &vy, &vz, &m, &id, &V);
		mat[i].fx = 0.;
		mat[i].fy = 0.;
		mat[i].fz = 0.;
	}
		printf("OK\n");



	if (f !=stdin) fclose(f);


	return mat;
}


int get_nbr_particles(const char * restrict fn)
{

	int nbr_part = 0;

	FILE *f;

	
	if ((f = fopen(fn, "r")) == NULL) 
	{
		printf("Could not open file");
		exit(1);
	}
	fscanf(f, "%d", &nbr_part);

	if (f !=stdin) fclose(f);

	return nbr_part;
}

particle_t getMinMax (particle_t * array, int nbr_particles) {
	int i;
	double minx = DBL_MAX;
	double maxx = DBL_MIN;
	double miny = DBL_MAX;
	double maxy = DBL_MIN;
	double minz = DBL_MAX;
	double maxz = DBL_MIN;
	double maxt, mint;
	particle_t tmp;
	
	for (i = 0; i<nbr_particles; i++){
		if (array[i].x < minx) minx = array[i].x;
		if (array[i].x > maxx) maxx = array[i].x;
		if (array[i].y < miny) miny = array[i].y;
		if (array[i].y > maxy) maxy = array[i].y;
		if (array[i].z < minz) minz = array[i].z;
		if (array[i].z > maxz) maxz = array[i].z;	
	}

	maxt = max(maxx,maxy);
	maxt = max(maxt,maxz);
	mint = min(minx,miny);
	mint = min(mint,minz);

	tmp.x = mint*SIZEOFSPACE;
	tmp.vx = maxt*SIZEOFSPACE;
	tmp.y = mint*SIZEOFSPACE;
	tmp.vy = maxt*SIZEOFSPACE;
	tmp.z = mint*SIZEOFSPACE;
	tmp.vz = maxt*SIZEOFSPACE;

	return tmp;

}



