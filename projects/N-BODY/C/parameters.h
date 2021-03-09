#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <stdbool.h> 

#define id(width,row,col) (width*col+row)
#define GRAV_CONSTANT 0.5
#define THETA 1.0
#define TIMESTEP 1.
#define SIZEOFSPACE 10.0
#define NBRITERATIONS 4


#define SW_DOWN 0
#define SE_DOWN 1
#define NW_DOWN 2
#define NE_DOWN 3
#define SW_UP 4
#define SE_UP 5
#define NW_UP 6
#define NE_UP 7


typedef struct particles particle_t ;
typedef struct nodes node;

struct particles{
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double fx;
	double fy;
	double fz;
	double m;
	int id;
	double V;
	node * node;
};



struct nodes{
	particle_t * particle;
	int sub_nbr_particles;
	int depth;

	double minx;
	double maxx;
	double miny;
	double maxy;
	double minz;
	double maxz;

	double mass;
	double centerx;
	double centery;
	double centerz;
	node * parent;
	node * children;
};

#endif /*PARAMETERS_H_*/
