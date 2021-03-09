#include "nbody_barneshut.h"

/*
Implementation of a barnes-hut algorithm for the N-Body problem.
*/
void nbodybarneshut (particle_t * array, int nbr_particles, int nbr_iterations) 
{
	int n;
	double step = TIMESTEP;
	node * root1;
	node * root2;
	node * root;
	particle_t tmp;

	printf("Creation of the tree ...");
	root1 = malloc(sizeof(node));	
	root2 = malloc(sizeof(node));	
	tmp = getMinMax(array, nbr_particles);
	init_tree(&tmp, root1);	
	init_tree(&tmp, root2);	
	printf("OK \n");

	printf("Construction of the tree from file ...");
	construct_bh_tree(array,root1, nbr_particles);
	printf("OK \n");
	printf("Init forces ...");
	compute_force_in_node(root1, root1);
	printf(" OK \n");
	for (n = 0 ; n  < nbr_iterations ; n++){
		printf("ITERATION %d \n",n+1);
		compute_force_in_node(root1, root1);
		compute_bh_force(root1);
		move_all_particles(root2, root1,step);
		root = root1;
		root1 = root2;
		root2 = root;
		clean_tree(root2);
	}

	printf("It remains %d particles in space \n",root1->sub_nbr_particles);	

	clean_tree(root1);
	clean_tree(root2);

	free(root1);
	free(root2);
}


/*

1. If the current node is an external node (and it is not body b), calculate the force exerced by the current node on b, and add this amount to b’s net force.
    
2. Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal node as a single body, and calculate the force it exerts on body b, and add this amount to b’s net force.

3. Otherwise, run the procedure recursively on each of the current node’s children.

Once the computation of the force applied to the particles is complete, the new position of the particles is computed, and a new tree corresponding to the new position is created. 

*/

/*
Move all the particles from node n to new_root
*/

void move_all_particles(node * new_root, node * n, double step) {
	int i;
	if(n->children != NULL){
		for (i = 0; i < 8; i++){
			move_all_particles(new_root, &n->children[i], step);
		}
	}else{
		particle_t * p = n->particle;
		move_particle(new_root, n, p,step);
	}
}

/*
Compute new position/velocity of the particle
*/

void move_particle(node * root, node * n, particle_t * p, double step) {
	double ax,ay,az;

	if ((p==NULL)||(n==NULL)) return;

	p->x += p->vx * step;	 
	p->y += p->vy * step;	 
	p->z += p->vz * step;	 
	ax = p->fx/p->m;
	ay = p->fy/p->m;
	az = p->fz/p->m;
	p->vx += ax*step;
	p->vy += ay*step;
	p->vz += az*step;
		
	if (! is_particle_out_of_scope(p,root)) {
		insert_particle(p,root);
	}else{
//		printf("\tparticle %d is out of scope. It will be destroyed at next iteration \n",p->id);
		n->particle = NULL;
	}
	
}

/*
Check if a particle is out of scope (lost body in space)
*/

bool is_particle_out_of_scope(particle_t * p, node * root){
	bool ret = false;
	if ((p->x < root->minx)||(p->y < root->miny)||(p->z < root->minz)) ret = true;
	if ((p->x > root->maxx)||(p->y > root->maxy)||(p->z > root->maxz)) ret = true;	
//	printf("\tmin abs : (%f : %f : %f)\t max abs : (%f : %f : %f)\n",root->minx,root->miny,root->minz, root->maxx,root->maxy,root->maxz);
//	printf("\tpar pos : (%f : %f : %f)\n",p->x,p->y,p->z);
	return ret;
}


/*
Clean tree root
*/
void clean_tree(node * root) {
	int i;
	if (root == NULL) {return;}

	if(root->children != NULL){
		for (i = 0; i < 8; i++){
			clean_tree(&root->children[i]);
		}
		free(root->children);
		root->children = NULL;
		root->sub_nbr_particles=0;
	}
/*
	}else{
		free(root->children);
		root->children = NULL;
	}
	free(root->children);
	root->children = NULL;
*/
}



/*
compute the forces on the BH tree
*/

void compute_bh_force(node * n) {
	int i;
	if(n->children != NULL){
		for (i = 0; i < 8; i++){
			compute_bh_force(&n->children[i]);
		}
	}else{
		particle_t * p = n->particle;
		compute_force_particle(n,p);
	}
}

/*
Compute force of node n on particle p
*/

void compute_force_particle(node * n, particle_t * p){
	int i;
	double diffx,diffy,diffz,distance;
	double size;
	if ((n==NULL)||(n->sub_nbr_particles==0)){ return;}

	if ((n->particle != NULL)&&(n->children==NULL)) {
		compute_force(p, n->centerx, n->centery,  n->centerz, n->mass) ;
	}
	else{
		size = n->maxx - n->minx;
		diffx = n->centerx - p->x;
		diffy = n->centery - p->y;
		diffz = n->centerz - p->z;
		distance = sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

//	The particle is far away. Use an approximation of the force
		if(size / distance < THETA) {
			compute_force(p, n->centerx, n->centery, n->centerz, n->mass);
		} else {

//      Otherwise, run the procedure recursively on each of the current node's children.
			for(i=0; i<8; i++) {
				compute_force_particle(&n->children[i], p);
			}
		}
	}
}



/*
Compute force 
*/

void compute_force(particle_t *p, double xpos, double ypos,  double zpos, double mass) {
	double xsep, ysep, zsep, dist_sq, gravity;

	xsep = xpos - p->x;
	ysep = ypos - p->y;
	zsep = zpos - p->z;
	dist_sq = max((xsep*xsep) + (ysep*ysep) + (zsep*zsep), 0.01);

	gravity = GRAV_CONSTANT*(p->m)*(mass)/ dist_sq / sqrt(dist_sq);

	p->fx += gravity*xsep;
	p->fy += gravity*ysep;
	p->fz += gravity*zsep;
}

/*
Compute all the forces in the particles
*/
void compute_force_in_node(node *n,node *root) {
	int i;
	if(n==NULL) return;

	if((n->particle != NULL)&&(n->children == NULL)) {
		particle_t*p = n->particle;
		p->fx = 0;
		p->fy = 0;
		p->fz = 0;
		compute_force_particle(root, p);
	}
	if(n->children != NULL) {
		for(i=0; i<8; i++) {
			compute_force_in_node(&n->children[i], root);
		}
	}
}




/*
Construction of the barnes-hut tree

Reminder:
Construction of the Barnes-Hut Tree in 2D
http://arborjs.org/docs/barnes-hut
*/

void construct_bh_tree(particle_t * array, node * root, int nbr_particles){
	int i;
	for (i=0;i < nbr_particles; i++){
		insert_particle(&array[i],root);
	}
}



/*
Add particle p in node n or one of its children
*/

void insert_particle(particle_t * p, node * n){
	int octrant ;
	double totalmass = 0.;
	double totalx = 0.;
	double totaly = 0.;
	double totalz = 0.;
	int i;
// there is no particle
	if ((n->sub_nbr_particles == 0)&&(n->children==NULL)) {
		n->particle = p;
		n->centerx = p->x;
		n->centery = p->y;
		n->centerz = p->z;
		n->mass = p->m;
		n->sub_nbr_particles++;
		p->node = n;
// There is already a particle
	}else{
		if (n->children==NULL){
			create_children(n);
			particle_t * particule_parent = n->particle;
// Insert the particle in the correct children
			octrant = get_octrant(particule_parent,n);
			n->particle = NULL;
			insert_particle(particule_parent,&n->children[octrant]);
		}
// insert the particle p
		octrant = get_octrant(p,n);
		insert_particle(p,&n->children[octrant]);

// Update mass and barycenter (sum of momentums / total mass)
		for(i=0; i<8; i++) {
			totalmass += n->children[i].mass;
			totalx += n->children[i].centerx*n->children[i].mass;
			totaly += n->children[i].centery*n->children[i].mass;
			totalz += n->children[i].centerz*n->children[i].mass;
		}
		n->mass = totalmass;
		n->centerx = totalx/totalmass;
		n->centery = totaly/totalmass;
		n->centerz = totalz/totalmass;
		p->node = n;
		n->sub_nbr_particles++;
	}
}

/*
create 8 children from 1 node
*/

void create_children(node * n){
	n->children = malloc(8*sizeof(node));

	double x12 = n->minx+(n->maxx-n->minx)/2.;
	double y12 = n->miny+(n->maxy-n->miny)/2.;
	double z12 = n->minz+(n->maxz-n->minz)/2.;

	init_node(&n->children[SW_DOWN], n, n->minx, x12, n->miny, y12, n->minz, z12 );
	init_node(&n->children[NW_DOWN], n, n->minx, x12, n->miny, y12, z12, n->maxz );

	init_node(&n->children[SE_DOWN], n, n->minx, x12, y12, n->maxy, n->minz, z12 );
	init_node(&n->children[NE_DOWN], n, n->minx, x12, y12, n->maxy, z12, n->maxz );

	init_node(&n->children[SW_UP], n, x12, n->maxx, n->miny, y12, n->minz, z12 );
	init_node(&n->children[NW_UP], n, x12, n->maxx, n->miny, y12, z12, n->maxz );

	init_node(&n->children[SE_UP], n, x12, n->maxx, y12, n->maxy, n->minz, z12 );
	init_node(&n->children[NE_UP], n, x12, n->maxx, y12, n->maxy, z12, n->maxz );	
}

/*
Init a node n attached to parent parent. 
*/

void init_node(node * n, node * parent,  double minx, double maxx, double miny, double maxy, double minz, double maxz ){
	n->parent=parent;
	n->children = NULL;
	n->minx = minx;
	n->maxx = maxx;
	n->miny = miny;
	n->maxy = maxy;
	n->minz = minz;
	n->maxz = maxz;
	n->depth = parent->depth + 1;
	n->particle = NULL;
	n->sub_nbr_particles = 0.;
	n->centerx = 0.;
	n->centery = 0.;
	n->centerz = 0.;
	n->mass = 0.;
}



/*
get the "octrant" where the particle resides (octrant is a generalization in 3D of a 2D quadrant)
*/

int get_octrant(particle_t * p, node * n){
	int octrant=-1;
	double xmin = n->minx;
	double xmax = n->maxx;
	double x_center = xmin+(xmax-xmin)/2;

	double ymin = n->miny;
	double ymax = n->maxy;
	double y_center = ymin+(ymax-ymin)/2;

	double zmin = n->minz;
	double zmax = n->maxz;
	double z_center = zmin+(zmax-zmin)/2;
	if (n==NULL) printf("ERROR: node is NULL \n");
	if (p==NULL) printf("ERROR: particle is NULL \n");

	// order : x -> y -> z
	if(p->x <= x_center) {
		if(p->y <= y_center) {
			if(p->z <= z_center) {
				octrant = SW_DOWN;
			}else{
				octrant = NW_DOWN;
			}
		} else {
			if(p->z <= z_center) {
				octrant = SE_DOWN;
			}else{
				octrant = NE_DOWN;
			}
		}
	} else {
		if(p->y <= y_center) {
			if(p->z <= z_center) {
				octrant = SW_UP;
			}else{
				octrant = NW_UP;
			}
		} else {
			if(p->z <= z_center) {
				octrant = SE_UP;
			}else{
				octrant = NE_UP;
			}
		}
	}
	return octrant;
}

/*
Init the tree

Remark :We use a particle struct to transfer min and max values from main
*/

void init_tree(particle_t * particle, node * root){
	root->minx = particle->x;
	root->maxx = particle->vx;
	root->miny = particle->y;
	root->maxy = particle->vy;
	root->minz = particle->z;
	root->maxz = particle->vz;
	root->particle = NULL;
	root->sub_nbr_particles = 0;
	root->parent = NULL;
	root->children = NULL;
	root->centerx = 0.;
	root->centery = 0.;
	root->centerz = 0.;
	root->mass = 0.;
	root->depth = 0;
}

/*
============================================
Utilities for testing
============================================
*/


/* print the tree */
void print_tree(node * root){
	node * tmp;
	int i;
	if (root->children!=NULL){
		for (i =0;i<8;i++){
			tmp = &root->children[i];
			print_tree(tmp);
		}
	}
	print_node(root);

}


/* 
print a node 
*/
void print_node(node * n){
	int d = n->depth;
	int i;
	for (i=0;i<d;i++){
		printf("\t");
	}
	printf("[level %d]",d);
	printf(" ([%f:%f:%f])",n->centerx, n->centery,n->centerz);
	printf(" Node ");
	printf(" M = %f", n->mass);
	printf(" has %d particles ", n->sub_nbr_particles);
	if (n->particle!=NULL){
		particle_t * p = n->particle;
		printf(". Particle ID = %d",p->id);
	}
	printf("\n");
}


/*
print a particle 
*/
void print_particle(particle_t * p){
	printf("[Particle %d]",p->id);
	printf(" position ([%f:%f:%f])",p->x, p->y, p->z);
	printf(" M = %f", p->m);
	printf("\n");
}



