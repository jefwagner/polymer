#include <math.h>

#include <stdio.h>

#include "polymer.h"
#include "rng.c"
#include "vec3.h"

static void disk( uint64_t *rng, double *o){
	do{
		o[0] = 2.0*rng_uniform(rng)-1.0;
		o[1] = 2.0*rng_uniform(rng)-1.0;
	}while(o[0]*o[0]+o[1]*o[1] > 1.0);
}

static void sphere( uint64_t *rng, double *o){
	do{
		o[0] = 2.0*rng_uniform(rng)-1.0;
		o[1] = 2.0*rng_uniform(rng)-1.0;
		o[2] = 2.0*rng_uniform(rng)-1.0;
	}while(o[0]*o[0]+o[1]*o[1]+o[2]*o[2] > 1.0);
}

static void shell( uint64_t *rng, double a, double b, double *o){
	int i;
	double r, scale;
	sphere(rng, o);
	r = sqrt( o[0]*o[0]+o[1]*o[1]+o[2]*o[2]);
	scale = cbrt( (b*b*b-a*a*a)+(a/r)*(a/r)*(a/r));
	for( i=0; i<3; i++){
		o[i] *= scale;
	}
}

static void donut( uint64_t *rng, double R, double a, double b, double *o){
	double u[2];
	double z, phi, rho;
	phi = TWOPI*rng_uniform(rng);
	if( R > a ){
		disk(rng, u);
		rho = sqrt( ((R+a)*(R+a) - (R-a)*(R-a))*0.5*(u[0]+1.0) + (R-a)*(R-a) );
	}else{
		double ll = R/a;
		do{
			disk(rng, u);
		}while(u[0] < -ll);
		rho = (R+a)*sqrt( (u[0]+ll)/(1.0+ll) );
	}
	z = b*u[1];

	o[0] = rho*cos(phi);
	o[1] = rho*sin(phi);
	o[2] = z;
}

static void ellipsoid( uint64 *rng, double a, double b, double *o){
	double u[3];
	sphere(rng, u);
	o[0] = a*u[0];
	o[1] = a*u[1];
	o[2] = b*u[2];
}

void rand1( uint64_t *rng, double r_min, double r_max, 
			double *p0, double *output){
	int i;
	shell( rng, r_min, r_max, output);
	for( i=0; i<3; i++){
		output[i] += p0[i];
	}
}

void rand2( uint64_t *rng, double r_min, double r_max,
		    double *p0, double *p1, double *output){
	int i;
	double r0;
	double tv[3], t0, dt;
	double nv[3], n0, nmax;
	double bv[3], tnb[3];

	double zhat[] = {0.0,0.0,1.0};
	double xhat[] = {1.0,0.0,0.0};
	double scale, tol = 1.0e-3;

	// Get the tangent vector between the points p0 and p1
	vec_dif( p1, p0, tv);
	// Distance to the center!
	t0 = 0.5*mag(tv);
	// Normalize the tangent vector
	scale = 0.5/t0;
	for( i=0; i<3; i++){
		tv[i] *= scale;
	}

	if( t0 > r_max ){
		// If the polymer is stretched out, put it right in the center.
		for( i=0; i<3; i++){
			output[i] = p0[i] + tv[i];
		}
	}else{
		// Otherwize find the normal and binormal vectors.
		cross( tv, zhat, nv);
		if( mag(nv) < tol ){
			cross( tv, xhat, nv);
		}
		scale = 1.0/mag(nv);
		for( i=0; i<3; i++){
			nv[i] *= scale;
		}
		cross( tv, nv, bv);
		scale = 1.0/mag(bv);

		// Base our work from between the two points
		for( i=0; i<3; i++){
			output[i] = p0[i] + t0*tv[i];
		}

		// Find the default radius and maximum normal extent
		r0 = 0.5*(r_min+r_max);
		nmax = sqrt(r_max*r_max - t0*t0);
		if( t0 > r0 ){
			// If distance is larger than the default use
			// and ellipsoidal distribution
			dt = r_max - t0;
			ellipsoid( rng, nmax, dt, tnb);
		}else{
			// Otherwise use a donut shaped distribution
			n0 = sqrt(r0*r0-t0*t0);
			dt = sqrt(r_max*r_max+t0*t0-r0*r0) - t0;
			donut( rng, n0, nmax-n0, dt, tnb);
		}
		// Then add the tangent, normal and binormal parts.
		output[0] += tv[0]*tnb[2] + nv[0]*tnb[1] + bv[0]*tnb[0];
		output[1] += tv[1]*tnb[2] + nv[1]*tnb[1] + bv[1]*tnb[0];
		output[2] += tv[2]*tnb[2] + nv[2]*tnb[1] + bv[2]*tnb[0];
	}
}

void rand3( uint64_t *rng, double r_min, double r_max,
			double *p0, double *p1, double *p2, double *output){
	output[0] = 0.0;
	output[1] = 0.0;
	output[2] = 0.0;
}

#define NEIGHBOR_MAX 20

int is_nan( double x){
	return isnan(x);
}

double delta_en( simulation *sim, uint16 monomer, double *pos){
	polymer *poly = sim->poly;
	uint16 num_b = poly->num_bonds[monomer];
	double r_co = sim->lim.r_co;
	// double en;
	double en_new = 0.0, en_old = 0.0;
	double *old_pos, *other_pos;
	uint16 i, bond, other;
	uint16 neighbors[NEIGHBOR_MAX];
	uint16 num_nei = 0;

	old_pos = &poly->positions[DIM*monomer];

	// Loop over the bonds
	for( i=0; i<num_b; i++){
		bond = poly->list_bonds[NB_MAX*monomer+i];
		other = linked_monomer(&poly->bonds[2*bond], monomer);
		other_pos = &poly->positions[DIM*other];
		// en = E_bond(pos, other_pos, sim->enp);
		// 	if( isnan(en) ){
		// 		en += 1.0;
		// 	}
		// en_new += en;
		en_new += E_bond( pos, other_pos, sim->enp);
		// en = E_bond(old_pos, other_pos, sim->enp);
		// 	if( isnan(en) ){
		// 		en += 1.0;
		// 	}
		// en_old += en;
		en_old += E_bond( old_pos, other_pos, sim->enp);
	}

	// Loop over close neighbors
	num_nei = 0;
	kd_ball_search( sim->kdtree, sim->kdtree->root, pos, 
					r_co, neighbors, &num_nei);
	for( i=0; i<num_nei; i++){
		if( neighbors[i] != monomer ){
			other_pos = &poly->positions[DIM*neighbors[i]];
			// en = E_pair( pos, other_pos, sim->enp);
			// if( isnan(en) ){
			// 	en += 1.0;
			// }
			// en_new += en;
			en_new += E_pair( pos, other_pos, sim->enp);
		}
	}
	num_nei = 0;
	kd_ball_search( sim->kdtree, sim->kdtree->root, old_pos, 
					r_co, neighbors, &num_nei);
	for( i=0; i<num_nei; i++){
		if( neighbors[i] != monomer ){
			other_pos = &poly->positions[DIM*neighbors[i]];
			// en = E_pair( old_pos, other_pos, sim->enp);
			// if( isnan(en) ){
			// 	en += 1.0;
			// }
			// en_old += en;
			en_old += E_pair( old_pos, other_pos, sim->enp);
		}
	}

	return en_new - en_old;
}

double energy( simulation *sim){
	polymer *poly = sim->poly;
	double r_co = sim->lim.r_co;
	uint16 i, j, mon0, mon1;
	double *p0, *p1;
	uint16 neighbors[NEIGHBOR_MAX];
	uint16 num_nei = 0;
	double en = 0.0;

	for( i=0; i<(sim->Nm-1); i++){
		mon0 = poly->bonds[2*i+0];
		mon1 = poly->bonds[2*i+1];
		p0 = &(poly->positions[DIM*mon0]);
		p1 = &(poly->positions[DIM*mon1]);
		en += E_bond( p0, p1, sim->enp);
	}

	for( i=0; i<sim->Nm; i++){
		p0 = &(poly->positions[DIM*i]);
		num_nei = 0;
		kd_ball_search( sim->kdtree, sim->kdtree->root, p0,
						r_co, neighbors, &num_nei);
		for( j=0; j<num_nei; j++){
			if( i != neighbors[j]){
				p1 = &(poly->positions[DIM*neighbors[j]]);
				if( p0[0] > p1[0] ){
					en += 0.5*E_pair( p0, p1, sim->enp);
				}
			}
		}
	}

	return en;
}

void mh_step( simulation *sim, uint16 monomer){
 	uint16 nb, b, *bond, li;
 	double *p0, *p1, new[3];
 	double en, a, r, kt=1.0;

 	b = sim->poly->list_bonds[NB_MAX*monomer+0];
 	bond = &sim->poly->bonds[2*b];
 	li = linked_monomer( bond, monomer);
 	p0 = &sim->poly->positions[DIM*li];
 	// Test for 1 or 2 monomers
 	nb = sim->poly->num_bonds[monomer];
 	if( nb == 1 ){
 		// If one monomer, use the rand1 distribution
 	  rand1( sim->rng, sim->lim.r_min, sim->lim.r_max, p0, new);
 	}else if( nb == 2){
 		// If two monomers, get the second monomer and use rand2
 	  b = sim->poly->list_bonds[NB_MAX*monomer+1];
 	  bond = &sim->poly->bonds[2*b];
 	  li = linked_monomer( bond, monomer);
 	  p1 = &sim->poly->positions[DIM*li];
 	  rand2( sim->rng, sim->lim.r_min, sim->lim.r_max, p0, p1, new);
 	}
 	// Get the change in energy
 	en = delta_en( sim, monomer, new);
 	// And do the Metropolis Hastings step
 	a = exp( -en/kt);
 	r = rng_uniform( sim->rng);
 	if( a>r ){
 	  // fprintf( stdout, "moving monommer %lu\n", monomer);

 	  kd_move( sim->kdtree, monomer, new);

 	  uint16 i;
 	  kd_node *t1;
 	  for( i=0; i<sim->Nm; i++){
 	  	t1 = kd_find( sim->kdtree, sim->kdtree->root, i);
 	  	if( t1 == NULL ){
 	  		fprintf( stdout, "lost monomer %lu !\n", i);
 	  	}
 	  }
 	}
}

void gibbs_step( simulation *sim){
	int i;
	for( i=0; i<sim->Nm; i++){
		mh_step( sim, i);
	}	
}

// void gibbs_step( simulation *sim){
//  	uint16 *order = (uint16 *) malloc( sim.Nm*sizeof(uint16));
//  	uint16 i, j, temp;
//  	uint16 nb, b, *bond, li;
//  	double *p0, *p1, new[3];
//  	double en, a, r, kt=1.0;
//  	// Do a Fisher-Yates shuffle to go through the polymer in random order
//  	for( i=0; i<sim.Nm; i++){
//  	  order[i] = i;
//  	}
//  	for( i=sim.Nm; i>1; i++){
//  	  j = rng_next(sim.rng)%i;
//  	  temp = order[i-1];
//  	  order[i-1] = order[j];
//  	  order[j] = order[i-1];
//  	}
//  	// Loop over the monomers
//  	for( i=0; i<sim.Nm; i++){
//  	  j = order[i];
//  	  // Get the first linked monomer
//  	  b = sim.poly->list_bonds[NB_MAX*j+0];
//  	  bond = &sim.poly->bonds[2*b];
//  	  li = linked_monomer( bond, j);
//  	  p0 = &sim.poly->positions[DIM*li];
//  	  // Test for 1 or 2 monomers
//  	  nb = sim.poly->num_bonds[j];
//  	  if( nb == 1 ){
//  	  	// If one monomer, use the rand1 distribution
//  	    rand1( sim, sim.lim.r_min, sim.lim.r_max, p0, new);
//  	  }else if( nb == 2){
//  	  	// If two monomers, get the second monomer and use rand2
//  	    b = sim.poly->list_bonds[NB_MAX*j+1];
//  	    bond = &sim.poly->bonds[2*b];
//  	    li = linked_monomer( bond, j);
//  	    p1 = &sim.poly->positions[DIM*li];
//  	    rand2( sim, sim.lim.r_min, sim.lim.r_max, p0, p1, new);
//  	  }
//  	  // Get the change in energy
//  	  en = delta_en( sim, j, new);
//  	  // And do the Metropolis Hastings step
//  	  a = exp( -en/kt);
//  	  r = rng_uniform( sim.rng);
//  	  if( a>r ){
//  	    kd_move( sim.kdtree, j, new);
//  	  }
//  	}
//  	free(order):
// }