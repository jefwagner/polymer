#include <stdlib.h>
#include <string.h>

//#include <stdio.h>

#include "polymer.h"
#include "vec3.h"
#include "rng.c"

polymer* polymer_malloc( uint16 Nm){
  polymer *poly = (polymer *) malloc( sizeof(polymer));
  poly->positions = (double *) malloc( DIM*Nm*sizeof(double));
  poly->num_bonds = (uint8 *) malloc( Nm*sizeof(uint8));
  poly->list_bonds = (uint16 *) malloc( NB_MAX*Nm*sizeof(uint16));
  poly->bonds = (uint16 *) malloc( 2*(Nm-1)*sizeof(uint16));
  if( poly->positions == NULL || 
      poly->num_bonds == NULL || 
      poly->list_bonds == NULL || 
      poly->bonds == NULL ){
    poly = NULL;
  }
  return poly;
}

void polymer_free(polymer *poly){
  free(poly->positions);
  free(poly->num_bonds);
  free(poly->list_bonds);
  free(poly->bonds);
  free(poly);
}

// Initialization for a linear polymer
#define OUTER_MAX 50
#define INNER_MAX 500
int polymer_linear(polymer *poly, uint16 Nm,
                   uint64_t *rng){
//  int count_outer = 0; 
  int count_inner = 0;
  int i, j;
  double dr, u, v, th, phi;
  double pos[3], dv[3], trial[3], *vp;
  bool monomer_overlap;

  /// Add the first monomer
  i=0;
  pos[0] = 0.0;
  pos[1] = 0.0;
  pos[2] = 0.0;
  memcpy(&(poly->positions[DIM*i]), pos, DIM*sizeof(double));
  poly->num_bonds[i] = 1;
  poly->list_bonds[NB_MAX*i+0] = i;
  // loop over the monomers in the middle of the polymer
  for(i=1; i<Nm; i++){
    // Find a point that doesn't overlap previous points
    count_inner = 0;
    do{
      // Generate a random point on the sphere
      dr = 0.8 + 0.4*rng_uniform(rng);
      u = rng_uniform(rng);
      v = rng_uniform(rng);
      th = acos(2.0*v-1.0);
      phi = TWOPI*u;
      dv[0] = dr*sin(th)*cos(phi);
      dv[1] = dr*sin(th)*sin(phi);
      dv[2] = dr*cos(th);
      // Check of it overlaps
      monomer_overlap = false;
      for( j=0; j<i; j++ ){
        vec_sum(pos,dv,trial);
        vp = &(poly->positions[DIM*j]);
        if( dist2(vp, trial) < 0.64 ){
          monomer_overlap = true;
          break;
        }
      }
      count_inner++;
      // repeat if it does
    }while( monomer_overlap && count_inner<=INNER_MAX);
    // Accept the position
    // if( count_inner > 5 ){
    //   fprintf( stdout, "%u . ", count_inner);
    // }
    memcpy(pos, trial, DIM*sizeof(double));
    // add the monomer
    memcpy(&(poly->positions[DIM*i]), pos, DIM*sizeof(double));
    poly->num_bonds[i] = 2;
    poly->list_bonds[NB_MAX*i+0] = i-1;
    poly->list_bonds[NB_MAX*i+1] = i;
    poly->bonds[2*(i-1)+0] = i-1;
    poly->bonds[2*(i-1)+1] = i;
  }
  // Fix the last monomer
  poly->num_bonds[Nm-1] = 1;

  return 0;
}
