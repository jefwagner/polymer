/*!*******************************************************************
 * Monte-Carlo polymer simulation
 * ==============================
 * \author Jef Wagner
 * \email wagnerj@union.edu
 * \date 2019-06-15
 *********************************************************************
 */
/*!
 * This file contains the main routine including
 *  - input : reading the parameter file
 *  - initialization : creating and initializing the structures
 *  - evaluation : running the metropolis-hastings algorithm
 *  - output : intermittantly writing the results to an HDF5 file
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "polymer.h"
#include "energy.h"
#include "rng.c"

#define INTERNAL_STEPS 10

void usage( const char* argv0);
int get_cl_args( int argc, const char **argv,
                 char *parameter_filename, int *n);

int main( int argc, const char **argv){
  int n, status;
  char parameter_filename[80];
  char output_filename[95];

  status = get_cl_args( argc, argv, parameter_filename, &n);
  if( status != 0){
    usage( argv[0]);
    return status;
  }

  /******************************************************
   * Simulation Initialization
   */
  simulation sim;
  en_params enp;
  sim.enp = &enp;

  status = read_parameter_file( &sim, parameter_filename, n);
  if( status != 0){
    usage( argv[0]);
    return status;
  }

  sprintf( output_filename, "%s_%03d.f5", parameter_filename, n);
  file_handle *output = open_output_file( &sim, output_filename, sim.num_steps+1);

  fprintf( stdout, "Read from parameter file : %s\n", parameter_filename);
  fprintf( stdout, " Nm = %d\n", sim.Nm);
  fprintf( stdout, " num_steps = %d\n", sim.num_steps);
  fprintf( stdout, " seed = %x\n", sim.seed);
  fprintf( stdout, " k = %f\n", enp.k);
  fprintf( stdout, " a = %f\n", enp.a);
  fprintf( stdout, " ep = %f\n", enp.ep);
  fprintf( stdout, " sig = %f\n", enp.sig);

  sim.lim.r_min = r_min_calc( sim.enp);
  sim.lim.r_max = r_max_calc( sim.enp);
  sim.lim.r_co = r_co_calc( sim.enp);
  fprintf( stdout, "Limits from parameters :\n");
  fprintf( stdout, " r_min = %f\n", sim.lim.r_min);
  fprintf( stdout, " r_max = %f\n", sim.lim.r_max);
  fprintf( stdout, " r_co  = %f\n", sim.lim.r_co);
  fprintf( stdout, "\n");

  rng_seed(sim.rng, sim.seed);
  sim.poly = polymer_malloc( sim.Nm);
  sim.kdtree = kd_tree_malloc( sim.Nm);
  polymer_linear( sim.poly, sim.Nm, sim.rng);
  kd_tree_fill( sim.kdtree, sim.Nm, sim.poly->positions);

  uint16 max_depth = (uint16) (log2( (double) sim.Nm) + sqrt( (double) sim.Nm) + 1.0);

  /****************************************************************
   * Write the initial data
   */
  write_parameters( output, &sim);
  double en = energy( &sim);
  write_polymer_state( output, &sim, en);

  /****************************************************************
   * Run the simulation
   */
  int i, j;
  for( i=0; i<sim.num_steps; i++){
    // INTERNAL_STEPS steps per gibbs step
    for( j=0; j<INTERNAL_STEPS; j++){
      gibbs_step( &sim);
    }
    // Check if tree is too unbalanced
    if( sim.kdtree->max_depth > max_depth ){
      kd_tree_fill( sim.kdtree, sim.Nm, sim.poly->positions);
    }
    en = energy( &sim);
    write_polymer_state( output, &sim, en);
  }

  /****************************************************************
   * Close the simulation
   */
  close_output_file( output);
  kd_tree_free( sim.kdtree);
  polymer_free( sim.poly);
  return status;
}

void usage( const char* argv0){
  fprintf( stdout, "%s parameter_filename [simulation_number]\n", argv0);
  fprintf( stdout, "See the README for full details\n\n");
}

int get_cl_args( int argc, const char **argv,
                 char *parameter_filename, int *n){
  int status;
  switch( argc){
    case 2:
      strcpy( parameter_filename, argv[1]);
      *n = 0;
      status = 0;
      break;
    case 3:
      strcpy( parameter_filename, argv[1]);
      *n = atoi( argv[2]);
      status = 0;
      break;
    case 1:
      fprintf( stderr, "Must provide parmeter filename.\n");
      status = 1;
      break;
    default:
      fprintf( stderr, "Too many command line arguments.\n");
      status = 1;
      break;
  }
  return status;
}

  // fprintf( stdout, "Testing RNG: \n");
  // fprintf( stdout, "RNG state: %lx : %lx \n", sim.rng.s[0], sim.rng.s[1]);
  // int i;
  // for(i=0; i<10; i++){
  //   fprintf( stdout, "  int: %lu  double: %f \n", 
  //            rng_next(&sim.rng), rng_uniform(&sim.rng));
  // }

  // fprintf( stdout, "Testing polymer linking: \n");
  // int b, j=0;
  // fprintf( stdout, "%d->", j);
  // i=0;
  // b = sim.poly->list_bonds[NB_MAX*j+0];
  // j = linked_monomer( &(sim.poly->bonds[2*b]), j);
  // fprintf( stdout, "%d->", j);
  // for( i=1; i<sim.Nm-1; i++){
  //   b = sim.poly->list_bonds[NB_MAX*j+1];
  //   j = linked_monomer( &(sim.poly->bonds[2*b]), j);
  //   fprintf( stdout, "%d->", j);
  //   if( i%16==15){
  //     fprintf( stdout, "\n");
  //   };
  // }
  // fprintf( stdout,"\n\n");
  //
  // fprintf( stdout, "Testing kd tree implementation\n");
  // fprintf( stdout, "  writing initial tree to file.h5\n");

  // fprintf( stdout, "  Testing max-min\n");
  // kd_node *parent = NULL;
  // kd_node *node = kd_max( sim.kdtree, sim.kdtree->root, 0, &parent);
  // fprintf( stdout, "   max_monomer : %u \n", node->monomer);
  // fprintf( stdout, "   max_parent : %u \n", parent->monomer);
  // double x_max = sim.poly->positions[DIM*node->monomer+0];
  // fprintf( stdout, "  x_max = %f \n\n", x_max);

  // close_output_file(output);

  // fprintf( stdout, "  Testing node removal\n");
  // node = sim.kdtree->root->left;
  // parent = sim.kdtree->root;
  // node = kd_remove( sim.kdtree, node, &parent);
  // fprintf( stdout, "   node_removed : %u\n", node->monomer);
  // fprintf( stdout, "   parent of removed : %u\n", parent->monomer);
  // fprintf( stdout, "   writing removed tree to removed.h5\n");
  // output = open_output_file( &sim, "removed.h5", 1);
  // kdarray = (uint16 *) malloc( 4*sim.Nm*sizeof(uint16));
  // kd_tree_to_array( sim.kdtree, kdarray);
  // write_kdtree_data(output, kdarray);
  // free( kdarray);
  // close_output_file(output);

  // fprintf( stdout, "  Testing adding node\n");
  // kd_add( sim.kdtree, sim.kdtree->root, node);
  // fprintf( stdout, "   writing removed tree to added.h5\n");
  // output = open_output_file( &sim, "added.h5", 1);
  // kdarray = (uint16 *) malloc( 4*sim.Nm*sizeof(uint16));
  // kd_tree_to_array( sim.kdtree, kdarray);
  // write_kdtree_data(output, kdarray);
  // free( kdarray);
  // close_output_file(output);

  // uint16 neighbors[512];
  // uint16 num_neigh = 0;
  // double *pos = &sim.poly->positions[DIM*(sim.Nm/2)];
  // kd_ball_search( sim.kdtree, sim.kdtree->root, 
  //                 pos, 5.0, 
  //                 neighbors, &num_neigh);
  // fprintf( stdout, "Testing kd_ball_search\n");
  // fprintf( stdout, " found %lu neighbors in ball\n", num_neigh);
  // int i;
  // for( i=0; i<num_neigh; i++){
  //   fprintf( stdout, "  %lu \n", neighbors[i]);
  // }

  // fprintf( stdout, "\nTesting Distributions: \n");
  // fprintf( stdout, " Testing rand1 : output file.h5/positions[2]\n");
  // double p0[3] = {0.0, 0.0, 0.0};
  // double new[3];
  // int i;
  // memcpy( &sim.poly->positions[DIM*0], p0, DIM*sizeof(double));
  // for( i=1; i<sim.Nm; i++){
  //   rand1( sim.rng, sim.lim.r_min, sim.lim.r_max, p0, new);
  //   memcpy( &sim.poly->positions[DIM*i], new, DIM*sizeof(double));
  // }
  // write_position_data(output, sim.poly->positions);

  // fprintf( stdout, " Testing rand2 : output file.h5/positions[3]\n");
  // sim.lim.r_min = 0.8;
  // sim.lim.r_max = 1.2;
  // double p1[3] = {1.0/ROOT3, 1.0/ROOT3, 1.0/ROOT3};
  // memcpy( &sim.poly->positions[DIM*0], p0, DIM*sizeof(double));
  // memcpy( &sim.poly->positions[DIM*1], p1, DIM*sizeof(double));
  // for( i=2; i<sim.Nm; i++){
  //   rand2( sim.rng, sim.lim.r_min, sim.lim.r_max, p0, p1, new);
  //   memcpy( &sim.poly->positions[DIM*i], new, DIM*sizeof(double));
  // }
  // write_position_data(output, sim.poly->positions);

  // fprintf( stdout, " Testing rand2 : output file.h5/positions[4]\n");
  // sim.lim.r_min = 0.7;
  // sim.lim.r_max = 1.3;
  // for(i=0; i<3; i++){
  //   p1[i] *= 1.8;
  // }
  // memcpy( &sim.poly->positions[DIM*0], p0, DIM*sizeof(double));
  // memcpy( &sim.poly->positions[DIM*1], p1, DIM*sizeof(double));
  // for( i=2; i<sim.Nm; i++){
  //   rand2( sim.rng, sim.lim.r_min, sim.lim.r_max, p0, p1, new);
  //   memcpy( &sim.poly->positions[DIM*i], new, DIM*sizeof(double));
  // }
  // write_position_data(output, sim.poly->positions);

  // fprintf( stdout, " Testing rand2 : output file.h5/positions[5]\n");
  // sim.lim.r_min = 0.7;
  // sim.lim.r_max = 1.3;
  // for(i=0; i<3; i++){
  //   p1[i] *= 2.1/1.8;
  // }
  // memcpy( &sim.poly->positions[DIM*0], p0, DIM*sizeof(double));
  // memcpy( &sim.poly->positions[DIM*1], p1, DIM*sizeof(double));
  // for( i=2; i<sim.Nm; i++){
  //   rand2( sim.rng, sim.lim.r_min, sim.lim.r_max, p0, p1, new);
  //   memcpy( &sim.poly->positions[DIM*i], new, DIM*sizeof(double));
  // }
  // write_position_data(output, sim.poly->positions);
